#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# main_cp_pipeline.sh
# Usage: sbatch main_cp_pipeline.sh SAMPLE_NAME REF_NAME
# cwd should include ref.fasta, ref.gb and raw reads in the following format: SAMPLE_NAME_R1_001.fastq.gz and SAMPLE_NAME_R2_001.fastq.gz


SAMPLE=$1
THREADS=16
REF_SEED="${2}.fasta"  # The fasta seed
REF_GB="${2}.gb"       # The GenBank ref for standardization
BIN_DIR="/home/mbata001/bioinformatics/cp_assembly"

# Find Reads
R1=$(ls ${SAMPLE}-*R1_001.fastq.gz ${SAMPLE}_R1_001.fastq.gz 2>/dev/null | head -n1)
R2=$(ls ${SAMPLE}-*R2_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz 2>/dev/null | head -n1)

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "? Could not find reads for $SAMPLE"
    exit 1
fi

#Create sample folder
mkdir -p "$SAMPLE"
cd "$SAMPLE"

echo "=== Starting Pipeline for $SAMPLE ==="

#load containers
enable_lmod
module load container_env
module load conda

################################################# STEP 1: fastp #############################################################
echo "Step 1: running fastp"
crun -p ~/envs/getorganelle_env conda run -n getorganelle fastp \
    -i ../"$R1" -I ../"$R2" \
    -o ${SAMPLE}_R1.trimmed.fastq.gz -O ${SAMPLE}_R2.trimmed.fastq.gz \
    -h ${SAMPLE}_fastp.html -j ${SAMPLE}_fastp.json \
    -q 15 -l 35 --thread $THREADS

TRIM1="${SAMPLE}_R1.trimmed.fastq.gz"
TRIM2="${SAMPLE}_R2.trimmed.fastq.gz"

############################################ STEP 2: GetOrganelle ############################################################
echo "Step 2: running GetOrganelle"
crun -p ~/envs/getorganelle_env conda run -n getorganelle get_organelle_from_reads.py \
    -1 "$TRIM1" -2 "$TRIM2" -s ../"$REF_SEED" -F embplant_pt \
    -t $THREADS -o getorganelle_output --overwrite

################################ STEP 3: Handle Assembly & Standardization ###################################################
echo "Step 3: checking if multiple scaffolds/assemblies for extension and/or standardization"
#nullglob makes it so that the list (array) doesn't include the query by default if no matches are found
shopt -s nullglob
COMPLETE_PATHS=(./getorganelle_output/*.complete.*.path_sequence.fasta)
shopt -u nullglob

FINAL_RAW_FASTA="${SAMPLE}_cpDNA_raw.fasta"

if [[ ${#COMPLETE_PATHS[@]} -ge 1 ]]; then
    echo "? GetOrganelle successful. Standardizing..."
    # If 2 isomers exist, std_modified.sh handles the 'race' to find the best one
    bash "${BIN_DIR}/std_modified.sh" -d "./getorganelle_output" -g "../$REF_GB" -o "$FINAL_RAW_FASTA" -p "$SAMPLE"
    # std_modified.sh outputs ${SAMPLE}_cpDNA_raw.fasta
else
    echo "? GetOrganelle failed to circularize. Attempting NOVOPlasty extension..."
    
    # Isolate longest scaffold as seed
    GO_SCAFFOLD=$(ls ./getorganelle_output/*.scaffolds.*.path_sequence.fasta 2>/dev/null | head -n1)
    if [[ -f "$GO_SCAFFOLD" ]]; then
        SEED_FILE="${SAMPLE}_seed.fasta"
        # calculate length of each scaffold
        awk '/^>/ { if(n) {print n"\t"h"\t"s}; h=$0; s=""; n=0; next; } {s=s""$0; n+=length($0)} END {if(n) print n"\t"h"\t"s}' "$GO_SCAFFOLD" \
        # save longest scaffold to SEED_FILE
        | sort -nr | head -n 1 | cut -f2,3 | tr '\t' '\n' > "$SEED_FILE"
        
        # Run NOVOPlasty script (which includes standardization)
        bash "${BIN_DIR}/novoplast_std.sh" "$SAMPLE" "$SEED_FILE" "../$REF_SEED" "$TRIM1" "$TRIM2"
        # novoplast_std.sh outputs ${SAMPLE}_cpDNA_raw.fasta
    else
        echo "?? No scaffolds found. Assembly failed."
        exit 1
    fi
fi

# Safety check: did we get a standardized file (${SAMPLE}_cpDNA_raw.fasta)?
if [[ ! -f "$FINAL_RAW_FASTA" ]]; then
    echo "?? Standardization failed to produce $FINAL_RAW_FASTA"
    exit 1
fi

##################################################### STEP 4: Mapping #######################################################
echo "Step 4: mapping trimmed reads to standardized genome"
bash "${BIN_DIR}/bwa_map.sh" "$SAMPLE" "$FINAL_RAW_FASTA" "$TRIM1" "$TRIM2" $THREADS

################################################ STEP 5: Pilon Polishing ####################################################
echo "Step 5: polishing genome"
# Provide java with more memory, if not it will crash. 48g might be an overkill, CHECK!!
export _JAVA_OPTIONS="-Xmx48g"
crun -p ~/envs/getorganelle_env conda run -n getorganelle pilon \
    --genome "$FINAL_RAW_FASTA" \
    --frags ${SAMPLE}_cpDNA.bam \
    --output ${SAMPLE}_cpDNA_polished \
    --changes --threads $THREADS

############################################### STEP 6: Re-map if needed ###################################################
echo "Step 6: re-mapping (if polishing did changes)"
POLISHED_FASTA="${SAMPLE}_cpDNA_polished.fasta"
if [[ -s ${SAMPLE}_cpDNA_polished.changes ]]; then
    echo "? Pilon made changes. Re-mapping..."
    bash "${BIN_DIR}/bwa_map.sh" "$SAMPLE" "$POLISHED_FASTA" "$TRIM1" "$TRIM2" $THREADS
else
    echo "Polishing didn't change anything, re-using previous map"
    cp "$FINAL_RAW_FASTA" "$POLISHED_FASTA"
fi

mv ${SAMPLE}_cpDNA.bam ${SAMPLE}_cpDNA_polished.bam
mv ${SAMPLE}_cpDNA.bam.bai ${SAMPLE}_cpDNA_polished.bam.bai

############################################### STEP 7: Coverage Stats ####################################################
echo "Step 7: getting stats from assembly"
crun -p ~/envs/getorganelle_env conda run -n getorganelle samtools coverage ${SAMPLE}_cpDNA_polished.bam \
    | tee ${SAMPLE}_coverage.txt

echo "=== Pipeline Finished for $SAMPLE ==="