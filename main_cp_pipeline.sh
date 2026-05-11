#!/bin/bash
# chloroplast_pipeline.sh
# Usage: ./chloroplast_pipeline.sh SAMPLE THREADS

SAMPLE=$1
THREADS=$2
REF_SEED="malva_cpgenome.fasta"

# Find R1/R2 files (handles -N, -NR, -acNR, etc.)
R1=$(ls ${SAMPLE}-*R1_001.fastq.gz ${SAMPLE}_R1_001.fastq.gz 2>/dev/null | head -n1)
R2=$(ls ${SAMPLE}-*R2_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz 2>/dev/null | head -n1)

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "? Could not find reads for $SAMPLE"
    exit 1
fi

# Make output directory
mkdir -p "$SAMPLE"
cd "$SAMPLE"

echo "=============================="
echo "Starting chloroplast pipeline"
echo "Sample: $SAMPLE"
echo "R1: $R1"
echo "R2: $R2"
echo "Reference seed: $REF_SEED"
echo "Threads: $THREADS"
echo "=============================="

##### Step 1: fastp #####
crun -p ~/envs/getorganelle_env conda run -n getorganelle fastp \
  -i ../"$R1" -I ../"$R2" \
  -o ${SAMPLE}_R1.trimmed.fastq.gz \
  -O ${SAMPLE}_R2.trimmed.fastq.gz \
  -h ${SAMPLE}_fastp.html \
  -j ${SAMPLE}_fastp.json \
  -q 15 -l 35 --thread $THREADS

echo "Step 1:fastp finished, starting Step 2: GetOrganelle"

##### Step 2: GetOrganelle #####
crun -p ~/envs/getorganelle_env conda run -n getorganelle get_organelle_from_reads.py \
  -1 ${SAMPLE}_R1.trimmed.fastq.gz \
  -2 ${SAMPLE}_R2.trimmed.fastq.gz \
  -s ../$REF_SEED \
  -F embplant_pt \
  -t $THREADS \
  -o getorganelle_output --overwrite \
  &> ${SAMPLE}_getorganelle.log

echo "Step 2:GetOrganelle finished, starting Step 3: Rename contigs"

##### Step 3: Rename contigs (use only complete assembly) ###################################



# Use arrays and nullglob to accurately count and store file paths
shopt -s nullglob
COMPLETE_PATHS=(./getorganelle_output/*.complete.*.path_sequence.fasta)
SCAFFOLD_PATHS=(./getorganelle_output/*.scaffolds.*.path_sequence.fasta)
shopt -u nullglob

ASSEMBLY_FILE=""

# Check if exactly two isomers were generated
if [[ ${#COMPLETE_PATHS[@]} -eq 2 ]]; then
    echo "? Two isomers found. Determining the standard orientation..."
    
    # Blast both against Malva and sum the 'plus' strand alignment lengths
    # Call blastn directly from its installation folder without activating the environment
    SCORE1=$(/home/mbata001/envs/miniconda3/envs/novoplasty_env/bin/blastn -query ../$REF_SEED -subject "${COMPLETE_PATHS[0]}" -perc_identity 90 -outfmt '6 sstrand length' | awk '$1=="plus" {sum+=$2} END {print sum+0}')
    SCORE2=$(/home/mbata001/envs/miniconda3/envs/novoplasty_env/bin/blastn -query ../$REF_SEED -subject "${COMPLETE_PATHS[1]}" -perc_identity 90 -outfmt '6 sstrand length' | awk '$1=="plus" {sum+=$2} END {print sum+0}')    
    # Pick the one with the highest forward-facing alignment
    if [ "$SCORE1" -gt "$SCORE2" ]; then
        echo "? Isomer 1 (${COMPLETE_PATHS[0]}) matches standard orientation."
        ASSEMBLY_FILE="${COMPLETE_PATHS[0]}"
    else
        echo "? Isomer 2 (${COMPLETE_PATHS[1]}) matches standard orientation."
        ASSEMBLY_FILE="${COMPLETE_PATHS[1]}"
    fi

# If only one complete assembly was found
elif [[ ${#COMPLETE_PATHS[@]} -eq 1 ]]; then
    ASSEMBLY_FILE="${COMPLETE_PATHS[0]}"
    echo "? Single complete assembly found."

# If no complete assembly, check scaffolds
elif [[ ${#SCAFFOLD_PATHS[@]} -ge 1 ]]; then
    SCAFFOLD_PATH="${SCAFFOLD_PATHS[0]}"
    # Count headers in the scaffold file
    SEQ_COUNT=$(grep -c "^>" "$SCAFFOLD_PATH")
    if [[ "$SEQ_COUNT" -eq 1 ]]; then
        ASSEMBLY_FILE="$SCAFFOLD_PATH"
        echo "? Found single scaffold in $SCAFFOLD_PATH. Proceeding..."
    else
        echo "? $SEQ_COUNT scaffolds found. Cannot proceed automatically."
    fi
fi


# Now check if we found a valid file to work with
if [[ -n "$ASSEMBLY_FILE" ]]; then
    # Step 3: Rename contigs
    awk '/^>/{print ">'"$SAMPLE"'_" ++i; next} {print}' "$ASSEMBLY_FILE" > ${SAMPLE}_cpDNA_raw.fasta
    cp ./getorganelle_output/get_org.log.txt get_org.log.txt
    rm -rf ./getorganelle_output 
else
    echo "? No complete assembly or single scaffold found for $SAMPLE. Quitting."
    exit 1
fi

###### Step 4: Map trimmed reads to raw assembly #####
../bwa_map.sh $SAMPLE ${SAMPLE}_cpDNA_raw.fasta ${SAMPLE}_R1.trimmed.fastq.gz ${SAMPLE}_R2.trimmed.fastq.gz $THREADS

echo "Step 4: Map trimmed reads to raw assembly finished, starting Step 5: pilon polishing"

##### Step 5: Pilon polishing #####
crun -p ~/envs/getorganelle_env conda run -n getorganelle pilon \
  --genome ${SAMPLE}_cpDNA_raw.fasta \
  --frags ${SAMPLE}_cpDNA.bam \
  --output ${SAMPLE}_cpDNA_polished \
  --changes \
  --threads $THREADS

echo "Step 5: pilon polishing finished, starting Step 6: Reuse or make new BAM"

##### Step 6: If Pilon made changes, remap. Otherwise reuse BAM #####
if [[ -s ${SAMPLE}_cpDNA_polished.changes ]]; then
    ../bwa_map.sh $SAMPLE ${SAMPLE}_cpDNA_polished.fasta ${SAMPLE}_R1.trimmed.fastq.gz ${SAMPLE}_R2.trimmed.fastq.gz $THREADS
else
    echo "? Pilon made no changes; reused previous BAM (no remapping)."
fi


mv ${SAMPLE}_cpDNA.bam ${SAMPLE}_cpDNA_polished.bam #even if pilon made changes, the bwa_map.sh file always outputs/rewrites {SAMPLE}_cpDNA.bam
mv ${SAMPLE}_cpDNA.bam.bai ${SAMPLE}_cpDNA_polished.bam.bai

echo "Step 6: Reuse or make new BAM finished, starting Step 6b: Ensure polished fasta exists"

##### Step 6b: Ensure polished fasta always exists #####
if [[ ! -f ${SAMPLE}_cpDNA_polished.fasta ]]; then
    cp ${SAMPLE}_cpDNA_raw.fasta ${SAMPLE}_cpDNA_polished.fasta
fi

echo "Step 6b: Ensure polished fasta exists finished, starting Step 7: Calculate coverage"

##### Step 7: Coverage #####
crun -p ~/envs/getorganelle_env conda run -n getorganelle samtools coverage ${SAMPLE}_cpDNA_polished.bam \
  | tee ${SAMPLE}_coverage.txt


echo "? Finished sample: $SAMPLE"
