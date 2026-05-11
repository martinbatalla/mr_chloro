#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

#novoplast_std.sh
#Usage: sbatch novoplast_std.sh SAMPLE_NAME /path/to/seed.fasta /path/to/malva.fasta /path/to/R1.fastq.gz /path/to/R2.fastq.gz

# $1 = sample name
# $2 = seed (getorganelle longest scaffold)
# $3 = ref sequence
# $4 = forward reads
# $5 = reverse reads

# $6 = insert size (figure an automatic way to get this from the json file). Now ignored and done automatically by the tool, until I have time to look into this.


# Load the module 
module load miniconda3 2>/dev/null || true

# This allows the shell to recognize the 'conda activate' command
source /home/mbata001/envs/miniconda3/etc/profile.d/conda.sh

# Activate the environment
conda activate novoplasty_env


# Set default seed to the original input
FINAL_SEED=$2

# Check orientation of the seed against the reference
SEED_STRAND=$(blastn -query $3 -subject $2 -perc_identity 80 -evalue 0.1 -outfmt '6 sstrand' | head -n 1)

if [ "$SEED_STRAND" == "minus" ]; then
    echo "Seed is inverted. Reverse-complementing before NOVOPlasty..."
    FINAL_SEED="./${1}_RC_seed.fasta"
    
    # Save the fasta header
    head -1 $2 > $FINAL_SEED
    
    # Grab the sequence, reverse it, complement it, and append it
    sed 1d $2 | tr -d '\n' | rev | tr 'ATCGatcg' 'TAGCtagc' >> $FINAL_SEED
fi


# Create config file

cat > config.txt <<EOF

Project:
-----------------------
Project name          = $1
Type                  = chloro
Genome Range          = 145000 - 170000
K-mer                 = 39
Max memory            = 64
Extended log          = 0
Save assembled reads  = no
Seed Input            = $FINAL_SEED
Extend seed directly  = yes
Reference sequence    = $3
Variance detection    = no
Chloroplast sequence  = 

Dataset 1:
-----------------------
Read Length           = 151
Insert size           = 75
Platform              = illumina
Single/Paired         = PE
Combined reads        = 
Forward reads         = $4
Reverse reads         = $5
Store Hash            =

Heteroplasmy:
-----------------------
MAF                   = 
HP exclude list       = 
PCR-free              = 

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = yes
Reduce ambigious N's  = 
Output path           = 
EOF

# Run novoplasty
NOVOPlasty.pl -c config.txt

# Run standardize_cpDNA.sh
std_output="./${1}_cpDNA_raw.fasta"
../standardize_cpDNA.sh -d . -g ../malva.gb -o $std_output -p $1

######################### EVENTUALLY COMBINE GETORGANELLE_ENV AND NOVOPLASTY_ENV SO ITS ALL IN ONE ENVIRONMENT #####################
conda deactivate

enable_lmod
module load container_env
module load conda

###### Step 4: Map trimmed reads to raw assembly #####
#Now check if we found a valid file to work with
if [[ -f "$std_output" ]]; then
    echo "Step 4: Map trimmed reads to raw assembly"
else
    echo "? No complete assembly or single circular scaffold found for $1. Quitting."
    exit 1
fi

#################### EVENTUALLY CHANGE THIS TO MAKE IT CUSTOMISABLE!! ############################
THREADS=16

../bwa_map.sh $1 $std_output $4 $5 $THREADS

echo "Step 4: Map trimmed reads to raw assembly finished, starting Step 5: pilon polishing"

##### Step 5: Pilon polishing #####
crun -p ~/envs/getorganelle_env conda run -n getorganelle pilon \
  --genome ${1}_cpDNA_raw.fasta \
  --frags ${1}_cpDNA.bam \
  --output ${1}_cpDNA_polished \
  --changes \
  --threads $THREADS

echo "Step 5: pilon polishing finished, starting Step 6: Reuse or make new BAM"

##### Step 6: If Pilon made changes, remap. Otherwise reuse BAM #####
if [[ -s ${1}_cpDNA_polished.changes ]]; then
    ../bwa_map.sh $1 ${1}_cpDNA_polished.fasta $4 $5 $THREADS
else
    echo "? Pilon made no changes; reused previous BAM (no remapping)."
fi


mv ${1}_cpDNA.bam ${1}_cpDNA_polished.bam #even if pilon made changes, the bwa_map.sh file always outputs/rewrites {1}_cpDNA.bam
mv ${1}_cpDNA.bam.bai ${1}_cpDNA_polished.bam.bai

echo "Step 6: Reuse or make new BAM finished, starting Step 6b: Ensure polished fasta exists"

##### Step 6b: Ensure polished fasta always exists #####
if [[ ! -f ${1}_cpDNA_polished.fasta ]]; then
    cp ${1}_cpDNA_raw.fasta ${1}_cpDNA_polished.fasta
fi

echo "Step 6b: Ensure polished fasta exists finished, starting Step 7: Calculate coverage"

##### Step 7: Coverage #####
crun -p ~/envs/getorganelle_env conda run -n getorganelle samtools coverage ${1}_cpDNA_polished.bam \
  | tee ${1}_coverage.txt


echo "? Finished sample: $1. Final, polished, standardized fasta saved as ${1}_cpDNA_polished.fasta"
