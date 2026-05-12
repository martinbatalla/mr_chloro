#!/bin/bash
# bwa_map.sh

module load container_env
module load conda

SAMPLE=$1
FASTA=$2
R1=$3
R2=$4
THREADS=$5

crun -p ~/envs/getorganelle_env /bin/bash <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate getorganelle

echo "Mapping reads for ${SAMPLE}..."

# Step 1: index
bwa-mem2 index "$FASTA"

# Step 2: map reads 
# Added -m 2G to samtools sort. 16 threads * 2G = 32G for sorting.
# This leaves plenty of room for bwa-mem2 (approx 15-20G) within your 64G limit.
# Reduce -m to 1G per thread (16G total for sorting)
# This leaves ~48G for bwa-mem2 and system overhead.
bwa-mem2 mem -t "$THREADS" "$FASTA" "$R1" "$R2" | \
samtools view -@ "$THREADS" -Sb - | \
samtools sort -@ "$THREADS" -m 1G -o "${SAMPLE}_cpDNA.bam" -

# Step 3: index BAM
if [[ -f "${SAMPLE}_cpDNA.bam" ]]; then
    samtools index "${SAMPLE}_cpDNA.bam"
    echo "Successfully created ${SAMPLE}_cpDNA.bam"
else
    echo "ERROR: Mapping failed."
    exit 1
fi
EOF