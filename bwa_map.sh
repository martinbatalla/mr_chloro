#!/bin/bash
# bwa_map.sh

# Assign arguments to variables for clarity
SAMPLE=$1
FASTA=$2
R1=$3
R2=$4
THREADS=$5

echo "=== Mapping reads for ${SAMPLE} ==="

# Step 1: Index the assembly
bwa-mem2 index "$FASTA"

# Step 2: Map reads and sort
# We use 1G per thread for sorting to stay safe within your 64G RAM limit
bwa-mem2 mem -t "$THREADS" "$FASTA" "$R1" "$R2" | \
samtools view -@ "$THREADS" -Sb - | \
samtools sort -@ "$THREADS" -m 1G -o "${SAMPLE}_cpDNA.bam" -

# Step 3: Index the resulting BAM
if [[ -f "${SAMPLE}_cpDNA.bam" ]]; then
    samtools index "${SAMPLE}_cpDNA.bam"
    echo "Successfully created ${SAMPLE}_cpDNA.bam"
else
    echo "ERROR: Mapping failed."
    exit 1
fi