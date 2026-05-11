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

# Step 1: index
bwa-mem2 index $FASTA

# Step 2: map reads
bwa-mem2 mem -t $THREADS $FASTA $R1 $R2 | \
samtools view -@ $THREADS -b | \
samtools sort -@ $THREADS -o ${SAMPLE}_cpDNA.bam

# Step 3: index BAM
samtools index ${SAMPLE}_cpDNA.bam
EOF
