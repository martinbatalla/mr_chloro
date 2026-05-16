#!/bin/bash
# all_main_pipeline.sh

REF=$1
if [[ -z "$REF" ]]; then echo "Usage: path/to/all_main_pipeline.sh <REF_NAME>"; exit 1; fi

# 1. Setup sample list
ls *_R1_001.fastq.gz | sed -E 's/([-_]?R1_.*)//' > CPsamples.txt
NUM_SAMPLES=$(wc -l < CPsamples.txt)
mkdir -p logs

# 2. Submit the Array Job and capture the Job ID
JOB_ID=$(sbatch <<EOF | grep -oP "\d+"
#!/bin/bash
#SBATCH --job-name=cp_array
#SBATCH --array=1-${NUM_SAMPLES}%20
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

SAMPLE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" CPsamples.txt)
ln -snf "\$(pwd)/logs/\${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}.out" "logs/\${SAMPLE}.log"
#Prevent all n jobs starting at the same time
sleep \$((RANDOM % 150))
bash ~/bioinformatics/cp_assembly/main_cp_pipeline.sh "\$SAMPLE" "$REF"
EOF
)

echo "Array job submitted: $JOB_ID"