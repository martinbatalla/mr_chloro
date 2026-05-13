#!/bin/bash
# all_main_pipeline.sh
# Submit chloroplast pipeline jobs for all samples automatically

REF=$1

# Make sure logs directory exists
mkdir -p logs

# Make a text file with a list of all the samples
ls *_R1_001.fastq.gz | sed -E 's/([-_]?R1.*)//' > CPsamples.txt
NUM_SAMPLES=$(wc -l < CPsamples.txt)

echo "Created CPsamples.txt for $NUM_SAMPLES samples"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=cp_array
#SBATCH --array=1-${NUM_SAMPLES}%20
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Extract the specific sample for THIS task
SAMPLE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" CPsamples.txt)

# Create a link so in order to have a named log file
ln -snf "$(pwd)/logs/\${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}.out" "logs/\${SAMPLE}.log"


mkdir -p logs

../main_cp_pipeline.sh \$SAMPLE $REF 16
EOF
