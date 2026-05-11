#!/bin/bash

# 1. Define the main working directory and reference
# Using $(pwd) ensures we use absolute paths
MAIN_DIR=$(pwd)
REF_SEQ="${MAIN_DIR}/malva_cpgenome.fasta"
SLURM_SCRIPT="${MAIN_DIR}/novoplast_std.sh"

# 2. Loop over every directory in the current folder
for sample_dir in */; do
    
    # Strip the trailing slash to get the clean sample name (e.g., "B124N")
    SAMPLE_NAME=${sample_dir%/}
    
    # Define the expected paths for this specific sample
    GO_DIR="${MAIN_DIR}/${SAMPLE_NAME}/getorganelle_output"
    R1_FILE="${MAIN_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_R1.trimmed.fastq.gz"
    R2_FILE="${MAIN_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_R2.trimmed.fastq.gz"
    
    # 3. Check if the getorganelle_output folder exists
    if [ -d "$GO_DIR" ]; then
        echo "========================================"
        echo "Processing Sample: $SAMPLE_NAME"
        
        # Define the path to the expected GetOrganelle fasta file
        GO_FASTA="${GO_DIR}/embplant_pt.K115.scaffolds.graph1.1.path_sequence.fasta"
        
        # Check if that specific fasta file actually generated successfully
        if [ -f "$GO_FASTA" ]; then

            # The "|| continue" is a safety net: if the cd fails for some reason, it skips to the next plant.
            cd "${MAIN_DIR}/${SAMPLE_NAME}" || continue
            
            # 4. Count the number of sequences
            SEQ_COUNT=$(grep -c "^>" "$GO_FASTA")
            echo "Found $SEQ_COUNT sequence(s) in the GetOrganelle output."
            
            # 5. Extract the longest sequence to use as the seed
            SEED_FILE="${MAIN_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_longest_seed.fasta"
            
            # The awk command flattens the FASTA, counts characters, sorts numerically, grabs the top 1, and reformats it back to FASTA
            awk '/^>/ { if(n) {print n"\t"h"\t"s}; h=$0; s=""; n=0; next; } {s=s""$0; n+=length($0)} END {if(n) print n"\t"h"\t"s}' "$GO_FASTA" \
            | sort -nr \
            | head -n 1 \
            | cut -f2,3 \
            | tr '\t' '\n' > "$SEED_FILE"
            
            echo "Longest sequence isolated and saved to $SEED_FILE"
            
            # 6. Submit the SLURM job
            echo "Submitting NOVOPlasty job to the cluster..."
            sbatch --job-name="$SAMPLE_NAME" --output="${SAMPLE_NAME}_novo.log" "$SLURM_SCRIPT" "$SAMPLE_NAME" "$SEED_FILE" "$REF_SEQ" "$R1_FILE" "$R2_FILE"

            #RETURN TO MAIN DIRECTORY FOR THE NEXT LOOP
            cd "$MAIN_DIR"

        else
            echo "Warning: getorganelle_output exists, but the FASTA file is missing."
        fi
    else
        # Optional: You can comment this out if you don't want it cluttering your terminal
        echo "Skipping $SAMPLE_NAME: No getorganelle_output folder found."
    fi
done
echo "========================================"
echo "Batch submission complete!"