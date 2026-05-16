#!/bin/bash

# Header for your report
echo -e "Sample\tStatus\tScaffold_Count\tMax_Scaffold_Size"
echo "----------------------------------------------------------------"

while read SAMPLE; do
    # 1. Check if circularized polished genome exists
    # Removed the backslashes from ${SAMPLE}
    if [[ -f "${SAMPLE}/${SAMPLE}_cpDNA_polished.fasta" ]]; then
        SIZE=$(grep -v '^>' "${SAMPLE}/${SAMPLE}_cpDNA_polished.fasta" | tr -d '\n' | wc -c)
        echo -e "${SAMPLE}\tCIRCULARIZED\t1\t${SIZE}"
    
    else
        # 2. If not, check GetOrganelle scaffold output
        # Removed backslashes from $SCAF_FILE and $SAMPLE
        SCAF_FILE=$(ls ${SAMPLE}/getorganelle_output/*.scaffolds.*.path_sequence.fasta 2>/dev/null | head -n1)
        
        if [[ -f "$SCAF_FILE" ]]; then
            COUNT=$(grep -c '^>' "$SCAF_FILE")
            # Calculate max scaffold size
            MAX_SIZE=$(awk '/^>/ {if(n) print n; n=0; next} {n+=length($0)} END {print n}' "$SCAF_FILE" | sort -nr | head -n1)
            echo -e "${SAMPLE}\tSCAFFOLDS\t${COUNT}\t${MAX_SIZE}"
        else
            echo -e "${SAMPLE}\tFAILED\t0\t0"
        fi
    fi
done < CPsamples.txt