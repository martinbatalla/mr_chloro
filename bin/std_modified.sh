#!/bin/bash
# std_modified.sh

set -eo pipefail
while getopts 'd:g:o:p:' options
do
  case $options in
    d) dir=$OPTARG ;;
    g) gb=$OPTARG ;;
    o) out=$OPTARG ;;
    p) prefix=$OPTARG ;;
  esac
done

standardize() {
  local fasta_in=$1
  local opt_num=$2
  local dir2="${dir}/_tmp_plastaumatic"
  
  echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tStandardizing Option ${opt_num}..."

  # 1. Clean FASTA of invalid residues
  mkdir -p "$dir2"
  local clean_fasta="${dir2}/clean_${opt_num}.fasta"
  head -n 1 "$fasta_in" > "$clean_fasta"
  sed 1d "$fasta_in" | tr -cd 'ATCGatcgNn' >> "$clean_fasta"

  samtools faidx "$clean_fasta"
  local seq_len=$(cut -f2 "${clean_fasta}.fai")

  # 2. Orientation Check
  local ref_fa="${dir2}/ref_internal.fa"
  echo ">ref" > "$ref_fa"
  awk '/ORIGIN/,/\/\//' "$gb" | grep -v "ORIGIN" | grep -v "\//" | tr -d "[0-9] " | tr -d '\n' >> "$ref_fa"
  
  local STRAND=$(blastn -query "$ref_fa" -subject "$clean_fasta" -perc_identity 80 -evalue 0.1 -outfmt '6 length sstrand' | sort -nr | head -n 1 | awk '{print $2}')
  
  local working_fa="${dir2}/oriented_${opt_num}.fa"
  if [ "$STRAND" == "minus" ]; then
    echo -e "Option ${opt_num} is inverted. Flipping..."
    head -1 "$clean_fasta" > "$working_fa"
    sed 1d "$clean_fasta" | tr -cd 'ATCGatcgNn' | rev | tr 'ATCGatcg' 'TAGCtagc' >> "$working_fa"
  else
    cp "$clean_fasta" "$working_fa"
  fi

  # 3. ROTATE TO LSC BOUNDARY USING SELF-BLAST
  # Extract the coordinates of the Inverted Repeats (ignoring the 100% self-to-self hit)
  blastn -task blastn -query "$working_fa" -subject "$working_fa" -perc_identity 99 -evalue 0.00001 -outfmt '6 qstart qend sstart send' | awk '$1!=$3 {print}' | head -1 > "${dir2}/_tmp.blast.out"
  
  if [ -s "${dir2}/_tmp.blast.out" ] && [[ $seq_len -le 250000 ]]; then 
    local coords=$(cat "${dir2}/_tmp.blast.out")
    
    # Sort the 4 boundaries to get C1, C2, C3, and C4 sequentially along the circular genome
    local C1=$(echo "$coords" | tr -s '[:space:]' '\n' | sort -n | sed -n '1p')
    local C2=$(echo "$coords" | tr -s '[:space:]' '\n' | sort -n | sed -n '2p')
    local C3=$(echo "$coords" | tr -s '[:space:]' '\n' | sort -n | sed -n '3p')
    local C4=$(echo "$coords" | tr -s '[:space:]' '\n' | sort -n | sed -n '4p')

    # Calculate gap sizes (Gap A is internal, Gap B wraps around the file edges)
    local gapA=$(( C3 - C2 ))
    local gapB=$(( seq_len - C4 + C1 ))

    # The larger gap is the Large Single Copy (LSC)
    local shift_pos=1
    if [ $gapA -gt $gapB ]; then
        shift_pos=$C2
    else
        shift_pos=$C4
    fi
    
    echo -e "IRs identified. Shifting start position to LSC boundary at base ${shift_pos}..."
    
    # Strip the header to create a continuous string of nucleotides
    sed 1d "$working_fa" | tr -d '\n' > "${dir2}/seq_only"
    
    # Slice and stitch
    echo ">${prefix}" > "${dir2}/out_${opt_num}.fa"
    if [ "$shift_pos" -eq 1 ]; then
        # Already at boundary, no cut needed
        cat "${dir2}/seq_only" | fold -w 80 >> "${dir2}/out_${opt_num}.fa"
    else
        local chunkA=$(cut -c ${shift_pos}-${seq_len} "${dir2}/seq_only")
        local chunkB=$(cut -c 1-$((shift_pos - 1)) "${dir2}/seq_only")
        echo "${chunkA}${chunkB}" | fold -w 80 >> "${dir2}/out_${opt_num}.fa"
    fi

  else 
    echo "No significant IR detected, outputting oriented genome."
    cat "$working_fa" | sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' | sed "s/^>.*$/>${prefix}/" | fold -w 80 > "${dir2}/out_${opt_num}.fa"
  fi
}

# --- MAIN RACE ---
mkdir -p ./_race_tmp
echo ">ref" > ./_race_tmp/ref.fa

# Extract reference sequence from GenBank for scoring orientation
awk '/ORIGIN/,/\/\//' "$gb" | grep -v "ORIGIN" | grep -v "\//" | tr -d "[0-9] " | tr -d '\n' >> ./_race_tmp/ref.fa

shopt -s nullglob

# Explicitly target only known assembly outputs
OPTIONS=()
# NOVOPlasty outputs:
OPTIONS+=(${dir}/Option_*_${prefix}.fasta)
OPTIONS+=(${dir}/Circularized_assembly_*_${prefix}.fasta)
OPTIONS+=(${dir}/Contigs_*_${prefix}.fasta)
# GetOrganelle outputs:
OPTIONS+=(${dir}/*.complete.*.path_sequence.fasta)
OPTIONS+=(${dir}/*.scaffolds.*.path_sequence.fasta)
# User-provided raw fastas (for Polish-Only mode):
OPTIONS+=(${dir}/*_raw.fasta)

shopt -u nullglob

# Safety check: exit if no valid tool outputs are found
if [ ${#OPTIONS[@]} -eq 0 ]; then
    echo "ERROR: No NOVOPlasty (Option_*) or GetOrganelle (*.complete.*) files found in $dir."
    exit 1
fi

BEST_SCORE=-1
WINNER=""
WINNING_NAME=""

for opt_path in "${OPTIONS[@]}"; do
    [ "$(basename "$opt_path")" == "$(basename "$out")" ] && continue
    
    # Extract name for the temp file
    num=$(basename "$opt_path" | sed 's/.fasta//' | sed 's/.path_sequence//')
    
    standardize "$opt_path" "$num"
    
    CURRENT_OUT="${dir}/_tmp_plastaumatic/out_${num}.fa"
    
    if [ -f "$CURRENT_OUT" ]; then
        # Score standardized output
        SCORE=$(blastn -query ./_race_tmp/ref.fa -subject "$CURRENT_OUT" -perc_identity 80 -outfmt '6 sstrand length' | awk '{if($1=="plus") sum+=$2; else sum+=($2*0.1)} END {print sum+0}')
        
        echo -e "Candidate: $num \t Weighted Score: $SCORE"
        
        # Comparison logic
        is_better=$(awk -v s="$SCORE" -v b="$BEST_SCORE" 'BEGIN {print (s >= b ? 1 : 0)}')
        
        if [ "$is_better" -eq 1 ]; then
            BEST_SCORE=$SCORE
            WINNING_NAME=$num
            cp "$CURRENT_OUT" "./_race_tmp/winner.fasta"
            WINNER="./_race_tmp/winner.fasta"
        fi
    else
        echo "?? Error: Standardization did not create $CURRENT_OUT"
    fi
done

if [ -n "$WINNER" ]; then
    mv "$WINNER" "$out"
    echo "--------------------------------------------------------"
    echo "SUCCESS: Best assembly ($WINNING_NAME) selected."
    echo "Finalized output saved to: $out"
    echo "--------------------------------------------------------"
else
    echo "ERROR: Standardization failed for all candidates."
    exit 1
fi

# Cleanup
rm -rf ${dir}/_tmp_plastaumatic ./_race_tmp