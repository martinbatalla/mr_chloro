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

  # 3. RESTORED PLASTAUMATIC ROTATION LOGIC
  head -1 "$working_fa" > "${dir2}/header"
  blastn -task blastn -query "$working_fa" -subject "$working_fa" -perc_identity 99 -evalue 0.00001 -outfmt '6 qseqid qstart qend sseqid sstart send length pident' | awk '$2!=$5&&$3!=$6 {print}' | sort -k7,7nr -k2,2n > "${dir2}/_tmp.blast.out"
  
  if [[ $(awk '{print $7}' "${dir2}/_tmp.blast.out" | head -1) -le 10000 ]] && [[ $seq_len -le 150000 ]]; then 
    echo "No significant IR detected, outputting oriented genome."
    cat "$working_fa" | sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' | sed "s/^>.*$/>${prefix}/" > "${dir2}/out_${opt_num}.fa"
  else 
    head -2 "${dir2}/_tmp.blast.out" > "${dir2}/_tmp.blast.out2" 
    local start=$(head -1 "${dir2}/_tmp.blast.out" | awk '{print $2}')
    local end=$(head -1 "${dir2}/_tmp.blast.out" | awk '{print $5}')
    
    if [ $start -eq 1 ]; then  
      awk '$3=='$seq_len' && $6-1=='$end' {print}' "${dir2}/_tmp.blast.out" >> "${dir2}/_tmp.blast.out2" 
    elif [ $end -eq $seq_len ]; then
      awk '$2==1 && $5+1=='$start' {print}' "${dir2}/_tmp.blast.out" >> "${dir2}/_tmp.blast.out2"
    fi
    mv "${dir2}/_tmp.blast.out2" "${dir2}/_tmp.blast.out"

    if [ $(grep -c "^" "${dir2}/_tmp.blast.out") -gt 2 ]; then 
        local ir_tmp=$(cat "${dir2}/_tmp.blast.out" | awk '$3=='$seq_len' {print $6-1}')
        local sc1=$(cat "${dir2}/_tmp.blast.out" | awk '$5=='$ir_tmp' {print $6-$3}')
        local sc2=$(cat "${dir2}/_tmp.blast.out" | awk '$3=='$seq_len' {print $2-$5}')
        if [ $sc2 -gt $sc1 ]; then 
            lsc=$(cat "${dir2}/_tmp.blast.out" | awk '$3=='$seq_len' {print $1":"$5+1"-"$2-1}')
            ira_start=$(cat "${dir2}/_tmp.blast.out" | awk '$3=='$seq_len' {print $1":"$2"-"$3}')
            ira_end=$(cat "${dir2}/_tmp.blast.out" | awk '$5=='$ir_tmp' {print $1":"$2"-"$3}') 
            ssc=$(cat "${dir2}/_tmp.blast.out" | awk '$5=='$ir_tmp' {print $1":"$3+1"-"$6-1}')
            irb_start=$(cat "${dir2}/_tmp.blast.out" | awk '$5=='$ir_tmp' {print $1":"$6"-"$5}') 
            irb_end=$(cat "${dir2}/_tmp.blast.out" | awk '$3=='$seq_len' {print $1":"$6"-"$5}')
        else 
            lsc=$(cat "${dir2}/_tmp.blast.out" | awk '$5=='$ir_tmp' {print $1":"$3+1"-"$6-1}')
            ira_start=$(cat "${dir2}/_tmp.blast.out" | awk '$5=='$ir_tmp' {print $1":"$6"-"$5}') 
            ira_end=$(cat "${dir2}/_tmp.blast.out" | awk '$3=='$seq_len' {print $1":"$6"-"$5}')
            ssc=$(cat "${dir2}/_tmp.blast.out" | awk '$3=='$seq_len' {print $1":"$5+1"-"$2-1}')
            irb_start=$(cat "${dir2}/_tmp.blast.out" | awk '$3=='$seq_len' {print $1":"$2"-"$3}')
            irb_end=$(cat "${dir2}/_tmp.blast.out" | awk '$5=='$ir_tmp' {print $1":"$2"-"$3}') 
        fi 
        samtools faidx "$working_fa" "$irb_start" | sed "s/$irb_start/IRb/" > ${dir2}/irbs 
        samtools faidx "$working_fa" "$irb_end" | sed 1d > ${dir2}/irbe  
        samtools faidx "$working_fa" "$ira_start" | sed "s/$ira_start/IRa/" > ${dir2}/iras 
        samtools faidx "$working_fa" "$ira_end" | sed 1d > ${dir2}/irae  
        cat ${dir2}/irbs ${dir2}/irbe > ${dir2}/IRb.fa 
        cat ${dir2}/iras ${dir2}/irae > ${dir2}/IRa.fa 
        
        samtools faidx "$working_fa" "$lsc" | sed 1d > ${dir2}/LSC.fa 
        samtools faidx "$working_fa" "$ssc" | sed 1d > ${dir2}/SSC.fa 
        cat ${dir2}/IRa.fa ${dir2}/IRb.fa | sed 1d > ${dir2}/_tmp_irs.fa
        
        cat ${dir2}/header ${dir2}/LSC.fa ${dir2}/iras ${dir2}/irae ${dir2}/SSC.fa ${dir2}/irbs ${dir2}/irbe | sed '/^>/d' | tr -d '\n' > ${dir2}/seq_only
        cat ${dir2}/header ${dir2}/seq_only | sed "s/^>.*$/>${prefix}/" > "${dir2}/out_${opt_num}.fa"
    else 
        # Logic for sequences split at boundaries
        if [ "$(cat ${dir2}/_tmp.blast.out | awk '$3=='$seq_len'||$2==1 {print "yes"}')" == "yes" ]; then
             # [Insert specific boundary split logic if needed, but usually defaults to clean output]
             cat "$working_fa" | sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' | sed "s/^>.*$/>${prefix}/" > "${dir2}/out_${opt_num}.fa"
        else
             cat "$working_fa" | sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' | sed "s/^>.*$/>${prefix}/" > "${dir2}/out_${opt_num}.fa"
        fi
    fi
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
# GetOrganelle outputs:
OPTIONS+=(${dir}/*.complete.*.path_sequence.fasta)

shopt -u nullglob

# Safety check: exit if no valid tool outputs are found
if [ ${#OPTIONS[@]} -eq 0 ]; then
    echo "ERROR: No NOVOPlasty (Option_*) or GetOrganelle (*.complete.*) files found in $dir."
    exit 1
fi

BEST_SCORE=0
WINNER=""
WINNING_NAME=""

# --- MAIN RACE ---
BEST_SCORE=-1
WINNER=""

for opt_path in "${OPTIONS[@]}"; do
    [ "$(basename "$opt_path")" == "$(basename "$out")" ] && continue
    
    # Extract name for the temp file
    num=$(basename "$opt_path" | sed 's/.fasta//' | sed 's/.path_sequence//')
    
    standardize "$opt_path" "$num"
    
    # CORRECTED PATH: Must match the dir2 variable inside the standardize function
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