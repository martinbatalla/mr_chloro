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

#Main scripts directory
BIN_DIR="/home/mbata001/bioinformatics/cp_assembly"

# Set default seed to the original input
FINAL_SEED=$2

# Check orientation of the seed against the reference
SEED_STRAND=$(blastn -query $3 -subject $2 -perc_identity 80 -evalue 0.1 -outfmt '6 length sstrand' | sort -nr | head -n 1 | awk '{print $2}')

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
K-mer                 = 75
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

STD_SCRIPT="${BIN_DIR}/std_modified.sh"
REF_GB="../malva.gb"

# Run standardize_cpDNA.sh
std_output="./${1}_cpDNA_raw.fasta"
bash "$STD_SCRIPT" -d . -g "$REF_GB" -o "$std_output" -p "$1"
