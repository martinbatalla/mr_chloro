# mr_assembly: Chloroplast Assembly Pipeline

A Nextflow pipeline designed for the automated assembly, polishing, and standardization of chloroplast genomes (optimized for *Nototriche* (Malvaceae) and related genera).

## Quick Start

You do not need to clone this repository manually. Nextflow will handle the download and execution automatically. 

**Run on a local machine (requires Docker):**
\`\`\`bash
nextflow run martinbatalla/mr_assembly -profile standard --input_dir /path/to/reads --ref_seed /path/to/seed.fasta --ref_gb /path/to/reference.gb
\`\`\`

**Run on an HPC cluster (requires Singularity & SLURM):**
\`\`\`bash
nextflow run martinbatalla/mr_assembly -profile hpc --input_dir /path/to/reads --ref_seed /path/to/seed.fasta --ref_gb /path/to/reference.gb
\`\`\`

## Dependencies

You do **not** need to install any bioinformatics tools (BWA, GetOrganelle, NOVOPlasty, etc.) locally. The pipeline is fully containerized. You only need:
* [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (version 23.10.0 or later)
* Docker (for local execution) OR Singularity (for HPC execution)

## Parameters

| Parameter | Description |
|-----------|-------------|
| `--input_dir` | Path to the directory containing paired-end fastq files (must end in `_R1_001.fastq.gz` and `_R2_001.fastq.gz`). |
| `--ref_seed` | Path to the `.fasta` seed file used for NOVOPlasty and GetOrganelle orientation. |
| `--ref_gb` | Path to the GenBank (`.gb`) reference file used for standardizing the starting position. |
| `--out_dir` | (Optional) Directory where results will be saved. Default: `results/`. |

## Pipeline Outputs

When the pipeline finishes, your output directory will contain a folder for each sample ID. Inside, you will find:
* **HTML/JSON quality reports:** Output of Fastp to check quality of trimming reads. Trimmed reads files are not saved to save space
* **`{sample_id}_cpDNA_polished.fasta`** Final, polished `.fasta` assembly
* **`${sample_id}_coverage.txt`** A `.txt` with basic stats of assembly, including the average depth of coverage


## Pipeline Architecture

1. **Quality Control:** Fastp
2. **Primary Assembly:** GetOrganelle 
3. **Secondary Assembly/Rescue:** NOVOPlasty (extends best scaffold produced by GetOrganelle)
4. **Standardization:** Custom Python/Bash logic to rotate assemblies to a uniform starting position and select correct isomer (if isomers present)
5. **Polishing:** BWA-MEM2 mapping followed by Pilon consensus correction
6. **Statistics** Samtools to get average depth of coverage of final, polished assembly