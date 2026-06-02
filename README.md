# MrChloro: Chloroplast Assembly Pipeline

A Nextflow pipeline designed for the automated assembly, polishing, and standardization of chloroplast genomes (optimized for *Nototriche* (Malvaceae) and related genera). The pipeline is based on the assembly of genomes by GetOrganelle or, if failed, the extension and completion of genomes from GetOrganelle's output using NOVOPlasty; genomes are then standardized and polished automatically.

Fully containerized tool; no need to pre-download or install any files, tools, or dependencies other than those which are typically already installed in an HPC cluster (or Docker, if running locally or via AWS) See **Dependencies** section.

## Pipeline Architecture

1. **Quality Control:** Fastp to trim adapters out of reads
2. **Primary Assembly:** GetOrganelle for initial assembly attempt
3. **Secondary Assembly/Rescue:** If no complete genome is outputted from GetOrganelle, NOVOPlasty for extension of best scaffold produced by GetOrganelle (automatically inverted if orientation does not match the standard)
4. **Standardization:** Custom Python/Bash logic to rotate assemblies to a uniform starting position and select correct isomer (if isomers present)
5. **Polishing:** BWA-MEM2 mapping of trimmed reads to standardized assembly followed by Pilon consensus correction
6. **Re-mapping:** If pilon produced changed, BWA-MEM2 mapping of reads to the new, polished genome
7. **Statistics** Samtools to get basic statistics (e.g.: average depth of coverage.) of final, polished assembly

## Quick Start

You do not need to clone this repository manually. Nextflow will handle the download and execution automatically. 

**Run on a local machine (requires Docker):**
```
nextflow run martinbatalla/mr_chloro -profile standard --input_dir /path/to/reads --ref_seed /path/to/seed.fasta --ref_gb /path/to/reference.gb
```

**Run on an HPC cluster (requires Singularity & SLURM):**
```
nextflow run martinbatalla/mr_chloro -profile hpc --input_dir /path/to/reads --ref_seed /path/to/seed.fasta --ref_gb /path/to/reference.gb
```

## Dependencies

You do **not** need to install any bioinformatics tools (BWA, GetOrganelle, NOVOPlasty, etc.) locally. The pipeline is fully containerized. You only need:
* [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (version 23.10.0 or later)
* Docker (for local execution) OR Singularity (for HPC execution)
 If running via an HPC cluster, Nextflow and dependencies are usually already installed, although they may need to be loaded by user (e.g.: `module load nextflow`)

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
* **`${sample_id}_cpDNA_raw.fasta`:** Raw, pre-polished assembly of genome as outputted by GetOrganelle or NOVOPlasty
* **`${sample_id}_cpDNA_polished.fasta`:** Final, polished `.fasta` assembly
* **`${sample_id}_coverage.txt`:** A `.txt` with basic stats of assembly, including the average depth of coverage
If NOVOPlasty was run, some additions outputs may include:
* **`${sample_id}_seed.fasta`:** Best scaffold output by GetOrganelle used as a seed to extend in NOVOPlasty
* **`Circularized_assembly_[1-9]_${sample_id}.fasta`:** Circularized genome outputted by NOVOPlasty
* **`Contigs_[1-9]_${sample_id}.fasta`:** Contigs outputted by NOVOPlasty
* **`Option_[1-9]_${sample_id}.fasta`:** Circulized genomes outputted by NOVOPlasty. When multiple options are present, these mostly display different orientations of the SSC (MrChloro selects the best one).


## Citation

**MrChloro** is currently pre-publication. If you use this pipeline in your research, please cite this GitHub repository directly:

> Batalla, M. I. (2026). *MrChloro*: A Nextflow pipeline for automated chloroplast genome assembly. GitHub repository. https://github.com/martinbatalla/mr_chloro

## Acknowledgements & Credits

This pipeline automates and standardizes the workflows of several open-source bioinformatics tools. If you use `MrChloro`, please ensure you also cite the primary tools it utilizes:

* **Fastp:** Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*.
* **GetOrganelle:** Jin, J.-J., et al. (2020). GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. *Genome Biology*.
* **NOVOPlasty:** Dierckxsens, N., et al. (2017). NOVOPlasty: de novo assembly of organelle genomes from whole genome data. *Nucleic Acids Research*.
* **BLAST+:** Camacho, C., et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*.
* **BWA-MEM2:** Vasimuddin, M., et al. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. *IEEE IPDPS*.
* **Pilon:** Walker, B. J., et al. (2014). Pilon: an integrated tool for comprehensive microbial variant detection and genome assembly improvement. *PLoS One*.
* **SAMtools:** Li, H., et al. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*.