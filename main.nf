#!/usr/bin/env nextflow

// This enables the modern "DSL2" syntax for Nextflow
nextflow.enable.dsl=2

log.info """\
    C H L O R O P L A S T   A S S E M B L Y   P I P E L I N E
    ========================================================
    Input Directory: ${params.input_dir}
    Reference Seed : ${params.ref_seed}
    Reference GB   : ${params.ref_gb}
    Output Dir     : ${params.out_dir}
    ========================================================
    """
    .stripIndent()

process FASTP {
    tag "${sample_id}" // Print sample name to the terminal while running
    container 'https://depot.galaxyproject.org/singularity/fastp%3A1.3.3--h43da1c4_0'

    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', pattern: "*.{html,json}"
    
    input:
    // Take a sample ID and a pair of reads from the input channel
    tuple val(sample_id), path(r1), path(r2)

    output:
    // Bundle the sample ID and the new trimmed files into an output channel
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.html", emit: html_report
    path "${sample_id}_fastp.json", emit: json_report

    script:
    """
    fastp \
        -i ${r1} -I ${r2} \
        -o ${sample_id}_R1.trimmed.fastq.gz -O ${sample_id}_R2.trimmed.fastq.gz \
        -h ${sample_id}_fastp.html -j ${sample_id}_fastp.json \
        -q 15 -l 35 --thread ${task.cpus} --disable_trim_poly_g
    """
}

workflow {
    // Capture all paired-end sequencing files matching pattern
    // flat: true ensures R1 and R2 are passed as an array [R1, R2]
    read_pairs_ch = Channel.fromFilePairs("${params.input_dir}/*_R{1,2}_001.fastq.gz", flat: true)

    // Feed the channel of raw reads directly into the FASTP process
    FASTP(read_pairs_ch)
}