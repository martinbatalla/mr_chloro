#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process FASTP {
    tag "${sample_id}" // Print sample name to the terminal while running
    container 'quay.io/biocontainers/fastp:1.3.3--h43da1c4_0'

    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'copy', pattern: "*.{html,json}"
    
    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
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

process GETORGANELLE {
    tag "${sample_id}"
    container 'quay.io/biocontainers/getorganelle:1.7.7.1--pyhdfd78af_0'
    
    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'copy', pattern: "*.fasta"

    input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2)
    path ref_seed
    path go_config_dir

    output:
    tuple val(sample_id), path("${sample_id}_getorganelle_output"), emit: assembly_dir
    tuple val(sample_id), path("${sample_id}_getorganelle_output/*.complete.*.path_sequence.fasta"), emit: complete_paths, optional: true
    tuple val(sample_id), path("${sample_id}_getorganelle_output/*.scaffolds.*.path_sequence.fasta"), emit: scaffold_paths, optional: true

    script:
    """
    get_organelle_from_reads.py \
        -1 ${trimmed_r1} \
        -2 ${trimmed_r2} \
        -s ${ref_seed} \
        -F embplant_pt \
        -t ${task.cpus} \
        -o ${sample_id}_getorganelle_output \
        --config-dir ${go_config_dir} \
        --overwrite
    """
}


process STANDARDIZE{
    tag "${sample_id}"
    container 'martinbatalla/mr_chloro:v2.1'
    
    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(status), path(getorganelle_fasta)
    path ref_gb

    output:
    tuple val(sample_id), val(status), path("${sample_id}_${status}_cpDNA_raw.fasta"), emit: raw_fasta

    script:
    """
    std_modified.sh \
    -d . \
    -g ${ref_gb} \
    -o "${sample_id}_${status}_cpDNA_raw.fasta" \
    -p ${sample_id}
    """  
}

process EXTRACT{
    tag "${sample_id}"
    container 'martinbatalla/mr_chloro:v2.1'   

    input:
    tuple val(sample_id),  path(scaffolds)

    output:
    tuple val(sample_id), path("${sample_id}_pre_seedfile.fasta"), emit: pre_seed_file

    script:
    """
    awk '/^>/ { if(n) {print n"\\t"h"\\t"s}; h=\$0; s=""; n=0; next; } {s=s""\$0; n+=length(\$0)} END {if(n) print n"\\t"h"\\t"s}' "${scaffolds}" \
    | sort -nr | head -n 1 | cut -f2,3 | tr '\\t' '\\n' > "${sample_id}_pre_seedfile.fasta"

    """

}

process BEST_FASTA{
    tag "${sample_id}"
    container 'martinbatalla/mr_chloro:v2.1'

    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(multi_scaffs)
    path ref_gb
    
    output:
    tuple val(sample_id), path("${sample_id}_seed.fasta"), emit: seed_file

    script:
    """
    std_modified.sh \
    -d . \
    -g ${ref_gb} \
    -o "${sample_id}_seed.fasta" \
    -p ${sample_id}
    """

}


process ORIENT{
    tag "${sample_id}"
    container  'quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7'

    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(pre_seed)
    path ref_seed

    output:
    tuple val(sample_id), path("${sample_id}_seed.fasta"), emit: seed_file

    script:
    """
    SEED_STRAND=\$(blastn -query "${ref_seed}" -subject "${pre_seed}" -perc_identity 80 -evalue 0.1 -outfmt '6 length sstrand' | sort -nr | head -n 1 | awk '{print \$2}')
  
    if [ "\$SEED_STRAND" == "minus" ]; then
    echo "Seed is inverted. Reverse-complementing before NOVOPlasty..."
    
    # Save the fasta header
    head -1 "${pre_seed}" > "${sample_id}_seed.fasta"
    
    # Grab the sequence, reverse it, complement it, and append it
    sed 1d "${pre_seed}" | tr -d '\\n' | awk '{ for(i=length(\$0); i>0; i--) printf "%s", substr(\$0, i, 1); print "" }' \
    | tr 'ATCGatcg' 'TAGCtagc' >> "${sample_id}_seed.fasta"
    else
        echo "Seed is in correct direction. Ready for NOVOPlasty..."
        cp "${pre_seed}" "${sample_id}_seed.fasta"
    fi
    """
}

process NOVOPLASTY{
    tag "${sample_id}"

    container 'quay.io/biocontainers/novoplasty:4.3.5--pl5321hdfd78af_0'

    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(seed_file), path(trimmed_r1), path(trimmed_r2)
    path ref_seed

    output:
    tuple val(sample_id), path("Option_*_${sample_id}.fasta"), emit: novo_isomers, optional: true
    tuple val(sample_id), path("Circularized_assembly_*_${sample_id}.fasta"), emit: novo_complete, optional: true
    tuple val(sample_id), path("Contigs_*_${sample_id}.fasta"), emit: novo_scaffolds, optional: true

    script:
    """
    cat > config.txt <<EOF

    Project:
    -----------------------
    Project name          = ${sample_id}
    Type                  = chloro
    Genome Range          = 145000 - 170000
    K-mer                 = 39
    Max memory            = 64
    Extended log          = 0
    Save assembled reads  = no
    Seed Input            = ${seed_file}
    Extend seed directly  = yes
    Reference sequence    = ${ref_seed}
    Variance detection    = no
    Chloroplast sequence  = 

    Dataset 1:
    -----------------------
    Read Length           = 151
    Insert size           = 75
    Platform              = illumina
    Single/Paired         = PE
    Combined reads        = 
    Forward reads         = ${trimmed_r1}
    Reverse reads         = ${trimmed_r2}
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

    """
}

process BWA_MAP{
    tag "${sample_id}"

    container 'martinbatalla/mr_chloro:v2.1'

    input:
    tuple val(sample_id), val(status), path(raw_fasta), path(trimmed_r1), path(trimmed_r2) 

    output:
    tuple val(sample_id), val(status), path("${sample_id}_${status}_cpDNA.bam"), path("${sample_id}_${status}_cpDNA.bam.bai"), emit: cp_bam

    script:
    """
    echo "=== Mapping reads for ${sample_id} (${status}) ==="

    # Step 1: Index the assembly
    bwa-mem2 index "${raw_fasta}"

    # Step 2: Map reads and sort
    # We use 1G per thread for sorting to stay safe within your 64G RAM limit
    bwa-mem2 mem -t ${task.cpus} ${raw_fasta} ${trimmed_r1} ${trimmed_r2} | \
    samtools view -@ ${task.cpus} -Sb - | samtools sort -@ ${task.cpus} -m 1G -o "${sample_id}_${status}_cpDNA.bam" -

    # Step 3: Index the resulting BAM
    if [[ -f "${sample_id}_${status}_cpDNA.bam" ]]; then
        samtools index "${sample_id}_${status}_cpDNA.bam"
        echo "Successfully created ${sample_id}_${status}_cpDNA.bam"
    else
        echo "ERROR: Mapping failed."
        exit 1
    fi

    """
}

process PILON{
    tag "${sample_id}"
    container 'quay.io/biocontainers/pilon:1.24--hdfd78af_0'

    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(sample_id), val(status), path(raw_fasta), path(cp_bam), path(cp_bai)

    output:
    tuple val(sample_id), val(status), path("${sample_id}_${status}_cpDNA_polished.fasta"), path("${sample_id}_${status}_cpDNA_polished.changes"), emit: polished_fasta

    script:
    """
    echo "Step 5: polishing genome (${status})"
    export _JAVA_OPTIONS="-Xmx48g"
    pilon \
        --genome ${raw_fasta} \
        --frags ${cp_bam} \
        --output ${sample_id}_${status}_cpDNA_polished \
        --changes --threads ${task.cpus}

    """
}

process REMAP{
    tag "${sample_id}"
    container 'martinbatalla/mr_chloro:v2.1'

    input:
    tuple val(sample_id), val(status), path(polished_fasta), path(pilon_changes), path(trimmed_r1), path(trimmed_r2)

    output:
    tuple val(sample_id), val(status), path("${sample_id}_${status}_polished.bam"), path("${sample_id}_${status}_polished.bam.bai"), emit: polished_bam

    script:
    """
    echo "=== Remapping reads for polished ${sample_id} (${status}) fasta ==="

    # Step 1: Index the assembly
    bwa-mem2 index "${polished_fasta}"

    # Step 2: Map reads and sort
    bwa-mem2 mem -t ${task.cpus} ${polished_fasta} ${trimmed_r1} ${trimmed_r2} | \
    samtools view -@ ${task.cpus} -Sb - | samtools sort -@ ${task.cpus} -m 1G -o "${sample_id}_${status}_polished.bam" -

    # Step 3: Index the resulting BAM
    if [[ -f "${sample_id}_${status}_polished.bam" ]]; then
        samtools index "${sample_id}_${status}_polished.bam"
        echo "Successfully created ${sample_id}_${status}_polished.bam"
    else
        echo "ERROR: Re-mapping failed."
        exit 1
    fi

    """
}

process STATS{
    tag "${sample_id}"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(status), path(final_bam), path(final_bai)

    output:
    tuple val(sample_id), val(status), path("${sample_id}_${status}_coverage.txt"), emit: coverage_stats

    script:
    """
    samtools coverage ${final_bam} | tee ${sample_id}_${status}_coverage.txt
    """
}


workflow {
    log.info """\
         MR. CHLORO: ASSEMBLE CHLOROPLAST GENOMES
    =================================================
    Input Directory: ${params.input_dir}
    Reference Seed : ${params.ref_seed}
    Reference GB   : ${params.ref_gb}
    Output Dir     : ${params.out_dir}
    =================================================
    """
    .stripIndent()


    read_pairs_ch = Channel.fromFilePairs("${params.input_dir}/*_R{1,2}_001.fastq.gz", flat: true)

    FASTP(read_pairs_ch)

    GETORGANELLE(FASTP.out.trimmed_reads, file(params.ref_seed), "${projectDir}/go_config")

    // Organize getorganelle output
    scaffold_branch_ch = GETORGANELLE.out.scaffold_paths
        .branch {
            multiple: [it[1]].flatten().size() > 1
            single:   [it[1]].flatten().size() == 1
        }

    EXTRACT(scaffold_branch_ch.single)

    ORIENT(EXTRACT.out.pre_seed_file, file(params.ref_seed))

    BEST_FASTA(scaffold_branch_ch.multiple, file(params.ref_gb))
    
    orient_fastp_ch = ORIENT.out.seed_file.mix(BEST_FASTA.out.seed_file).join(FASTP.out.trimmed_reads)
    NOVOPLASTY(orient_fastp_ch, file(params.ref_seed))

    // Map all outputs into triplets: [sample_id, "status", fasta_path]
    def complete_assemblies = GETORGANELLE.out.complete_paths
        .mix(NOVOPLASTY.out.novo_complete, NOVOPLASTY.out.novo_isomers)
        .map { sample_id, fasta -> tuple(sample_id, "complete", fasta) }

    def draft_assemblies = GETORGANELLE.out.scaffold_paths
        .mix(NOVOPLASTY.out.novo_scaffolds)
        .map { sample_id, fasta -> tuple(sample_id, "draft", fasta) }

    // Mix all outputs together into one channel
    def all_tagged_assemblies = complete_assemblies.mix(draft_assemblies)

    // Group by sample_id, and only keep the "complete" assembly if it exists
    standardize_ch = all_tagged_assemblies
        .groupTuple(by: 0)
        .map { sample_id, statuses, fastas ->
            def complete_idx = statuses.indexOf("complete")
            if (complete_idx != -1) {
                // Return the complete genome
                return tuple(sample_id, "complete", fastas[complete_idx])
            } else {
                // If no complete genome exists, pass the draft scaffold forward
                return tuple(sample_id, "draft", fastas[0])
            }
        }

    STANDARDIZE(standardize_ch, file(params.ref_gb))

    // Update the joins to account for the new status string in the tuple
    bwa_channel = STANDARDIZE.out.raw_fasta.join(FASTP.out.trimmed_reads)
    BWA_MAP(bwa_channel)

    // Join on both sample_id and status (indices 0 and 1) to be safe
    pilon_channel = STANDARDIZE.out.raw_fasta.join(BWA_MAP.out.cp_bam, by: [0, 1])
    PILON(pilon_channel)

    // Check if pilon changed anything (index 3 is the .changes file)
    remap_ch = PILON.out.polished_fasta
        .branch {
            no_change: it[3].size() == 0
            change:   it[3].size() > 0
        }
    
    remap_channel = remap_ch.change.join(FASTP.out.trimmed_reads)
    REMAP(remap_channel)

    no_change_bams = remap_ch.no_change
        .join(BWA_MAP.out.cp_bam, by: [0, 1])
        .map { sample_id, status, polished_fasta, pilon_changes, orig_bam, orig_bai -> 
            tuple(sample_id, status, orig_bam, orig_bai) 
        }

    final_stats_ch = REMAP.out.polished_bam.mix(no_change_bams)

    STATS(final_stats_ch)
}