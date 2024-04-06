#!/usr/bin/env nextflow

params.sample_sheet = "/lustre/scratch126/casm/team154pc/at31/kidney/data/lcm_patient_manifest.tsv"
params.out_dir = "out/"

/*
 * Merge and index bams
 */
process merge_bams {
    tag "${meta.sample_id}"
    label "long10gb"
    publishDir "${params.out_dir}/merged_normal_bams/", mode:"copy"

    input:
        tuple val(meta), path(bams)

    output:
        tuple val(meta), path("${meta.sample_id}.bam"), path("${meta.sample_id}.bam.bai")
    
    script:
        """
        module load samtools
        samtools merge -f \
            --threads ${task.cpus} \
            ${meta.sample_id}.bam \
            $bams
        samtools index ${meta.sample_id}.bam
        """
}

/*
 * Define the workflow
 */
workflow {
    
    // get metadata + bam paths
    Channel.fromPath(params.sample_sheet, checkIfExists: true)
    | splitCsv(header: true)
    | map { row ->
        meta = [donor_id: row.donor_id, sample_id: row.sample_id, 
                biopsy_id: row.biopsy_id, sample_tumour_status: row.sample_tumour_status]
        [meta, 
        file(row.biopsy_bam_file, checkIfExists: true)] 
    } 
    | set { ch_input }

    // get normal bams for merge
    ch_input 
    | branch { meta, biopsy_bam_file -> 
        tumour: meta.sample_tumour_status == "Tumour"
            return tuple( meta, biopsy_bam_file )
        normal: meta.sample_tumour_status == "Normal"
            return tuple( meta, biopsy_bam_file )}
    | set { ch_tumour_statuses }

    ch_tumour_statuses.normal 
    | map { meta, biopsy_bam_file ->
        meta = meta.subMap('donor_id', 'sample_id') 
        [meta, biopsy_bam_file] }
    | groupTuple 
    | map { meta, biopsy_bam_file -> [meta, biopsy_bam_file.flatten()] }
    | merge_bams
}