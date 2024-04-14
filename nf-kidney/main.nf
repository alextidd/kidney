#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// merge and index normal bams (for pseudobulk), create new sample ID (*b denoting normal)
process merge_normal_bams {
    tag "${meta.donor_id}"
    label "long10gb"
    publishDir "${params.out_dir}/merged_normal_bams/", mode:"copy"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.donor_id}_merged.bam"), path("${meta.donor_id}_merged.bam.bai")
    
    script:
    """
    module load samtools
    samtools merge -f \
        --threads ${task.cpus} \
        ${meta.donor_id}_merged.bam \
        ${bams}
    samtools index ${meta.donor_id}_merged.bam
    """
}

// define the workflow
workflow {

    // print help message, supply typical command line usage for the pipeline
    if (params.help) {
    log.info paramsHelp("nextflow run nf-chemo-trees --sample_sheet sample_sheet.csv")
    exit 0
    }

    // print summary of supplied parameters
    log.info paramsSummaryLog(workflow)
    
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
        meta = meta.subMap('donor_id') 
        [meta, biopsy_bam_file] }
    | groupTuple 
    | map { meta, biopsy_bam_file -> [meta, biopsy_bam_file.flatten()] }
    | merge_normal_bams

}