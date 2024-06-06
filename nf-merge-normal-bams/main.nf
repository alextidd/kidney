#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// merge and index normal bams (for pseudobulk), create new sample ID (*b denoting normal)
process merge_normal_bams {
    tag "${meta.donor_id}"
    label "week16core10gb"
    publishDir "${params.out_dir}/merged_normal_bams/", mode:"copy"
    conda '/nfs/users/nfs_a/at31/miniforge3/envs/samtools'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.donor_id}_merged.bam"), path("${meta.donor_id}_merged.bam.bai")
    
    script:
    """
    samtools merge -f \
        --threads ${task.cpus} \
        ${meta.donor_id}_merged.bam \
        ${bams}
    samtools index ${meta.donor_id}_merged.bam
    """
}

process nextflow {
    tag "${meta.donor_id}"
    label "week16core10gb"
    publishDir "${params.out_dir}/${meta.donor_id}/", mode:"copy"

    input:
    tuple val(meta), path(query_bam), path(normal_bam)

    script:
    """
    module purge
    module add bcftools
    module add R/3.6.1 # or choose the one under which you have all your packages installed
    module add samtools
    #module add singularity # I disable this to make nextflow use my own NanoSeq version. If not, it takes the info from the nextflow file and downloads all that is needed automatically
    module add nextflow

    # Activate access to IRODs:
    iinit 

    cd /lustre/scratch125/casm/team268im/fa8/117/CONTAMINATION_TEST_48390
    export PATH=/lustre/scratch125/casm/team268im/fa8/117/TWINSUK/PLATE14/DEBUG_INDEL_ERROR/NanoSeq-hotfix-3.5.5/bin:$PATH
    nextflow run ./NanoSeq-develop/Nextflow/NanoSeq_main.nf  \
        --jobs 200 -qs 20000 -profile lsf \
        --remap false \
        --grch37 true \
        --dsa_d 2 \
        --cov_Q 15 --var_b 0 --var_n 3 --var_z 12 -resume --study 6050  --sample_sheet plate.csv \
        --var_v 0.01 --var_a 50 --var_d 2 \
        --var_r 144 --var_x 8 \
        --indel_rb 2 \
        --var_q 40 \
        --outDir .
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

    // 

}