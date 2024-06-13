#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// get corrected mutation burdens
process get_mut_burden_corrected {
  tag "${meta.id}"

  input:
  tuple val(meta),
        path(mut_burden),
        path(trint_subs_obs_corrected)
  
  output:
  tuple val(meta),
        path("${meta.id}.mut_burden_corrected.tsv"),
        path(trint_subs_obs_corrected)
  
  script:
  """
  cat ${mut_burden} |
  awk -F"\t" -v OFS="\t" '{if(\$1=="corrected") { print "${meta.id}",\$0}}' \
  > ${meta.id}.mut_burden_corrected.tsv
  """
}

process make_sigprofiler_matrix {
  
}

// define the workflow
workflow {

    // print summary of supplied parameters
    log.info paramsSummaryLog(workflow)
    
    // get metadata + bam paths
    Channel.fromPath(params.sample_sheet)
    | splitCsv(header: true)
    | map { row ->
        meta = [id: row.id, nanoseq_dir: row.nanoseq_dir]
        [meta, 
        file(row.nanoseq_dir + "/" + row.id + ".mut_burden.tsv", checkIfExists: true),
        file(row.nanoseq_dir + "/" + row.id + ".trint_subs_obs_corrected.tsv", checkIfExists: true)] 
    } 
    | set { ch_input }

    // get mutation burdens
    ch_input 
    | get_mut_burden_corrected

}