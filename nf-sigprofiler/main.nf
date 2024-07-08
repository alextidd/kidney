#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// get corrected mutation burdens
process get_mut_burden_corrected {
  tag "${meta.id}"
  label "normal"

  input:
  tuple val(meta),
        path(mut_burden),
        path(trint_subs_obs_corrected)
  
  output:
  tuple path("${meta.id}.mut_burden_corrected.tsv"),
        path(trint_subs_obs_corrected)
  
  script:
  """
  cat ${mut_burden} |
  awk -F"\t" -v OFS="\t" '{if(\$1=="corrected") { print "${meta.id}",\$0}}' \
  > ${meta.id}.mut_burden_corrected.tsv
  """
}

process make_sigprofiler_matrix {
  label "normal"

  input:
  tuple path(mut_burden_corrected),
        path(trint_subs_obs_corrected)

  output:
  tuple path(mut_burden_corrected),
        path("sigprofiler_matrix.tsv")

  script:
  """
  sigprofiler_matrix_maker_v2.py ./ sigprofiler_matrix.tsv
  """
}

process correct_wgs {
  label "normal"
  publishDir "${params.out_dir}/", mode: "copy", 
    pattern: "sigprofiler_matrix_wgs.tsv"

  input:
  tuple path(mut_burden_corrected),
        path(sigprofiler_matrix)
  
  output:
  path("sigprofiler_matrix_wgs.tsv")
  
  script:
  """
  cat *.mut_burden_corrected.tsv > mut_burden_corrected.tsv
  WGS_correction.R \
    mut_burden_corrected.tsv \
    sigprofiler_matrix.tsv \
    sigprofiler_matrix_wgs.tsv
  """
}

process run_sigprofiler {
  publishDir "${params.out_dir}/", mode: "copy"
  label "week100gb"
  errorStrategy = 'retry'
  maxRetries = 5

  input:
  path(sigprofiler_matrix_wgs)

  output:
  path("SBS96/*")
  path("Seeds.txt")

  script:
  """
  #!/usr/bin/env python
  from SigProfilerExtractor import sigpro as sig
  if __name__ == "__main__":
    sig.sigProfilerExtractor(
      input_type = "matrix", 
      output = "./", 
      input_data = "${sigprofiler_matrix_wgs}", 
      reference_genome = "${params.genome_build}",
      minimum_signatures = 1,
      maximum_signatures = 5,
      nmf_replicates = 100,
      cpu = -1)
  """
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
    | collect( flat: false )
    | map { it.transpose() }
    | make_sigprofiler_matrix
    | correct_wgs
    | run_sigprofiler
}