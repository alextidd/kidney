#!/bin/bash
# cd /lustre/scratch126/casm/team154pc/at31/kidney; ~/bin/jsub lsf -q week -n nanoseq -m 2g -l log "bash src/02_run_nanoseq.sh" | bsub

# load singularity
module load singularity

# exclude samples that have already run to completion from sample sheet
ls out/nanoseq/outNextflow/*/post/*.indel.vcf.gz | cut -d/ -f4 | 
grep -v -f - out/nanoseq/sample_sheet.csv \
> out/nanoseq/sample_sheet.csv.tmp

# run nanoseq pipeline
# increase maximum normal VAF (var_v) because query variants will be in the matched normal (due to pseudobulking)
# increase coverage (var_z) because matched normal is huge (was 12 before)
nextflow run ./NanoSeq_develop/Nextflow/NanoSeq_main.nf  \
  --jobs 200 -qs 300 -profile lsf_singularity \
  -w work/ \
  --sample_sheet out/nanoseq/sample_sheet.csv.tmp \
  --remap false \
  --grch37 true  \
  --dsa_d 2 \
  --cov_Q 15 \
  --var_b 0 \
  --var_n 3  \
  --var_a 50 \
  --var_d 2 \
  --var_r 144 \
  --var_x 8 \
  --indel_rb 2 \
  --var_q 40 \
  --var_v 0.05 \
  --var_z 25 \
  --outDir out/nanoseq/ \
  -c config/resources.config \
  -resume \
  -N at31@sanger.ac.uk
