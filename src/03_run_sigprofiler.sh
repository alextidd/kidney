#!/bin/bash

# dirs
wd=/lustre/scratch126/casm/team154pc/at31/kidney
sp_dir=out/sigprofiler/
mkdir -p $sp_dir

# generate sample sheet for completed nanoseq outputs
(
  echo "id,nanoseq_dir"
  while read -r nanoseq_dir ; do
    id=$(echo $nanoseq_dir | cut -d/ -f4)
    echo "$id,$wd/$nanoseq_dir" 
  done < <(ls -d out/nanoseq/outNextflow/*/post/)
) | cat > $sp_dir/sample_sheet.csv

# run nf-sigprofiler
nextflow run nf-sigprofiler \
    --sample_sheet out/sigprofiler/sample_sheet.csv \
    --out_dir out/sigprofiler/ \
    -w nf-sigprofiler/work/