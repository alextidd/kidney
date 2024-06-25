#!/bin/bash
# cd /lustre/scratch126/casm/team154pc/at31/kidney; bsub -q week -M2000 -R "span[hosts=1] select[mem>2000] rusage[mem=2000]" -J sigprofiler -o log/sigprofiler.%J.out -e log/sigprofiler.%J.err "bash src/03_run_sigprofiler.sh"

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

# modules
module load python/3.12.3 
module load sigprofiler/1.1.23-GRCh37

# run nf-sigprofiler
nextflow run nf-sigprofiler \
    --sample_sheet out/sigprofiler/sample_sheet.csv \
    --out_dir out/sigprofiler/ \
    --genome_build GRCh37 \
    -w nf-sigprofiler/work/ \
    -resume \
    -N at31@sanger.ac.uk