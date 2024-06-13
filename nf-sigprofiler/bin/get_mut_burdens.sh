#!/usr/bin/env bash

# Create file of mutation numbers and burdens from all samples in nanoseq project
# Usage: get_mut_burdens.sh <directory_with_summary_dirs> <output_file_name>

# Grab variables
input_dir="$1"
output_file="$2"

# write header
echo "sample	muts	total	burden	burden_lci	burden_uci" > $output_file

# grab trinuc corrected burdens from summary files and write to output
grep "corrected" summary*/results.mut_burden.tsv >> $output_file

# remove path from sample names
sed -i 's!/results.mut_burden.tsv:corrected!!g' $output_file
sed -i 's!summary_!!g' $output_file