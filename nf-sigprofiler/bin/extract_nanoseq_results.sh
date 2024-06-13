#!/usr/bin/env bash

# Input a canapps project number to get the zipped summary files for nanoseq data, unzip, and rename (else they are all called 'summary' and will overwrite)
# Usage: extract_nanoseq_results.sh <canapps_project_number> <output_directory>

# Grab variables
project_number="$1"
output_dir="$2"

# Make output directory and subdirectories for the data, and for the observed and corrected mutation counts
mkdir "$output_dir"
cd "$output_dir"
mkdir all_trint_subs_obs_corrected

# Copy data from irods
cp /nfs/irods-cgp-*/intproj/"$project_number"/sample/PD*/*.v1.summary.tar.gz ./


for file in *.gz; 
do
	tar -xf "$file" # unzip file
	mv summary summary_"$(cut -d'.' -f1 <<<"$file")" # rename so they don't overwrite eachother
	cp summary_"$(cut -d'.' -f1 <<<"$file")"/results.trint_subs_obs_corrected.tsv  all_trint_subs_obs_corrected/ # copy mutation file to subdirectory
	mv all_trint_subs_obs_corrected/results.trint_subs_obs_corrected.tsv all_trint_subs_obs_corrected/"$(cut -d'.' -f1 <<<"$file")"_results.trint_subs_obs_corrected.tsv # rename so they don't overwrite eachother
done


#Insert R script here to make matrix when done