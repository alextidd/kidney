# Makes sigprofiler-style mutational matrix from nanoseq "trint_subs_obs_corrected.tsv" output files
# Usage: python ./sigprofiler_matrix_maker.py <directory_containing_input_files> <output_path>
# Use within conda env "general" (or your own, that has pandas installed) 

## Import required packages ##
from sys import argv, exit # For parsing CLI options and exiting cleanly
from glob import glob # For grabbing all input files from directory
import pandas as pd # For handling dataframes
import os

## Check number of input arguments ##
if len(argv) != 3:
    print("Incorrect number of arguments supplied")
    print()
    print("Usage: python ./sigprofiler_matrix_maker.py <directory_containing_input_files> <output_path>")
    exit(1)

## Get CLI paths ##
input_directory = argv[1] #input_directory='/lustre/scratch119/casm/team267ms/yw2/kidney/Nanoseq_summary/all_trint_subs_obs_corrected'
output_file = argv[2]

## Initialise mutation matrix with first column ##
print("Initialising matrix...")
sbs96_classes = ['A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T', 'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T', 'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T', 'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T', 'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T', 'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T', 'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T', 'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T'] # Hardcoded SBS96 mutation classes

matrix = pd.DataFrame(sbs96_classes, columns=["MutationType"])

## Grab list of files from specified directory ##
print("Getting input files...")
files = sorted(glob(f"{input_directory}/*trint_subs_obs_corrected.tsv")) # Returns a sorted list of all .tsv files in input directory
# glob(os.path.join(input_directory,'/*.trint_subs_obs_corrected.tsv'))
# glob(input_directory + '/*.trint_subs_obs_corrected.tsv')

## Iterate through file list and add contents to mutation matrix ##
print("Generating mutation matrix")
for file in files: # For each file...
    sample_name = file.replace(f"{input_directory}/", "")
    sample_name = sample_name.replace(".trint_subs_obs_corrected.tsv", "") # store the sample name (i.e. no common text)
    print(f"Parsing sample {sample_name}")
    sample_data = pd.read_table(file).reset_index() # store the sample data (we need to reset index otherwise it won't merge with the mutation matrix properly)
    sample_data = sample_data["trint_onto_genome"].round().astype("int32") # grab the column we need, round it to nearest integer and store as int32
    matrix[sample_name] = sample_data # store our nicely formatted column in our full mutation matrix with sample_name as header
    print(f"Sample {sample_name} added to matrix")

## Write output file ##
print(f"All files parsed. Writing matrix to {output_file}")
matrix.to_csv(output_file, sep="\t", index=False)
print(f"Matrix written to {output_file}")

## Exit cleanly ##
exit(0)
