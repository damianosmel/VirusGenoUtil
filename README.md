# VirusGenoUtil
This is a project to contribute utilities for a viral genomic data explorer.

## Dependencies
* Python 3.7.6
* BioPython 1.74

## Tested functionalities
1. Convert nucleic-acid sequence ORF multiple sequence alignemnts, provived by virulign aligner, saved as multi-fasta alignments (MFA) to variants per target sequence

2. Convert amino-acid sequence ORF multiple sequence alignemnts, provived by virulign aligner, saved as multi-fasta alignments (MFA) to variants per target sequence

## How to use
0. install virulign (on a server if you need to run for many sequences)
1. Depending on the sequence **type**:
* for amino-acid sequence:
  check the example [virulign_msa_aa_run.sh](bash/virulign_msa_aa_run.sh) used for running virulign for all ncbi sars cov2 as of May 15, 
* change this bash appropriately and run it for your own paths

2. put all output MFA fasta (one per ORF) to an alignments folder

3. go to [main.py](code/main.py) and change appropriately for the input and output paths,
* **input** should be the path of such alignment folder
* please change the refernce ncbi id
* create a file that contains your email (will be needed for programmatically accessing the Entrez service and fetch the genbank file of the reference sequence)
* **output** will be a folder containing a variant csv file for all ORFs of one target sequence, for example see [output](test/output/test_variants/)
