from MFA2CSV import MFA2CSV
### ### ### ### ### ### ### ### ### ### ### ####
# virus reference + target genomes -> variants #
### ### ### ### ### ### ### ### ### ### ### ####

###
# install virulign
###

###
# run the provided bash for your specified (target) sequences
###

###
# convert MFA --> variants tabular file per target
###

###                    ###
###     SARS-Cov2      ###
###                    ###
# example for S (Spike) protein and 2 target sequences
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/data"
alignments_folder = "alignments"
xmls_folder = "xmls"
out_path = "/home/damian/Documents/L3S/projects/sars_cov2/variants"

mfa_file = "S_mfa_nt.fasta"
orf_xml_file = "S.xml"
ncbi_ref_id = "NC_045512.2"
print("=== ~ ===")
MFA2CSV_S_protein_test = MFA2CSV(data_path, alignments_folder, xmls_folder, ncbi_ref_id, out_path)
MFA2CSV_S_protein_test.run(orf_xml_file, mfa_file)
print("=== ~ ===")