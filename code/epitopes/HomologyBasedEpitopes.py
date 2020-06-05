from os.path import join, exists
from os import scandir
from code.utils import is_fasta_file_extension

from pandas import read_csv

class HomologyBasedEpitopes:
	"""
	Class to extract epitopes for subject virus sequence
	based on homology searching of experimental identified immunodominant regions

	Reference: "Defining Immunodominant Regions within SARS-CoV Genome"
	from the paper: Grifoni, Alba, et al. "A sequence homology and bioinformatic approach can
	predict candidate targets for immune responses to SARS-CoV-2." Cell host & microbe (2020).

	Example:
	input: BLAST file of immunodominant regions of sars cov1 to sars cov2
	output: file containing the aligned sequences, response frequency and identity
	"""

	def __init__(self, ncbi_ids, data_path, subject_proteins_folder, blast_folder, out_path):
		"""
		HomologyBasedEpitopes constructor

		Parameters
		----------
		ncbi_ids : dict of str: str
			dictionary to keep the NCBI ids for target organism and its close-related
		data_path : str
			input data root
		subject_proteins_folder : str
			subject virus proteins folder
		blast_folder : str
			blast folder
		out_path : str
			output data root
		"""
		self.ncbi_ids = ncbi_ids
		self.data_path = data_path
		self.subject_proteins_folder = join(data_path, subject_proteins_folder)
		self.blast_folder = join(data_path, blast_folder)
		self.out_path = out_path

	@staticmethod
	def process_relative_protein_info(protein_info):
		"""
		Process close relative virus protein info tag to get
		protein id, start end of immunodominant region, max RF (response frequency)

		Parameters
		----------
		protein_info : str

		Returns
		-------
		protein_id : str
			relative virus protein id
		protein_start : str
			start of immunodominant region
		protein_end : str
			end of immunodominant region
		protein_max_RF : str
			immunodominant region maximum RF
		"""
		protein_info_parts = protein_info.split(",")
		protein_id = protein_info_parts[0].split("_")[0]
		protein_start_end = protein_info_parts[1].split("=")[1].split("-")
		protein_start,protein_end = protein_start_end[0],protein_start_end[1]
		protein_max_RF = protein_info_parts[2].split("=")[1]

		return protein_id,protein_start,protein_end,protein_max_RF

	def blast_out2csv(self,blast_out_file):
		"""
		Convert a blast alignment (output format = 7, commented tabular) to csv file

		Parameters
		----------
		blast_out_file : str
			blast output filename

		Returns
		-------
		None
		"""
		print("Convert blast output file to homology based epitopes csv")
		epitopes_name = "homology_based_epitopes.csv"
		if not exists(join(self.out_path, epitopes_name)):
			print("Create file: {}".format(epitopes_name))
			with open(join(self.out_path, epitopes_name), "w") as epitopes_out:
				epitopes_out.write("relative_organism,relative_prot_id,relative_epi_start,relative_epi_end,relative_epi_sequence,relative_max_RF,target_organism,target_prot_id,target_epi_start,target_epi_end,target_epi_sequence,blast_identity,prediction_method\n")

		with open(join(self.out_path,epitopes_name),"a") as epitopes_out:
			print("Update file: {}".format(epitopes_name))
			blast_alignments_df = read_csv(join(self.blast_folder,blast_out_file),sep="\t",header=0)
			for index, blast_alignment in blast_alignments_df.iterrows():
				relative_prot_id,relative_epi_start,relative_epi_end,relative_max_RF = HomologyBasedEpitopes.process_relative_protein_info(blast_alignment["qseqid"])
				target_prot_id = blast_alignment["sseqid"].split("|")[3]
				target_start, target_end = str(blast_alignment["sstart"]), str(blast_alignment["send"])
				blast_identity = str(blast_alignment["pident"])
				print("Write homology based epitope")
				epitope_row = ",".join([self.ncbi_ids["relative_organism"],relative_prot_id,relative_epi_start,relative_epi_end,"unknown_seq",relative_max_RF,self.ncbi_ids["target_organism"],target_prot_id,target_start,target_end,"unknown_seq",blast_identity,"homology"])
				epitopes_out.write(epitope_row+"\n")

	def find_epitopes(self):
		"""
		Find homology based epitopes

		Parameters
		----------

		Returns
		-------
		None
		"""
		with scandir(self.blast_folder) as blast_dir:
			for blast_content in blast_dir:
				if blast_content.name[-3:] == "tsv" and blast_content.is_file():
					print("Process blast output file: {}".format(blast_content.name))
					self.blast_out2csv(blast_content.name)
