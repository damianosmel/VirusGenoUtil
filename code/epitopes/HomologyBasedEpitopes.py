from os.path import join, exists
from os import scandir
from code.utils import is_fasta_file_extension
from Bio import SearchIO


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

	def __init__(self, data_path, subject_proteins_folder, blast_folder, out_path):
		"""
		HomologyBasedEpitopes constructor

		Parameters
		----------
		data_path : str
			input data root
		subject_proteins_folder : str
			subject virus proteins folder
		blast_folder : str
			blast folder
		out_path : str
			output data root
		"""
		self.data_path = data_path
		self.subject_proteins_folder = join(data_path, subject_proteins_folder)
		self.blast_folder = join(data_path, blast_folder)
		self.out_path = out_path

	def blast_out2csv(self,blast_out_file):
		"""
		Convert a blast alignment (output format = 7, commented tabular) to csv file

		Credits: https://www.biostars.org/p/321102/

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
				epitopes_out.write("sars_cov_seq,sars_cov_max_RF,sars_cov2_seq,sars_cov2_protein_id, sars_cov2_alignment_pos,identity\n")
		with open(join(self.out_path,epitopes_name),"a") as epitopes_out:
			print("Update file: {}".format(epitopes_name))
			blast_alignments = SearchIO.parse(join(self.blast_folder,blast_out_file), 'blast-tab', comments=True)
			for blast_alignment in blast_alignments:
				for hsp in blast_alignment.hsps:
					print(hsp)
					print("---")
	def find_epitopes(self):
		"""
		Find homology based epitopes

		Parameters
		----------

		Returns
		-------

		"""
		with scandir(self.blast_folder) as blast_dir:
			for blast_content in blast_dir:
				print(blast_content.name)
				if blast_content.name[-3:] == "txt" and blast_content.is_file():
					print("Process blast output file: {}".format(blast_content.name))
					self.blast_out2csv(blast_content.name)
