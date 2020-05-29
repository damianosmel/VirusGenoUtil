from os.path import join
from code.epitopes.Protein import Protein
from os import scandir
from code.utils import is_fasta_file_extension


class Immunodominance:
	"""
	Class to calculate the immunodominance regions per amino acid position
	for a virus with experimentally identified epitopes based on the IEDB immunome browser

	Reference: "Defining Immunodominant Regions within SARS-CoV Genome"
	from the paper: Grifoni, Alba, et al. "A sequence homology and bioinformatic approach can
	predict candidate targets for immune responses to SARS-CoV-2." Cell host & microbe (2020).

	Example:
	input: immunobrowser csv file
	output: sliding average of lower bound of 95% confidence interval (CI) of response frequency (RF)
	"""

	def __init__(self, data_path, proteins_folder, exp_epitopes_folder, out_path):
		"""
		Immunodominance constructor

		Parameters
		----------
		data_path : str
			input data root
		proteins_folder : str
			proteins folder
		exp_epitopes_folder : str
			experimental epitopes folder
		out_path : str
			output data root

		Returns
		-------
		None
		"""

		# self.protein_rec = protein_rec
		self.data_path = data_path
		self.proteins_path = join(data_path,proteins_folder)
		self.exp_epitopes_path = join(data_path, exp_epitopes_folder)
		# self.protein_immunobrowser_csv = join(data_path,protein_immunobrowser_csv)
		self.out_path = out_path
		# self.protein_rec.letter_annotations["lower_bound_RF"] = self.protein_rec.get_protein_length() * [0.0]
		self.proteins = {}

	def load_proteins(self):
		"""
		Load all proteins of virus with experimental epitopes

		Parameters
		----------

		Returns
		-------
		None
		"""
		with scandir(self.proteins_path) as proteins_dir:
			for proteins_content in proteins_dir:
				print(proteins_content.name)
				if is_fasta_file_extension(proteins_content.name) and proteins_content.is_file():
					protein = Protein(join(self.proteins_path, proteins_content.name))
					protein.load()
					if protein.get_short_id() not in self.proteins:
						print("Saving protein with id: {}".format(protein.get_short_id()))
						protein_rec = protein.get_record()
						self.proteins[protein.get_short_id()] = protein_rec
		print(self.proteins)

	def sum_per_position(self):
		"""
		Sum RF per amino acid position using the start-end positions and lower bound values

		Parameters
		----------

		Returns
		-------
		"""

	def slide_window_RF(self, window_size, save_plot):
		"""
		Slide window to define contiguous regions for lower bound CI of RF
		Reference: Supplementary of Grifoni et al.

		Parameters
		----------
		window_size : int
			sliding window size
		save_plot : bool
			create and save plot of sliding window values

		Returns
		-------
		"""

	def set_RF_zero(self, protein_rec):
		protein_rec.letter_annotations["lower_bound_RF"] = len(protein_rec.seq) * [0.0]
		return protein_rec

	def process_all_proteins(self, window_size, save_plot):
		"""
		Process Immunodominance regions for each loaded protein

		Parameters
		----------
		window_size : int
			sliding window size
		save_plot : bool
			create and save plot of sliding window values

		Returns
		-------
		"""

		self.load_proteins()
		for protein_id,protein_rec in self.proteins.items():
			self.proteins[protein_id] = self.set_RF_zero(protein_rec)
			#print(self.proteins[protein_id].letter_annotations["lower_bound_RF"][0:10])