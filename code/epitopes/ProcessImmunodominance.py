from os.path import join
from code.epitopes.Protein import Protein
from os import scandir
from code.utils import is_fasta_file_extension


class ProcessImmunodominance:
	"""
	Class to calculate the immunodominance regions per amino acid position
	for a virus with experimentally identified epitopes based on the IEDB immunome browser

	Example:
	input: immunobrowser csv file
	output: sliding average of lower bound of 95% confidence interval (CI) of response frequency (RF)
	"""

	def __init__(self, data_path, exp_epitopes_folder, out_path):
		"""
		ProcessImmunodominance constructor

		Parameters
		----------
		data_path : str
			input data root
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
		with scandir(self.exp_epitopes_path) as exp_epitopes_dir:
			for epitopes_content in exp_epitopes_dir:
				if is_fasta_file_extension(epitopes_content.name) and epitopes_content.is_file():
					protein_rec = Protein(join(self.exp_epitopes_path, epitopes_content.name))
					if protein_rec.id not in self.proteins:
						self.proteins[protein_rec.id] = protein_rec

	def sum_per_position(self):
		"""
		Sum RF per amino acid position using the csv range and lower bound values
		:return:
		"""
