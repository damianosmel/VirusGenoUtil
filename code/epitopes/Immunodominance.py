from os.path import join
from code.epitopes.Protein import Protein
from os import scandir
from code.utils import is_fasta_file_extension
from pandas import read_csv,Series
import numpy as np
import matplotlib.pyplot as plt


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
		self.data_path = data_path
		self.proteins_path = join(data_path, proteins_folder)
		self.exp_epitopes_path = join(data_path, exp_epitopes_folder)
		self.out_path = out_path
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
					protein_id = proteins_content.name.split(".")[0]
					if protein_id not in self.proteins:
						print("Saving protein with id: {}".format(protein_id))
						protein_rec = protein.get_record()
						self.proteins[protein_id] = protein_rec
		print(self.proteins)

	def load_immunome_csv(self, protein_id):
		"""
		Load immunome csv

		Parameters
		----------

		Returns
		-------
		Pandas.DataFrame
			loaded immunome dataframe
		"""
		print("Load immunome csv")
		immunome_name = "RF_" + protein_id + ".csv"
		immunome_df = read_csv(join(self.exp_epitopes_path, immunome_name), sep=",", header=0)
		return immunome_df

	def sum_per_position(self, immunome, protein_id):
		"""
		Sum RF per amino acid position using the start-end positions and lower bound values

		Parameters
		----------
		immunome : Pandas.DataFrame
			immunome pandas dataframe
		protein_id : str
			id of protein currently processed

		Returns
		-------
		None
		"""
		print("Sum lower bound RF per epitope and amino-acid position")

		epitope_position = int(immunome["position"]) - 1  # 0-index
		epitope_lower_bound = float(immunome["lowerbound"])
		self.proteins[protein_id].letter_annotations["lower_bound_rf"][epitope_position] += epitope_lower_bound

	def compute_sliding_avg(self, protein_id, window_size, save_plot):
		"""
		Average over sliding window to define contiguous regions for lower bound CI of RF

		Reference: Supplementary of Grifoni et al.

		Credits: https://stackoverflow.com/questions/13728392/moving-average-or-running-mean

		Parameters
		----------
		protein_id : str
			id of protein to compute average sliding window
		window_size : int
			sliding window size
		save_plot : bool
			create and save plot of sliding window values

		Returns
		-------
		list of float
			sliding_average
		"""
		print("Compute sliding average RF")
		lower_bound_rf = np.asarray(self.proteins[protein_id].letter_annotations["lower_bound_rf"])
		#sliding_average = np.convolve(lower_bound_rf, np.ones((window_size,)) / float(window_size), mode='valid')
		sliding_average = Series(lower_bound_rf).rolling(window=window_size).mean().iloc[window_size - 1:].values
		assert len(self.proteins[protein_id].letter_annotations["lower_bound_rf"]) - sliding_average.shape[
			0] > 0, "AssertionError: lower bound RF values should be larger than sliding average values"
		sliding_average = np.concatenate((np.array(
			self.proteins[protein_id].letter_annotations["lower_bound_rf"][0:window_size - 1]), sliding_average),
			axis=0)
		if save_plot:
			fig = plt.figure()
			plt.plot(sliding_average,color='k')
			plt.axhline(y=0.3, color='k', linestyle='--')
			plt.ylim(0, 1)
			plt.xlim(0, len(self.proteins[protein_id].seq))
			plt.xlabel("Amino acid position")
			plt.ylabel("Lower bound\nsliding average")
			plt.title(self.proteins[protein_id].id)
			avg_plot_name = "sliding_avg_" + protein_id + ".png"
			fig.savefig(join(self.out_path, avg_plot_name), bbox_inches='tight', dpi=600)
			print("Saved {}".format(join(self.out_path, avg_plot_name)))
		return sliding_average.tolist()

	def set_RF_zero(self, protein_rec):
		"""
		Set protein sequence feature RF to 0
		across all its amino acid

		Parameters
		----------
		protein_rec :

		Returns
		-------
		Bio.SeqIO.SeqRecord
			protein with set lower
		"""
		print("Init lower bound RF to 0")
		protein_rec.letter_annotations["lower_bound_rf"] = len(protein_rec.seq) * [0.0]
		protein_rec.letter_annotations["sliding_avg_rf"] = len(protein_rec.seq) * [0.0]
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
		dict of str: Bio.SeqIO.SeqRecord
			processed proteins
		"""
		print("RF=0 -> load(RF,immunome.csv) -> slide window")
		self.load_proteins()
		for protein_id, protein_rec in self.proteins.items():
			print("--- ---")
			print("Process protein: {}".format(protein_id))
			self.proteins[protein_id] = self.set_RF_zero(protein_rec)
			# read immuno csv
			immunome_df = self.load_immunome_csv(protein_id)
			# sum lower bound RF per position
			immunome_df.apply(self.sum_per_position, protein_id=protein_id, axis=1)
			self.proteins[protein_id].letter_annotations["sliding_avg_rf"] = self.compute_sliding_avg(protein_id,
			                                                                                          window_size,
			                                                                                          save_plot)

		return self.proteins
