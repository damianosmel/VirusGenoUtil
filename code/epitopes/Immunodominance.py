from os.path import join
from code.epitopes.Protein import Protein
from os import scandir
from code.utils import is_fasta_file_extension
from pandas import read_csv
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
		self.proteins_path = join(data_path,proteins_folder)
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

	def load_immunome_csv(self,protein_id):
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
		# get 0-index of epitope start position
		epitope_start = int(immunome["Mapped Start Position"]) - 1
		if epitope_start < 0:
			epitope_start = 0
		epitope_end = int(immunome["Mapped End Position"])
		epitope_lower_bound = float(immunome["Lower Bound of 95% CI"])
		for amino_acid_pos in range(epitope_start,epitope_end):
			self.proteins[protein_id].letter_annotations["lower_bound_rf"][amino_acid_pos] += epitope_lower_bound
			if epitope_lower_bound > 0:
				self.proteins[protein_id].letter_annotations["found_epitopes_count"][amino_acid_pos] += 1.0

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
		"""
		print("Compute sliding average RF")
		lower_bound_rf = np.asarray(self.proteins[protein_id].letter_annotations["sliding_avg_rf"])
		average_sliding_window = np.convolve(lower_bound_rf, np.ones((window_size,)) / window_size, mode='valid')
		if save_plot:
			fig = plt.figure()
			plt.plot(average_sliding_window)
			plt.axhline(y=0.3, color='k', linestyle='--')
			plt.xlabel("Amino acid position")
			plt.ylabel("Lower bound\nsliding average")
			plt.title(self.proteins[protein_id].id)
			avg_plot_name = "3sliding_avg_" + protein_id + ".png"
			fig.savefig(join(self.out_path, avg_plot_name), bbox_inches='tight', dpi=600)
			print("Saved {}".format(join(self.out_path, avg_plot_name)))

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
		protein_rec.letter_annotations["found_epitopes_count"] = len(protein_rec.seq) * [0.0]
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
		"""
		print("RF=0 -> load(RF,immunome.csv) -> slide window")
		self.load_proteins()
		for protein_id,protein_rec in self.proteins.items():
			print("Process protein: {}".format(protein_id))
			self.proteins[protein_id] = self.set_RF_zero(protein_rec)
			# read immuno csv
			immunome_df = self.load_immunome_csv(protein_id)
			# sum lower bound RF per position
			immunome_df.apply(self.sum_per_position,protein_id=protein_id,axis=1)
			# self.proteins[protein_id].letter_annotations["sliding_avg_rf"] = [lower_bound_sum/epitopes_count for lower_bound_sum, epitopes_count in zip(self.proteins[protein_id].letter_annotations["lower_bound_rf"],self.proteins[protein_id].letter_annotations["found_epitopes_count"])]

			for aa_pos in range(len(protein_rec.seq)):
				if self.proteins[protein_id].letter_annotations["found_epitopes_count"][aa_pos] > 0:
					pos_avg = self.proteins[protein_id].letter_annotations["lower_bound_rf"][aa_pos] / self.proteins[protein_id].letter_annotations["found_epitopes_count"][aa_pos]
				else:
					pos_avg = self.proteins[protein_id].letter_annotations["lower_bound_rf"][aa_pos]
				self.proteins[protein_id].letter_annotations["sliding_avg_rf"][aa_pos] = pos_avg

			self.compute_sliding_avg(protein_id,window_size,save_plot)