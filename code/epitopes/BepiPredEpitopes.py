class BepiPredEpitopes:
	"""
	Class to extract epitopes found by BepiPred 2.0
	Following parameters from the reference paper

	Reference: "Defining Immunodominant Regions within SARS-CoV Genome"
	from the paper: Grifoni, Alba, et al. "A sequence homology and bioinformatic approach can
	predict candidate targets for immune responses to SARS-CoV-2." Cell host & microbe (2020).

	Example:
	input: bepipred output tsv file of continous epitopes for the structural proteins of SARS-CoV-2
	(please update) output: file containing the aligned sequences, response frequency and identity
	"""

	def __init__(self, data_path, bepipred_folder, out_path):
		"""
		BepiPredEpitopes constructor

		Parameters
		----------
		data_path : str
			input data root
		bepipred_folder : str
			bepipred folder
		out_path : str
			output data root

		Returns
		-------
		None
		"""
		self.data_path = data_path
		self.bepipred_folder = bepipred_folder
		self.out_path = out_path