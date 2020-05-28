from Bio import SeqIO


class Protein:
	"""
	Class to keep all needed info for Protein objects

	Example:
	protein --> id, organism, amino acid sequence
	"""

	def __init__(self, protein_file):
		"""
		Protein constructor

		Parameters
		----------
		protein_file : str
			protein file

		Returns
		-------
		None
		"""
		self.protein_record = SeqIO.parse(protein_file, "fasta")
		print("Protein record for {} was created.".format(self.protein_record.id))

	def get_protein_id(self):
		"""
		Get protein id

		Parameters
		----------

		Returns
		-------
		str
			protein id
		"""
		return str(self.protein_record.id)

	def get_protein_length(self):
		"""
		Get protein length

		Parameters
		----------

		Returns
		-------
		int
			protein length
		"""
		return len(self.protein_record)
