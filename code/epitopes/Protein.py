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
		self.protein_file = protein_file
		self.record = None
		self.short_id = None

	def load(self):
		"""
		Load protein

		Parameters
		----------

		Returns
		-------
		None
		"""
		self.record = SeqIO.parse(self.protein_file, "fasta").__next__()
		print("Protein record for {} was created.".format(self.record.id))

	def get_record(self):
		"""
		Get protein record

		Parameters
		----------

		Returns
		-------
		Bio.SeqIO.SeqRecord
			loaded protein sequence record
		"""
		return self.record

	def get_short_id(self):
		"""
		Get protein short id

		Parameters
		----------

		Returns
		-------
		str
			protein short id
		"""
		return self.record.id.split("|")[1]

	def get_name(self):
		"""
		Get protein name from NCBI identifier
		e.g. gi|1796318598|ref|YP_009724390.1|:1-1273 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]

		Returns
		-------
		str
			extract protein name
		"""
		return self.record.description.split(self.record.name)[1].lstrip()

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
