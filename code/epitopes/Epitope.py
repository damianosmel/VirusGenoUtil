import itertools


class Epitope:
	"""
	Epitope to store the most important information of an epitope
	"""
	# credits: https://stackoverflow.com/questions/1045344/how-do-you-create-an-incremental-id-in-a-python-class/54318273#54318273
	new_id = itertools.count()

	def __init__(self, virus_taxid, protein_ncbi_id, host_taxid, cell_type, hla_restriction, response_frequency,
	             region_seq, region_start, region_stop, is_imported, external_links, prediction_process):
		"""
		Epitope construstor

		Parameters
		----------
		virus_taxid  : str
			virus NCBI taxonomy id
		protein_ncbi_id : str
			protein (antigen) NCBI id
		host_taxid : str
			host NCBI taxonomy
		cell_type : str
			cell type
		hla_restriction : str
			HLA restriction (for T-cells)
		response_frequency : float
			response frequency of epitope
		region_seq : str
			epitope region sequence
		region_start : int
			epitope region start
		region_stop : int
			epitope region stop
		is_imported : bool
			epitope information is imported (True), otherwise is predicted (False)
		external_links : list of str
			IEDB reference links
		prediction_process : str
			epitope identification process
		"""
		self.id = next(Epitope.new_id)
		self.virus_taxid = virus_taxid
		self.protein_ncbi_id = protein_ncbi_id
		self.host_taxid = host_taxid
		self.cell_type = cell_type
		if cell_type == "B cell":
			self.hla_restriction = "not applicable"
		else:
			self.hla_restriction = hla_restriction
		self.response_frequency = response_frequency
		self.seq = region_seq
		self.region_start = region_start
		self.region_stop = region_stop
		self.is_imported = is_imported
		self.prediction_process = prediction_process
		self.external_links = external_links
		self.is_linear = True  # self.check_linearity()
		self.epitope_fragments = []
		if not self.is_linear:
			# self.fragment()
			pass

	def check_linearity(self):
		pass

	def fragment(self):
		# create epitope fragment for each start-stop
		# ask to get all epitope fragments
		pass

	# add each epitope to the fragments

	def get_all_attributes(self):
		"""
		Return all the attributes of the epitope

		Returns
		-------
		dict of str
			all epitope attributes returned in a dictionary
		"""
		return {"epitope_id": self.id,
		        "virus_taxid": self.virus_taxid,
		        "protein_ncbi_id": self.protein_ncbi_id,
		        "host_taxid": self.host_taxid,
		        "cell_type": self.cell_type,
		        "hla_restriction": self.hla_restriction,
		        "response_frequency": self.response_frequency,
		        "region_seq": self.seq,
		        "region_start": self.region_start,
		        "region_stop": self.region_stop,
		        "is_imported": self.is_imported,
		        "external_links": self.external_links,
		        "prediction_process": self.prediction_process,
		        "fragments": self.epitope_fragments,
		        "is_linear": self.is_linear
		        }
