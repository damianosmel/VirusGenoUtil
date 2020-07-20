import itertools


class EpitopeObject:
	"""
	Epitope to store the most important information of an epitope
	"""
	new_id = itertools.count().next

	def __init__(self, virus_taxid, protein_ncbi_id, host_taxid, cell_type, hla_restriction, response_frequency,
	             region_seq, region_start, region_stop, is_imported, prediction_process):
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
		prediction_process : str
			epitope identification process
		"""
		self.id = EpitopeObject.new_id
		self.virus_taxid = virus_taxid
		self.protein_ncbi_id = protein_ncbi_id
		self.host_taxid = host_taxid.split("_")[1]
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
		self.is_linear = self.check_linearity()
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
