from code.utils import is_fasta_file_extension
from code.epitopes.Protein import Protein
from code.epitopes.Epitope import Epitope
from code.epitopes.EpitopeFragment import EpitopeFragment
from os.path import join, splitext, isfile, exists
from os import scandir
from pandas import read_csv, unique
from math import isnan, sqrt


class IEDBEpitopes:
	"""
	Class to extract B and T-cells epitopes for
	each virus taxon id found in virus_protein folder

	Input epitopes per cell type are downloaded from http://www.iedb.org/database_export_v3.php
	IEDB paper reference:
	Vita R, Mahajan S, Overton JA, Dhanda SK, Martini S, Cantrell JR, Wheeler DK, Sette A, Peters B.
	The Immune Epitope Database (IEDB): 2018 update. Nucleic Acids Res. 2019 Jan 8;47(D1):D339-D343.
	doi: 10.1093/nar/gky1006. PMID: 30357391
	"""

	def __init__(self, data_path, cell_epitopes_folder, viruses_folder, host_taxon_id, host_name, assay_type,
	             output_path):
		"""
		IEDBEpitopes constructor

		Parameters
		----------
		data_path : str
			input data root
		cell_epitopes_folder : str
			B and T cell IEDB epitopes folder
		viruses_folder : str
			viruses protein folder
		host_taxon_id : str
			NCBI taxon id of host organism whose cells used for the experimental identification of epitopes
		host_name : str
			name of host organism whose cells used for the experimental identification of epitopes
		assay_type : str
			IEDB experimental assay type
		output_path : str
			output data root

		Returns
		-------
		None
		"""
		self.data_path = data_path
		self.cell_epitopes_path = join(data_path, cell_epitopes_folder)
		self.viruses_path = join(data_path, viruses_folder)

		self.tcell_iedb_assays, self.bcell_iedb_assays = None, None
		self.host_taxon_id = host_taxon_id.split("_")[1]
		self.host_name = str(host_name)
		self.assay_type = assay_type
		self.url_prefix = "http://purl.obolibrary.org/obo/NCBITaxon_"
		self.output_path = output_path
		self.current_virus_proteins = {}
		self.current_virus_taxon_id = None
		self.current_virus_epitopes = []  # all Epitope objects for current virus

		self.epitope_id = 0
		self.epitope_fragment_id = 0
		self.load_iedb_csvs()

	def epitope_fragments2csv(self):

		pass

	def virus_epitopes2csv(self):
		print("Save current virus epitopes to csv")

		epitopes_name = "imported_iedb_epitopes.tsv"
		if not exists(join(self.output_path, epitopes_name)):
			print("Create file: {}".format(epitopes_name))
			with open(join(self.output_path, epitopes_name), "w") as epitopes_out:
				epitopes_out.write(
					"epitope_id\tvirus_taxid\thost_taxid\tprotein_ncbi_id\tcell_type\thla_restriction\tresponse_frequency\tepitope_sequence\tepitope_start\tepitope_stop\tis_imported\texternal_links\tprediction_process\tis_linear\n")

		with open(join(self.output_path, epitopes_name), "a") as epitopes_out:
			print("Update file: {}".format(epitopes_name))
			for epitope in self.current_virus_epitopes:
				print("Write IEDB imported epitope")
				epi_attributes = epitope.get_all_attributes()
				epitope_row = "\t".join(
					[str(epi_attributes["epitope_id"]), epi_attributes["virus_taxid"], epi_attributes["host_taxid"],
					 epi_attributes["protein_ncbi_id"], epi_attributes["cell_type"], epi_attributes["hla_restriction"],
					 str(epi_attributes["response_frequency"]), epi_attributes["region_seq"],
					 str(epi_attributes["region_start"]), str(epi_attributes["region_stop"]),
					 str(epi_attributes["is_imported"]), ",".join(epi_attributes["external_links"]),
					 str(epi_attributes["prediction_process"]), str(epi_attributes["is_linear"])])
				epitopes_out.write(epitope_row + "\n")
				if epi_attributes["is_linear"]:  # write epitope
					pass

	def load_iedb_csvs(self):
		"""
		Load IEDB B and T cells assays csvs

		Returns
		-------
		None
		"""
		assert isfile(join(self.cell_epitopes_path,
		                   "tcell_full_v3.csv")), "AssertionError: IEDB Tcell assays csv was not found in {}".format(
			self.cell_epitopes_path)
		self.tcell_iedb_assays = read_csv(join(self.cell_epitopes_path, "tcell_full_v3.csv"), sep=",", header=1)
		assert isfile(join(self.cell_epitopes_path,
		                   "bcell_full_v3_1000.csv")), "AssertionError: IEDB Bcell assays csv was not found in {}".format(
			self.cell_epitopes_path)
		self.bcell_iedb_assays = read_csv(join(self.cell_epitopes_path, "bcell_full_v3_1000.csv"), sep=",", header=1)


	def subset_iedb_by_host_assay_type(self):
		"""
		Subset IEDB records by host taxon id and assay type

		Returns
		-------
		None
		"""
		print("Get working subset of iedb assays")
		if self.assay_type == "positive":
			print("Select positive assays")
			self.tcell_iedb_assays = self.tcell_iedb_assays.loc[
				self.tcell_iedb_assays["Qualitative Measure"].str.contains(self.assay_type, case=False)]
			# print(self.tcell_iedb_assays.head())
			self.bcell_iedb_assays = self.bcell_iedb_assays.loc[
				self.bcell_iedb_assays["Qualitative Measure"].str.contains(self.assay_type, case=False)]

		print("Select epitopes experimental identified in host id: {} = {}".format(self.host_taxon_id, self.host_name))
		host_taxid = self.url_prefix + self.host_taxon_id
		self.tcell_iedb_assays = self.tcell_iedb_assays.loc[self.tcell_iedb_assays["Host IRI"] == host_taxid]
		self.bcell_iedb_assays = self.bcell_iedb_assays.loc[self.bcell_iedb_assays["Host IRI"] == host_taxid]
		print("B cell subset number of non-unique epitopes: {}".format(self.bcell_iedb_assays.shape[0]))
		print("T cell subset number of non-unique epitopes: {}".format(self.tcell_iedb_assays.shape[0]))

		self.tcell_iedb_assays.to_csv(join(self.cell_epitopes_path, "tcell_positive_host.csv"), index=False)

	def subset_Bcells_by_epi_type(self):
		"""
		Subset B cells to keep only B-cells linear and discontinuous epitopes

		Returns
		-------
		None
		"""
		print("Select only linear or discontinuous epitopes")
		self.bcell_iedb_assays = self.bcell_iedb_assays.loc[(self.bcell_iedb_assays["Object Type"] == "Discontinuous peptide") | (self.bcell_iedb_assays["Object Type"] == "Linear peptide")]
		print("B cell subset number of non-unique epitopes: {}".format(self.bcell_iedb_assays.shape[0]))

	def subset_iedb_by_virus_id(self):
		"""
		Subset iedb record to get the ones related to virus id

		Returns
		-------
		Pandas.DataFrame, Pandas.DataFrame
			B-cell IEDB assays related only to virus id, T-cell IEDB assay only to virus id
		"""
		print("Get iedb only for taxon id={}".format(self.current_virus_taxon_id))
		tcell_iedb_virus = self.tcell_iedb_assays[
			self.tcell_iedb_assays["Organism IRI"].str.split("NCBITaxon_").str[-1] == self.current_virus_taxon_id]
		bcell_iedb_virus = self.bcell_iedb_assays[
			self.bcell_iedb_assays["Organism IRI"].str.split("NCBITaxon_").str[-1] == self.current_virus_taxon_id]
		print("B cell subset for virus contains {} non unique epitopes".format(bcell_iedb_virus.shape[0]))
		print("T cell subset for virus contains {} non unique epitopes".format(tcell_iedb_virus.shape[0]))
		return bcell_iedb_virus, tcell_iedb_virus

	def process_all_viruses(self):
		"""
		Process and get IEDB epitopes for each protein of each virus

		Returns
		-------
		None
		"""
		print("Process all viruses found in {}".format(self.viruses_path))
		self.subset_iedb_by_host_assay_type()
		self.subset_Bcells_by_epi_type()
		with scandir(self.viruses_path) as viruses_dir:
			for content in viruses_dir:
				if content.is_dir() and "taxon_" in content.name:
					self.process_virus_proteins(content)
					self.virus_epitopes2csv()
		print("=== ~ ===")

	def load_virus_proteins(self, virus_proteins_path):
		"""
		Load virus proteins from argument path

		Parameters
		----------
		virus_proteins_path : str
			virus proteins path

		Returns
		-------
		None
		"""
		print("Load virus proteins")
		self.current_virus_proteins = {}
		with scandir(virus_proteins_path) as proteins_dir:
			for proteins_content in proteins_dir:
				# print(proteins_content.name)
				if is_fasta_file_extension(proteins_content.name) and proteins_content.is_file():
					protein = Protein(join(virus_proteins_path, proteins_content.name))
					protein.load()
					protein_id = splitext(proteins_content.name)[0]
					if protein_id not in self.current_virus_proteins:
						# print("Saving protein with id: {}".format(protein_id))
						# protein_rec = protein.get_record()
						self.current_virus_proteins[protein_id] = protein

		# print(self.current_virus_proteins)
		print("====")

	def process_virus_proteins(self, virus_proteins_path):
		"""
		Process information for IEDB epitopes for each protein of the virus

		Parameters
		----------
		virus_proteins_path : str
			virus protein folder path

		Returns
		-------
		None
		"""
		self.current_virus_taxon_id = splitext(virus_proteins_path.name)[0].split("_")[1]
		print("=== Virus ===")
		print("Process proteins of virus with taxon id: {}".format(self.current_virus_taxon_id))
		self.load_virus_proteins(virus_proteins_path)
		bcells_current_virus, tcells_current_virus = self.subset_iedb_by_virus_id()
		self.current_virus_epitopes = []  # clear current virus epitopes
		# print(bcells_current_virus.head())
		# print("--")
		# print(tcells_current_virus.head())
		tcells_current_virus.to_csv(join(self.cell_epitopes_path, "tcell_virus.csv"), index=False)
		# self.process_Bcells(bcells_current_virus)
		self.process_Tcells(tcells_current_virus)

	def process_Bcells(self, bcells_current_virus):
		"""
		Process B-cells assays for current virus

		Parameters
		----------
		bcells_current_virus : Pandas.DataFrame
			B-cells assays for current virus

		Returns
		-------

		"""
		print("Process B-cells")
		for protein_id, protein in self.current_virus_proteins.items():
			print("---")
			print("Process protein: {}".format(protein_id))
			pass
		print("====")

	def map_epitope2allele(self, idbe_assay):
		"""
		Map epitopes to allele information

		Parameters
		----------
		idbe_assay : Pandas.DataFrame
			dataframe created from csv of IDBE assay tab
		save_epitopes : bool
			Save retrieved epitopes sequences (True), don't save (False)

		Returns
		-------
		dict of str : str
			dictionary mapping unique epitopes to their allele info
		"""
		print("Map unique epitope sequence to allele information")
		uniq_epitopes2allele = {}
		unique_epitopes = unique(idbe_assay["Description"])

		for uniq_epi in unique_epitopes:
			if uniq_epi not in uniq_epitopes2allele:
				all_allele_names = []
				for _, row in idbe_assay.loc[idbe_assay["Description"] == uniq_epi, ["Allele Name"]].iterrows():
					allele = str(row["Allele Name"])
					if "HLA" not in allele:
						allele = "unknown"
					else:
						allele = allele.replace(" ", "_")
					if allele not in all_allele_names:
						all_allele_names.append(allele)

				if "unknown" in all_allele_names and len(all_allele_names) > 1:
					all_allele_names.remove("unknown")
				uniq_epitopes2allele[uniq_epi] = ",".join(all_allele_names)

		return uniq_epitopes2allele

	def find_epitope_regions(self, protein_rec, epitopes):
		"""
		Find the epitopes regions (start-stop) in the protein record
		Credits: Chapter 20.1.8 Biopython (http://biopython.org/DIST/docs/tutorial/Tutorial.html)
		Parameters
		----------
		protein_rec : Bio.SeqIO.SeqRecord
			protein record
		epitopes : list of str
			epitope sequences

		Returns
		-------
			list of list of int
			sorted list of start end position of epitopes regions
		"""
		print("Map unique epitope sequence to start end positions")
		epi_regions = []
		for epi in epitopes:
			start = protein_rec.seq.find(epi)
			if start != -1:
				end = start + len(epi) - 1
				epi_regions.append([start, end])
		epi_regions.sort(key=lambda x: x[0])  # sort ascending by starting position
		return epi_regions

	def find_epitope_external_links(self, idbe_assay, uniq_epitopes):
		"""
		Find external links for each unique epitope

		Parameters
		----------
		idbe_assay : Pandas.DataFrame
			dataframe created from csv of IDBE assay tab
		uniq_epitopes : list of str
			list of unique epitopes for which the external links will be extracted

		Returns
		-------
		dict of str : list of str
			dictionary with key the unique epitope and the value the list of external links
		"""
		print("Extract external links for each unique epitope")
		epi2external_links = {}
		for epi in uniq_epitopes:
			epi2external_links[epi] = []
			for _, row in idbe_assay.loc[idbe_assay["Description"] == epi, ["Reference IRI"]].iterrows():
				if str(row["Reference IRI"]) not in epi2external_links[epi]:
					epi2external_links[epi].append(str(row["Reference IRI"]))
		return epi2external_links

	def calculate_RF_score(self, idbe_assay, uniq_epitopes):
		"""
		Calculate RF score for T cell assays, using:
		RF = (r-sqrt(r))/t,
		where r is the # positive responding assays
		and t is the # of the total tested assays

		Parameters
		----------
		idbe_assay : Pandas.DataFrame
			dataframe created from csv of IDBE assay tab
		uniq_epitopes : list of str
			list of unique epitope for which the RF score will be calculated

		Returns
		-------
		dict of str : float
			dictionary with key the unique epitope and the value the RF score
		"""
		print("Map unique epitope sequence to RF score")
		rf_scores = {}
		for epi in uniq_epitopes:
			# sum r and t for all assays testing the current epitopte

			t, r = 0.0, 0.0
			for _, row in idbe_assay.loc[idbe_assay["Description"] == epi, ["Number of Subjects Tested",
			                                                                "Number of Subjects Responded"]].iterrows():
				if isnan(row["Number of Subjects Tested"]):
					t = t + 1.0
				else:
					t = t + row["Number of Subjects Tested"]
				if isnan(row["Number of Subjects Responded"]):
					r = r + 1.0
				else:
					r = r + row["Number of Subjects Responded"]
			if t == 0:
				rf_score = 0.0
			else:
				rf_score = (r - sqrt(r)) / t
			rf_scores[epi] = rf_score
		return rf_scores

	def process_Tcells(self, tcells_current_virus):
		"""
		Process T-cells assays for current virus

		Parameters
		----------
		tcells_current_virus : Pandas.DataFrame
			T-cells assays for current virus

		Returns
		-------

		"""
		print("Process T-cells")
		for protein_id, protein in self.current_virus_proteins.items():
			print("---")
			print("Process protein name: {}".format(protein.get_name()))
			tcells_current_protein = tcells_current_virus.loc[
				tcells_current_virus["Antigen Name"] == protein.get_name()]
			print("Non unique epitopes for protein = {}".format(tcells_current_protein.shape[0]))

			if tcells_current_protein.shape[0] == 0:
				print("Skip protein as no epitopes are found in IEDB")
			else:
				# map epitope to allele
				epitope2allele = self.map_epitope2allele(tcells_current_protein)

				# save unique epitopes into fasta file
				protein_record = protein.get_record()

				# find epitope regions
				epi_regions = self.find_epitope_regions(protein_record, list(epitope2allele.keys()))

				# calculate RF score
				epi2rf_score = self.calculate_RF_score(tcells_current_protein, list(epitope2allele.keys()))

				# extract external links per unique epitope
				epi2external_links = self.find_epitope_external_links(tcells_current_protein,
				                                                      list(epitope2allele.keys()))

				# create an Epitope object for each identified IEDB epitope
				host_taxon_id = self.host_taxon_id
				is_imported = True
				prediction_process = "IEDB_import"

				for i, region in enumerate(epi_regions):
					reg_start, reg_end = region[0], region[-1] + 1
					epi_frag = protein_record.seq[reg_start:reg_end]
					reg_start = reg_start + 1  # increaase by one to convert 0-index to 1-index
					assert epi_frag in epi2rf_score, "Epitope with sequence {} does not have calculated RF score".format(
						epi_frag)
					external_links = epi2external_links[epi_frag]
					rf_score = float("{:.4f}".format(epi2rf_score[epi_frag]))
					if rf_score > 0.0:
						epitope = Epitope(self.current_virus_taxon_id, protein_id, host_taxon_id, "T cell",
						                  epitope2allele[epi_frag], rf_score, str(epi_frag), reg_start, reg_end,
						                  is_imported,
						                  external_links, prediction_process)
						self.current_virus_epitopes.append(epitope)
				all_attributes = self.current_virus_epitopes[-1].get_all_attributes()
				print(all_attributes)
		print("====")
