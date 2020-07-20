from code.utils import is_fasta_file_extension
from code.epitopes.Protein import Protein
from os.path import join, splitext, isfile
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
		self.host_taxon_id = "http://purl.obolibrary.org/obo/NCBITaxon_" + host_taxon_id.split("_")[1]
		self.host_name = str(host_name)
		self.assay_type = assay_type

		self.output_path = output_path
		self.current_virus_proteins = {}
		self.current_virus_taxon_id = None
		self.current_epitopes = {}  # Epitope

		self.epitope_id = 0
		self.epitope_fragment_id = 0
		self.load_iedb_csvs()

	def epitopes2csv(self):
		pass

	def epitope_fragments2csv(self):
		pass

	def virus_epitopes2csv(self):
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
			self.tcell_iedb_assays = self.tcell_iedb_assays[
				self.tcell_iedb_assays["Qualitative Measure"].str.contains(self.assay_type, case=False)]
			# print(self.tcell_iedb_assays.head())
			self.bcell_iedb_assays = self.bcell_iedb_assays[
				self.bcell_iedb_assays["Qualitative Measure"].str.contains(self.assay_type, case=False)]
		print(self.host_taxon_id)
		print("Select epitopes experimental identified in host id: {} = {}".format(self.host_taxon_id.split("_")[1],
		                                                                           self.host_name))
		self.tcell_iedb_assays = self.tcell_iedb_assays[self.tcell_iedb_assays["Host IRI"] == self.host_taxon_id]
		self.bcell_iedb_assays = self.bcell_iedb_assays[self.bcell_iedb_assays["Host IRI"] == self.host_taxon_id]
		print(self.bcell_iedb_assays.head())
		print("--")
		print(self.tcell_iedb_assays.head())

	def subset_iedb_by_virus_id(self):
		"""
		Subset iedb record to get the ones related to virus id

		Returns
		-------
		Pandas.DataFrame, Pandas.DataFrame
			B-cell IEDB assays related only to virus id, T-cell IEDB assay only to virus id
		"""
		# self.current_virus_taxon_id = 694009
		print("Get iedb only for taxon id={}".format(self.current_virus_taxon_id))
		url_taxon_id = "http://purl.obolibrary.org/obo/NCBITaxon_" + str(self.current_virus_taxon_id)
		return self.bcell_iedb_assays[self.bcell_iedb_assays["Organism IRI"] == url_taxon_id], self.tcell_iedb_assays[
			self.tcell_iedb_assays["Organism IRI"] == url_taxon_id]

	def process_all_viruses(self):
		"""
		Process and get IEDB epitopes for each protein of each virus

		Returns
		-------
		None
		"""
		print("Process all viruses found in {}".format(self.viruses_path))
		self.subset_iedb_by_host_assay_type()
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
						print("Saving protein with id: {}".format(protein_id))
						# protein_rec = protein.get_record()
						self.current_virus_proteins[protein_id] = protein

		print(self.current_virus_proteins)
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

		print(bcells_current_virus.head())
		print("--")
		print(tcells_current_virus.head())

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
			tcells_current_protein = tcells_current_virus.loc[tcells_current_virus["Antigen Name"] == protein_id]

			# map epitope to allele
			epitope2allele = self.map_epitope2allele(tcells_current_protein)

			# save unique epitopes into fasta file
			protein_record = self.proteins[protein_id].get_record()

			# find epitope regions
			epi_regions = self.find_epitope_regions(protein_record, list(epitope2allele.keys()))

			# calculate RF score
			epi2rf_score = self.calculate_RF_score(tcells_current_protein, list(epitope2allele.keys()))
		print("====")
