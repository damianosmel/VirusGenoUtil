from code.utils import is_fasta_file_extension
from code.epitopes.Protein import Protein
from os.path import join, splitext, isfile
from os import scandir
from pandas import read_csv

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

	def __init__(self, data_path, cell_epitopes_folder, viruses_folder, host_taxon_id, host_name, assay_type, output_path):
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
		self.host_taxon_id, self.host_name = host_taxon_id, host_name
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
		assert isfile(join(self.cell_epitopes_path, "tcell_full_v3.csv")), "AssertionError: IEDB Tcell assays csv was not found in {}".format(self.cell_epitopes_path)
		self.tcell_iedb_assays = read_csv(join(self.cell_epitopes_path, "tcell_full_v3.csv"), sep=",", header=1)
		assert isfile(join(self.cell_epitopes_path,"bcell_full_v3.csv")), "AssertionError: IEDB Bcell assays csv was not found in {}".format(self.cell_epitopes_path)
		self.bcell_iedb_assays = read_csv(join(self.cell_epitopes_path, "bcell_full_v3.csv"), sep=",", header=1)

	def subset_iedb_assays_by_host_assay(self):
		"""

		Returns
		-------

		"""
		print("Get working subset of iedb assays")
		if self.assay_type == "positive":
			self.tcell_iedb_assays = self.tcell_iedb_assays[self.tcell_iedb_assays["Qualitative Measure"].str.contains(self.assay_type,case=False)]
			print(self.tcell_iedb_assays.head())

	def process_all_viruses(self):
		"""
		Process and get IEDB epitopes for each protein of each virus

		Returns
		-------

		"""
		print("Process all viruses found in {}".format(self.viruses_path))
		self.subset_iedb_assays()
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
				print(proteins_content.name)
				if is_fasta_file_extension(proteins_content.name) and proteins_content.is_file():
					protein = Protein(join(virus_proteins_path, proteins_content.name))
					protein.load()
					protein_id = splitext(proteins_content.name)[0]
					if protein_id not in self.current_virus_proteins:
						print("Saving protein with id: {}".format(protein_id))
						protein_rec = protein.get_record()
						self.current_virus_proteins[protein_id] = protein_rec
		print(self.current_virus_proteins)

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
		self.current_virus_taxon_id = splitext(virus_proteins_path.name)[0]
		print("=== Virus ===")
		print("Process proteins of virus with taxon id: {}".format(self.current_virus_taxon_id))
		self.load_virus_proteins(virus_proteins_path)
		self.process_Bcells()
		self.process_Tcells()

	def process_Bcells(self):
		print("Process B-cells")
		for protein_id, protein_rec in self.current_virus_proteins.items():
			print("---")
			print("Process protein: {}".format(protein_id))
			pass

	def process_Tcells(self):
		print("Process T-cells")
		for protein_id, protein_rec in self.current_virus_proteins.items():
			print("---")
			print("Process protein: {}".format(protein_id))
			subset_virus_protein = self.tcell_iedb_assays.loc[self.tcell_iedb_assays["Antigen IRI"].split("/")[-1] == protein_id and self.tcell_iedb_assays["Organism IRI"].split("_")[-1]==self.current_virus_taxon_id]
			for _, row in subset_virus_protein.iterrows():
				print(row)
