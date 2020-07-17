from code.utils import is_fasta_file_extension
from code.epitopes.Protein import Protein
from os.path import join, splitext
from os import scandir


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

	def __init__(self, data_path, cell_epitopes_folder, viruses_folder, host_taxon_id, host_name, output_path):
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
		output_path : str
			output data root

		Returns
		-------
		None
		"""
		self.data_path = data_path
		self.cell_epitopes_path = join(data_path, cell_epitopes_folder)
		self.viruses_path = join(data_path, viruses_folder)
		self.host_taxon_id, self.host_name = host_taxon_id, host_name
		self.output_path = output_path
		self.current_virus_proteins = {}
		self.current_virus_taxon = None
		self.current_epitopes = {}  # EpitopeObject

	def virus_epitopes2csv(self):
		pass

	def process_all_viruses(self):
		"""
		Process and get IEDB epitopes for each protein of each virus

		Parameters
		----------

		Returns
		-------

		"""
		print("Process all viruses found in {}".format(self.viruses_path))
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

		"""
		self.current_virus_taxon = splitext(virus_proteins_path.name)[0]
		print("=== Virus ===")
		print("Process proteins of virus with taxon id: {}".format(self.current_virus_taxon))
		self.load_virus_proteins(virus_proteins_path)

	def process_Bcells(self):
		print("Process B-cells")
		for protein_id, protein_rec in self.proteins.items():
			print("--- ---")
			print("Process protein: {}".format(protein_id))
			pass

	def process_Tcells(self):
		print("Process T-cells")
		for protein_id, protein_rec in self.proteins.items():
			print("--- ---")
			print("Process protein: {}".format(protein_id))
			pass
