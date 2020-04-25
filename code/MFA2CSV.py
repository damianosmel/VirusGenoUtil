from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import xml.etree.ElementTree as ET
from os.path import join
from utils import create_dir,is_fasta_file_extension


class MFA2CSV:
	"""
	Class to convert a multi-fasta alignment file (MFA)
	to variants csv per target sequence

	Example:
		pos:  012345678
		ref:  ATGGAGAGC
	target_1: ATGGT--GC

	then an example of produced variant csv for target1 will be:
	ref_id   |   target_id   |   ref_location  | ref | alt | type
	ncbi_ref | ncbi_target   |        4        |  A  |  T  | mismatch
	ncbi_ref | ncbi_target   |        5        |  G  |  -  | deletion
	"""

	def __init__(self, data_path, alignments_folder, xmls_folder, ncbi_ref_id, out_path):
		"""
		MFA2CSV constructor

		Parameters
		----------
		data_path : str
			absolute path where data lie
		alignments_folder : str
			alignments folder name
		xmls_folder : str
			xmls folder name
		ncbi_ref_id : str
			NCBI reference id
		out_path : str
			absolute path where output will be placed

		Returns
		-------
		None
		"""
		self.data_path = data_path
		self.out_path = out_path
		self.alignments_path = join(data_path,alignments_folder)
		self.xmls_path = join(data_path,xmls_folder)
		self.ref_id = ncbi_ref_id
		self.ref_record = None
		create_dir(out_path)

	def parse_ref(self, xml_file):
		"""
		Parse reference sequence from xml file (created by virusalign)

		Parameters
		----------
		xml_file : str
			xml file name

		Returns
		-------
		None
		"""
		print("Parsing ORF reference sequence, please wait..")
		root = ET.parse(join(self.xmls_path, xml_file)).getroot()
		assert "referenceSequence" in root.attrib.keys(), "AssertionError: {} does not contain reference sequence".format(xml_file)
		self.ref_seq = Seq(root.attrib["referenceSequence"])
		ref_name, ref_description = None, None
		if "name" in root.attrib.keys():
			ref_name = root.attrib["name"]
		if "description" in root.attrib.keys():
			ref_description = root.attrib["description"]

		self.ref_record = SeqRecord(Seq(root.attrib["referenceSequence"]),
		                   id=self.ref_id, name=ref_name,
		                   description=ref_description)

		print("Parsed reference record:\n {}".format(self.ref_record))
		print("---")

	def run(self, xml_file, mfa_file):
		"""
		Convert MFA file to variants csv per target sequence

		Parameters
		----------
		xml_file : str
			XML file containing the sequence of the ORF element
		mfa_file : str
			MFA file containing the multiple-fasta alignment

		Returns
		-------
		None
		"""
		# parse ref
		self.parse_ref(xml_file)
		# parse multi-fasta alignment
		assert is_fasta_file_extension(mfa_file), "AssertionError: Currently supporting only multi-fasta alignment files."
		print("Parsing MAF, please wait..")
		multiple_alignments = AlignIO.read(join(self.alignments_path, mfa_file), "fasta")
		for alignment in multiple_alignments:
			print(alignment)
		print("---")