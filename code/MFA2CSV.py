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
		print("Parse ORF reference sequence")
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

	def alignment2csv(self, alignment):
		"""
		Convert pairwise alignment to csv
		Credits: following loosy the VCF format tutorial at:
		https://faculty.washington.edu/browning/beagle/intro-to-vcf.html

		Parameters
		----------
		alignment : Bio.Align.MultipleSeqAlignment
			aligned target sequence

		Returns
		-------
		None
		"""
		print("Convert alignment of target id {} to variants csv".format(alignment.id))
		variants_name = self.ref_record.name + "_" + alignment.id + "_variants.csv"
		with open(join(self.out_path,variants_name),"w") as variants_out:
			variants_out.write("POS,REF_ID,TARGET_ID,REF,ALT\n")
			for ref_pos,ref_base in enumerate(self.ref_record.seq):
				if ref_base != alignment.seq[ref_pos]:  # variant detected
					print("Write variant")
					variant_row = ",".join([str(ref_pos+1),self.ref_record.id,alignment.id,ref_base,alignment.seq[ref_pos]])
					variants_out.write(variant_row + "\n")

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
		print("Parse MAF")
		multiple_alignment = AlignIO.read(join(self.alignments_path, mfa_file), "fasta")
		for alignment in multiple_alignment:
			self.alignment2csv(alignment)
		print("---")