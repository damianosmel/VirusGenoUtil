from os import makedirs


def create_dir(base_path):
	makedirs(base_path, exist_ok=True)


def is_fasta_file_extension(file_name):
	if file_name[-4:] == ".fna":
		return True
	elif file_name[-4:] == ".faa":
		return True
	elif file_name[-6:] == ".fasta":
		return True
	elif file_name[-6:] == ".fastq":
		return True
	else:
		return False