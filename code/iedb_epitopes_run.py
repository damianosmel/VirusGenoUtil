from code.utils import create_dir
from code.epitopes.IEDBEpitopes import IEDBEpitopes

###                                                          ###
###                  Import IEDB epitopes                    ###
###                                                          ###

###                                                          ###
###     input = IEDB Bcells.csv + Tcells.csv                 ###
###     output = epitope csv for all taxon ids               ###
###                                                          ###

# input
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/iedb_data"
cells_epitopes_folder = "cell_epitopes"
virus_proteins_folder = "viruses_proteins"
### Parameters August 2020 ###
# host_taxon, host_name = "taxon_9606", "Homo Sapiens"
# assay_type = "positive"
host_taxon, host_name = "taxa_all", None
assay_type = "all"  # all = positive and negative

# output
output_path = "/home/damian/Documents/L3S/projects/sars_cov2/IEDB_epitopes_ViruSurf_sept2020"
create_dir(output_path)

IEDBEpitopes = IEDBEpitopes(data_path, cells_epitopes_folder, virus_proteins_folder,
                            host_taxon, host_name, assay_type, output_path)
IEDBEpitopes.process_all_viruses()
