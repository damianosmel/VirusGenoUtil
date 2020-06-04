from code.epitopes.Immunodominance import Immunodominance
from code.epitopes.HomologyBasedEpitopes import HomologyBasedEpitopes

from code.utils import create_dir

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
# 1)sars cov1 proteins + immunome browser -> lower bound sliding window per position#
# 2)sars cov2 proteins -BLAST-> sars cov1 -> similarity per position                #
# 3)lower bound + similarity per position -> experimental epitopes for sars cov2    #
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #


###                                          ###
### B cells epitopes sars-cov1 --> sars-cov2 ###
###                                          ###

# input
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/data/exp_epitopes_input"
query_proteins_folder = "sars_cov1_proteins"
exp_epitopes_folder = "Bcells"
window_size = 10
save_plot = False  # True
sliding_avg_cutoff = 0.3
save_immunodominant_reg = True
# output
out_path = "/home/damian/Documents/L3S/projects/sars_cov2/exp_epitopes/Bcells"
create_dir(out_path)

# print("==== ====")
# Immunodominance_Bcells = Immunodominance(data_path, query_proteins_folder, exp_epitopes_folder, out_path)
# processed_proteins = Immunodominance_Bcells.process_all_proteins(window_size, save_plot, sliding_avg_cutoff,
#                                                                  save_immunodominant_reg)

print("==== ====")
subject_proteins_folder = "sars_cov2_proteins"
blast_folder = "blast"
HomologyBasedEpitopes_Bcells = HomologyBasedEpitopes(data_path, subject_proteins_folder, blast_folder,out_path)
HomologyBasedEpitopes_Bcells.find_epitopes()
print("=== * ===")
print("== *** ==")
