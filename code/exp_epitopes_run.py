from code.epitopes.Protein import Protein
from code.epitopes.Immunodominance import Immunodominance
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
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/data/sars_cov1"
proteins_folder = "proteins"
exp_epitopes_folder = "Bcells"
window_size = 10
save_plot = False  # True
sliding_avg_cutoff = 0.3
save_immunodominant_reg = True
# output
out_path = "/home/damian/Documents/L3S/projects/sars_cov2/exp_epitopes/Bcells"
create_dir(out_path)

print("==== ====")
Immunodominance_Bcells = Immunodominance(data_path, proteins_folder, exp_epitopes_folder, out_path)
processed_proteins = Immunodominance_Bcells.process_all_proteins(window_size, save_plot, sliding_avg_cutoff,
                                                                 save_immunodominant_reg)
print("=== * ===")
print("== *** ==")