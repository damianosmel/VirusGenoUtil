from code.epitopes.Protein import Protein
from code.epitopes.ProcessImmunodominance import ProcessImmunodominance
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
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/data"
exp_epitopes_folder = "sars_cov1"

# output
out_path = "/home/damian/Documents/L3S/projects/sars_cov2_Bcells"

processImmuno_Bcells = ProcessImmunodominance(data_path, exp_epitopes_folder, out_path)
