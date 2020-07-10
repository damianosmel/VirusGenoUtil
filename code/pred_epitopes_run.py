from code.utils import create_dir

###                                          ###
### Bepipred 2.0 output SARS-CoV-2 --> csv   ###
###                                          ###

#input

data_path = "/home/damian/Documents/L3S/projects/sars_cov2/sars_cov2_data/pred_epitopes_input"
pred_epitopes_folder = "bepipred_out"

# output
out_path = "/home/damian/Documents/L3S/projects/sars_cov2/pred_epitopes/bepipred"
create_dir(out_path)