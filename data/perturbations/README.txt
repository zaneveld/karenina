# parse_perturbation_file

Only take in 1 file using the option --pert_file_path
File can have 1 or more perturbations
(no warnings for wrong format or data as of now)

Each line of entry is converted to a dictionary
Each dictionary is stored in a perturbations list
The perturbations list is put into a list of treatments

EXAMPLE Format:
params	lambda	0.005	update_mode	replace	axes	x	y	z
