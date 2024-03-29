# 2023.11.30 CellBender Script
# Assumes you're using raw_feature_bc_matrix.h5 files output from CellRanger Count (v7.2.0 as of 2023.11.30)
# IO https://cellbender.readthedocs.io/en/latest/usage/index.html

# conda create -n cellbender python=3.7


# Import modules
import os
import re
# import subprocess

# Define global variables

cellranger_folder = "/data/CRGGH/heustonef/hpapdata/cellranger_scRNA/" # Required: folder containing cellranger output data
cellbender_file = "cellbender.swarm" # Optional: None, or name of file to write cellbender commands. Useful if submitting as swmarm file
gpu_run = False
flags = "--cpu-threads 4" # Optional: CellBender flags to include in command
ignore_folders = "SC_RNA_COUNTER_CS"

if cellbender_file:
	cellbender_file = os.path.join(cellranger_folder, cellbender_file)
if gpu_run is True:
		flags = ''.join((" --cuda ", flags))

# Create list of raw_feature_bc_matrix.h5 files
sample_list = []
for root, _, files in os.walk(cellranger_folder):
	if ignore_folders in root:
		continue
	else:
		for file in files:
			if re.match("raw_feature_bc_matrix.h5", file):
				sample_list.append(os.path.join(root, file))
				print(''.join(("Found ", os.path.join(root, file), "\n")))

# Run CellBender

for cellbender_input in sample_list:
	export_folder = os.path.dirname(cellbender_input)
	cellbender_output = os.path.basename(cellbender_input).replace("raw", "cb")
	cellbender_cmd = ''.join(( 
		"cellbender remove-background ",
		flags,
		" --input ",
		cellbender_input,
		" --output ",
		os.path.join(export_folder, cellbender_output), 
		'; \\\n'
	))
	print(cellbender_cmd)
	if cellbender_file:
		cd_cmd = ''.join(("cd ", export_folder, '; \\\n'))
		with open(cellbender_file, 'a') as cellbender_out:
			cellbender_out.write(cd_cmd)
			cellbender_out.write(cellbender_cmd)
			cellbender_out.write('ptrepack --complevel 5 cb_feature_bc_matrix_filtered.h5:/matrix cb-seurat_feature_bc_matrix_filtered.h5:/matrix')
			cellbender_out.write('\n\n')
	elif cellbender_file == None:
		print("running cellbender on files")
		# subprocess.run(cellbender_cmd, capture_output=True)
	
# note: to process through seurat, need to reformat using PyTables
#  ptrepack --complevel 5 cb_feature_bc_matrix_filtered.h5:/matrix cb-seurat_feature_bc_matrix_filtered.h5:/matrix	


print("Done")

