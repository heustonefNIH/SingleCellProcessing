# 2025.12.11 CellBender Script
# Assumes you're using raw_feature_bc_matrix.h5 files output from CellRanger Count (v10.0.0 as of 2025.12.11)
# IO https://cellbender.readthedocs.io/en/latest/usage/index.html


# Import modules
import os
from collections import defaultdict

# Define global variables

#swarm -f cellbender.swarm -g 64 -t 8 --time=24:00:00 --merge-output --module cellbender --sbatch "--mail-type=BEGIN,END,FAIL"
cellranger_folder = "/data/CRGGH/heustonef/huMuscle/aadm_x18/" # Required: folder containing cellranger output data
cellbender_file = "cellbender.swarm" # Optional: None, or name of file to write cellbender commands. Useful if submitting as swmarm file
gpu_run = "--cuda"
flags = "--cpu-threads $SLURM_CPUS_PER_TASK" # Optional: CellBender flags to include in command

swarm_statement=f"""
#swarm -f {cellbender_file} -g 64 -t 8 --time=24:00:00 --partition=gpu --gres=gpu:XX:X --merge-output --module cellbender --sbatch "--mail-type=BEGIN,END,FAIL"
"""



if cellbender_file:
	cellbender_file = os.path.join(cellranger_folder, cellbender_file)
	with open(cellbender_file, 'w') as cellbender:
		cellbender.write(swarm_statement)
if gpu_run is True:
		flags = ''.join((" --cuda ", flags))

# Create list of raw_feature_bc_matrix.h5 files
sample_list = defaultdict(dict)
for sampleID in os.listdir(cellranger_folder):
    sample_path = os.path.join(cellranger_folder, sampleID)
    if not os.path.isdir(sample_path):
        continue  # skip non-directories

    outs_dir = os.path.join(sample_path, "outs")
    if not os.path.isdir(outs_dir):
        continue  # skip samples without outs

    # Now look directly in the outs folder
    for fname in os.listdir(outs_dir):
        if fname == "raw_feature_bc_matrix.h5":
            matrix_path = os.path.join(outs_dir, fname)
            sample_list[sampleID][matrix_path] = os.path.join(
                outs_dir, "cb_feature_bc_matrix.h5"
            )

for sampleID, mapping in sorted(sample_list.items()):
	for paths, outfile in mapping.items():
		cellbender_cmd=f"""#Sample {sampleID}
cd {paths}; \\
cellbender remove-background {flags} \\
{gpu_run} \\
--input {paths} \\
--output {outfile}; \\
ptrepack --complevel 5 cb_feature_bc_matrix_filtered.h5:/matrix cb-seurat_feature_bc_matrix_filtered.h5:/matrix

"""
	with open (cellbender_file, 'a') as cellbender:
		cellbender.write(cellbender_cmd)

print("Done")

