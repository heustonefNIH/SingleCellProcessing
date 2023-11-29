# 2023.09.22
# Code is meant to loop through folders and generate the cellranger.swarm file to be run on Biowulf HPC


# Import libraries

import os
import re

#Define variables
fastq_dir = "/Users/heustonef/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/scRNA_ATAC-Seq/huMuscle/testFolder"
swarmfile_name = 'huMuscle-RNA.swarm'
match_pattern = '^EH\d{3}'
ignore_folders = "AACT5JKM5"
skippedsamplesfile_name = 'SkippedSamples.txt'
sequencedata_type = 'rna'

if re.search('atac', sequencedata_type, re.IGNORECASE):
	cellranger_module = 'cellranger-atac'
	ref_genome_cmd = "--reference=/fdb/cellranger-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"

if re.search('rna', sequencedata_type, re.IGNORECASE):
	cellranger_module = 'cellranger'
	ref_genome_cmd = "--transcriptome=$CELLRANGER_REF/refdata-gex-GRCh38-2020-A"


os.chdir(fastq_dir)
match_pattern_compiled = re.compile(match_pattern)

# Create swarm file
with open(swarmfile_name, 'w') as swarmfile:
	swarmfile.write(statement :=' '.join(('#swarm -f', 
										  swarmfile_name, 
										  ' -g 64 -t 12 --time=48:00:00 --merge-output --module', 
										  cellranger_module, 
										  '--sbatch \"--mail-type=BEGIN,END,FAIL\"\n\n')))
with open(skippedsamplesfile_name, 'w') as skippedsamples:
	skippedsamples.write('The following HPAP sample IDs were not written to the swarm file:\n')

# Get a list of samples to run
sample_list = []
for sampleDir, subdirectories, files in os.walk(fastq_dir):
	for filename in files:
		if re.search(match_pattern, filename):
			sample_list.append(sampleDir)    #And add it to the lit of sample IDs to read
			break

# GLoop through directory and collect *fastq.gz files
for sampleDir in sample_list: # for each sample id
	fastq_list = [] #Start a new list for each sample ID
	for _, _, fastq_files in os.walk(sampleDir):
		fastq_list.extend([x for x in fastq_files if x.endswith('fastq.gz')]) # Now all the fastq files in sampleDir dir are in a list
		# Get batch IDs
		if len(fastq_list) > 0:
			#Removing section on batch identification for now, since it's not relevant to getting this running--See "MakeSwarmFile_cellranger.py" for original code
			sample_set = set([re.match('.*(?=_S\d+_L\d+)', x).group(0) for x in fastq_list])
			# Write swarm file for sample
			for sample in sample_set:
				with open(swarmfile_name, 'a') as swarmfile:
					swarmfile.write(''.join(('FASTQ_PATH=', sampleDir, '; \\\n')))
					swarmfile.write('ulimit -u 10240 -n 16384; \\\n')
					swarmfile.write(''.join((cellranger_module, ' count --id=', sample, ' \\\n')))
					swarmfile.write(''.join(('\t ', ref_genome_cmd, ' \\\n')))
					swarmfile.write(''.join(('\t --fastqs=\"$FASTQ_PATH\" \\\n')))
					swarmfile.write(''.join(('\t --sample=', sample, ' \\\n')))
					swarmfile.write(''.join(('\t --expect-cells=REPLACE \\\n'))) # specifically to manually insert expected cell number for troubleshooting run
					swarmfile.write(''.join(('\t --localcores=$SLURM_CPUS_PER_TASK \\\n')))
					swarmfile.write(''.join(('\t --localmem=34 \\\n')))
					swarmfile.write(''.join(('\t --jobmode=slurm \\\n')))
					swarmfile.write(''.join(('\t --maxjobs=10')))
					swarmfile.write('\n\n')
			else:
				print(sample)
				print("Fastq file name format in", sample, "is not interpretable by current version of MakeSwarmFile_cellranger-HU.py")
				with open(skippedsamplesfile_name, 'a') as skippedsamples:
					skippedsamples.write(''.join((sample, '\n')))

if os.path.isfile(swarmfile_name):
	print(''.join(("Created ", swarmfile_name)))
else:
	print("Swarmfile creation failed.")





