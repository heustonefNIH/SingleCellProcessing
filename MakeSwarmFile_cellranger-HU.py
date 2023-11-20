# 2023.09.22
# Code is meant to loop through folders and generate the cellranger.swarm file to be run on Biowulf HPC


# Import libraries

import os
import re

#Define variables

fastq_dir = "/data/CRGGH/heustonef/huMuscle/fastq/"
swarmfile_name = 'huMuscle.swarm'
skippedsamplesfile_name = 'SkippedSamples.txt'
sequencedata_type = 'rna'

if re.search('atac', sequencedata_type, re.IGNORECASE):
    cellranger_module = 'cellranger-atac'
    ref_genome_cmd = "--reference=/fdb/cellranger-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"

if re.search('rna', sequencedata_type, re.IGNORECASE):
    cellranger_module = 'cellranger'
    ref_genome_cmd = "--transcriptome=$CELLRANGER_REF/refdata-gex-GRCh38-2020-A"


os.chdir(fastq_dir)


# Create swarm file
with open(swarmfile_name, 'w') as swarmfile:
    swarmfile.write(statement :=' '.join(('#swarm -f', swarmfile_name, ' -g 64 -t 12 --time=48:00:00 --merge-output --module', cellranger_module, '--sbatch \"--mail-type=BEGIN,END,FAIL\"\n\n')))

with open(skippedsamplesfile_name, 'w') as skippedsamples:
    skippedsamples.write('The following HPAP sample IDs were not written to the swarm file:\n')

# Get a list of samples to run
sample_list = []
for sampleDir in os.listdir(hpap_dir): # for each sample ID
    if re.search(r'^EH-\d{3}$', sampleDir) is not None: # Make sure it follows the HPAP nomenclature
        sample_list.append(sampleDir)    #And add it to the lit of sample IDs to read

# GLoop through directory and collect *fastq.gz files
for sampleDir in sample_list: # for each sample id
    fastq_list = [] #Start a new list for each sample ID
    batches = set()
    for sampleID, subdirs, files in os.walk(sampleDir):
        if re.search(sequencedata_type, sampleID, re.IGNORECASE) is not None:
            if sampleID.endswith('fastq') in sampleID:
                for file in (file for file in files if not file.startswith('.')): #for each fastq file
                    if file.endswith('.fastq.gz') and re.search('L\d{3}', file):    # Check it's nomenclature and fix it if necessary
                        # print('renaming file', file)
                        # origfile = file
                        # if not re.search('_S\d+_L\d+', file) and re.search('\.L\d{3}\.', file): # if file is in "<SampleID>.L00X" and does not contain "_S\d+_L\d+" format
                        #     file = re.sub(r'\.(L\d{3})\.', r'_S1_\1_', file)
                        # if not file.endswith('001.fastq.gz') and 'fastq-data.fastq.gz' in file: # if file does not end as "001.fastq.gz"
                        #     file = re.sub(r'fastq-data.fastq.gz', '001.fastq.gz', file)
                        # if file.endswith('001.fastq.gz') and 'fastq-data.fastq.gz' in file:
                        #     file = re.sub(r'fastq-data.fastq.gz', '.fastq.gz', file)
                        # if re.search(r'_10xscRNA_', file):
                        #     file = re.sub(r'_10xscRNA_', '_', file)
                        # os.rename(os.path.join(sampleID, origfile), os.path.join(sampleID, file))
                        fastq_list.append(file) # Now all the fastq files in sampleDir dir are in a list
            # Get batch IDs
                    if len(fastq_list) > 0:
                        if ([m := re.search('.*_S\d+_L\d+', x) for x in fastq_list] and m is not None):
                            batchID = [re.search('EH-?\d{3}(_.*?)(?=_S\d+_L\d+)', x).group(1) for x in fastq_list]
                        else:
                            print("Could not find batch nomenclature in sampleDir", sampleDir)
                        batches.update(batchID)
        # Make swarm file
                if len(batches) > 0:
                    sampleID = sampleID.replace(" ", "\ ")
                    for batch_id in batches:
                        sample = ",".join(set([re.match('.*(?=_S\d+_L\d+)', x).group(0) for x in fastq_list if batch_id in x]))
                    # Write swarm file for sample
                        with open(swarmfile_name, 'a') as swarmfile:
                            swarmfile.write(''.join(('FASTQ_PATH=', sampleID, '; \\\n')))
                            swarmfile.write('ulimit -u 10240 -n 16384; \\\n')
                            swarmfile.write(''.join((cellranger_module, ' count --id=', ''.join((sampleDir, batch_id)), ' \\\n')))
                            swarmfile.write(''.join(('\t ', ref_genome_cmd, ' \\\n')))
                            swarmfile.write(''.join(('\t --fastqs=\"$FASTQ_PATH\" \\\n')))
                            swarmfile.write(''.join(('\t --sample=', sample, ' \\\n')))
                            swarmfile.write(''.join(('\t --localcores=$SLURM_CPUS_PER_TASK \\\n')))
                            swarmfile.write(''.join(('\t --localmem=34 \\\n')))
                            swarmfile.write(''.join(('\t --jobmode=slurm \\\n')))
                            swarmfile.write(''.join(('\t --maxjobs=10')))
                            swarmfile.write('\n\n')
                else:
                    print(sampleID)
                    print("Fastq file name format in", sampleDir, "is not interpretable by current version of MakeSwarmFile_cellranger.py")
                    with open(skippedsamplesfile_name, 'a') as skippedsamples:
                        skippedsamples.write(''.join((sampleDir, '\n')))

if os.path.isfile(swarmfile_name):
    print(''.join(("Created ", swarmfile_name)))
else:
    print("Swarmfile creation failed.")





