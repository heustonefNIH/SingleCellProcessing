# 2025.12.05
# Code is meant to generate the cellranger.swarm file to be run on Biowulf


# Import libraries

import os
import re
from collections import defaultdict
from fx_cellranger_command import call_cellranger_command

## REQUIRED UPDATES
# Currenlty only works with properly labeled files and files with the fastq-data.fastq.gz naming convention.

#Define variables

fastq_dir = "/data/CRGGH/heustonef/huMuscle/fastq_files"
swarmfile_name = 'test.txt'
skippedsamplesfile_name = 'SkippedSamples.txt'
sequencedata_type = 'rna'



if re.search('atac', sequencedata_type, re.IGNORECASE):
    cellranger_module = 'cellranger-atac'
    ref_genome_cmd = "--reference=/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2024-A"

if re.search('rna', sequencedata_type, re.IGNORECASE):
    cellranger_module = 'cellranger'
    ref_genome_cmd = "--transcriptome=$CELLRANGER_REF/refdata-gex-GRCh38-2024-A"

# Create swarm file
with open(swarmfile_name, 'w') as swarmfile:
    swarmfile.write(statement :=' '.join(('#swarm -f', swarmfile_name, ' -g 64 -t 12 --time=48:00:00 --merge-output --module', cellranger_module, '--sbatch \"--mail-type=BEGIN,END,FAIL\"\n\n')))

# Generate dynamic filename matching pattern
fastq_pattern = re.compile(
    r"""^(?P<sampleID>.+?) # match sample id
    _S(?P<chipsample>\d{1,2}) # match illumina sample ID in standard nomenclature
    _.*_ # match all the lane stuff
    (?P<readID>R1|R2|I1|I2) #match read type
    .* # match anythign between lane and the file type
    \.fastq\.gz$ # make sure it's a fastq file
    """,
    re.X
)

# Create fastq list

samples = defaultdict(dict)
for dirpath, dirnames, filenames in os.walk(fastq_dir):
    for fname in filenames:
        if not fname.endswith(".fastq.gz"):
            continue
        m = fastq_pattern.match(fname)
        if not m:
            continue
        sampleID = m.group("sampleID")
        readID = m.group("readID")
        full_path = os.path.join(dirpath, fname)

        samples[sampleID][readID] = full_path


# filter for complete sets
required_reads = {"R1", "R2", "I1", "I2"}
complete_samples = {
    sampleID: paths
    for sampleID, paths in samples.items()
    if required_reads.issubset(paths.keys())
    
}

# Print results

for sampleID, paths in sorted(complete_samples.items()): 
    read_path = next(iter(paths.values()))
    sample_path = os.path.dirname(read_path)
# Write swarm file for sample
    with open(swarmfile_name, 'a') as swarmfile:
        swarmfile.write(call_cellranger_command(sampleID, sample_path, ref_genome_cmd, cellranger_module))
missing_samples = set(samples) - set(complete_samples)
if len(missing_samples) > 0:
    with open(skippedsamplesfile_name, 'w') as skippedsamples:
        skippedsamples.write('The following sample IDs were not written to the swarm file:\n')
    for incomplete_sample in missing_samples:
        print("Missing entries for", incomplete_sample, "is not interpretable by current version of MakeSwarmFile_cellranger.py")
        with open(skippedsamplesfile_name, 'a') as skippedsamples:
            skippedsamples.write(''.join((incomplete_sample, '\n')))

if os.path.isfile(swarmfile_name):
    print(''.join(("Created ", swarmfile_name)))
else:
    print("Swarmfile creation failed.")





