# 2022.12.20
# This will replicate the folder structure of each sample and transfer only necessary files


# Import libraries

import os
import re
from pathlib import Path
import shutil


sc_dir = "/data/CRGGH/heustonef/hpapdata/cellranger_scRNA/"
target_dir = "/data/CRGGH/heustonef/hpapdata/cellranger_scRNA/scRNA_transfer"
target_files = ['web_summary.html', 'metrics_summary.csv']

search_term = 'EH00'
sequencedata_type = "rna"


os.chdir(sc_dir)

# Compile search_term
search_term = re.compile(search_term)

# Get a list of samples to run
sample_list = []
for sampleDir in os.listdir(sc_dir): # for each sample ID
    if re.search(search_term, sampleDir) is not None: 
        print(sampleDir)
        with open(os.path.join(sampleDir, "_cmdline"), 'r') as cmdline:
            if sequencedata_type is not None:
                cellranger_cmd = cmdline.readline().split(' ')[0]
                if re.search('cellranger-atac', cellranger_cmd):
                    sequencedata_type = 'atac'
                    filteredData_dir = "filtered_peak_bc_matrix"
                    target_files = target_files + ['filtered_peak_bc_matrix.h5', 'barcodes.tsv', 'matrix.mtx', 'peaks.bed']
                    print("Reading ATAC files")
                elif re.search('cellranger$', cellranger_cmd):
                    sequencedata_type = 'rna'
                    filteredData_dir = "filtered_feature_bc_matrix"
                    target_files = target_files + ['filtered_feature_bc_matrix.h5', 'barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz']
                    print("Reading RNA files")
                else:
                    print("Could not determine sequencedata_type")
        target_path = os.path.join(target_dir, sampleDir, "outs", filteredData_dir)
        Path(target_path).mkdir(parents = True, exist_ok = True)
        # Copy run information
        for file in (file for file in os.listdir(sampleDir) if file.startswith('_')):
            shutil.copy2(os.path.join(sampleDir, file), os.path.join(target_dir, sampleDir, file))
        # Copy target files
        for dirpath, subdirs, files in os.walk(sampleDir):
            for file in files:
                if re.search(file, target_files) and  dirpath in target_path:
                    shutil.copy2(os.path.join(dirpath, file), os.path.join(target_dir, dirpath, file))
                    if re.search(file, "web_summary.html") and not file == "web_summary.html":
                        os.rename(os.path.join(target_dir, dirpath, file), os.path.join(target_dir, dirpath,''.join((sampleDir, '-', file))))
                    if re.search(file, "metrics_summary.csv") and not file == "metrics_summary.csv":
                        os.rename(os.path.join(target_dir, dirpath, file), os.path.join(target_dir, dirpath,''.join((sampleDir, '-', file))))


