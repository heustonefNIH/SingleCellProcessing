# 2022.12.20
# This will replicate the folder structure of each sample and transfer only necessary files


# Import libraries

import os
import re
from pathlib import Path
import shutil


sc_dir = "/data/CRGGH/heustonef/huMuscle"
target_dir = "/data/CRGGH/heustonef/huMuscle/transfer_folder"
additional_files = ['web_summary.html', 'metrics_summary.csv']

search_term = 'EH0'

os.chdir(sc_dir)

# Compile search_term
search_term = re.compile(search_term)

# Get a list of samples to run
sample_list = []
for sampleDir in os.listdir(sc_dir): # for each sample ID
    if re.search(search_term, sampleDir) is not None: 
        print(sampleDir)
        filteredData_dir = re.compile('^filtered.+matrix$')
        filteredData_dir = list(filter(filteredData_dir.match, os.listdir(sampleDir)))
        if len(filteredData_dir) != 1:
            print(''.join(("Multiple matches for filtered...matrix dir found in ", sampleDir, "... skipping")))
        else:
            filteredData_dir = filteredData_dir[0]
            target_patterns = ['h5$', '^barcodes', '^features', '^matrix'] + additional_files
            target_files = re.compile('|'.join(target_patterns))
            #target_path = os.path.join(target_dir, sampleDir, "outs", filteredData_dir)
            target_path = os.path.join(target_dir, sampleDir, filteredData_dir)
            Path(target_path).mkdir(parents = True, exist_ok = True)
            # Copy run information
            for file in (file for file in os.listdir(sampleDir) if file.startswith('_')):
                shutil.copy2(os.path.join(sampleDir, file), os.path.join(target_dir, sampleDir, file))
            # Copy target files
            for dirpath, subdirs, files in os.walk(sampleDir):
                for file in files:
                    if target_files.match(file) and  dirpath in target_path:
                        shutil.copy2(os.path.join(dirpath, file), os.path.join(target_dir, dirpath, file))
                        if re.search(file, "web_summary.html") and not file == "web_summary.html":
                            os.rename(os.path.join(target_dir, dirpath, file), os.path.join(target_dir, dirpath,''.join((sampleDir, '-', file))))
                        if re.search(file, "metrics_summary.csv") and not file == "metrics_summary.csv":
                            os.rename(os.path.join(target_dir, dirpath, file), os.path.join(target_dir, dirpath,''.join((sampleDir, '-', file))))


