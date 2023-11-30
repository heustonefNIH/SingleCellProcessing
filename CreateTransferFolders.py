# 2022.12.20
# This will replicate the folder structure of each sample and transfer only necessary files


# Import libraries

import os
import re
from pathlib import Path
import shutil


sc_dir = "/Users/heustonef/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/SingleCellMetaAnalysis/GitRepositories/SingleCellProcessing/testFolder"
transfer_dir = "/Users/heustonef/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/SingleCellMetaAnalysis/GitRepositories/SingleCellProcessing/transfer_folder"
additional_files = ['web_summary.html', 'metrics_summary.csv']
req_outs_folder = True
search_term = 'EH0'
ignore_folders = ["raw_feature_bc_matrix", "analysis"]

os.chdir(sc_dir)
target_patterns = ['h5$', '^barcodes', '^features', '^matrix', '^_'] + additional_files # '^_' copies run information

# Compile search_terms
search_term = re.compile(search_term)
target_files = re.compile('|'.join(target_patterns))
ignore_folders = re.compile('|'.join(ignore_folders))

# Get a list of samples to run
for sampleDir in os.listdir(sc_dir): # for each sample ID
    if req_outs_folder == True and not os.path.exists(os.path.join(sampleDir, "outs")):
        continue # If we require an "outs" folder and none exists, skip this round of "for sampleDir" and move to the next one
    if re.search(search_term, sampleDir) is not None: 
        print(sampleDir)
        filteredData_dir = re.compile('^filtered.+matrix$')
        filteredData_dir = list(filter(filteredData_dir.match, [os.path.basename(x[0]) for x in os.walk(sampleDir)]))
        if len(filteredData_dir) != 1:
            print(''.join(("Multiple matches for filtered...matrix dir found in ", sampleDir, "... skipping")))
        else:
            sampleDir_transfer = os.path.join(transfer_dir, sampleDir)
            Path(sampleDir_transfer).mkdir(parents = True, exist_ok = True)
            # Copy target files
            # right way to do this is to make a list of all files to be transfered, then transfer everything in the list
            transfer_files = []
            for root, subdirs, files in os.walk(sampleDir):
                if re.search(ignore_folders, root):
                    continue
                for file in files:
                    if target_files.search(file) and not file.startswith("\."):
                        transfer_files.append(os.path.join(root, file))
        print("finished making file list")
        # Recreate folder structure in transfer folder
        for file in transfer_files:
            dir_structure = list(set([os.path.dirname(x) for x in transfer_files]))
            for file_path in dir_structure:
                file_path = os.path.join(transfer_dir, file_path)
                Path(file_path).mkdir(parents=True, exist_ok=True)
        for file in transfer_files:
            shutil.copy2(os.path.join(sc_dir, file), os.path.join(transfer_dir, file))
            if re.match("web_summary\.html|metrics_summary\.csv", os.path.basename(file)):
                os.rename(os.path.join(transfer_dir, file), os.path.join(transfer_dir, os.path.dirname(file), ''.join((sampleDir, '-', os.path.basename(file)))))
    print("finished first for loop")


