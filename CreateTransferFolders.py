# 2025.12.12
# This will copy files from the outs folder only


# Import libraries

import os
import re
from pathlib import Path
import shutil
from collections import defaultdict


sc_dir = "./testFolder/"
transfer_dir = "./summary_transfer"
transfer_files = ['web_summary', 'metrics_summary', 'cb_']
rename_files = ['web_summary', 'metrics_summary']
# req_outs_folder = True
search_term = 'scrna'
ignore_folders = ["raw_feature_bc_matrix", "analysis", " SC_RNA_COUNTER_CS"]

# Compile search_terms
search_term = re.compile(search_term)
target_pattern = re.compile('|'.join(transfer_files))
rename_pattern = re.compile('|'.join(rename_files))
ignore_folders = re.compile('|'.join(ignore_folders))

# Create the transfer folder
Path(os.path.join(transfer_dir, 'web_summaries')).mkdir(parents = True, exist_ok= True)

# start the list of things to transfer
sample_list=defaultdict(dict)

for sampleID in os.listdir(sc_dir):
    if os.path.isdir(os.path.join(sc_dir, sampleID, "outs")):
        current_path = os.path.join(sc_dir, sampleID, "outs")
        for fname in os.listdir(current_path):
            if target_pattern.search(fname):
                print(f'copying {fname}')
                target_path=os.path.join(transfer_dir, sampleID, "outs")
                Path(target_path).mkdir(parents = True, exist_ok=True)
                shutil.copy2(os.path.join(current_path, fname), os.path.join(target_path, fname))
                # rename test
                if rename_pattern.search(fname):
                    new_fname='_'.join((sampleID, fname))
                    os.rename(os.path.join(target_path, fname), os.path.join(target_path, new_fname))
                    if 'web_summary' in new_fname:
                        os.rename(os.path.join(target_path, new_fname), os.path.join(transfer_dir, 'web_summaries', new_fname))


