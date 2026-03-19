import os
import re
import logging

def hpap_rename(fname, dirpath, target_pattern, logfile, dry_run):
    if not target_pattern.match(fname):
        # apply renaming rules
        new_name = fname
        new_name = re.sub(r'\.(L\d{3})\.', r'_S1_\1_', new_name)   # does nothing if no match
        new_name = re.sub('_fastq-data', '', new_name)
        new_name = re.sub('_10xscRNA_', '_', new_name)
        new_name = re.sub(r'HPAP(\d{3})', r'HPAP-\1', new_name)
        if new_name == fname:
            return None
        
        src = os.path.join(dirpath, fname)
        dst = os.path.join(dirpath, new_name)

        logger = logging.getLogger(__name__)
        logger.info(f"{'Dry run:' if dry_run else 'Renaming'} {os.path.basename(src)} -> {os.path.basename(dst)}")

    if not dry_run:
        os.rename(src, dst)

    return new_name

