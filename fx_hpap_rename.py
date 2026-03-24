import os
import re
import logging

logger = logging.getLogger(__name__)

def hpap_rename(fname, dirpath, run_mode):

    new_name = fname
    # apply renaming rules
    new_name = re.sub('_fastq-data', '', new_name)
    new_name = re.sub('_10xscRNA_', '_', new_name)
    new_name = re.sub(r'HPAP(\d{3})', r'HPAP-\1', new_name)
    new_name = re.sub(r'\.(L\d{3})\.', r'_\1_', new_name) 
    new_name = re.sub(r'\.(S\d+_L\d{3})', r'_\1', new_name) 
    new_name = re.sub(r'_(R\d|I\d)\.fastq\.gz$', r'_\1_001.fastq.gz', new_name)

    src = os.path.join(dirpath, fname)
    dst = os.path.join(dirpath, new_name)

    #Safety checks for renaming
    if not new_name:
        logger.error("Renaming %s creates an invalid filename(%r)", fname, new_name)
        return None
    if new_name in {".", ".."}:
        logger.error("Renaming %s creates an invalid filename(%r)", fname, new_name)
        return None
    
    if not new_name == fname:
        logger.info("Attempting to rename %s -> %s", src, dst)

    if src == dst:
        logger.debug("Source and destination are identical: %s", src)
        return None

    if os.path.exists(dst):
        logger.error("Cannot rename %s to %s because destination already exists.", src, dst)
        return None

    if run_mode:
        logger.info("DRY RUN: Would rename %s -> %s", src, dst)
        return new_name

    os.rename(src, dst)
    logger.info("Renamed %s -> %s", src, dst)
    return new_name

