# 2025.12.11 CellBender Script
# Assumes you're using raw_feature_bc_matrix.h5 files output from CellRanger Count (v10.0.0 as of 2025.12.11)


# Import libraries
import os
import logging
import re

from collections import defaultdict
from fx_logging_utils import setup_logging
from cellbender_argparse import get_args

def main():
    args = get_args()
    cellranger_folder = args.cellranger_folder
    sample_id_format = args.sample_id_format
    cellbender_file = args.cellbender_file
    cuda = args.cuda
    gpu_partition = args.gpu_partition
    flags = args.flags

    logfile = args.logfile
    debug = args.debug

    #Set up logging
    setup_logging(logfile, debug=debug)
    logger = logging.getLogger(__name__)
    logger.info("Starting CellBender Swarm File Generation", 
                extra={'console_only': True}
    )
    logger.info(f"""Running cellbender_swarm with the following parameters:
    cellranger_folder: {cellranger_folder}
    sample_id_format: {sample_id_format}
    cellbender_file: {cellbender_file}
    cuda: {cuda}
    gpu_partition: {gpu_partition}
    flags: {flags}
    debug: {debug}
    logfile: {logfile}
    """
    )

    if cellbender_file:
        cellbender_file = os.path.join(cellranger_folder, cellbender_file)
    if gpu_partition:
        gres_flag = f"--gres={gpu_partition}"
        if cuda: #add --cuda flag if cuda is True and gpu_partition is specified
            gres_flag = f"--cuda \\
                  {gpu_partition}"
    if sample_id_format:
        sample_id_regex = re.compile(sample_id_format)
    
    swarm_statement=f'swarm -f {cellbender_file} -g 32 
    --time=8:00:00 --partition=gpu {gres_flag} -t 8 
    --merge-output --module cellbender 
    --sbatch "--mail-type=BEGIN,END,FAIL"\n\n'

    # Create list of raw_feature_bc_matrix.h5 files
    sample_list = defaultdict(dict)
    for sampleID in os.listdir(cellranger_folder):
        if sample_id_regex and not sample_id_regex.search(sampleID):
            logger.debug(f"Skipping {sampleID} - does not match sample ID format")
            continue
        
        sample_path = os.path.abspath(sampleID)
        outs_dir = os.path.join(sample_path, "outs")
        
        if not os.path.isdir(sample_path) and not os.path.isdir(outs_dir):
            continue  # skip non-directories and directories without outs folder

        # Now look directly in the outs folder
        for fname in os.listdir(outs_dir):
            if fname == "raw_feature_bc_matrix.h5":
                matrix_path = os.path.join(outs_dir, fname)
                sample_list[sampleID][matrix_path] = os.path.join(
                    outs_dir, "cb_feature_bc_matrix.h5"
                )
                logger.debug(f"Found matrix file for sample {sampleID}: {matrix_path}")
                
    #only start swarm file if there's at least one sample to process
    if len(sample_list) == 0:
        logger.error("No samples found matching the specified sample ID format.")
        return
    else:
        with open (cellbender_file, 'a') as cellbender:
            cellbender.write(swarm_statement)
    # Write swarm file
    for sampleID, matrixFile in sorted(sample_list.items()):
        for paths, outfile in matrixFile.items():
            logger.info(f"Processing sample {sampleID} with matrix file {paths}")
            cellbender_cmd=f"""#Sample {sampleID}
    cd {paths}; \\
    cellbender remove-background {flags} \\
    {gres_flag} \\
    --input {paths} \\
    --output {outfile}; \\
    ptrepack --complevel 5 cb_feature_bc_matrix_filtered.h5:/matrix cb_seurat_feature_bc_matrix_filtered.h5:/matrix

    """
        with open (cellbender_file, 'a') as cellbender:
            cellbender.write(cellbender_cmd)

    print("Done")

if __name__ == "__main__":
    main()