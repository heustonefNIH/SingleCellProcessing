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
    logger.info(f"""Running cellbender_swarm with the following parameters:
    cellranger_folder: {cellranger_folder}
    sample_id_format: {sample_id_format}
    cellbender_file: {cellbender_file}
    Cellbender flags:
    cuda: {cuda}
    gpu_partition: {gpu_partition}
    {flags}
    
    debug: {debug}
    logfile: {logfile}
    """, 
    extra = {'console_only': True}
    )

    if gpu_partition:
        gres_flag = f"--gres={gpu_partition}"
        if cuda: 
            cuda = "--cuda"
            gres_flag = "".join(["--cuda \\", "\n--gres=", gpu_partition])
    if sample_id_format:
        sample_id_regex = re.compile(sample_id_format)
    else:
        sample_id_regex = None
    logger.debug(f"Using sample ID regex: {sample_id_regex}")
    
    swarm_statement=(
        f'#swarm -f {cellbender_file} -g 64 '
        f'--time=24:00:00 --gres={gpu_partition} -t 8 '
        '--merge-output --module cellbender '
        '--sbatch "--mail-type=BEGIN,END,FAIL"\n\n')

    # Create list of raw_feature_bc_matrix.h5 files
    sample_list = defaultdict(dict)
    for sampleID in os.listdir(cellranger_folder):
        if not sample_id_regex.search(sampleID):
            logger.debug(f"Skipping {sampleID} - does not match sample ID format")
            continue
        if not os.path.isdir(os.path.join(cellranger_folder, sampleID, "outs")):
            logger.debug(f"Skipping {sampleID} - no outs directory")
            continue  # skip non-directories and directories without outs folder
        sample_path = os.path.join(cellranger_folder, sampleID)
        outs_dir = os.path.join(sample_path, "outs")
        
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
        logger.info(f"Found {len(sample_list)} samples to process. Generating {cellbender_file}...", 
                    extra={'console_only': True})
        with open (cellbender_file, 'w') as cellbender:
            cellbender.write(swarm_statement)
        # Write swarm file
        for sampleID, matrixFile in sorted(sample_list.items()):
            for matrix_path, outfile in matrixFile.items():
                cd_path = os.path.dirname(matrix_path)
                logger.info(f"Processing sample {sampleID} with matrix file {matrix_path}")
                cellbender_cmd=(
                    f"#Sample {sampleID}\n"
                    f"cd {cd_path}; \\\n"
                    f"cellbender remove-background \\\n"
                    f"--input {matrix_path} \\\n"
                    f"--output {outfile} \\\n"
                    f"{cuda} \\\n"
                    f"{flags}; \\\n"
                    f"ptrepack --complevel 5 cb_feature_bc_matrix_filtered.h5:/matrix cb_seurat_feature_bc_matrix_filtered.h5:/matrix\n"
                    f"\n\n"
                    )
            with open (cellbender_file, 'a') as cellbender:
                cellbender.write(cellbender_cmd)

    if os.path.isfile(cellbender_file):
        logger.info(''.join(("Created ", cellbender_file)), 
                    extra={'console_only': True}
        )
    else:
        logger.error(''.join(("Failed to create ", cellbender_file)), 
                    extra={'console_only': True}
        )
if __name__ == "__main__":
    main()