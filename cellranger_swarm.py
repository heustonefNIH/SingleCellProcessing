# 2025.12.05
# Code is meant to generate the cellranger.swarm file to be run on Biowulf


# Import libraries

import os
import re
import logging

from collections import defaultdict
from fx_logging_utils import setup_logging
from fx_cellranger_command import call_cellranger_command
from fx_hpap_rename import hpap_rename
from cellranger_argparse import get_args


def main():
    args = get_args()
    fastq_dir = args.fastq_dir
    sample_id_format = args.sample_id_format
    logfile = args.logfile
    sequencedata_type = args.data_type 
    swarmfile_name = args.swarmfile_name

    dry_run = args.dry_run
    allow_loose_match = args.allow_loose_match
    path_restrictions = args.path_restrictions

    #Set up logging
    setup_logging(logfile, debug=args.debug)
    logger = logging.getLogger(__name__)
    logger.info("Starting Cellranger_swarm", 
                extra={'console_only': True}
    )

    #define data type
    if sequencedata_type=='atac':
        cellranger_module = 'cellranger-atac'
        ref_genome_cmd = "--reference=/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2024-A"
        logger.info("set atac variables")
    elif sequencedata_type=='rna':
        cellranger_module = 'cellranger'
        ref_genome_cmd = "--transcriptome=$CELLRANGER_REF/refdata-gex-GRCh38-2024-A"
        logger.info("set rna variables")
    else:
        logger.error("Error: data_type must be 'rna' or 'atac'")
        return

    logger.info(f"""Running Cellranger_swarm.py with the following parameters:
        Fastq directory: {fastq_dir}
        Sample ID format: {sample_id_format}
        Logfile: {logfile}
        Data type: {sequencedata_type}
        Swarmfile name: {swarmfile_name}
        Dry run: {dry_run}
        Allow loose match: {allow_loose_match}
        Path restrictions: {path_restrictions}\n""",
        extra={'console_only': True}
    )
    # Generate dynamic filename matching pattern
    sample_id_core = sample_id_format.lstrip("^")
    internal_tracking = r"(?:_[^_]+)*" if allow_loose_match else ""

    fastq_file_prefix = re.compile(rf"^(?:{sample_id_core}).*\.fastq\.gz$")
    fastq_file_pattern = re.compile(
    rf"""^(?P<sampleID>{sample_id_core})
        (?P<internalTracking>{internal_tracking})
        _S(?P<chipsample>\d{{1,2}})
        _(?P<laneID>L\d{{3}})
        _(?P<readID>R1|R2|I1|I2)
        .*\.fastq\.gz$
    """,
    re.X,
    )      
    
    # Adjust sample_id_format --allow_loose_match
    if allow_loose_match:
        sample_id_pattern = sample_id_core
    else:
        sample_id_pattern = rf"{sample_id_core}{internal_tracking}"
    
    logger.debug("Using sample_id_pattern: %s", sample_id_pattern)
 
    # Create fastq list
    samples = defaultdict(dict)

    for dirpath, _, filenames in os.walk(fastq_dir):
        for fname in filenames:
            if path_restrictions and not re.search(path_restrictions, dirpath):
                continue
            if not fastq_file_prefix.match(fname):
                continue
            m = fastq_file_pattern.match(fname)
            if not m and args.allow_renaming:
                # apply renaming rules
                new_name = hpap_rename(
                    fname,
                    dirpath,
                    fastq_file_pattern,
                    dry_run=dry_run
                )
                if new_name:
                    fname = new_name
                    m = fastq_file_pattern.match(fname)
            if m:
                logger.debug("Matched sampleID: %s, %s", m.group("sampleID"), fname)
            else:
                logger.debug("No match")
                logger.warning(f"Skipping {fname} after renaming attempt; still does not match pattern.")
                continue
            sampleID = m.group("sampleID")
            readID = m.group("readID")
            full_path = os.path.join(dirpath, fname)

            samples[sampleID][readID] = full_path


    # filter for complete sets
    required_reads = {"R1", "R2", "I1"}
    complete_samples = {
        sampleID: paths
        for sampleID, paths in samples.items()
        if required_reads.issubset(paths.keys())
    }
    logger.info(f"{len(complete_samples)} complete samples found.", 
                extra={'console_only': True}
    )

    # Create swarm file
    if len(complete_samples) > 0:
         with open(swarmfile_name, 'w') as swarmfile:
                header = (
                    f"#swarm -f {swarmfile_name} -g 64 -t 12 --time=48:00:00 "
                    f'--merge-output --module {cellranger_module} --sbatch "--mail-type=BEGIN,END,FAIL"\n\n'
                )
                swarmfile.write(header)

    # Print results
    for sampleID, paths in sorted(complete_samples.items()): 
        read_path = next(iter(paths.values()))
        sample_path = os.path.dirname(read_path)
    # Write swarm file for sample
        with open(swarmfile_name, 'a') as swarmfile:
            swarmfile.write(call_cellranger_command(sampleID, sample_path, ref_genome_cmd, cellranger_module))

    # Log missing samples
    missing_samples = set(samples) - set(complete_samples)
    if len(missing_samples) > 0:
        logger.warning('The following sample IDs were not written to the swarm file:')
        for incomplete_sample in missing_samples:
            logger.warning(incomplete_sample)

    if os.path.isfile(swarmfile_name):
        logger.info(''.join(("Created ", swarmfile_name)), 
                    extra={'console_only': True}
        )
    else:
        logger.warning("Swarmfile creation failed.")



if __name__ == "__main__":
    main()


