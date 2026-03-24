# 2025.12.05
# Code is meant to generate the cellranger.swarm file to be run on Biowulf


# Import libraries

import os
import re
import logging

from collections import defaultdict
from fx_logging_utils import setup_logging
from fx_cellranger_command import call_cellranger_command, arc_library_csv
from fx_hpap_rename import hpap_rename
from cellranger_argparse import get_args


def main():
    args = get_args()
    fastq_dir = args.fastq_dir
    sample_id_format = args.sample_id_format
    logfile = args.logfile
    data_type = args.data_type 
    swarmfile_name = args.swarmfile_name
    allow_loose_match = args.allow_loose_match
    path_restrictions = args.path_restrictions

    dry_run = args.dry_run
    debug = args.debug
    run_mode = dry_run or debug


    #Set up logging
    setup_logging(logfile, debug=debug)
    logger = logging.getLogger(__name__)
    logger.info("Starting Cellranger_swarm", 
                extra={'console_only': True}
    )

    #define data type
    if data_type=='atac':
        cellranger_module = 'cellranger-atac'
        ref_genome_cmd = "--reference=/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2024-A"
        logger.info("set atac variables")
    elif data_type=='rna':
        cellranger_module = 'cellranger'
        ref_genome_cmd = "--transcriptome=$CELLRANGER_REF/refdata-gex-GRCh38-2024-A"
        logger.info("set rna variables")
    elif data_type=='multi':
        cellranger_module = 'cellranger-arc'
        ref_genome_cmd = "--reference=/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2024-A"
        logger.info("set multi variables")
    else:
        logger.error("Error: data_type must be 'rna', 'atac', or 'multi'")
        return

    logger.info(f"""Running Cellranger_swarm.py with the following parameters:
        Fastq directory: {fastq_dir}
        Sample ID format: {sample_id_format}
        Logfile: {logfile}
        Data type: {data_type}
        Swarmfile name: {swarmfile_name}
        Dry run: {dry_run}
        Allow loose match: {allow_loose_match}
        Path restrictions: {path_restrictions}\n""",
        extra={'console_only': True}
    )
    # Generate dynamic filename matching pattern
    internal_tracking = r"(?:_[^_]+)*" if allow_loose_match else ""

    # fastq_file_prefix = re.compile(rf"^(?:{sample_id_core}).*\.fastq\.gz$")
    fastq_file_pattern = re.compile(
    rf"""^(?P<sampleID>{sample_id_format})
        (?P<internalTracking>{internal_tracking})
        _S(?P<chipsample>\d{{1,2}})
        _(?P<laneID>L\d{{3}})
        _(?P<readID>R1|R2|I1|I2)
        _001\.fastq\.gz$
    """,
    re.X,
    )
    if path_restrictions:
        logger.info(f"Applying path restrictions: {path_restrictions}")
        path_restrictions = re.compile(path_restrictions)     
    
    logger.debug("Using sample_id_pattern: %s", sample_id_format)
 
    # Create fastq list
    samples = defaultdict(lambda: defaultdict(dict))

    for dirpath, _, filenames in os.walk(fastq_dir):
        for fname in filenames:
            if path_restrictions and not re.search(path_restrictions, dirpath):
                continue
            if not fname.endswith('.fastq.gz'):
                continue
            m = fastq_file_pattern.match(fname)
            if not m and args.allow_renaming:
                # apply renaming rules
                new_name = hpap_rename(
                    fname,
                    dirpath,
                    run_mode=run_mode
                )
                if new_name:
                    fname = new_name
                    m = fastq_file_pattern.match(fname)
            if m:
                logger.debug("Matched sampleID: %s, %s", m.group("sampleID"), fname)
                sampleID = m.group("sampleID")
                internal_tracking = m.group("internalTracking")
                readID = m.group("readID")
                full_path = os.path.join(dirpath, fname)
                samples[sampleID][internal_tracking][readID] = full_path
            else:
                logger.debug("No match")
                logger.warning(f"Skipping {fname} after renaming attempt; still does not match pattern.")
                continue



    # filter for complete sets
    required_reads = {"R1", "R2", "I1"}
    complete_samples = {}
    for sampleID, library_id in samples.items():
        complete_libraries={}
        for tracking, reads in library_id.items():
            if required_reads.issubset(reads.keys()):
                complete_libraries[tracking] = reads
            else:
                logger.warning(
                    "incomplete library for %s_%s: found %s",
                    sampleID, tracking, list(reads.keys())
                )
        if complete_libraries:
            logger.debug("Complete libraries for sample %s: %s", sampleID, list(complete_libraries.keys()))
            complete_samples[sampleID] = complete_libraries

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

        # Generate dict of results
        for sampleID, library_id in samples.items():
            sample_names = [
                f"{sampleID}{track}" if track else sampleID
                for track in library_id.keys()
            ]
            sample_names = sorted(sample_names)
            library_tracker = next(iter(library_id.values()))
            read_path = next(iter(library_tracker.values()))
            sample_path = os.path.dirname(read_path)
            
            # Generate sample argument based on data_type
            if data_type == 'rna':
                sample_arg = ",".join(sample_names)
                sample_arg = f"--sample={sample_arg}"
            elif data_type == 'multi':
                logger.info("Generating library CSV for multiome sample %s", sampleID)
                csv_path = arc_library_csv(
                    sampleID = sampleID, 
                    GEX_id = sampleID+"_GEX", 
                    ATAC_id = sampleID+"_ATAC", 
                    sample_path = sample_path
                )
                sample_arg = f"--library={csv_path}"
        # Write swarm file for sample
            with open(swarmfile_name, 'a') as swarmfile:
                swarmfile.write(call_cellranger_command(
                    output_ID = sampleID, 
                    sample_path = sample_path, 
                    ref_genome_cmd = ref_genome_cmd, 
                    sample_args = sample_arg,
                    cellranger_module = cellranger_module))

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


