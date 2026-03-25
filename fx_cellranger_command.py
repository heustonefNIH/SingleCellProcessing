import shlex
import os
import logging

logger = logging.getLogger(__name__)




def call_cellranger_command(output_ID, sample_path, ref_genome_cmd, sample_args, cellranger_module, data_type):
    if data_type == 'rna':
        sample_path = shlex.quote(sample_path)
        path_variables = f"""FASTQ_PATH={sample_path}; \\"""
        sample_args = f"""--fastqs="$FASTQ_PATH" \\
--samples={sample_args}"""
    elif data_type in ['atac', 'multi']:
        library_path = os.path.join(sample_path, f"{output_ID}_library.csv")
        library_path = shlex.quote(library_path)
        sample_path = shlex.quote(sample_path)
        path_variables = f"""LIBRARY_PATH={library_path}; \\"""
        sample_args = f"""--libraries=\"$LIBRARY_PATH\""""
    else:
        raise ValueError("Unsupported data_type: %s" % data_type)
    
    return f"""# SAMPLE {output_ID} 
{path_variables}
ulimit -u 10240 -n 16384; \\
{cellranger_module} count --id={output_ID} \\
{ref_genome_cmd} \\
{sample_args} \\
--create-bam=false \\
--localcores=$SLURM_CPUS_PER_TASK \\
--localmem=62 \\
--maxjobs=10

"""

def arc_library_csv(sampleID, sample_names, sample_path, gex_identifier = None, atac_identifier = None):
    sample_path = os.path.abspath(sample_path)
    csv_path = os.path.join(sample_path, f"{sampleID}_library.csv")
    
    sample_set = set(sample_names)
    if gex_identifier or atac_identifier:
        gex_identifier = gex_identifier or []
        atac_identifier = atac_identifier or []

        gex_samples = sorted(
            s for s in sample_set
            if any(token in s for token in gex_identifier)
        )
        atac_samples = sorted(
            s for s in sample_set
            if any(token in s for token in atac_identifier)
        )
    else:
        # fallback to original behavior
        gex_samples = sorted(
            s for s in sample_set
            if "GEX" in s
        )
        atac_samples = sorted(
            s for s in sample_set
            if "ATAC" in s
        )
    if not gex_samples:
        logger.warning(f"No GEX sample names found for {sampleID}: {sample_names}")
    if not atac_samples:
        logger.warning(f"No ATAC sample names found for {sampleID}: {sample_names}")

    with open(csv_path, 'w') as csv_file:
        csv_file.write("fastqs,sample_id,library_id\n")
        for sample in gex_samples:
            csv_file.write(f"{sample_path},{sample},Gene Expression\n")
        for sample in atac_samples:
            csv_file.write(f"{sample_path},{sample},Chromatin Accessibility\n")
    return csv_path