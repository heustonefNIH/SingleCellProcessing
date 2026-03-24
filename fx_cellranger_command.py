import shlex
import os
import logging

logger = logging.getLogger(__name__)

def call_cellranger_command(output_ID, sample_path, ref_genome_cmd, sample_args, cellranger_module):
	sample_path = shlex.quote(sample_path)
	return f"""# SAMPLE {output_ID} 
FASTQ_PATH={sample_path}; \\
ulimit -u 10240 -n 16384; \\
{cellranger_module} count --id={output_ID} \\
{ref_genome_cmd} \\
--fastqs="$FASTQ_PATH" \\
--sample={sample_args} \\
--create-bam=false \\
--localcores=$SLURM_CPUS_PER_TASK \\
--localmem=62 \\
--maxjobs=10

"""

def arc_library_csv(sampleID, GEX_id, ATAC_id, sample_path):
	csv_content = f"""fastqs,sample,library_type,
{sample_path},{GEX_id},Gene Expression,
{sample_path},{ATAC_id},ATAC,
"""
	csv_path = os.path.join(sample_path, f"{sampleID}_library.csv")
	with open(csv_path, 'w') as csv_file:
		csv_file.write(csv_content)
	return csv_path