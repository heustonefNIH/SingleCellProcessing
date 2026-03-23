import shlex

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