import shlex

def call_cellranger_command(sampleID, sample_path, ref_genome_cmd, cellranger_module):
	sample_path = shlex.quote(sample_path)
	return f"""# SAMPLE {sampleID} 
FASTQ_PATH={sample_path}; \\
ulimit -u 10240 -n 16384; \\
{cellranger_module} count --id={sampleID} \\
{ref_genome_cmd} \\
--fastqs="$FASTQ_PATH" \\
--sample={sampleID} \\
--create-bam=false \\
--localcores=$SLURM_CPUS_PER_TASK \\
--localmem=62 \\
--maxjobs=10

"""

