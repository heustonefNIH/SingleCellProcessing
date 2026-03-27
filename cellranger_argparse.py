import argparse

def get_args():
	"""
    Build and parse command-line arguments for the Cellranger_swarm.py script.
    Returns:
        argparse.Namespace
    """
	parser = argparse.ArgumentParser(
		description="Generate cellranger swarm file from fastq directory"	
	)
	parser.add_argument(
		'-f', '--fastq_dir',
		required=True,
		help='Directory containing fastq files'
	)
	parser.add_argument(
		'-id', '--sample_id_format',
		default='HPAP-\d{3}',
		metavar="REGEX",
		help="Format of the sample IDs (e.g., '^HPAP-?\d{3}')"
	)
	parser.add_argument(
		'--logfile',
		default='cellranger_swarm.log',
		help='Log file to record renaming events and missing sample information'
	)
	parser.add_argument(
		'-d', '--data_type',
		choices=['rna', 'atac', 'multi'],
		default='rna',
		help='Type of sequencing data (rna, atac, or multi)'
	)	
	parser.add_argument(
		'-s', '--swarmfile_name',
		default='cellranger_v10.0.0.swarm',
		help='Name of the swarm file to be created'
	)
	parser.add_argument(
		'--allow_renaming',
		action='store_true',
		help='Run fx_hpap_rename on files that do not match the expected pattern'
	)
	parser.add_argument(
		'--dry_run',
		action='store_true',
		help='Perform a dry run without actually writing the swarm file'
	)
	parser.add_argument(
		'--allow_loose_match',
		action='store_true',
		help='Allow wildcards at the end of --id regex'
	)
	parser.add_argument(
		'--path_restrictions',
		help='Require specific or regex path structure for fastq files (e.g., Single Cell RNA-Seq or Single-cell Multiome (ATAC+RNA))'
	)
	parser.add_argument(
		'--gex_identifier',
		help='Custom identifier for GEX samples in multiome datasets (default: "GEX"). Required if data_type is "multi"'
	)
	parser.add_argument(
		'--atac_identifier',
		help='Custom identifier for ATAC samples in multiome datasets (default: "ATAC"). Required if data_type is "multi"'
	)
	parser.add_argument(
		'--ref_genome',
		default='refdata-gex-GRCh38-2024-A',
		help='Reference genome to use for cellranger in path (e.g., "refdata-gex-GRCh38-2024-A").'
	)
	parser.add_argument(
		'--debug',
		action='store_true',
		help='Enable debug logging'
	)
	return parser.parse_args()