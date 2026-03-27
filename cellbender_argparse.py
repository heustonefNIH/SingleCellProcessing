import argparse

def get_args():
	"""
    Build and parse command-line arguments for the cellbender.py script.
    Returns:
        argparse.Namespace
    """
	parser = argparse.ArgumentParser(
		description="Generate cellbender swarm file from cellranger output directory"	
	)
	parser.add_argument(
		'-f', '--cellranger_folder',
		required=True,
		help='Directory containing cellranger output data'
	)
	parser.add_argument(
		'-id', '--sample_id_format',
		default=None,
		metavar="REGEX",
		help="Format of the sample IDs (e.g., '^HPAP-?\d{3}'). Can be single regex or list"
	)
	parser.add_argument(
		'-o', '--cellbender_file',
		default='cellbender.swarm',
		help='Name of the swarm file to be created. If None, commands will be printed to console instead of written to file.'
	)
	parser.add_argument(
		'--cuda',
		action='store_false',
		help='Include --cuda flag in cellbender command for GPU acceleration (Default: True)'
	)
	parser.add_argument(
		'--gpu_partition',
		default='gpu:v100x:1',
		help='Slurm gpu partition passed to --gres flag (default: "gpu:v100x:1")'
	)
	parser.add_argument(
		'--flags',
		default='--cpu-threads $SLURM_CPUS_PER_TASK',
		help='Additional CellBender flags to include in command (default: "--cpu-threads $SLURM_CPUS_PER_TASK")'
	)
	parser.add_argument(
		'--logfile',
		default='cellbender_swarm.log',
		help='Log file to record events and errors during swarm file generation'
	)
	parser.add_argument(
		'--debug',
		action='store_true',
		help='Enable debug logging'
	)
	return parser.parse_args()
