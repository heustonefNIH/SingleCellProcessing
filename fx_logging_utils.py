import logging
import sys

def setup_logging(logfile, debug = False):
	logger = logging.getLogger()
	logger.setLevel(logging.DEBUG)

	#Clean heandlers
	for handler in logger.handlers:
		logger.removeHandler(handler)
	
	# Create file handler
	file_handler = logging.FileHandler(logfile, mode='w')
	file_handler.setLevel(logging.DEBUG)
	file_format = logging.Formatter(
		'%(asctime)s - %(levelname)s - %(message)s',
		datefmt='%Y-%m-%d %H:%M:%S'
		)
	file_handler.setFormatter(file_format)
	logger.addHandler(file_handler)

	# Create console handler
	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.setLevel(logging.DEBUG if debug else logging.INFO)
	console_format = logging.Formatter('%(levelname)s: %(message)s')
	console_handler.setFormatter(console_format)
	logger.addHandler(console_handler)

	logger.debug("Logging setup complete. Log file: %s", logfile)