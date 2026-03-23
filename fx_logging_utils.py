import logging
import sys

class ConsoleOnlyFilter(logging.Filter):
	def filter(self, record):
		return getattr(record, 'console_only', False)

def setup_logging(logfile, debug = False):
	root = logging.getLogger()
	root.setLevel(logging.DEBUG)

	#Clean handlers
	for handler in root.handlers:
		root.removeHandler(handler)
	
	# Create file handler
	file_handler = logging.FileHandler(logfile, mode='w')
	file_handler.setLevel(logging.DEBUG)
	file_format = logging.Formatter(
		'%(asctime)s - %(levelname)s - %(message)s',
		datefmt='%Y-%m-%d %H:%M:%S'
		)
	file_handler.setFormatter(file_format)
	root.addHandler(file_handler)

	# Create console handler
	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.setLevel(logging.INFO)
	console_handler.addFilter(ConsoleOnlyFilter())
	console_format = logging.Formatter('%(levelname)s: %(message)s')
	console_handler.setFormatter(console_format)
	root.addHandler(console_handler)

	root.debug("Logging setup complete. Log file: %s", logfile)