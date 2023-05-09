from configparser import ConfigParser


#
from reader import Reader
#
# Load configuration
configuration = ConfigParser()
configuration.read('configuration.ini')
#
# Create a reader
aReader = Reader(configuration)
aReader.run()

print('> done')