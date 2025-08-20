"""Command line scripts for spherapy.

Contains:
 - create_credentials (spherapy-create-credentials)

"""

from spherapy import credentials


def create_credentials():
	"""Helper script to call create credentials function."""
	credentials.createCredentials()

def clear_credentials():
	"""Helper script to call clear credentials function."""
	credentials.clearCredentials()
