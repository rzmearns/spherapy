"""Functions to fetch and store spacetrack credentials on the system."""

import configparser
import getpass

import keyring
import keyring.errors

import spherapy

# used to store/fetch the username in keyring (abuse of API)
USERNAME_KEY = "spherapy_username"

def fetchConfigCredentials(config:configparser.ConfigParser) -> dict[str, str|None]:
	"""Fetch spacetrack credentials from config file.

	Args:
		config: ConfigParser object containing config fields

	Returns:
		dict containing username and password
	"""
	username = config['credentials'].get('SpacetrackUser')
	if username == 'None':
		username = None

	password = config['credentials'].get('SpacetrackPasswd')
	if password == 'None': 		# noqa: S105, password not hardcoded, just parsing empty config file
		password = None

	return {'user':username, 'passwd':password}

def fetchKeyringCredentials() -> dict[str,str|None]:
	"""Fetch spacetrack credentials from system keyring.

	Returns:
		dict containing username and password

	Raises:
		ValueError: If no system keyring exists
	"""
	username = _fetchUser()
	if username is None:
		return {'user':None, 'passwd':None}
	password = _fetchPass(username)
	return {'user':username, 'passwd':password}

def _fetchUser() -> str|None:
	try:
		username = keyring.get_password(spherapy.service_name, USERNAME_KEY)
	except keyring.errors.NoKeyringError:
		raise ValueError('No Keyring exists on this machine: '
							'did you forget to set env variable SPHERAPY_CONFIG_DIR '
							'if using on a headless machine?') \
							from keyring.errors.NoKeyringError
	return username

def _fetchPass(username:str) -> str|None:
	try:
		password = keyring.get_password(spherapy.service_name, username)
	except keyring.errors.NoKeyringError:
		raise ValueError('No Keyring exists on this machine: '
							'did you forget to set env variable SPHERAPY_CONFIG_DIR') \
							from keyring.errors.NoKeyringError
	return password

def _reloadCredentials(config_parser:configparser.ConfigParser|None=None):
	if config_parser is None:
		spherapy.spacetrack_credentials = fetchKeyringCredentials()
	else:
		spherapy.spacetrack_credentials = fetchConfigCredentials(config_parser)

def storeCredentials(user:None|str=None, passwd:None|str=None) -> bool:
	"""Store spacetrack credentials in system keyring.

	Can set either username or password, but can't set a password without a username

	Args:
		user: [optional] username string
		passwd: [optional] password string

	Returns:
		True if successfully stored

	Raises:
		KeyError: raised if trying to store a password without a username
	"""
	if user is not None:
		keyring.set_password(spherapy.service_name, USERNAME_KEY, user)
	else:
		user = _fetchUser()
		if user is None:
			raise KeyError("Can't set a Spacetrack password without a username: "
							"no existing username")
	if passwd is not None:
		keyring.set_password(spherapy.service_name, user, passwd)

	# check stored credentials match input
	stored_user = _fetchUser()
	stored_pass = _fetchPass(stored_user) if stored_user is not None else None
	if (stored_user == user and stored_pass == passwd):
		_reloadCredentials()
		return True

	return False

def clearCredentials() -> None:
	"""Clear the stored spacetrack credentials.

	Raises:
		KeyError: If no credentials are currently stored.
	"""
	user = _fetchUser()
	if user is not None:
		keyring.delete_password(spherapy.service_name, user)
	else:
		raise KeyError("No stored USER - no credentials to delete.")
	keyring.delete_password(spherapy.service_name, USERNAME_KEY)

	_reloadCredentials()

def createCredentials():
	"""Script helper function to create and store credentials.

	Called by command line script. Requires user input.
	"""
	user = input('Please enter your spacetrack username:')
	passwd = getpass.getpass(prompt='Please enter your spacetrack password:')
	print("These credentials are being saved in your system keyring") 	# noqa: T201 print needed
	if not storeCredentials(user=user, passwd=passwd):
		print("Could not save credentials in system keyring") 	# noqa: T201 print needed
