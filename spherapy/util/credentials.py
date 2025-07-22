import spherapy
import configparser
import keyring
import keyring.errors
import getpass

# used to store/fetch the username in keyring (abuse of API)
USERNAME_KEY = "spherapy_username"

def fetchConfigCredentials(config:configparser.ConfigParser) -> dict[str, str|None]:
	username = config['credentials'].get('SpacetrackUser')
	if username == 'None':
		username = None

	password = config['credentials'].get('SpacetrackPasswd')
	if password == 'None':
		password = None

	return {'user':username, 'passwd':password}

def fetchKeyringCredentials() -> dict[str,str|None]:
	username = _fetchUser()
	if username is None:
		return {'user':None, 'passwd':None}
	password = _fetchPass(username)
	return {'user':username, 'passwd':password}

def _fetchUser() -> str|None:
	try:
		username = keyring.get_password(spherapy.service_name, USERNAME_KEY)
	except keyring.errors.NoKeyringError:
		raise ValueError('No Keyring exists on this machine: did you forget to set env variable SPHERAPY_CONFIG_DIR')
	return username

def _fetchPass(username:str) -> str|None:
	try:
		password = keyring.get_password(spherapy.service_name, username)
	except keyring.errors.NoKeyringError:
		raise ValueError('No Keyring exists on this machine: did you forget to set env variable SPHERAPY_CONFIG_DIR')
	return password

def storeCredentials(user:None|str=None, passwd:None|str=None) -> bool:
	if user is not None:
		keyring.set_password(spherapy.service_name, USERNAME_KEY, user)
	else:
		user = _fetchUser()
		if user is None:
			raise KeyError("Can't set a Spacetrack password without a username: no existing username")
	if passwd is not None:
		keyring.set_password(spherapy.service_name, user, passwd)

	# check stored credentials match input
	stored_user = _fetchUser()
	if stored_user is not None:
		stored_pass = _fetchPass(stored_user)
	else:
		stored_pass = None
	if stored_user == user and stored_pass == passwd:
		return True

	return False

def createCredentials():
	user = input('Please enter your spacetrack username:')
	passwd = getpass.getpass(prompt='Please enter your spacetrack password:')
	print("These credentials are being saved in your system keyring")
	if not storeCredentials(user=user, passwd=passwd):
		print("Could not save credentials in system keyring")
