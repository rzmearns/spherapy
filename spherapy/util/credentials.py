import keyring
import keyring.errors
import getpass

# used to store/fetch the username in keyring (abuse of API)
USERNAME_KEY = "spherapy_username"
# program name
service_id = "spherapy"

try:
	keyring.get_keyring()
	method = 'keyring'
except keyring.errors.NoKeyringError:
	method = 'fallback'

def fetchCredentials() -> dict[str,str|None]:
	username = _fetchUser()
	if username is None:
		return {'user':None, 'passwd':None}
	password = _fetchPass(username)
	return {'user':username, 'passwd':password}

def _fetchUser() -> str|None:
	if method=='keyring':
		username = keyring.get_password(service_id, USERNAME_KEY)
	else:
		username = None
	return username

def _fetchPass(username:str) -> str|None:
	if method=='keyring':
		password = keyring.get_password(service_id, username)
	else:
		password = None
	return password

def storeCredentials(user:None|str=None, passwd:None|str=None) -> bool:
	if user is not None:
		keyring.set_password(service_id, USERNAME_KEY, user)
	else:
		user = _fetchUser()
		if user is None:
			raise KeyError("Can't set a Spacetrack password without a username: no existing username")
	if passwd is not None:
		keyring.set_password(service_id, user, passwd)

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
	print(f"These credentials are being saved in your system keyring")
	if not storeCredentials(user=user, passwd=passwd):
		print(f"Could not save credentials in system keyring")
