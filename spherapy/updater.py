import spherapy
import spherapy.util.spacetrack as spacetrack
import spherapy.util.spacetrack as celestrak
import pickle

def updateTLEs(sat_id_list: list[int]) -> list[str]:

	if spacetrack.doCredentialsExist():
		modified_list = spacetrack.updateTLEs(sat_id_list)
	else:
		modified_list = celestrak.updateTLEs(sat_id_list)

	return modified_list

def createCredentials():
	user = input('Please enter your spacetrack username:')
	passwd = input('Please enter your spacetrack password:')
	print(f"These credentials are being saved in {spherapy.spacetrack_cred_path}/.credentials")
	creds = {'user':user,'passwd':passwd}
	with open(f"{spherapy.spacetrack_cred_path.absolute()}/.credentials",'wb') as fp:
		pickle.dump(creds,fp)

def getTLEFilePaths(sat_id_list:list[int]) -> list[str]:
	if spacetrack.doCredentialsExist():
		return [ spacetrack.getTLEFilePath(sat_id) for sat_id in sat_id_list ]
	else:
		return [ celestrak.getTLEFilePath(sat_id) for sat_id in sat_id_list ]