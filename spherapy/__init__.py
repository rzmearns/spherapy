import pickle
import configparser
import spherapy.util.paths_u as paths_u

config = configparser.ConfigParser()
config.read('spherapy.conf')

tle_path = paths_u.sanitisePath(config['paths'].get('TLE_Path'))
spacetrack_cred_path = paths_u.sanitisePath(config['paths'].get('SpaceTrackCredential_Path'))


try:
	with open(f"{spacetrack_cred_path.absolute()}/.credentials",'rb') as fp:
		spacetrack_credentials = pickle.load(fp)
except Exception as e:
	print(e)
	spacetrack_credentials = {'user':None, 'passwd':None}