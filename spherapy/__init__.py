import pickle
import configparser
from pathlib import Path

pkg_dir = Path(__file__).parents[1]
config = configparser.ConfigParser()
config.read(f"{pkg_dir.absolute()}/spherapy.conf")

tle_path = Path(f"{pkg_dir}/{config['paths'].get('TLE_Path')}"
spacetrack_cred_path = Path(f"{pkg_dir}/{config['paths'].get('SpaceTrackCredential_Path')}"


try:
	with open(f"{spacetrack_cred_path.absolute()}/.credentials",'rb') as fp:
		spacetrack_credentials = pickle.load(fp)
except Exception as e:
	print(e)
	spacetrack_credentials = {'user':None, 'passwd':None}