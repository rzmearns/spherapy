import configparser
from pathlib import Path
from spherapy.util import credentials

pkg_dir = Path(__file__).parents[1]
config = configparser.ConfigParser()
config.read(f"{pkg_dir.absolute()}/spherapy.conf")

tle_path = Path(f"{pkg_dir}/{config['paths'].get('TLE_Path')}")

spacetrack_credentials = credentials.fetchCredentials()