import configparser
from importlib import metadata, resources
import os
import pathlib

import platformdirs

from spherapy.util import credentials


def _creatPackagedTLEListing() -> None|dict[int,pathlib.Path]:
	packaged_tles = {}
	package_file_listing = metadata.files('spherapy')
	package_dir = resources.files('spherapy')
	if package_file_listing is not None and len(package_file_listing) > 0:
		for path in package_file_listing:
			if 'TLEs' in path.parts:
				try:
					tle_id = int(path.stem)
				except ValueError:
					print("Can't import packaged TLEs, skipping...") #noqa: T201
					return None
				traversable = package_dir.joinpath(f"{path.relative_to('spherapy')}")
				resource_path = pathlib.Path(f'{traversable}')
				packaged_tles[tle_id] = resource_path
	return packaged_tles

service_name = "spherapy"
service_author = "MSL"

try:
	__version__ = metadata.version(service_name)
except metadata.PackageNotFoundError:
	print(f'Version of {service_name} is uknown') #noqa: T201
	__version__ = "x.x.x"

config_dir_str = os.getenv('SPHERAPY_CONFIG_DIR')
if config_dir_str is None:
	use_config_file = False
	config_dir = pathlib.Path(platformdirs.user_config_dir(appname=service_name,
															version=__version__,
															ensure_exists=True))
else:
	use_config_file = True
	config_dir = pathlib.Path(config_dir_str)

# will load a config even if not using
config = configparser.ConfigParser()
read_files = config.read(config_dir.joinpath('spherapy.conf'))

if use_config_file and len(read_files) == 0:
	raise configparser.Error(f"Couldn't find spherapy.conf at {config_dir}")

# Locate dir to store TLE files
if not use_config_file:
	# use default user directories
	_data_dir = pathlib.Path(platformdirs.user_data_dir(appname=service_name, ensure_exists=True))
	tle_dir = _data_dir.joinpath('TLEs')
else:
	# use the config file
	_TLE_dir_str = config['paths'].get('TLE_path')

	if _TLE_dir_str is None:
		raise ValueError("Config file has not TLE_path")

	_temp_dir = pathlib.Path(_TLE_dir_str)
	if _TLE_dir_str == '':
		# TLE dir not set in config file
		_data_dir = pathlib.Path(platformdirs.user_data_dir(appname=service_name,
															ensure_exists=True))
		tle_dir = _data_dir.joinpath('TLEs')
	elif _temp_dir.is_absolute():
		tle_dir = _temp_dir
	else: 						# noqa: PLR5501
		# _temp_dir is relative
		if config_dir.is_dir():
			tle_dir = config_dir.joinpath(_TLE_dir_str)
		else:
			tle_dir = config_dir.parent.joinpath(_TLE_dir_str)

tle_dir.mkdir(parents=True, exist_ok=True)

# Load credentials

spacetrack_credentials:dict[str,str|None] = {'user':None, 'passwd':None}
# load spacetrack_credentials from the relevant source
if not use_config_file:
	credentials._reloadCredentials() 			# noqa: SLF001
else:
	credentials._reloadCredentials(config) 		# noqa: SLF001


packaged_TLEs = _creatPackagedTLEListing() 	#noqa: N816
