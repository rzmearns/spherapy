"""Functions to update TLE of a satcat ID, and fetch where the file is stored."""
import datetime as dt
import pathlib

import spherapy
from spherapy.util import epoch_u, spacetrack
import spherapy.util.spacetrack as celestrak


def updateTLEs(sat_id_list: list[int]) -> list[int]:
	"""Fetch most recent TLE for provided list of satcat IDs.

	Fetch most recent TLE for provided list of satcat IDs,
	If credentials are stored, use Spacetrack as TLE source, this will fetch all historical
		TLEs for that satcat ID (from most recent stored TLE up to current)
	If no credentials are stored, use Celestrak as TLE source, this will fetch only the most
		recent TLE

	Args:
		sat_id_list: list of satcat ids to fetch

	Returns:
		list of satcat ids successfully fetched
	"""
	if spacetrack.doCredentialsExist():
		modified_list = spacetrack.updateTLEs(sat_id_list)
	else:
		modified_list = celestrak.updateTLEs(sat_id_list)

	return modified_list

def getTLEFilePaths(sat_id_list:list[int], use_packaged:bool=False) -> list[pathlib.Path]:
	"""Fetch list of paths to TLE files.

	Fetch paths to the file storing the list of TLEs for each provided satcat ID
	If credentials are stored, return path containing all historical TLEs (fetched from Spacetrack)
	If no credentials are stored, return path to .temptle containing only most recent TLE
		(from Celestrak)
	You can force the use of TLE data packaged with spherapy by setting use_packaged=True
		This data will be out of date, and should only be used in examples.

	Args:
		sat_id_list: list of satcat ids to find paths
		use_packaged: use the default TLEs packaged with spherapy, these will be out of date.

	Returns:
		list of paths
	"""
	if use_packaged and spherapy.packaged_TLEs is None:
		raise ValueError('There are no TLEs packaged with spherapy')
	if use_packaged and spherapy.packaged_TLEs is not None:
		try:
			path_list = []
			attempted_sat_id = 0
			for sat_id in sat_id_list:
				attempted_sat_id = sat_id
				path_list.append(spherapy.packaged_TLEs[sat_id])
		except KeyError as e:
			raise KeyError(f'TLE for {attempted_sat_id} is not packaged with spherapy.') from e
		else:
			return path_list
	elif spacetrack.doCredentialsExist():
		return [ spacetrack.getTLEFilePath(sat_id) for sat_id in sat_id_list ]
	return [ celestrak.getTLEFilePath(sat_id) for sat_id in sat_id_list ]

def getStoredEpochLimits(sat_id_list:list[int], use_packaged:bool=False) \
							-> dict[int, None|tuple[dt.datetime, dt.datetime|None]]:
	"""Returns limiting epochs for the stored TLEs for each sat in sat_id_list.

	[description]

	Args:
		sat_id_list (list[int]): [description]
		use_packaged (bool): [description] (default: `False`)

	Returns:
		list[tuple[dt.datetime, dt.datetime]]: [description]
	"""
	terminating_epochs = {}
	if use_packaged and spherapy.packaged_TLEs is not None:
		for sat_id in sat_id_list:
			packaged_path = spherapy.packaged_TLEs[sat_id]
			terminating_epochs[sat_id] = epoch_u.getStoredEpochs(packaged_path)
	elif spacetrack.doCredentialsExist():
		for sat_id in sat_id_list:
			terminating_epochs[sat_id] = spacetrack.getStoredEpochs(sat_id)
	else:
		for sat_id in sat_id_list:
			terminating_epochs[sat_id] = celestrak.getStoredEpochs(sat_id)

	return terminating_epochs
