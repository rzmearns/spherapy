"""Functions to update TLE of a satcat ID, and fetch where the file is stored."""

import pathlib

from spherapy.util import spacetrack
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

def getTLEFilePaths(sat_id_list:list[int]) -> list[pathlib.Path]:
	"""Fetch list of paths to TLE files.

	Fetch paths to the file storing the list of TLEs for each provided satcat ID
	If credentials are stored, return path containing all historical TLEs (fetched from Spacetrack)
	If no credentials are stored, return path to .temptle containing only most recent TLE
		(from Celestrak)

	Args:
		sat_id_list: list of satcat ids to find paths

	Returns:
		list of paths
	"""
	if spacetrack.doCredentialsExist():
		return [ spacetrack.getTLEFilePath(sat_id) for sat_id in sat_id_list ]
	return [ celestrak.getTLEFilePath(sat_id) for sat_id in sat_id_list ]
