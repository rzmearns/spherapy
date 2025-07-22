import spherapy
import spherapy.util.spacetrack as spacetrack
import spherapy.util.spacetrack as celestrak
import pickle
import pathlib

def updateTLEs(sat_id_list: list[int]) -> list[str]:

	if spacetrack.doCredentialsExist():
		modified_list = spacetrack.updateTLEs(sat_id_list)
	else:
		modified_list = celestrak.updateTLEs(sat_id_list)

	return modified_list

def getTLEFilePaths(sat_id_list:list[int]) -> list[pathlib.Path]:
	if spacetrack.doCredentialsExist():
		return [ spacetrack.getTLEFilePath(sat_id) for sat_id in sat_id_list ]
	else:
		return [ celestrak.getTLEFilePath(sat_id) for sat_id in sat_id_list ]