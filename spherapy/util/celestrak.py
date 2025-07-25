"""Functions to fetch TLEs from Celestrak.

Attributes:
	MAX_RETRIES: number of times to try and reach Celestrak.
	TIMEOUT: timeout for connection and request servicing by Celestrak.
"""

import datetime as dt
import logging
import pathlib

import requests

import spherapy
from spherapy.util import epoch_u

MAX_RETRIES=3
TIMEOUT=10

logger = logging.getLogger(__name__)

def updateTLEs(sat_id_list:list[int]) -> list[int]:
	"""Fetch most recent TLE for satcat IDs from celestrak.

	Fetch most recent TLE for provided list of satcat IDs, and store in file.
	Will try MAX_RETRIES before raising a ValueError

	Args:
		sat_id_list: list of satcat ids to fetch

	Returns:
		list: list of satcat ids successfully fetched

	Raises:
		TimeoutError
	"""
	logger.info("Using CELESTRAK to update TLEs")
	modified_list = []
	for sat_id in sat_id_list:
		url = f'https://celestrak.org/NORAD/elements/gp.php?CATNR={sat_id}'
		retry_num = 0
		fetch_successful = False
		while retry_num < MAX_RETRIES:
			retry_num += 1
			r = requests.get(url, timeout=TIMEOUT)
			if r.status_code == requests.codes.success:
				fetch_successful = True
				tle_file = getTLEFilePath(sat_id)
				dat_list = r.text.split('\r\n')
				with tle_file.open('w') as fp:
					fp.write(f'0 {dat_list[0].rstrip()}\n')
					fp.write(f'{dat_list[1].rstrip()}\n')
					fp.write(f'{dat_list[2].rstrip()}')
					modified_list.append(sat_id)
				break

		if retry_num == MAX_RETRIES and fetch_successful:
			logger.error('Could not fetch celestrak information for sat_id: %s', sat_id)
			raise TimeoutError(f'Could not fetch celestrak information for sat_id: {sat_id}')

	return modified_list

def getTLEFilePath(sat_id:int) -> pathlib.Path:
	"""Gives path to file where celestrak TLE is stored.

	Celestrak TLEs are stored in {satcadID}.temptle

	Args:
		sat_id: satcat ID

	Returns:
		path to file
	"""
	return spherapy.tle_dir.joinpath(f'{sat_id}.temptle')

def getStoredEpochs(sat_id:int) -> None|tuple[dt.datetime, dt.datetime|None]:
	"""Return the start and end epoch for {sat_id}.temptle .

	Args:
		sat_id: satcat id to check

	Returns:
		(first epoch datetime, last epoch datetime)
		None if no spacetrack tle stored for sat_id
	"""
	tle_path = getTLEFilePath(sat_id)
	return epoch_u.getStoredEpochs(tle_path)
