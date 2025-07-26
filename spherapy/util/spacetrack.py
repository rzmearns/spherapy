"""Functions to fetch TLEs from Spacetrack.

Attributes:
	MAX_RETRIES: number of times to try and reach Spacetrack.
"""

import datetime as dt
import logging
import pathlib

from typing import TypedDict
from typing_extensions import NotRequired

import spacetrack as sp

import spherapy
from spherapy.util import elements_u, epoch_u

MAX_RETRIES=3

logger = logging.getLogger(__name__)

class _TLEGetter:
	"""Container around python spacetrack library."""
	def __init__(self, sat_id_list:list[int], user:None|str=None, passwd:None|str=None):
		# initialise the client
		if user is None:
			self.username = spherapy.spacetrack_credentials['user']
		else:
			self.username = user
		if passwd is None:
			self.password = spherapy.spacetrack_credentials['passwd']
		else:
			self.password = passwd
		if self.username is None or self.password is None:
			raise InvalidCredentialsError('No Spacetrack Credentials have been entered')

		self.modified_ids = []
		try:
			self.stc = sp.SpaceTrackClient(self.username, self.password)
			for sat_id in sat_id_list:
				logger.info("Trying to update TLEs for %s", sat_id)
				if not self.checkTLEFileExists(sat_id) or self.getNumPastTLEs(sat_id) == 0:
					res = self.fetchAll(sat_id)
					if res is not None:
						self.modified_ids.append(res)
				else:
					res = self.fetchLatest(sat_id)
					if res is not None:
						self.modified_ids.append(res)

		except sp.AuthenticationError:
			raise InvalidCredentialsError('Username and password are incorrect!') from sp.AuthenticationError 	#noqa: E501

	def getModifiedIDs(self) -> list[int]:
		return self.modified_ids

	def checkTLEFileExists(self, sat_id:int) -> bool:
		return getTLEFilePath(sat_id).exists()

	def fetchAll(self, sat_id:int) -> int|None:
		"""Fetch all TLEs for sat_id.

		Args:
			sat_id: satcat id to fetch

		Returns:
			sat_id if succesfully fetched
			None if couldn't fetch
		"""
		attempt_number = 0
		while attempt_number < MAX_RETRIES:
			try:
				opts = self._calcRequestAllOptions(sat_id)
				logger.info("Requesting TLEs from spacetrack with following options: %s", opts)
				resp_line_iterator = self.stc.tle(**opts)
				resp_lines = list(resp_line_iterator)
				tle_dict_list = elements_u.dictify3LEs(resp_lines)
				self._writeTLEsToNewFile(sat_id, tle_dict_list)
				break
			except TimeoutError:
				attempt_number += 1

		if attempt_number == MAX_RETRIES:
			logger.error("Could not fetch All TLEs for sat %s: failed %s times.",
							sat_id, attempt_number)
			return None

		return sat_id

	def getNumPastTLEs(self, sat_id:int) -> int:
		"""Retrieve number of TLEs stored for sat_id."""
		with getTLEFilePath(sat_id).open('r') as fp:
			lines = fp.readlines()
		return int(len(lines)/3)

	def fetchLatest(self, sat_id:int) -> int|None:
		"""Fetch the latest TLEs for sat_id, and append them to the TLE file.

		Args:
			sat_id: satcat id to fetch

		Returns:
			sat_id if succesfully fetched
			None if couldn't fetch
		"""
		attempt_number = 0
		while attempt_number < MAX_RETRIES:
			try:
				_, days_since_last_epoch = self._findLocalLastEpoch(sat_id)
				if days_since_last_epoch > 0.0:
					opts = self._calcRequestPartialOptions(sat_id)
					logger.info("Requesting TLEs from spacetrack with following options: %s", opts)
					resp_line_iterator = self.stc.tle(**opts)
					resp_lines = list(resp_line_iterator)
					tle_dict_list = elements_u.dictify3LEs(resp_lines)
					self._writeTLEsToFile(sat_id, tle_dict_list)
				break
			except TimeoutError:
				attempt_number += 1
		if attempt_number == MAX_RETRIES:
			logger.error("Could not fetch the latest TLEs for sat %s: failed %s times.",
							sat_id, attempt_number)
			raise TimeoutError(f"Could not fetch the latest TLEs for sat {sat_id}: "
								f"failed {attempt_number} times.")

		return sat_id

	def _findLocalLastEpoch(self, sat_id:int) -> tuple[str, float]:
		"""Find the last TLE epoch in a TLE file."""
		# get penultimate and ultimate epochs
		with getTLEFilePath(sat_id).open('r') as fp:
			lines = fp.readlines()
		while lines[-1] == '':
			lines = lines[:-1]
		last_epoch_line = lines[-2]
		last_tle_epoch = last_epoch_line.split()[3]

		last_epoch_datetime = epoch_u.epoch2datetime(last_tle_epoch)
		delta = dt.datetime.now(tz=dt.timezone.utc) - last_epoch_datetime

		return last_tle_epoch, delta.days

	def _calcRequestAllOptions(self, sat_id:int) -> "_RequestOptions":
		"""Format spacetrack library request options. Request all TLEs."""
		request_options:_RequestOptions = {
				'norad_cat_id': sat_id,
				'orderby': 'epoch asc',
				'limit': 500000,
				'format': '3le',
				'iter_lines': True}

		return request_options

	def _calcRequestPartialOptions(self, sat_id:int) -> "_RequestOptions":
		"""Format spacetrack library request options. Request all TLEs since last epoch."""
		_, days_since_last_epoch = self._findLocalLastEpoch(sat_id)

		request_options:_RequestOptions = {
				'norad_cat_id': sat_id,
				'orderby': 'epoch asc',
				'epoch': f'>now-{days_since_last_epoch+1}',
				'limit': 500000,
				'format': '3le',
				'iter_lines': True}

		return request_options

	def _writeTLEsToFile(self, sat_id:int,
								tle_dict_list:list[dict[int, elements_u.ElementsLineDict]]):
		"""Append TLEs to TLE file."""
		last_epoch, _ = self._findLocalLastEpoch(sat_id)
		first_new_idx = None
		for tle_idx, tle_dict in enumerate(tle_dict_list):
			epoch = tle_dict[1]['fields'][3]
			if epoch_u.epochLaterThan(epoch, last_epoch):
				first_new_idx = tle_idx
				break
		if first_new_idx is not None:
			with getTLEFilePath(sat_id).open('a') as fp:
				for tle_dict in tle_dict_list[first_new_idx:]:
					fp.write('\n')
					fp.write(elements_u.stringify3LEDict(tle_dict))

	def _writeTLEsToNewFile(self, sat_id:int,
									tle_dict_list:list[dict[int, elements_u.ElementsLineDict]]):
		"""New TLE, write entire content to new file."""
		with getTLEFilePath(sat_id).open('w') as fp:
			# don't want to begin or end file with newline
			fp.write(elements_u.stringify3LEDict(tle_dict_list[0]))
			for tle_dict in tle_dict_list[1:]:
				fp.write('\n')
				fp.write(elements_u.stringify3LEDict(tle_dict))

class InvalidCredentialsError(Exception):
	"""Error indicating invalid credentials used to access spacetrack."""
	def __init__(self, message:str): # noqa: D107
		super().__init__(message)

class _RequestOptions(TypedDict):
	norad_cat_id: int
	orderby: str
	epoch: NotRequired[str]
	limit: int
	format: str
	iter_lines: bool

def updateTLEs(sat_id_list:list[int], user:None|str=None, passwd:None|str=None) -> list[int]:
	"""Fetch most recent TLE for satcat IDs from spacetrack.

	Fetch most recent TLEs for provided list of satcat IDs, and append to file.
		If no TLE file yet exists, will download all historical TLEs
	Will try MAX_RETRIES before raising a TimeoutError

	Args:
		sat_id_list: list of satcat ids to fetch
		user: [Optional] overriding spacetrack username to use
		passwd: [Optional] overriding spacetrack password to use

	Returns:
		list:list of satcat ids successfully fetched

	Raises:
		TimeoutError
	"""
	logger.info("Using SPACETRACK to update TLEs")
	g = _TLEGetter(sat_id_list,user=user,passwd=passwd)

	return g.getModifiedIDs()

def getTLEFilePath(sat_id:int) -> pathlib.Path:
	"""Gives path to file where spacetrack TLE is stored.

	Spacetrack TLEs are stored in {satcadID}.tle

	Args:
		sat_id: satcat ID

	Returns:
		path to file
	"""
	return spherapy.tle_dir.joinpath(f'{sat_id}.tle')

def getStoredEpochs(sat_id:int) -> None|tuple[dt.datetime, dt.datetime|None]:
	"""Return the start and end epoch for {sat_id}.tle .

	Args:
		sat_id: satcat id to check

	Returns:
		(first epoch datetime, last epoch datetime)
		None if no spacetrack tle stored for sat_id
	"""
	tle_path = getTLEFilePath(sat_id)

	return epoch_u.getStoredEpochs(tle_path)


def doCredentialsExist() -> bool:
	"""Checks if spacetrack credentials have been loaded into Spherapy.

	Returns:
		bool:
	"""
	user_stored = False
	passwd_stored = False
	if spherapy.spacetrack_credentials['user'] is not None:
		user_stored = True
	if spherapy.spacetrack_credentials['passwd'] is not None:
		passwd_stored = True

	return (user_stored and passwd_stored)
