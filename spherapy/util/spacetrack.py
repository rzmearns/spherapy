import pathlib
import logging
import spacetrack as sp
import spherapy
import datetime as dt
import spherapy.util.epoch_u as epoch_u
import sys
import os
from progressbar import progressbar
import json

MAX_RETRIES=3
TIMEOUT=60

logger = logging.getLogger(__name__)

class TLEGetter:
	def __init__(self, sat_id_list:list[int], user:str=None, passwd:str=None):
		# initialise the client
		if user == None:
			self.username = spherapy.spacetrack_credentials['user']
		else:
			self.username = user
		if passwd == None:			
			self.password = spherapy.spacetrack_credentials['passwd']
		else:
			self.password = passwd
		if self.username is None or self.password is None:
			raise InvalidCredentials('No Spacetrack Credentials have been entered')		
		
		self.modified_ids = []
		try:
			self.stc = sp.SpaceTrackClient(self.username, self.password)
			ii = 0
			for sat_id in progressbar(sat_id_list):
				pc = ii/len(sat_id_list)*100
				bar_str = int(pc)*'='
				space_str = (100-int(pc))*'  '
				print(f'Loading {pc:.2f}% ({ii} of {len(sat_id_list)}) |{bar_str}{space_str}|\r')

				print(f"{sat_id=}")

				if not self.checkTLEFileExists(sat_id) or self.getNumPastTLEs(sat_id) == 0:
					res = self.fetchAll(sat_id)
					if res is not None:
						self.modified_ids.append(res)
				else:
					res = self.fetchLatest(sat_id)
					if res is not None:
						self.modified_ids.append(res)
				ii+=1
		except sp.AuthenticationError:
			raise InvalidCredentials('Username and password are incorrect!')

	def getModifiedIDs(self) -> list[int]:
		return self.modified_ids

	def checkTLEFileExists(self, sat_id:int) -> bool:
		return getTLEFilePath(sat_id).exists()

	def fetchAll(self, sat_id:int) -> int|None:
		attempt_number = 0
		while attempt_number < MAX_RETRIES:
			try:
				opts = self._calcRequestAllOptions(sat_id)
				resp_block = self.stc.tle(*opts)
				self._writeResponseToNewFile(sat_id, resp_block)
				break
			except TimeoutError:
				attempt_number += 1

		if attempt_number == MAX_RETRIES:
			logger.error("Could not fetch All TLEs for sat %s: failed %s times.",sat_id, attempt_number)
			return None
		
		return sat_id

	def getNumPastTLEs(self, sat_id:int) -> int:
		with open(getTLEFilePath(sat_id), 'r') as fp:
			lines = fp.readlines()
		return int(len(lines)/3)

	def fetchLatest(self, sat_id:int) -> int|None:
		attempt_number = 0
		while attempt_number < MAX_RETRIES:
			try:
				_, days_since_last_epoch = self._findLocalLastEpoch(sat_id)
				if days_since_last_epoch > 0.0:
					opts = self._calcRequestPartialOptions(sat_id)
					resp_block = self.stc.tle(*opts)
					self._writeResponseToFile(sat_id, resp_block)
				break
			except TimeoutError as e:
				print(e)
				attempt_number += 1
		if attempt_number == MAX_RETRIES:
			print(f"Could not fetch the latest TLEs for sat {sat_id}: failed {attempt_number} times.", file=sys.stderr)
			return None

		return sat_id

	def _findLocalLastEpoch(self, sat_id:int) -> tuple[float, float]:
		"""Find the last TLE epoch in a TLE file"""
		# get penultimate and ultimate epochs
		with getTLEFilePath(sat_id).open('r') as fp:
			lines = fp.readlines()
		while lines[-1] == '':
			lines = lines[:-1]
		last_epoch_line = lines[-2]
		# pe_line = lines[-5]
		# pe_datetime = epoch_u.epoch2datetime(float(pe_line.split()[3]))
		try:
			last_epoch = float(last_epoch_line.split()[3])
		except IndexError:
			print(last_epoch_line)
			raise IndexError

		last_epoch_datetime = epoch_u.epoch2datetime(last_epoch)
		delta = dt.datetime.now(tz=dt.timezone.utc) - last_epoch_datetime

		return last_epoch, delta.days

	def _calcRequestAllOptions(self, sat_id:int) -> dict[str, str|int]:
		"""Format spacetrack library request options. Request all TLEs"""

		request_options = {}
		request_options['norad_cat_id'] = sat_id
		request_options['order_by'] = 'epoch_asc'
		request_options['limit'] = 500000
		request_options['format'] = '3le'

		return request_options

	def _calcRequestPartialOptions(self, sat_id:int) -> dict[str, str|int]:
		"""Format spacetrack library request options. Request all TLEs since last epoch"""
		_, days_since_last_epoch = self._findLocalLastEpoch(sat_id)

		request_options = {}
		request_options['norad_cat_id'] = sat_id
		request_options['order_by'] = 'epoch_asc'
		request_options['epoch'] = f'>now-{days_since_last_epoch+1}'
		request_options['limit'] = 500000
		request_options['format'] = '3le'

		return request_options

	def _writeResponseToFile(self, sat_id:int, resp:str):
		"""Append response to TLE file"""
		last_epoch, _ = self._findLocalLastEpoch(sat_id)
		resp = resp[:-1]
		resp_lines = resp.split('\n')
		next_idx = 0
		for ii, line in enumerate(resp_lines):
			# Try block for debugging, only sometimes failing, trying to investigate.
			try:
				next_idx = ii
				if line[0] == '1' and float(line.split()[3])>last_epoch:
					break
			except IndexError:
				print(f'{last_epoch=}')
				print(f'{ii=}')
				print(f'{line=}')
				print(f'{line.split()}')
				print(f'{resp_lines}')

		if resp_lines[0] != '':
			with open(getTLEFilePath(sat_id), 'a') as fp:
				fp.write('\n')
				fp.write('\n'.join(resp_lines[next_idx-1:]))

	def _writeResponseToNewFile(self, sat_id:int, resp:str):
		"""New TLE, write entire content to new file"""
		with open(getTLEFilePath(sat_id), 'w') as fp:
			if resp[-1] == '\n':
				resp = resp[:-1]
			fp.write(resp)

class InvalidCredentials(Exception):
	def __init__(self, message):
		super().__init__(message)
		return
	
def updateTLEs(sat_id_list:list[int], user:str=None, passwd:str=None) -> list[int]:
	print(f"Using SPACETRACK to update TLEs")
	g = TLEGetter(sat_id_list,user=user,passwd=passwd)
	
	return g.getModifiedIDs()

def getTLEFilePath(sat_id:int) -> pathlib.Path:
	return spherapy.tle_dir.joinpath(f'{sat_id}.tle')

def doCredentialsExist() -> bool:
	user_stored = False
	passwd_stored = False
	if spherapy.spacetrack_credentials['user'] is not None:
		user_stored = True
	if spherapy.spacetrack_credentials['passwd'] is not None:
		passwd_stored = True

	if user_stored and passwd_stored:
		return True
	
	return False