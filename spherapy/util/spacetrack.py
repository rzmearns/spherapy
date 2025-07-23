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
				logger.info("Requesting TLEs from spacetrack with following options: %s", opts)
				resp_line_iterator = self.stc.tle(**opts)
				resp_lines = list(resp_line_iterator)
				tle_dict_list = self._dictify3LEs(resp_lines)
				self._writeTLEsToNewFile(sat_id, tle_dict_list)
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
				print(f'{days_since_last_epoch=}')
				if days_since_last_epoch > 0.0:
					opts = self._calcRequestPartialOptions(sat_id)
					logger.info("Requesting TLEs from spacetrack with following options: %s", opts)
					resp_line_iterator = self.stc.tle(**opts)
					resp_lines = list(resp_line_iterator)
					tle_dict_list = self._dictify3LEs(resp_lines)
					self._writeTLEsToFile(sat_id, tle_dict_list)
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
		request_options['orderby'] = 'epoch asc'
		request_options['limit'] = 500000
		request_options['format'] = '3le'
		request_options['iter_lines'] = True

		return request_options

	def _calcRequestPartialOptions(self, sat_id:int) -> dict[str, str|int]:
		"""Format spacetrack library request options. Request all TLEs since last epoch"""
		_, days_since_last_epoch = self._findLocalLastEpoch(sat_id)

		request_options = {}
		request_options['norad_cat_id'] = sat_id
		request_options['orderby'] = 'epoch asc'
		request_options['epoch'] = f'>now-{days_since_last_epoch+1}'
		request_options['limit'] = 500000
		request_options['format'] = '3le'
		request_options['iter_lines'] = True

		return request_options

	def _writeTLEsToFile(self, sat_id:int, tle_dict_list:list[dict[int,dict[str,list[str]|str]]]):
		"""Append response to TLE file"""
		last_epoch, _ = self._findLocalLastEpoch(sat_id)
		first_new_idx = None
		for tle_idx, tle_dict in enumerate(tle_dict_list):
			try:
				epoch = float(tle_dict[1]['fields'][3])
			except ValueError as e:
				print(f'ValueError:{e}:{tle_idx}:{tle_dict}')
			if epoch > last_epoch:
				first_new_idx = tle_idx
				break
		if first_new_idx is not None:
			with getTLEFilePath(sat_id).open('a') as fp:
				for tle_dict in tle_dict_list[first_new_idx:]:
					fp.write('\n')
					fp.write(self._stringify3LEDict(tle_dict))

	def _writeTLEsToNewFile(self, sat_id:int, tle_dict_list:list[dict[int,dict[str,list[str]|str]]]):
		"""New TLE, write entire content to new file"""
		with open(getTLEFilePath(sat_id), 'w') as fp:\
			# don't want to begin or end file with newline
			fp.write(self._stringify3LEDict(tle_dict_list[0]))
			for tle_dict in tle_dict_list[1:]:
				fp.write('\n')
				fp.write(self._stringify3LEDict(tle_dict))

	def _dictify3LEs(self, lines:list[str]) -> list[dict[int,dict[str,list[str]|str]]]:
		if len(lines) == 0:
			raise ValueError("No data")
		if len(lines)%3 != 0:
			# If not obvious what lines relate to eachother, can't make assumption -> abort.
			raise ValueError('Incomplete TLEs present, aborting')

		list_3les = []
		tle = {0:{'fields':[], 'line_str':''},
				1:{'fields':[], 'line_str':''},
				2:{'fields':[], 'line_str':''}}
		for ii, line in enumerate(lines):
			fields = line.split()
			# insert fields into and store original line in tle dict
			tle[int(fields[0])]['fields'] = fields
			tle[int(fields[0])]['line_str'] = line
			if int(fields[0]) == 2:
				# store parsed tle
				list_3les.append(tle.copy())
				# empty dict
				tle[0] = {'fields':[], 'line_str':''}
				tle[1] = {'fields':[], 'line_str':''}
				tle[2] = {'fields':[], 'line_str':''}

		return list_3les

	def _stringify3LEDict(self, tle_dict:dict[int,dict[str,list[str]|str]]) -> str:
		lines = [line_dict['line_str'] for line_dict in tle_dict.values()]
		tle_str = '\n'.join(lines)
		return tle_str

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