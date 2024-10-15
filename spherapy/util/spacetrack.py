import spacetrack as sp
import spherapy
import datetime as dt
import spherapy.util.epoch_u as epoch_u
import sys
import os
from progressbar import progressbar
import json

MAX_RETRIES=3

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
		return os.path.exists(getTLEFilePath(sat_id))

	def fetchAll(self, sat_id:int) -> str|None:
		retries = 0
		while retries < MAX_RETRIES:
			try:
				res_str = self.stc.tle(norad_cat_id=sat_id, orderby='epoch asc', limit=500000, format='3le')
				with open(getTLEFilePath(sat_id), 'w') as fp:
					if res_str[-1] == '\n':
						res_str = res_str[:-1]
					fp.write(res_str)
				break
			except TimeoutError as e:
				retries += 1
		if retries == MAX_RETRIES:
			print(f"Could not fetch All TLEs for sat {sat_id}: failed {retries} times.", file=sys.stderr)
			return None
		
		return sat_id
		

	def getNumPastTLEs(self, sat_id:int) -> int:
		with open(getTLEFilePath(sat_id), 'r') as fp:
			lines = fp.readlines()
		return int(len(lines)/3)

	def fetchLatest(self, sat_id:int) -> str|None:
		retries = 0
		while retries < MAX_RETRIES:
			try:
				# get penultimate and ultimate epochs
				with open(getTLEFilePath(sat_id), 'r') as fp:
					lines = fp.readlines()
				while lines[-1] == '':
					lines = lines[:-1]
				le_line = lines[-2]
				# pe_line = lines[-5]
				# pe_datetime = epoch_u.epoch2datetime(float(pe_line.split()[3]))
				try:
					le = float(le_line.split()[3])
				except IndexError:
					print(le_line)
					raise IndexError
				le_datetime = epoch_u.epoch2datetime(le)
				delta = dt.datetime.now(tz=dt.timezone.utc) - le_datetime
				if delta.days != 0:
					res_str = self.stc.tle(norad_cat_id=sat_id, orderby='epoch asc', epoch=f'>now-{delta.days+1}', limit=500000, format='3le')
					res_str = res_str[:-1]
					res_lines = res_str.split('\n')
					for ii, line in enumerate(res_lines):
						# Try block for debugging, only sometimes failing, trying to investigate.
						try:
							if line[0] == '1' and float(line.split()[3])>le:
								break
						except IndexError:
							print(f'{le=}')
							print(f'{ii=}')
							print(f'{line=}')
							print(f'{line.split()}')
							print(f'{res_lines}')
					next_index = ii
					if res_lines[0] != '':
						with open(getTLEFilePath(sat_id), 'a') as fp:
							fp.write('\n')
							fp.write('\n'.join(res_lines[next_index-1:]))
				break
			except TimeoutError as e:
				print(e)
				retries += 1
		if retries == MAX_RETRIES:
			print(f"Could not fetch the latest TLEs for sat {sat_id}: failed {retries} times.", file=sys.stderr)
			return None

		return sat_id

class InvalidCredentials(Exception):
	def __init__(self, message):
		super().__init__(message)
		return
	
def updateTLEs(sat_id_list:list[int], user:str=None, passwd:str=None) -> list[int]:
	print(f"Using SPACETRACK to update TLEs")
	g = TLEGetter(sat_id_list,user=user,passwd=passwd)
	
	return g.getModifiedIDs()

def getTLEFilePath(sat_id:int) -> str:
	return f'{spherapy.tle_path.absolute()}/{sat_id}.tle'

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