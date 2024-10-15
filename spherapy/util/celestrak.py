import spherapy
import spacetrack
import requests

MAX_RETRIES=3

def updateTLEs(sat_id_list:list[int], user:str=None, passwd:str=None) -> list[int]:
	modified_list = []
	for sat_id in sat_id_list:
		url = f'https://celestrak.org/NORAD/elements/gp.php?CATNR={sat_id}'
		for ii in range(MAX_RETRIES):
			r = requests.get(url)
			if r.status_code == 200:
				break
		
		if ii == MAX_RETRIES-1:
			raise ValueError(f'Could not fetch celestrak information for sat_id: {sat_id}')
		
		fname = getTLEFilePath(sat_id)
		dat_list = r.text.split('\r\n')
		with open(fname, 'w') as fp:
			fp.write(f'0 {dat_list[0].rstrip()}\n')
			fp.write(f'{dat_list[1].rstrip()}\n')
			fp.write(f'{dat_list[2].rstrip()}')
			modified_list.append(sat_id)
		
	return modified_list

def getTLEFilePath(sat_id:int) -> str:
	return f'{spherapy.tle_path.absolute()}/{sat_id}.temptle'