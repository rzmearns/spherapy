"""Utility functions for converting different epochs.

Attributes:
	GMST_epoch: [description]
"""

import datetime as dt
import pathlib

import numpy as np

from spherapy.util import elements_u

GMST_epoch = dt.datetime(2000, 1, 1, 12, 0, 0)


def epoch2datetime(epoch_str: str) -> dt.datetime:
	"""Converts a fractional epoch string to a datetime object.

	Args:
	epoch_str:	fractional year epoch

	Returns:
		equivalent datetime
	"""
	if not isinstance(epoch_str, str):
		epoch_str = str(epoch_str)

	year = int(epoch_str[:2])
	if year < 50:  # noqa: PLR2004
		year += 2000
	else:
		year += 1900

	fractional_day_of_year = float(epoch_str[2:])

	base = dt.datetime(year, 1, 1, tzinfo=dt.timezone.utc)
	return base + dt.timedelta(days=fractional_day_of_year) - dt.timedelta(days=1)

def epochEarlierThan(epoch_a:str, epoch_b:str) -> bool:
	"""Check if epoch A is earlier than epoch B.

	Args:
		epoch_a: TLE epoch string A
		epoch_b: TLE epoch string B

	Returns:
		True if epoch A is earlier than epoch B
	"""
	datetime_a = epoch2datetime(epoch_a)
	datetime_b = epoch2datetime(epoch_b)
	return datetime_a < datetime_b

def epochLaterThan(epoch_a:str, epoch_b:str) -> bool:
	"""Check if epoch A is later than epoch B.

	Args:
		epoch_a: TLE epoch string A
		epoch_b: TLE epoch string B

	Returns:
		True if epoch A is later than epoch B
	"""
	datetime_a = epoch2datetime(epoch_a)
	datetime_b = epoch2datetime(epoch_b)
	return datetime_a > datetime_b

def datetime2TLEepoch(date: dt.datetime) -> str:
	"""Converts a datetime to a TLE epoch string.

	Args:
		date: Datetime object

	Returns:
		TLE epoch string, with fractional seconds
	"""
	tzinfo = date.tzinfo
	year_str = str(date.year)[-2:]
	day_str = str(date.timetuple().tm_yday).zfill(3)
	fraction_str = str(
		(date - dt.datetime(date.year, date.month, date.day, tzinfo=tzinfo)).total_seconds()
		/ dt.timedelta(days=1).total_seconds()
	)[1:]

	return year_str + day_str + fraction_str


def datetime2sgp4epoch(date: dt.datetime) -> float:
	"""Converts a datetime to an sgp4 epoch.

	Args:
		date: Datetime object

	Returns:
		SGP4 epoch, with fractional seconds
	"""
	tzinfo = date.tzinfo
	sgp4start = dt.datetime(1949, 12, 31, 0, 0, 0, tzinfo=tzinfo)
	delta = date - sgp4start
	return delta.days + delta.seconds / 86400


def findClosestDatetimeIndices(test_arr:np.ndarray[tuple[int], np.dtype[np.datetime64]],
								source_arr:np.ndarray[tuple[int], np.dtype[np.datetime64]]) \
								-> np.ndarray[tuple[int], np.dtype[np.int64]]:
	"""Find the index of the closest datetime in source arr for each datetime in test_arr.

	Both search_arr and source_arr must be sorted.

	Args:
		test_arr: Mx1 array of datetimes, will find closest time for each element in this array
		source_arr: Nx1: array of datetimes to compare to

	Returns:
		Mx1 array of indices of source_arr
	"""
	# convert to unix timestamps first to save memory
	test_arr_fl = np.vectorize(lambda x:x.timestamp())(test_arr)
	source_arr_fl = np.vectorize(lambda x:x.timestamp())(source_arr)
	idx_arr = np.zeros(test_arr_fl.shape, dtype=np.int64)
	for ii in range(len(test_arr_fl)):
		diff = np.abs(source_arr_fl-test_arr_fl[ii])
		idx_arr[ii] = np.argmin(diff,axis=0)

	return idx_arr

def getStoredEpochs(tle_path:pathlib.Path) -> None|tuple[dt.datetime, dt.datetime|None]:
	"""Return the start and end epoch for tle_path.

	Args:
		tle_path: tle file

	Returns:
		(first epoch datetime, last epoch datetime)
		None if no spacetrack tle stored for sat_id
	"""
	if not tle_path.exists():
		return None

	with tle_path.open('r') as fp:
		lines = fp.readlines()

	first_tle_line_1 = elements_u.split3LELineIntoFields(lines[1])
	last_tle_line_1 = elements_u.split3LELineIntoFields(lines[-2])

	first_epoch_dt = epoch2datetime(first_tle_line_1['fields'][3])
	last_epoch_dt = epoch2datetime(last_tle_line_1['fields'][3])

	return (first_epoch_dt, last_epoch_dt)
