"""Class for series of timestamps.

This module provides:
- TimeSpan: a series of timestamps to be used by an orbit object
"""
import datetime as dt
import logging

from typing import cast
from typing_extensions import Self

from astropy.time import Time as astropyTime
from dateutil.relativedelta import relativedelta
import numpy as np
from skyfield.api import load
import skyfield.timelib

logger = logging.getLogger(__name__)


class TimeSpan:
	"""A series of timestamps.

	Attributes:
		start: The first timestamp
		end: The last timestamp
		time_step: The difference between timestamps in seconds
						Can be None if irregular steps
		time_period: The difference between end and start in seconds

	"""
	def __init__(self, t0:dt.datetime, timestep:str='1S', timeperiod:str='10S'):
		"""Creates a series of timestamps in UTC.

		Difference between each timestamp = timestep
		Total duration = greatest integer number of timesteps less than timeperiod
			If timeperiod is an integer multiple of timestep,
			then TimeSpan[-1] - TimeSpan[0] = timeperiod
			If timeperiod is NOT an integer multiple of timestep,
			then TimeSpan[-1] = int(timeperiod/timestep) (Note: int cast, not round)

		Does not account for Leap seconds, similar to all Posix compliant UTC based time
			representations. see: https://numpy.org/doc/stable/reference/arrays.datetime.html#datetime64-shortcomings
			for equivalent shortcomings.
		Always contains at least two timestamps

		Parameters
		t0: datetime defining the start of the TimeSpan.
				If timezone naive, assumed to be in UTC
				If timezone aware, will be converted to UTC
		timestep: String describing the time step of the time span.
					The string is constructed as an integer or float, followed by a time unit:
					(d)ays, (H)ours, (M)inutes, (S)econds, (mS) milliseconds, (uS) microseconds
					(the default is '1S')
		timeperiod:	String describing the time period of the time span.
					The string is constructed as an integer or float, followed by a time unit:
					(d)ays, (H)ours, (M)inutes, (S)econds, (mS) milliseconds, (uS) microseconds
					(the default is '1d')

		Raises:
		------
		ValueError
		"""
		# convert and/or add timezone info
		if t0.tzinfo is not None:
			t0 = t0.astimezone(dt.timezone.utc)
		else:
			t0 = t0.replace(tzinfo=dt.timezone.utc)

		self.start:dt.datetime = t0
		self.time_step:None|dt.timedelta = self._parseTimestep(timestep)
		self.end:dt.datetime = self._parseTimeperiod(t0, timeperiod)
		self.time_period = self.end - self.start

		logger.info("Creating TimeSpan between %s and %s with %s seconds timestep",
						self.start, self.end, self.time_step.total_seconds())


		if self.start >= self.end:
			logger.error("Timeperiod: %s results in an end date earlier than the start date",
							timeperiod)
			raise ValueError(f"Invalid timeperiod value: {timeperiod}")

		if self.time_period < self.time_step:
			logger.error("Timeperiod: %s is smaller than the timestep %s.", timeperiod, timestep)
			raise ValueError(f"Invalid timeperiod value, too small: {timeperiod}")

		# Calculate step numbers and correct for non integer steps.
		n_down = int(np.floor(self.time_period / self.time_step))
		if n_down != self.time_period / self.time_step:
			logger.info("Rounding to previous integer number of timesteps")
			self.end = n_down * self.time_step + self.start
			self.time_period = self.end - self.start

			logger.info("Adjusting TimeSpan to %s -> %s with %s seconds timestep",
							self.start, self.end, self.time_step.total_seconds())


		self._timearr = np.arange(self.start.replace(tzinfo=None),
									self.end.replace(tzinfo=None)+self.time_step,
									self.time_step).astype(dt.datetime)
		self._timearr = np.vectorize(lambda x: x.replace(tzinfo=dt.timezone.utc))(self._timearr)
		self._skyfield_timespan = load.timescale()


	def __hash__(self) -> int:
		"""Returns a hash of the timespan."""
		return hash((self.start, self.end, self.time_step))

	# Make it callable and return the data for that entry
	def __call__(self) -> np.ndarray[tuple[int], np.dtype[np.datetime64]]:
		"""Returns the internal _timearr when the TimeSpan is called."""
		return self._timearr

	def __getitem__(self, idx:None|int|np.integer|tuple|list|np.ndarray|slice=None) \
						-> None|dt.datetime|np.ndarray[tuple[int], np.dtype[np.datetime64]]: 	# noqa: PLR0911
		"""Returns an index or slice of the TimeSpan as an array of datetime objects."""
		if idx is None:
			return self._timearr
		if isinstance(idx, int|np.integer):
			return self._timearr[idx]
		if isinstance(idx, tuple):
			if len(idx) == 2: 								#noqa: PLR2004
				return self._timearr[idx[0]:idx[1]]
			return self._timearr[idx[0]:idx[1]:idx[2]]
		if isinstance(idx, list):
			return self._timearr[[idx]]
		if isinstance(idx, np.ndarray):
			return self._timearr[idx]
		if isinstance(idx,slice):
			return self._timearr[idx]
		raise TypeError('index is unknown type')

	def __eq__(self, other:object) -> bool:
		"""Checks if other TimeSpan is equal to self."""
		if not isinstance(other, TimeSpan):
			return NotImplemented

		if len(self) != len(other):
			return False

		return bool(np.all(self.asDatetime() == other.asDatetime()))

	def __add__(self, other:Self) -> 'TimeSpan':
		"""Apeends other TimeSpan to this one, does not order timesteps."""
		self_copy = TimeSpan(self.start)

		self_copy.start = self.start
		self_copy.time_step = None

		self_copy.end = other.end
		self_copy.time_period = other.end - self_copy.start

		self_copy._timearr = np.hstack((self._timearr, other._timearr))

		return self_copy

	def asAstropy(self, idx:None|int=None, scale:str='utc') -> astropyTime:
		"""Return ndarray of TimeSpan as astropy.time objects.

		Args:
			idx: timestamp index to return inside an astropyTime object
					if no index supplied, returns all timestamps
			scale: astropy time scale, can be one of
					('tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1', 'utc')

		Returns:
		-------
		ndarray
		"""
		if idx is None:
			return astropyTime(self._timearr, scale=scale)
		return astropyTime(self._timearr[idx], scale=scale)

	def asDatetime(self, idx:None|int=None) \
		-> dt.datetime|np.ndarray[tuple[int], np.dtype[np.datetime64]]:
		"""Return ndarray of TimeSpan as datetime objects.

		Args:
			idx: timestamp index to return as a datetime
					If no index supplied, returns whole array

		Returns:
		-------
		ndarray
		"""
		if idx is None:
			return self._timearr
		return self._timearr[idx]

	def asSkyfield(self, idx:int) -> skyfield.timelib.Time:
		"""Return TimeSpan element as Skyfield Time object.

		Args:
			idx: timestamp index to return as skyfield time object

		Returns:
		-------
		Skyfield Time
		"""
		datetime = self._timearr[idx]
		datetime = datetime.replace(tzinfo=dt.timezone.utc)
		return self._skyfield_timespan.from_datetime(datetime)

	def asText(self, idx:int) -> str:
		"""Returns a text representation of a particular timestamp.

		Timestamp will be formatted as YYYY-mm-dd HH:MM:SS

		Args:
			idx: timestamp index to format

		Returns:
			str:
		"""
		return self._timearr[idx].strftime("%Y-%m-%d %H:%M:%S")

	def secondsSinceStart(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
		"""Return ndarray with the seconds of all timesteps since the beginning.

		Returns:
			array of seconds since start for each timestamp
		"""
		diff = self._timearr - self.start
		days = np.vectorize(lambda x: x.days)(diff)
		secs = np.vectorize(lambda x: x.seconds)(diff)
		usecs = np.vectorize(lambda x: x.microseconds)(diff)
		return days * 86400 + secs + usecs * 1e-6

	def getClosest(self, t_search:dt.datetime) -> tuple[dt.datetime, int]:
		"""Find the closest time in a TimeSpan.

		Parameters
		----------
		t_search : datetime to search for in TimeSpan
				If timezone naive, assumed to be in UTC
				If timezone aware, will be converted to UTCtime to find

		Returns:
		-------
		datetime, int
			Closest datetime in TimeSpan, index of closest date in TimeSpan
		"""
		# convert and/or add timezone info
		if t_search.tzinfo is not None:
			t_search = t_search.astimezone(dt.timezone.utc)
		else:
			t_search = t_search.replace(tzinfo=dt.timezone.utc)
		diff = self._timearr - t_search
		out = np.abs(np.vectorize(lambda x: x.total_seconds())(diff))
		res_index = int(np.argmin(out))
		res_datetime = self.asDatetime(res_index)

		# res_datetime must be a dt.datetime since res_index is an int, so cast
		res_datetime = cast("dt.datetime", res_datetime)

		return res_datetime, res_index

	def areTimesWithin(self, t_search:dt.datetime|np.ndarray[tuple[int], np.dtype[np.datetime64]])\
						-> np.ndarray[tuple[int],np.dtype[np.bool_]]:
		"""Find if the provided times are within the timespan.

		Args:
			t_search: times to check if within timespan
				If timezone naive, assumed to be in UTC
				If timezone aware, will be converted to UTCtime to find

		Returns:
			ndarray of bools, True if within timespan
		"""
		if isinstance(t_search, dt.datetime):
			if t_search.tzinfo is not None:
				t_search = t_search.astimezone(dt.timezone.utc)
			else:
				t_search = t_search.replace(tzinfo=dt.timezone.utc)

		elif isinstance(t_search, np.ndarray):
			if t_search[0].tzinfo is not None:
				t_search = np.vectorize(lambda x:x.astimezone(dt.timezone.utc))(t_search)
			else:
				t_search = np.vectorize(lambda x:x.replace(tzinfo=dt.timezone.utc))(t_search)

		else:
			raise TypeError('t_search is of an unrecognised type,'
							' should be a datetime object or ndarray')

		return np.logical_and(t_search>np.asarray(self.start), t_search<np.asarray(self.end))

	def getFractionalIndices(self, t_search:dt.datetime|np.ndarray[tuple[int],
											np.dtype[np.datetime64]])\
							-> np.ndarray[tuple[int],np.dtype[np.float64]]:
		"""Find the fractional indices of timespan at wich t_search should be inserted.

		Find the indices in the original timespan at which each value of t_search should be
		inserted to maintain the sorted order.
		The integer part of the index indicates the value immediately prior to the value in
		t_search, while the fractional part represents the point between the two adjacent indices
		in the timespan at which the t_search value falls.
		For example (using integers rather than datetime objects:
			timespan = [0, 1, 5, 6, 10]
			t_search = [0.5, 2, 3, 8]
			timespan.getFractionalIndices(t_search) = [0.5, 1.25, 1.5, 3.5]
		All values of t_search must be within the timespan, otherwise the output is undefined.

		Args:
			t_search: times to locate within timespan
				If timezone naive, assumed to be in UTC
				If timezone aware, will be converted to UTCtime to find

		Returns:
			ndarray of fractional indices
		"""
		if isinstance(t_search, dt.datetime):
			if t_search.tzinfo is not None:
				t_search = t_search.astimezone(dt.timezone.utc)
			else:
				t_search = t_search.replace(tzinfo=dt.timezone.utc)
			t_search = np.asarray(t_search)
		elif isinstance(t_search, np.ndarray):
			if t_search[0].tzinfo is not None:
				t_search = np.vectorize(lambda x:x.astimezone(dt.timezone.utc))(t_search)
			else:
				t_search = np.vectorize(lambda x:x.replace(tzinfo=dt.timezone.utc))(t_search)
		else:
			raise TypeError('t_search is of an unrecognised type, '
							'should be a datetime object or ndarray')

		post_idxs = np.searchsorted(self._timearr, t_search) 		# type: ignore [arg-type]
		pre_idxs = post_idxs-1
		return (t_search - self._timearr[pre_idxs])\
				/(self._timearr[post_idxs]-self._timearr[pre_idxs]) + pre_idxs


	def _parseTimeperiod(self, t0:dt.datetime, timeperiod:str) -> dt.datetime:
		last_idx = 0
		for idx, letter in enumerate(timeperiod):
			last_idx = idx
			if not (letter.isdigit() or letter == '.'):
				break

		val = float(timeperiod[0:last_idx])
		unit = timeperiod[last_idx:]
		values = np.zeros(6)
		if unit == 'd':
			values[5] = val
		elif unit == 'H':
			values[4] = val
		elif unit == 'M':
			values[3] = val
		elif unit == 'S':
			values[2] = val
		elif unit == 'mS':
			values[0] = 1000 * val
		elif unit == 'uS':
			values[0] = val
		else:
			logger.error("Invalid timeperiod unit: Valid units are (y)ears, (m)onths, (W)eeks,"
						" (d)ays, (H)ours, (M)inutes, (S)econds, (mS) milliseconds,"
						" (uS) microseconds")
			raise ValueError(f"Invalid timeperiod unit:{unit}")

		return t0 + relativedelta(days=values[5], hours=values[4], minutes=values[3],
								seconds=values[2], microseconds=values[0])

	def _parseTimestep(self, timestep:str) -> dt.timedelta:
		last_idx = 0
		for idx, letter in enumerate(timestep):
			last_idx = idx
			if not (letter.isdigit() or letter == '.'):
				break

		val = float(timestep[0:last_idx])
		unit = timestep[last_idx:]
		values = np.zeros(6)
		if unit == 'd':
			values[5] = val
		elif unit == 'H':
			values[4] = val
		elif unit == 'M':
			values[3] = val
		elif unit == 'S':
			values[2] = val
		elif unit == 'mS':
			values[1] = val
		elif unit == 'uS':
			values[0] = val
		else:
			logger.error("Invalid timestep unit: Valid units are (d)ays, (H)ours," \
						" (M)inutes, (S)econds, (mS) milliseconds, (uS) microseconds")
			raise ValueError(f"Invalid timestep unit:{unit}")

		return dt.timedelta(days=values[5], seconds=values[2],
							microseconds=values[0], milliseconds=values[1],
							minutes=values[3], hours=values[4])

	def __len__(self) -> int:
		"""Returns the length of the TimeSpan."""
		return len(self._timearr)


	def cherryPickFromIndices(self, idxs:int|tuple|slice):
		"""Adjust TimeSpan to only contain the indices specified by idxs.

		Args:
			idxs: numpy style indexing of TimeSpan
		"""
		self._timearr = self._timearr[idxs]
		self.start = self._timearr[0]
		self.end = self._timearr[-1]
		self.time_step = None
		self.time_period = self.end - self.start

	@classmethod
	def fromDatetime(cls, dt_arr:np.ndarray[tuple[int], np.dtype[np.datetime64]],
							timezone:dt.timezone=dt.timezone.utc) -> 'TimeSpan':
		"""Create a TimeSpan from an array of datetime objects.

		Args:
			dt_arr: 1D array of datetime objects
			timezone: timezone to apply to each element of datetime array.
				Default: dt.timezone.utc
		"""
		for ii in range(len(dt_arr)):
			dt_arr[ii] = dt_arr[ii].replace(tzinfo=timezone)
		start = dt_arr[0]
		end = dt_arr[-1]
		tperiod = (end-start).total_seconds()
		tstep = tperiod/len(dt_arr)
		t = cls(start,f'{tstep}S',f'{tperiod}S')
		t._timearr = dt_arr.copy()
		t.start = start
		t.end = end
		t.time_step = None
		t.time_period = t.end-t.start

		return t
