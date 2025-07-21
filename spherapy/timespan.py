import datetime as dt
from dateutil.relativedelta import relativedelta
import logging
import numpy as np
from astropy.time import Time as astropyTime
from skyfield.api import load

import spherapy.util.exceptions as exceptions

logger = logging.getLogger(__name__)


class TimeSpan(object):

	def __init__(self, t0, timestep='1S', timeperiod='10S'):
		"""Creates a series of timestamps in UTC.
		Difference between each timestamp = timestep
		Total duration = integer number of timesteps closest to timeperiod
			If timeperiod is an integer multiple of timestep,
			then TimeSpan[-1] - TimeSpan[0] = timeperiod
			If timeperiod is NOT an integer multiple of timestep,

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

		Raises
		------
		ValueError
		"""
		# convert and/or add timezone info
		if t0.tzinfo is not None:
			t0 = t0.astimezone(dt.timezone.utc)
		else:
			t0 = t0.replace(tzinfo=dt.timezone.utc)

		self.start = t0
		self.time_step = self._parseTimestep(timestep)
		self.end = self._parseTimeperiod(t0, timeperiod)
		self.time_period = self.end - self.start	

		logger.info("Creating TimeSpan between {} and {} with {} seconds timestep"
					.format(self.start, self.end, self.time_step.total_seconds()))


		if self.start >= self.end:
			logger.error("Timeperiod: {} results in an end date earlier than the start date".
							format(timeperiod))
			raise ValueError("Invalid timeperiod value: {}".format(timeperiod))

		if self.time_period < self.time_step:
			logger.error("Timeperiod: {} is smaller than the timestep {}.".
							format(timeperiod, timestep))
			raise ValueError("Invalid timeperiod value, too small: {}".format(timeperiod))

		# Calculate step numbers and correct for non integer steps.
		n_down = int(np.floor(self.time_period / self.time_step))
		if n_down != self.time_period / self.time_step:
			logger.info("Rounding to previous integer number of timesteps")
			self.end = n_down * self.time_step + self.start
			self.time_period = self.end - self.start
			
			logger.info("Adjusting TimeSpan to {} -> {} with {} seconds timestep".
						format(self.start, self.end, self.time_step.total_seconds()))


		self._timearr = np.arange(self.start.replace(tzinfo=None), self.end.replace(tzinfo=None)+self.time_step, self.time_step).astype(dt.datetime)
		self._timearr = np.vectorize(lambda x: x.replace(tzinfo=dt.timezone.utc))(self._timearr)
		self._skyfield_timespan = load.timescale()


    # Make it callable and return the data for that entry
	def __call__(self, idx=None):
		return self._timearr

	def __getitem__(self, idx=None):
		if idx == None:
			return self._timearr
		elif isinstance(idx, int) or isinstance(idx, np.integer):
			return self._timearr[idx]
		elif isinstance(idx, tuple):
			if len(idx) == 2:
				return self._timearr[idx[0]:idx[1]]
			else:
				return self._timearr[idx[0]:idx[1]:idx[2]]
		elif isinstance(idx, list):
			return self._timearr[[idx]]
		elif isinstance(idx,slice):
			return self._timearr[idx]
		else:
			raise TypeError('index is unknown type')

	def __eq__(self, other):
		if not isinstance(other, TimeSpan):
			return NotImplemented

		return np.all(self.asDatetime() == other.asDatetime())

	def __add__(self, other):
		self_copy = TimeSpan(self.start)

		self_copy.start = self.start
		self_copy.init_timestep_str = ''
		self_copy.init_timeperiod_str = ''
		self_copy.time_step = None

		self_copy.end = other.end
		self_copy.time_period = other.end - self_copy.start

		self_copy._timearr = np.hstack((self._timearr, other._timearr))

		return self_copy

	def asAstropy(self, *args, scale='utc'):
		"""
		Return ndarray of TimeSpan as astropy.time objects
		
		Returns
		-------
		ndarray
		"""
		# scale = kwargs.get('scale')
		if len(args) == 1 and isinstance(args[0], int):
			return astropyTime(self._timearr, scale=scale)[args[0]]
		else:
			return astropyTime(self._timearr, scale=scale)

	def asDatetime(self, *args):
		"""
		Return ndarray of TimeSpan as datetime objects	
		
		Returns
		-------
		ndarray
		"""
		if len(args) == 1 and isinstance(args[0], int):
			return self._timearr[args[0]]
		else:
			return self._timearr

	def asSkyfield(self, *args):
		"""
		Return TimeSpan element as Skyfield Time object
		
		Returns
		-------
		Skyfield Time
		"""
		if len(args) == 1 and isinstance(args[0], int):
			datetime = self._timearr[args[0]]
			datetime = datetime.replace(tzinfo=dt.timezone.utc)
			return self._skyfield_timespan.from_datetime(datetime)
		else:
			raise IndexError

	def asText(self, *args):
		if len(args) == 1 and isinstance(args[0], int):
			return self._timearr[args[0]].strftime("%Y-%m-%d %H:%M:%S")
		else:
			raise IndexError

	def secondsSinceStart(self):
		'''
		Return ndarray with the seconds of all timesteps since the beginning.
		'''
		diff = self._timearr - self.start
		days = np.vectorize(lambda x: x.days)(diff)
		secs = np.vectorize(lambda x: x.seconds)(diff)
		usecs = np.vectorize(lambda x: x.microseconds)(diff)
		return days * 86400 + secs + usecs * 1e-6

	def getClosest(self, t_search):
		"""Find the closest time in a TimeSpan
				
		Parameters
		----------
		t_search : {datetime}
			time to find
		
		Returns
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
		return self.asDatetime(res_index), res_index

	def _parseTimeperiod(self, t0, timeperiod):
		for index, letter in enumerate(timeperiod):
			if not (letter.isdigit() or letter == '.'):
				break

		val = float(timeperiod[0:index])
		unit = timeperiod[index:]
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
						+ " (d)ays, (H)ours, (M)inutes, (S)econds, (mS) milliseconds, (uS) microseconds")
			raise ValueError("Invalid timeperiod unit:{}".format(unit))

		return t0 + relativedelta(days=values[5], hours=values[4], minutes=values[3],
								seconds=values[2], microseconds=values[0])	# noqa: E126

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
			logger.error("Invalid timestep unit: Valid units are (d)ays, (H)ours,"
						+ "(M)inutes, (S)econds, (mS) milliseconds, (uS) microseconds")
			raise ValueError("Invalid timestep unit:{}".format(unit))

		return dt.timedelta(days=values[5], seconds=values[2], 
							microseconds=values[0], milliseconds=values[1],
							minutes=values[3], hours=values[4])

	def __len__(self):
		return len(self._timearr)


	def cherryPickFromIndices(self, idxs):
		self._timearr = self._timearr[idxs]
		self.start = self._timearr[0]
		self.end = self._timearr[-1]
		self.init_timeperiod_str = None
		self.init_timestep_str = None
		self.time_step = None
		self.time_period = self.end - self.start

	@classmethod
	def fromDatetime(cls, dt_arr:np.ndarray[dt.datetime], timezone=dt.timezone.utc):
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
		t.init_timeperiod_str = None
		t.init_timestep_str = None
		t.time_step = None
		t.time_period = t.end-t.start

		return t
