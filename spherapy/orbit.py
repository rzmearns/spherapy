
import numpy as np
import sys
import logging
from progressbar import progressbar
import pickle

from astropy import units as u

from hapsira.bodies import Earth, Mars, Sun, Moon
from hapsira.twobody import Orbit as hapsiraOrbit
from hapsira.ephem import Ephem

from sgp4.api import Satrec, WGS72
from astropy.time import Time as astropyTime
import datetime as dt
# Don't need to import all of skyfield just for EarthSatellite loading.
# TODO: figure out the lightest import to use
from skyfield.api import load, EarthSatellite

import spherapy.util.list_u as list_u
import spherapy.util.epoch_u as epoch_u
import spherapy.util.orbital_u as orbit_u
import spherapy.util.exceptions as exceptions
import spherapy.util.constants as consts

logger = logging.getLogger(__name__)

# TODO:
# functional call to get position at time


class Orbit(object):
	'''
	Contains timestepped array of orbital position and velocity data for a satellite, coordinate system depending on the central body.
	Contains timestepped array of sun position.
	Timing data taken from a TimeSpan object.
		
	Attributes
	----------
	timespan: {satplot.TimeSpan}
		Timespan over which orbit is to be simulated
		
	pos: (N, 3) np.array
		Cartesian coordinates of the position of the satellite at each time in a TimeSpan (units: km)
		Coordinate system depending on the central body.
		
	vel: (N, 3) np.array
		Cartesian coordinates of the velocity of the satellite at each time in a TimeSpan (units: m/s)
		Coordinate system depending on the central body.

	sun: (N, 3) np.array
		Vector of the sun's position at each timestep (GCRS; units: km)
		
	period: float
		Period of the satellite's orbit (units: s)
		
	period_steps: float
		Number of timesteps per orbital period of the satellite.
		Useful to check the timestep is appropriate for this orbit.
		Recommended: 10 < period_steps < 100 or so? Should we write warnings if outside this range?
	
	'''
	def __init__(self, *args, **kwargs):
		'''
		The constructor should never be called directly.
		Use one of:
			Orbit.from_tle()
			Orbit.from_tle_orbital_param()
			Orbit.from_orbital_param()
			Orbit.from_list_of_positions()
		'''
		# Should always be called from a class method, however,
		# If no self.gen_type, spit an error here
		self.gen_type = kwargs.get('type')
		if self.gen_type is None:
			logger.error("Orbit() should not be called directly.")
			raise ValueError("Orbit() should not be called directly.")

		if len(args) == 0:
			logger.error("Orbit() should not be called directly.")
			raise ValueError("Orbit() should not be called directly.")			

		self.calc_astrobodies = kwargs.get('astrobodies')

		self.timespan = args[0]
		
		# CONSTRUCTOR HANDLING
		if self.gen_type == 'TLE':
			# List of skyfield EarthSatellites, one for each TLE in the file
			sat_list = args[1]
			tle_dates = [sat.epoch.utc_datetime() for sat in sat_list]
			ts = load.timescale(builtin=True)
			if self.timespan.start < tle_dates[0] - dt.timedelta(days=14):
				logger.error("Timespan begins before provided TLEs (+14 days)")
				print("Timespan begins before provided TLEs (+14 days)", file=sys.stderr)
				# raise exceptions.OutOfRange("Timespan begins before provided TLEs (+14 days)")
			elif self.timespan.start > tle_dates[-1] + dt.timedelta(days=14):
				logger.error("Timespan begins after provided TLEs (+14 days)")
				print("Timespan begins after provided TLEs (+14 days)", file=sys.stderr)
				# raise exceptions.OutOfRange("Timespan begins after provided TLEs (+14 days)")
			elif self.timespan.end > tle_dates[-1] + dt.timedelta(days=14):
				logger.error("Timespan ends after provided TLEs (+14 days)")
				print("Timespan ends after provided TLEs (+14 days)", file=sys.stderr)
			# 	raise exceptions.OutOfRange("Timespan ends after provided TLEs (+14 days)")

			# Find closest listed TLE to start date
			start_datetime, start_index = list_u.get_closest(tle_dates, self.timespan.start)
			end_datetime, end_index = list_u.get_closest(tle_dates, self.timespan.end)
			if start_index == -1 or start_index == len(tle_dates):
				valid_sats = [sat_list[-1]]
				valid_tle_epoch_dates = [tle_dates[-1]]
			elif end_index == -1 or end_index == len(tle_dates):
				valid_sats = sat_list[start_index:]
				valid_tle_epoch_dates = tle_dates[start_index:]
			else:
				valid_sats = sat_list[start_index:end_index+1]
				valid_tle_epoch_dates = tle_dates[start_index:end_index+1]

			timesteps = self.timespan.asDatetime()
			
			if len(valid_tle_epoch_dates) > 1:
				valid_epoch_middates = np.diff(valid_tle_epoch_dates)/2 + np.asarray(valid_tle_epoch_dates[:-1])
							
				for ii, date in enumerate(valid_epoch_middates):
					valid_epoch_middates[ii] = date.replace(tzinfo=self.timespan.timezone)

				trans_indices = []
				for date in valid_epoch_middates:
					trans_indices.append(np.argmin(np.abs(np.vectorize(dt.timedelta.total_seconds)(timesteps-date))))

				sat_rec = valid_sats[0].at(ts.utc(timesteps[:trans_indices[0]]))
				pos = sat_rec.position.km
				vel = sat_rec.velocity.km_per_s * 1000
				self.pos = pos.T
				self.vel = vel.T
				TLE_epochs = np.tile(valid_tle_epoch_dates[0],trans_indices[0])
				self.TLE_epochs = TLE_epochs

				for ii in range(len(trans_indices)-1):
					sat_rec = valid_sats[ii+1].at(ts.utc(timesteps[trans_indices[ii]:trans_indices[ii+1]]))
					pos = sat_rec.position.km
					vel = sat_rec.velocity.km_per_s * 1000
					self.pos = np.vstack((self.pos, pos.T))
					self.vel = np.vstack((self.vel, vel.T))
					TLE_epochs = np.tile(valid_tle_epoch_dates[ii+1],trans_indices[ii+1]-trans_indices[ii])
					self.TLE_epochs = np.concatenate((self.TLE_epochs,TLE_epochs))

				sat_rec = valid_sats[-1].at(ts.utc(timesteps[trans_indices[-1]:]))
				pos = sat_rec.position.km
				vel = sat_rec.velocity.km_per_s * 1000
				self.pos = np.vstack((self.pos, pos.T))
				self.vel = np.vstack((self.vel, vel.T))
				TLE_epochs = np.tile(valid_tle_epoch_dates[ii+1],len(timesteps) - trans_indices[-1])
				self.TLE_epochs = np.concatenate((self.TLE_epochs,TLE_epochs))

			else:
				sat_rec = valid_sats[0].at(ts.utc(timesteps))
				pos = sat_rec.position.km
				vel = sat_rec.velocity.km_per_s * 1000
				self.pos = pos.T
				self.vel = vel.T
				self.TLE_epochs = np.tile(valid_tle_epoch_dates[0],len(timesteps))
		
			self.period = 2 * np.pi / valid_sats[-1].model.no_kozai * 60
			self.period_steps = int(self.period / self.timespan.time_step.total_seconds())
			self.name = valid_sats[-1].name
			try:
				self.sat = valid_sats[-1]
			except:
				# TODO: figure out what the exact error is and cause of when valid_sats doesn't exist
				pass

		elif self.gen_type == 'FAKE_TLE':
			satrec = args[1]
			a = kwargs.get('a')
			ts = load.timescale(builtin=True)
			sat = EarthSatellite.from_satrec(satrec, ts)
			# self.sat = sat

			self.pos = np.empty((self.timespan.num_steps, 3))
			self.vel = np.empty((self.timespan.num_steps, 3))

			for ii, timestep in enumerate(self.timespan.asDatetime()):	
				self.pos[ii, :] = sat.at(ts.utc(timestep.replace(tzinfo=self.timespan.timezone))).position.km 
				self.vel[ii, :] = sat.at(ts.utc(timestep.replace(tzinfo=self.timespan.timezone))).velocity.km_per_s * 1000

			self.period = 2 * np.pi / sat.model.no_kozai * 60
			self.period_steps = int(self.period / self.timespan.time_step.total_seconds())
			self.names = [kwargs['name']]


		elif self.gen_type == 'ANALYTICAL':
			body = kwargs.get('body')

			a = kwargs.get('a') * u.km
			ecc = kwargs.get('ecc') * u.one
			inc = kwargs.get('inc') * u.deg
			raan = kwargs.get('raan') * u.deg
			argp = kwargs.get('argp') * u.deg
			mean_nu = kwargs.get('mean_nu') * u.deg
			
			logger.info("Creating analytical orbit")
			self.orb = hapsiraOrbit.from_classical(body, a, ecc, inc, raan, argp, mean_nu)

			logger.info("Creating ephemeris for orbit, using timespan")
			self.ephem = Ephem.from_orbit(self.orb, self.timespan.asAstropy())

			self.pos = np.asarray(self.ephem.rv()[0])
			self.vel = np.asarray(self.ephem.rv()[1]) * 1000

			self.period = self.orb.period.unit.in_units('s') * self.orb.period.value
			self.period_steps = int(self.period / self.timespan.time_step.total_seconds())
			self.names = [kwargs['name']]

		elif self.gen_type == 'POS_LIST':
			self.pos = args[1]
			# Assume linear motion between each position at each timestep; 
			# Then assume it stops at the last timestep.
			self.vel = self.pos[1:] - self.pos[:-1]
			self.vel = np.concatenate((self.vel, np.array([[0, 0, 0]])))
			self.period = 1
			self.period_steps = 1
			self.names = ['Sat from position list']
			# Note that this doesn't define a period.
			logger.warning("Warning: When generating satellite orbit from list of positions, `period` and `period_steps` will not be defined.")
		else: 
			logger.error("Invalid orbit generation option {}. Valid options are TLE, FAKE_TLE, POS_LIST, and ANALYITCAL.".format(self.gen_type))
			raise ValueError("Invalid orbit generation option {}.".format(self.gen_type))

		if self.calc_astrobodies:
			logger.info('Creating ephemeris for sun using timespan')
			# Timescale for sun position calculation should use TDB, not UTC
			# The resultant difference is likely very small
			ephem_sun = Ephem.from_body(Sun, astropyTime(self.timespan.asAstropy(scale='tdb')), attractor=Earth)
			self.sun = np.asarray(ephem_sun.rv()[0].to(u.km))
			ephem_moon = Ephem.from_body(Moon, astropyTime(self.timespan.asAstropy(scale='tdb')), attractor=Earth)
			self.moon = np.asarray(ephem_moon.rv()[0].to(u.km))
		else:
			self.sun = None
			self.moon = None

	@classmethod	
	def fromListOfPositions(cls, timespan, positions, astrobodies=True):
		"""Create an orbit by explicitly specifying the position of the 
		satellite at each point in time. Useful for simplified test cases; but
		may lead to unphysical orbits.
		
		Parameters
		----------
		timespan : {satplot.TimeSpan}
			Timespan over which orbit is to be simulated
			
		positions: (N,3) np.array
			Position of the satellite in GCRS at each point in time.
			
		Returns
		-------
		satplot.Orbit
		"""
		return cls(timespan, positions, type='POS_LIST', astrobodies=astrobodies)
	
	@classmethod	
	def fromTLE(cls, timespan, tle_path, fp=sys.stdout, astrobodies=True):
		"""Create an orbit from an existing TLE or a list of historical TLEs
				
		Parameters
		----------
		timespan : {satplot.TimeSpan}
			Timespan over which orbit is to be simulated
		tle_path : {path to file containing TLE(s)}
			A plain text file containing the TLE or list of TLE(s) with each TLE line on a new line in the file.
		
		Returns
		-------
		satplot.Orbit
		"""
		sat_list = load.tle_file(tle_path)
		return cls(timespan, sat_list, type='TLE', astrobodies=astrobodies)

	@classmethod	
	def multiFromTLE(cls, timespan, tle_path, fp=sys.stdout, astrobodies=True):
		"""Create an orbit from an existing TLE or a list of historical TLEs
				
		Parameters
		----------
		timespan : {satplot.TimeSpan}
			Timespan over which orbit is to be simulated
		tle_path : {path to file containing TLE(s)}
			A plain text file containing the TLE or list of TLE(s) with each TLE line on a new line in the file.
		
		Returns
		-------
		satplot.Orbit
		"""
		with open(tle_path) as fp:
			lines = fp.readlines()
		num_sats = int(len(lines)/3)
		orbit_list = []
		for ii in progressbar(range(num_sats)):
			pc = ii/num_sats*100
			bar_str = int(pc)*'='
			space_str = (100-int(pc))*'  '
			print(f'Loading {pc:.2f}% ({ii} of {num_sats}) |{bar_str}{space_str}|\r')
			with open('temp.tle','w') as fp:
				for jj in range(3):
					fp.write(lines[ii*3+jj])
			sat_list = load.tle_file('temp.tle')
			orbit_list.append(cls(timespan, sat_list, type='TLE', astrobodies=astrobodies))
		return orbit_list


	@classmethod
	def fromTLEOrbitalParam(cls, timespan, a=6978, ecc=0, inc=0, raan=0, argp=0, mean_nu=0, name='Fake TLE', astrobodies=True):
		"""Create an orbit from orbital parameters, propagated using sgp4.
		
		Orbits created using this class method will respect gravity corrections such as J4, 
		allowing for semi-analytical sun-synchronous orbits.
		
		Parameters
		----------
		timespan : {satplot.TimeSpan}
			Timespan over which orbit is to be simulated
		a : {float}, optional
			semi-major axis of the orbit in km (the default is 6978 ~ 600km 
			above the earth.)
		ecc : {float}, optional
			dimensionless number 0 < ecc < 1 (the default is 0, which is a 
			circular orbit)
		inc : {float}, optional
			inclination of orbit in degrees (the default is 0, which represents 
			an orbit around the Earth's equator)
		raan : {float}, optional
			right-ascension of the ascending node (the default is 0)
		argp : {float}, optional
			argument of the perigee in degrees (the default is 0, which 
			represents an orbit with its semimajor axis in the plane of the 
			Earth's equator)
		mean_nu : {float}, optional
			mean anomaly in degrees (the default is 0, which represents an orbit
			that is beginning at periapsis)
		
		Returns
		-------
		satplot.Orbit
		"""
		if a < consts.R_EARTH + consts.EARTH_MIN_ALT:
			logger.error("Semimajor axis, {}, is too close to Earth".format(a))
			raise exceptions.OutOfRange("Semimajor axis, {}, is too close to Earth".format(a))
		if ecc > 1 or ecc < 0:
			logger.error("Eccentricity, {}, is non circular or eliptical".format(ecc))
			raise exceptions.OutOfRange("Eccentricity, {}, is non circular or eliptical".format(ecc))
		if inc > 180 or inc < -180:
			logger.error("Inclination, {}, is out of range, should be -180 < inc < 180".format(inc))
			raise exceptions.OutOfRange("Inclination, {}, is out of range, should be -180 < inc < 180".format(inc))
		if raan > 360 or raan < 0:
			logger.error("RAAN, {}, is out of range, should be 0 < inc < 360".format(inc))
			raise exceptions.OutOfRange("RAAN, {}, is out of range, should be 0 < inc < 360".format(inc))
		if argp > 360 or argp < 0:
			logger.error("Argument of periapsis, {}, is out of range, should be 0 < argp < 360".format(inc))
			raise exceptions.OutOfRange("Argument of periapsis, {}, is out of range, should be 0 < argp < 360".format(inc))
		if mean_nu > 360 or mean_nu < 0:
			logger.error("Mean anomaly, {}, is out of range, should be 0 < mean_nu < 360".format(inc))
			raise exceptions.OutOfRange("Mean anomaly, {}, is out of range, should be 0 < mean_nu < 360".format(inc))

		t0_epoch = epoch_u.datetime2sgp4epoch(timespan.start)
		satrec = Satrec()
		mean_motion = orbit_u.calc_meanmotion(a * 1e3)

		satrec.sgp4init(WGS72,		  # gravity model  # noqa: E128
								'i',  # mode  # noqa: E128
								1,  # satnum  # noqa: E128
								t0_epoch,  # mode  # noqa: E128
								0,  # bstar [/earth radii] # noqa: E128
								0,  # ndot [revs/day] # noqa: E128
								0,  # nddot [revs/day^3]  # noqa: E128
								ecc,  # ecc  # noqa: E128
								np.deg2rad(argp),  # arg perigee [radians] # noqa: E128
								np.deg2rad(inc),  # inclination [radians] # noqa: E128
								np.deg2rad(mean_nu),  # mean anomaly [radians] # noqa: E128
								mean_motion * 60,  # mean motion [rad/min]  # noqa: E128
								np.deg2rad(raan))  # raan [radians] # noqa: E126, E128

		return cls(timespan, satrec, a=a, type='FAKE_TLE', astrobodies=astrobodies)
		
	@classmethod
	def fromOrbitalParam(cls, timespan, body='Earth', a=6978, ecc=0, inc=0, raan=0, argp=0, mean_nu=0, name='Analytical', astrobodies=True):
		"""Create an orbit from orbital parameters
		
		Parameters
		----------
		timespan : {satplot.TimeSpan}
			Timespan over which orbit is to be simulated
		body : {str}, optional
			Central body for the satellite: ['Earth', 'Moon', 'Mars', 'Sun'] (the default is 'Earth')
		a : {float}, optional
			semi-major axis of the orbit in km (the default is 6978 ~ 600km above the earth.)
		ecc : {float}, optional
			dimensionless number 0 < ecc < 1 (the default is 0, which is a circular orbit)
		inc : {float}, optional
			inclination of orbit in degrees (the default is 0, which represents 
			an orbit around the Earth's equator)
		raan : {float}, optional
			right-ascension of the ascending node (the default is 0)
		argp : {float}, optional
			argument of the perigee in degrees (the default is 0, which 
			represents an orbit with its semimajor axis in the plane of the 
			Earth's equator)
		mean_nu : {float}, optional
			mean anomaly in degrees (the default is 0, which represents an orbit
			that is beginning at periapsis)
		
		Returns
		-------
		satplot.Orbit
		"""

		if body.upper() == 'EARTH':
			central_body = Earth
			min_a = consts.R_EARTH + consts.EARTH_MIN_ALT
		elif body.upper() == 'SUN':
			central_body = Sun
			min_a = consts.R_SUN + consts.SUN_MIN_ALT
		elif body.upper() == 'MARS':
			central_body = Mars
			min_a = consts.R_MARS + consts.MARS_MIN_ALT
		elif body.upper() == 'MOON':
			central_body = Moon
			min_a = consts.R_MOON + consts.MOON_MIN_ALT
		else:
			logger.error("Invalid central body {}. Valid options are Earth, Sun, Mars, Moon".format(body))
			raise ValueError("Invalid central body {}.".format(body))

		if a < min_a:
			logger.error("Semimajor axis, {}, is too close to the central body, {}".format(a, body.upper()))
			raise exceptions.OutOfRange("Semimajor axis, {}, is too close to the central body, {}".format(a, body.upper()))
		if ecc > 1 or ecc < 0:
			logger.error("Eccentricity, {}, is non circular or eliptical".format(ecc))
			raise exceptions.OutOfRange("Eccentricity, {}, is non circular or eliptical".format(ecc))
		if inc > 180 or inc < -180:
			logger.error("Inclination, {}, is out of range, should be -180 < inc < 180".format(inc))
			raise exceptions.OutOfRange("Inclination, {}, is out of range, should be -180 < inc < 180".format(inc))
		if raan > 360 or raan < 0:
			logger.error("RAAN, {}, is out of range, should be 0 < inc < 360".format(inc))
			raise exceptions.OutOfRange("RAAN, {}, is out of range, should be 0 < inc < 360".format(inc))
		if argp > 360 or argp < 0:
			logger.error("Argument of periapsis, {}, is out of range, should be 0 < argp < 360".format(inc))
			raise exceptions.OutOfRange("Argument of periapsis, {}, is out of range, should be 0 < argp < 360".format(inc))
		if mean_nu > 360 or mean_nu < 0:
			logger.error("Mean anomaly, {}, is out of range, should be 0 < mean_nu < 360".format(inc))
			raise exceptions.OutOfRange("Mean anomaly, {}, is out of range, should be 0 < mean_nu < 360".format(inc))

		return cls(timespan, body=central_body, a=a, ecc=ecc, inc=inc, raan=raan, argp=argp, mean_nu=mean_nu, type='ANALYTICAL', astrobodies=astrobodies)

	@classmethod
	def load(cls, file):
		with open(file, 'rb') as fp:
			logger.info(f'Loading orbit file from {file}')
			return pickle.load(fp)


	def getPosition(self, time):
		'''Return the position at the specified index or (closest) time
				
		Parameters
		----------
		time : {int, datetime.datetime, astropy.Time}
			The index of the timespan, or a particular time to fetch the position
			If a particular time is specified, the nearest timestep within the timespan will be returned
		
		Returns
		-------
		(3,) ndarray
			Orbital position (km)
		
		Raises
		------
		ValueError
			Raises a ValueError if 'time' is not an integer or a time type
		'''
		if isinstance(time, int):
			return self.pos[time, :]
		elif isinstance(time, dt.datetime):
			closest_time, index = list_u.get_closest(self.timespan.asDatetime(), time)
			if abs((closest_time - time).total_seconds()) > 60:
				logger.warning("{}, is more than 60s from the nearest timestep, choose a smaller timestep".format(time))
			return self.pos[index, :]
		elif isinstance(time, astropyTime):
			closest_time, index = list_u.get_closest(self.timespan.asDatetime(), time.to_datetime())
			if abs((closest_time - time.to_datetime()).total_seconds()) > 60:
				logger.warning("{}, is more than 60s from the nearest timestep, choose a smaller timestep".format(time))
			return self.pos[index, :]
		else:
			logger.error("{} cannot be used to index an orbital position".format(time))
			raise ValueError("{} cannot be used to index an orbital position".format(time))

	def getVelocity(self, time):
		'''Return the position at the specified index or (closest) time
				
		Parameters
		----------
		time : {int, datetime.datetime, astropy.Time}
			The index of the timespan, or a particular time to fetch the velocity
			If a particular time is specified, the nearest timestep within the timespan will be returned
		
		Returns
		-------
		(3,) ndarray
			Orbital velocity (m/s)
		
		Raises
		------
		ValueError
			Raises a ValueError if 'time' is not an integer or a time type
		'''
		if isinstance(time, int):
			return self.vel[time, :]
		elif isinstance(time, dt.datetime):
			closest_time, index = list_u.get_closest(self.timespan.asDatetime(), time)
			if abs((closest_time - time).total_seconds()) > 60:
				logger.warning("{}, is more than 60s from the nearest timestep, choose a smaller timestep".format(time))
			return self.vel[index, :]
		elif isinstance(time, astropyTime):
			closest_time, index = list_u.get_closest(self.timespan.asDatetime(), time.to_datetime())
			if abs((closest_time - time.to_datetime()).total_seconds()) > 60:
				logger.warning("{}, is more than 60s from the nearest timestep, choose a smaller timestep".format(time))

			return self.vel[index, :]
		else:
			logger.error("{} cannot be used to index an orbital velocity".format(time))
			raise ValueError("{} cannot be used to index an orbital velocity".format(time))

	# def save(self):
	# 	if self.gen_type == 'TLE':
	# 		file_name = f'orbit_TLE_{self.sat.model.satnum}_{self.sat.jdsatepoch}.pickle'

	# 	elif self.gen_type == 'ANALYTICAL':			
	# 		start = self.timespan.start.strftime('%s')
	# 		step = self.timespan.time_step.days * 86400 + self.timespan.time_step.seconds
	# 		period = self.timespan.time_period.days * 86400 + self.timespan.time_period.seconds
	# 		file_name = f'orbit_ANLT_{start}_{step}_{period}_{self.orb.epoch}_{self.orb.a}_{self.orb.ecc}_{self.orb.inc}_{self.orb.raan}_{self.orb.argp}_{self.orb.nu}.pickle'

	# 	elif self.gen_type == 'FAKE_TLE':
	# 		file_name = f'orbit_FAKETLE_{self.sat.model.satnum}_{self.sat.jdsatepoch}_{self.sat.ecco}_{self.sat.inclo}_{self.sat.nodeo}_{self.sat.argpo}_{self.sat.mo}.pickle'
		
	# 	elif self.gen_type == 'FILE':
	# 		file_name = f'orbit_genfrom_{self.filename}.pickle'

	# 	elif self.gen_type == 'POS_LIST':
	# 		file_name = 'pos_list'

	# 	with open(f'data/orbits/{file_name}', 'wb') as fp:
	# 		pickle.dump(self, fp)

	# 	return file_name
