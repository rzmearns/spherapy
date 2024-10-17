
import numpy as np
import sys
import logging
from progressbar import progressbar
import pickle

from astropy import units as astropy_units

from hapsira.bodies import Earth, Mars, Sun, Moon
from hapsira.twobody import Orbit as hapsiraOrbit
from hapsira.ephem import Ephem

from sgp4.api import Satrec, WGS72
from astropy.time import Time as astropyTime
import datetime as dt
# Don't need to import all of skyfield just for EarthSatellite loading.
# TODO: figure out the lightest import to use
from skyfield.api import load, EarthSatellite, wgs84

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
		# If no gen_type, spit an error here
		raise_exceptions = kwargs.get('raise_exceptions')
		gen_type = kwargs.get('type')
		if gen_type is None:
			logger.error("Orbit() should not be called directly.")
			raise ValueError("Orbit() should not be called directly.")

		if len(args) == 0:
			logger.error("Orbit() should not be called directly.")
			raise ValueError("Orbit() should not be called directly.")			

		calc_astrobodies = kwargs.get('astrobodies')

		timespan = args[0]
		
		# CONSTRUCTOR HANDLING
		if gen_type == 'TLE':
			# List of skyfield EarthSatellites, one for each TLE in the file
			skyfld_earthsats = args[1]
			tle_dates = [sat.epoch.utc_datetime() for sat in skyfld_earthsats]
			
			# check timespan is valid for provided TLEs
			if timespan.start < tle_dates[0] - dt.timedelta(days=14):
				logger.error("Timespan begins before provided TLEs (+14 days)")
				print("Timespan begins before provided TLEs (+14 days)", file=sys.stderr)
				if raise_exceptions:
					raise exceptions.OutOfRange("Timespan begins before provided TLEs (+14 days)")
			elif timespan.start > tle_dates[-1] + dt.timedelta(days=14):
				logger.error("Timespan begins after provided TLEs (+14 days)")
				print("Timespan begins after provided TLEs (+14 days)", file=sys.stderr)
				if raise_exceptions:
					raise exceptions.OutOfRange("Timespan begins after provided TLEs (+14 days)")
			elif timespan.end > tle_dates[-1] + dt.timedelta(days=14):
				logger.error("Timespan ends after provided TLEs (+14 days)")
				print("Timespan ends after provided TLEs (+14 days)", file=sys.stderr)
				if raise_exceptions:
					raise exceptions.OutOfRange("Timespan ends after provided TLEs (+14 days)")

			data_dict = self._propagateFromTLE(timespan, tle_dates, skyfld_earthsats, calc_astrobodies)
		

		elif gen_type == 'FAKE_TLE':
			a = kwargs.get('a')
			ecc = kwargs.get('ecc')
			inc = kwargs.get('inc')
			raan = kwargs.get('raan')
			argp = kwargs.get('argp')
			mean_nu = kwargs.get('mean_nu')

			data_dict = self._propagateAnalytical(timespan, a, ecc, inc, raan, argp, mean_nu)
			data_dict['name'] = kwargs['name']


		elif gen_type == 'ANALYTICAL':
			body = kwargs.get('body')
			a = kwargs.get('a') * astropy_units.km
			ecc = kwargs.get('ecc') * astropy_units.one
			inc = kwargs.get('inc') * astropy_units.rad
			raan = kwargs.get('raan') * astropy_units.rad
			argp = kwargs.get('argp') * astropy_units.rad
			mean_nu = kwargs.get('mean_nu') * astropy_units.rad
						
			data_dict = self._genAnalytical(timespan, a, ecc, inc, raan, argp, mean_nu)
			data_dict['name'] = kwargs['name']


		elif gen_type == 'POS_LIST':
			data = self._dfltDataDict()
			
			data['timespan'] = timespan
			data['pos'] = args[1]
			# Assume linear motion between each position at each timestep; 
			# Then assume it stops at the last timestep.
			vel = self.pos[1:] - self.pos[:-1]
			data['vel'] = np.concatenate((self.vel, np.array([[0, 0, 0]])))
			data['period'] = None
			data['period_steps'] = 1
			data['name'] = ['Sat from position list']
			# Note that this doesn't define a period.
			logger.warning("Warning: When generating satellite orbit from list of positions, `period` and `period_steps` will not be defined.")			
			data['gen_type'] = 'position list'

		else: 
			logger.error("Invalid orbit generation option {}. Valid options are TLE, FAKE_TLE, POS_LIST, and ANALYITCAL.".format(gen_type))
			raise ValueError("Invalid orbit generation option {}.".format(gen_type))


		if calc_astrobodies:
			logger.info('Creating ephemeris for Sun using timespan')
			# Timescale for sun position calculation should use TDB, not UTC
			# The resultant difference is likely very small
			ephem_sun = Ephem.from_body(Sun, astropyTime(timespan.asAstropy(scale='tdb')), attractor=Earth)
			data_dict['sun_pos'] = np.asarray(ephem_sun.rv()[0].to(astropy_units.km))

			logger.info('Creating ephemeris for Moon using timespan')
			ephem_moon = Ephem.from_body(Moon, astropyTime(timespan.asAstropy(scale='tdb')), attractor=Earth)
			data_dict['moon_pos'] = np.asarray(ephem_moon.rv()[0].to(astropy_units.km))

		self._attributiseDataDict(data_dict)

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
		skyfld_earth_sats = load.tle_file(tle_path)
		return cls(timespan, skyfld_earth_sats, type='TLE', astrobodies=astrobodies)

	@classmethod
	def fromPropagatedOrbitalParam(cls, timespan, a=6978, ecc=0, inc=0, raan=0, argp=0, mean_nu=0, name='Fake TLE', astrobodies=True):		
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


		return cls(timespan, body=Earth, a=a, ecc=ecc, inc=inc, raan=raan, argp=argp, mean_nu=mean_nu, type='FAKE_TLE', astrobodies=astrobodies)
		
	@classmethod
	def fromAnalyticalOrbitalParam(cls, timespan, body='Earth', a=6978, ecc=0, inc=0, raan=0, argp=0, mean_nu=0, name='Analytical', astrobodies=True):
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
			closest_time, index = list_u.get_closest(self.timespan[:], time)
			if abs((closest_time - time).total_seconds()) > 60:
				logger.warning("{}, is more than 60s from the nearest timestep, choose a smaller timestep".format(time))
			return self.pos[index, :]
		elif isinstance(time, astropyTime):
			closest_time, index = list_u.get_closest(self.timespan[:], time.to_datetime())
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
			closest_time, index = list_u.get_closest(self.timespan[:], time)
			if abs((closest_time - time).total_seconds()) > 60:
				logger.warning("{}, is more than 60s from the nearest timestep, choose a smaller timestep".format(time))
			return self.vel[index, :]
		elif isinstance(time, astropyTime):
			closest_time, index = list_u.get_closest(self.timespan[:], time.to_datetime())
			if abs((closest_time - time.to_datetime()).total_seconds()) > 60:
				logger.warning("{}, is more than 60s from the nearest timestep, choose a smaller timestep".format(time))

			return self.vel[index, :]
		else:
			logger.error("{} cannot be used to index an orbital velocity".format(time))
			raise ValueError("{} cannot be used to index an orbital velocity".format(time))

	def _dfltDataDict(self):
		d = {}
		# metadata
		d['name'] = None 		#str
		d['satcat_id'] = None 	#int
		d['gen_type'] = None 	#str

		# time data
		d['timespan'] = None	#spherapy timespan
		d['TLE_epochs'] = None 	#ndarray(n,) datetimes

		# pos data
		d['pos'] = None 		#ndarray(3,n)
		d['vel'] = None 		#ndarray(3,n)
		d['lat'] = None 		#ndarray(n,)
		d['lon'] = None 		#ndarray(n,)
		d['sun_pos'] = None 	#ndarray(3,n)
		d['moon_pos'] = None 	#ndarray(3,n)

		# orbit_data
		d['central_body'] = None#str
		d['period'] = None 		#float
		d['period_steps'] = None#int
		d['semi_major'] = None 	#ndarray(n,) [km]
		d['ecc'] = None 		#ndarray(n,) [radians]
		d['inc'] = None 		#ndarray(n,) [radians]
		d['raan'] = None 		#ndarray(n,) [radians]
		d['argp'] = None 		#ndarray(n,) [radians]

		return d

	def _attributiseDataDict(self, d):
		for k,v in d.items():
			setattr(self,k,v)

	def _propagateFromTLE(self, timespan, tle_dates, skyfld_earthsats, gen_astrobodies=True) -> dict:
			
			data = self._dfltDataDict()
			skyfld_ts = load.timescale(builtin=True)

			# Find closest listed TLE to start date
			start_datetime, start_index = list_u.get_closest(tle_dates, timespan.start)
			end_datetime, end_index = list_u.get_closest(tle_dates, timespan.end)

			# timespan starts after last TLE epoch
			if start_index == -1 or start_index == len(tle_dates):
				tspan_skyfld_earthsats = [skyfld_earthsats[-1]]
				tspan_tle_epochs = [tle_dates[-1]]

			# timespan starts before last TLE epoch, but ends after last TLE epoch
			elif end_index == -1 or end_index == len(tle_dates):
				tspan_skyfld_earthsats = skyfld_earthsats[start_index:]
				tspan_tle_epochs = tle_dates[start_index:]
			
			# timespan lies completely within available TLE epochs
			else:
				tspan_skyfld_earthsats = skyfld_earthsats[start_index:end_index+1]
				tspan_tle_epochs = tle_dates[start_index:end_index+1]

			timesteps = timespan[:]
			
			# tspan stretches across more than 1 TLE
			# 	-> calculate median date for each TLE pair
			# TODO: should we be back propagating for a TLE?
			# 	-> find closest tspan indices to these median dates
			# 	-> perform propagation using new TLE when required
			if len(tspan_tle_epochs) > 1:
				tspan_epoch_middates = np.diff(tspan_tle_epochs)/2 + np.asarray(tspan_tle_epochs[:-1])
							
				for ii, date in enumerate(tspan_epoch_middates):
					tspan_epoch_middates[ii] = date.replace(tzinfo=timespan.timezone)

				trans_indices = []
				for date in tspan_epoch_middates:
					trans_indices.append(np.argmin(np.abs(np.vectorize(dt.timedelta.total_seconds)(timesteps-date))))

				# skyfield doesn't appear to have a way to concatenate 'Geocentric' objecs
				sat_rec = tspan_skyfld_earthsats[0].at(skyfld_ts.utc(timesteps[:trans_indices[0]]))
				pos = sat_rec.position.km.T
				vel = sat_rec.velocity.km_per_s.T * 1000
				l, l2 = wgs84.latlon_of(sat_rec)
				lat = l.degrees
				lon = l2.degrees
				ecc = np.tile(sat_rec.model.ecco,trans_indices[0])
				inc = np.tile(sat_rec.model.inclo,trans_indices[0])
				semi_major = np.tile(sat_rec.model.a * consts.R_EARTH,trans_indices[0])
				raan = np.tile(sat_rec.model.nodeo,trans_indices[0])
				argp = np.tile(sat_rec.model.argpo,trans_indices[0])
				TLE_epochs = np.tile(tspan_tle_epochs[0],trans_indices[0])



				for ii in range(len(trans_indices)-1):
					sat_rec = tspan_skyfld_earthsats[ii+1].at(skyfld_ts.utc(timesteps[trans_indices[ii]:trans_indices[ii+1]]))
					pos = np.vstack((pos, sat_rec.position.km.T))
					vel = np.vstack((vel, sat_rec.velocity.km_per_s.T * 1000))
					l, l2 = wgs84.latlon_of(sat_rec)
					lat = np.concatenate((lat,l.degrees))
					lon = np.concatenate((lon,l2.degrees))
					ecc = np.concatenate((ecc,
										np.tile(sat_rec.model.ecco,trans_indices[ii+1])))
					inc = np.concatenate((inc,
										np.tile(sat_rec.model.inclo,trans_indices[ii+1])))
					semi_major = np.concatenate((semi_major,
										np.tile(sat_rec.model.a * consts.R_EARTH,trans_indices[ii+1])))
					raan = np.concatenate((raan,
										np.tile(sat_rec.model.nodeo,trans_indices[ii+1])))
					argp = np.concatenate((argp,
										np.tile(sat_rec.model.argpo,trans_indices[ii+1])))

					TLE_epochs = np.concatenate((TLE_epochs,
								  				np.tile(tspan_tle_epochs[ii+1],trans_indices[ii+1]-trans_indices[ii])))



				sat_rec = tspan_skyfld_earthsats[-1].at(skyfld_ts.utc(timesteps[trans_indices[-1]:]))
				pos = np.vstack((pos, sat_rec.position.km.T))
				vel = np.vstack((vel, sat_rec.velocity.km_per_s.T * 1000))
				l, l2 = wgs84.latlon_of(sat_rec)
				lat = np.concatenate((lat,l.degrees))
				lon = np.concatenate((lon,l2.degrees))
				ecc = np.concatenate((ecc,
									np.tile(sat_rec.model.ecco,trans_indices[-1])))
				inc = np.concatenate((inc,
									np.tile(sat_rec.model.inclo,trans_indices[-1])))
				semi_major = np.concatenate((semi_major,
									np.tile(sat_rec.model.a * consts.R_EARTH,trans_indices[-1])))
				raan = np.concatenate((raan,
									np.tile(sat_rec.model.nodeo,trans_indices[-1])))
				argp = np.concatenate((argp,
									np.tile(sat_rec.model.argpo,trans_indices[-1])))
				TLE_epochs = np.concatenate((TLE_epochs, 
									  			np.tile(tspan_tle_epochs[ii+1],len(timesteps) - trans_indices[-1])))


			# tspan is only single TLE
			else:
				sat_rec = tspan_skyfld_earthsats[0].at(skyfld_ts.utc(timesteps))
				pos = sat_rec.position.km.T
				vel = sat_rec.velocity.km_per_s.T * 1000
				l, l2 = wgs84.latlon_of(sat_rec)
				lat = l.degrees
				lon = l2.degrees
				ecc = np.tile(sat_rec.model.ecco, len(timesteps))
				inc = np.tile(sat_rec.model.inclo, len(timesteps))
				semi_major = np.tile(sat_rec.model.a * consts.R_EARTH, len(timesteps))
				raan = np.tile(sat_rec.model.nodeo, len(timesteps))
				argp = np.tile(sat_rec.model.argpo, len(timesteps))
				TLE_epochs = np.tile(tspan_tle_epochs[0],len(timesteps))

			data['timespan'] = timespan
			data['gen_type'] = 'propagated from TLE'
			data['central_body'] = 'Earth'
			data['pos'] = pos
			data['vel'] = vel
			data['lat'] = lat
			data['lon'] = lon
			data['ecc'] = ecc
			data['inc'] = inc
			data['semi_major'] = semi_major
			data['raan'] = raan
			data['argp'] = argp
			data['TLE_epochs'] = TLE_epochs
			data['period'] = 2 * np.pi / tspan_skyfld_earthsats[-1].model.no_kozai * 60
			data['period_steps'] = int(data['period'] / timespan.time_step.total_seconds())
			data['name'] = tspan_skyfld_earthsats[-1].name

			try:
				self.sat = tspan_skyfld_earthsats[-1]
			except:
				# TODO: figure out what the exact error is and cause of when tspan_skyfld_earthsats doesn't exist
				pass

			return data

	def _genAnalytical(self, timespan, body, a, ecc, inc, raan, argp, mean_nu):
		
		data = self._dfltDataDict()
		logger.info("Creating analytical orbit")
		orb = hapsiraOrbit.from_classical(body, a, ecc, inc, raan, argp, mean_nu)		
		logger.info("Creating ephemeris for orbit, using timespan")
		ephem = Ephem.from_orbit(orb, timespan.asAstropy())

		pos = np.asarray(ephem.rv()[0])
		vel = np.asarray(ephem.rv()[1]) * 1000

		period = orb.period.unit.in_units('s') * orb.period.value
		period_steps = int(period / timespan.time_step.total_seconds())

		data['timespan'] = timespan
		data['gen_type'] = 'propagated from TLE'
		data['central_body'] = body
		data['pos'] = pos
		data['vel'] = vel
		data['lat'] = None
		data['lon'] = None
		data['ecc'] = ecc.value
		data['inc'] = inc.value
		data['semi_major'] = a.value
		data['raan'] = raan.value
		data['argp'] = argp.value
		data['TLE_epochs'] = None
		data['period'] = period
		data['period_steps'] = period_steps

		return data

	def _propagateAnalytical(self, timespan, body, a, ecc, inc, raan, argp, mean_nu):
		
		data = self._dfltDataDict()
		skyfld_ts = load.timescale(builtin=True)

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
						argp,  # arg perigee [radians] # noqa: E128
						inc,  # inclination [radians] # noqa: E128
						mean_nu,  # mean anomaly [radians] # noqa: E128
						mean_motion * 60,  # mean motion [rad/min]  # noqa: E128
						raan)  # raan [radians] # noqa: E126, E128

		skyfld_earthsat = EarthSatellite.from_satrec(satrec, skyfld_ts)

		pos = np.empty((timespan.num_steps, 3))
		vel = np.empty((timespan.num_steps, 3))

		for ii, timestep in enumerate(timespan[:]):	
			pos[ii, :] = skyfld_earthsat.at(skyfld_ts.utc(timestep.replace(tzinfo=timespan.timezone))).position.km 
			vel[ii, :] = skyfld_earthsat.at(skyfld_ts.utc(timestep.replace(tzinfo=timespan.timezone))).velocity.km_per_s * 1000

		period = 2 * np.pi / skyfld_earthsat.model.no_kozai * 60
		period_steps = int(period / timespan.time_step.total_seconds())
			
		data['timespan'] = timespan
		data['gen_type'] = 'propagated from orbital param'
		data['central_body'] = 'Earth'
		data['pos'] = pos
		data['vel'] = vel
		data['lat'] = None
		data['lon'] = None
		data['ecc'] = ecc
		data['inc'] = inc
		data['semi_major'] = a
		data['raan'] = raan
		data['argp'] = argp
		data['TLE_epochs'] = t0_epoch
		data['period'] = period
		data['period_steps'] = period_steps

		return data

	def serialise(self):
		raise NotImplementedError