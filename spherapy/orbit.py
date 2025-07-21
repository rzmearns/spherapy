
import numpy as np
import sys
import logging
from progressbar import progressbar
import pathlib
import pickle

from astropy import units as astropy_units

from hapsira.bodies import Earth, Mars, Sun, Moon
from hapsira.twobody import Orbit as hapsiraOrbit
from hapsira.ephem import Ephem

from sgp4.api import Satrec, WGS72
from skyfield.framelib import itrs

from astropy.time import Time as astropyTime
import datetime as dt
# Don't need to import all of skyfield just for EarthSatellite loading.
# TODO: figure out the lightest import to use
from skyfield.api import load, EarthSatellite, wgs84

from spherapy.timespan import TimeSpan

import spherapy.util.list_u as list_u
import spherapy.util.epoch_u as epoch_u
import spherapy.util.orbital_u as orbit_u
import spherapy.util.exceptions as exceptions
import spherapy.util.constants as consts

import typing
from typing import Any, TypedDict

import types

logger = logging.getLogger(__name__)


class OrbitAttrDict(TypedDict):
	name: None|str
	satcat_id: None|int
	gen_type: None|str
	timespan: None|TimeSpan
	TLE_epochs: None|np.ndarray[tuple[int], np.dtype[dt.datetime]]
	pos: None|np.ndarray[tuple[int,int], np.dtype[np.float_]]
	pos_ecef: None|np.ndarray[tuple[int,int], np.dtype[np.float_]]
	vel_ecef: None|np.ndarray[tuple[int,int], np.dtype[np.float_]]
	vel: None|np.ndarray[tuple[int,int], np.dtype[np.float_]]
	lat: None|np.ndarray[tuple[int], np.dtype[np.float_]]
	lon: None|np.ndarray[tuple[int], np.dtype[np.float_]]
	sun_pos: None|np.ndarray[tuple[int,int], np.dtype[np.float_]]
	moon_pos: None|np.ndarray[tuple[int,int], np.dtype[np.float_]]
	alt: None|np.ndarray[tuple[int], np.dtype[np.float_]]
	eclipse: None|np.ndarray[tuple[int], np.dtype[np.bool_]]
	central_body: None|str
	period: None|float
	period_steps: None|int
	semi_major: None|np.ndarray[tuple[int], np.dtype[np.float_]]
	ecc: None|np.ndarray[tuple[int], np.dtype[np.float_]]
	inc: None|np.ndarray[tuple[int], np.dtype[np.float_]]
	raan: None|np.ndarray[tuple[int], np.dtype[np.float_]]
	argp: None|np.ndarray[tuple[int], np.dtype[np.float_]]

def createEmptyOrbitAttrDict() -> OrbitAttrDict:
	return OrbitAttrDict({
	'name': None,
	'satcat_id':None,
	'gen_type':None,
	'timespan':None,
	'TLE_epochs':None,
	'pos':None,
	'pos_ecef':None,
	'vel_ecef':None,
	'vel':None,
	'lat':None,
	'lon':None,
	'sun_pos':None,
	'moon_pos':None,
	'alt':None,
	'eclipse':None,
	'central_body':None,
	'period':None,
	'period_steps':None,
	'semi_major':None,
	'ecc':None,
	'inc':None,
	'raan':None,
	'argp':None,
	})

def validateTypedDict(obj:Any, typ: Any) -> bool:
	for attr_name, spec_attr_type in typ.__annotations__.items():
		try:
			attr_value = obj.get(attr_name, None)
		except KeyError:
			# Check for missing keys
			raise KeyError('Typed Dict is missing key: %s', attr_name)
		if type(spec_attr_type) is types.UnionType:
			# check union of types
			good_type = False
			for sub_type in typing.get_args(spec_attr_type):
				# choose validator
				if sub_type in (None, types.NoneType):
					validator = NoneValidator
				else:
					validator = genericValidator
				if validator(attr_value, sub_type):
					good_type = True
					break
			if not good_type:
				raise TypeError(f'{attr_name} is type {type(attr_value)}, should be {spec_attr_type}')
		elif not genericValidator(attr_value, spec_attr_type):
			# check all other types
			raise TypeError(f'{attr_name} is type {type(attr_value)}, should be {spec_attr_type}')
			return False
	return True

def NoneValidator(obj:Any, check_type:Any) -> bool:
	if type(obj) is types.NoneType:
		return True
	return False

def genericValidator(obj:Any, check_type:Any) -> bool:
	if type(check_type) is types.GenericAlias:
		if isinstance(obj, check_type.__origin__):
			return True
	else:
		if isinstance(obj, check_type):
			return True
	return False

class Orbit(object):
	'''
	Contains timestepped array of orbital position and velocity data for a satellite, coordinate system depending on the central body.
	Contains timestepped array of sun position.
	Timing data taken from a TimeSpan object.

	Attributes
	----------
	timespan: TimeSpan
		Timespan over which orbit is to be simulated

	pos: (N, 3) np.array
		Cartesian coordinates of the position of the satellite at each time in a TimeSpan (units: km)
		Coordinate system depending on the central body.

	vel: (N, 3) np.array
		Cartesian coordinates of the velocity of the satellite at each time in a TimeSpan (units: m/s)
		Coordinate system depending on the central body.

	sun_pos: (N, 3) np.array
		Vector of the sun's position at each timestep (GCRS; units: km)

	period: float
		Period of the satellite's orbit (units: s)

	period_steps: float
		Number of timesteps per orbital period of the satellite.
		Useful to check the timestep is appropriate for this orbit.
		Recommended: 10 < period_steps < 100 or so? Should we write warnings if outside this range?

	'''

	def __init__(self, data:OrbitAttrDict, calc_astrobodies:bool=False):
		'''
		The constructor should never be called directly.
		Use one of:
			Orbit.fromTLE()
			Orbit.fromListOfPositions()
			Orbit.fromPropagatedOrbitalParam()
			Orbit.fromAnalyticalOrbitalParam()
		'''
		# Should always be called from a class method, however,
		# If no gen_type, spit an error here

		if not validateTypedDict(data, OrbitAttrDict):
			logger.error("Orbit() should not be called directly, use one of the fromXX constructors.")
			raise ValueError("Orbit() should not be called directly, use one of the fromXX constructors")

		self.name = data['name']
		self.satcat_id = data['satcat_id']
		self.gen_type = data['gen_type']

		#time data
		self.timespan = data['timespan']
		self.TLE_epochs = data['TLE_epochs']

		#pos data
		self.pos = data['pos']
		self.pos_ecef = data['pos_ecef']
		self.vel_ecef = data['vel_ecef']
		self.vel = data['vel']
		self.lat = data['lat']
		self.lon = data['lon']
		self.alt = data['alt']

		# These should be None from the constructor
		self.sun_pos = data['sun_pos']
		self.moon_pos = data['moon_pos']
		self.eclipse = data['eclipse']

		#orbit_data
		self.central_body = data['central_body']
		self.period = data['period']
		self.period_steps = data['period_steps']
		self.semi_major = data['semi_major']
		self.ecc = data['ecc']
		self.inc = data['inc']
		self.raan = data['raan']
		self.argp = data['argp']

		# Check required fields are not empty
		if self.timespan is None:
			raise AttributeError('Error in creating Orbit object, no Timespan assigned')
		if self.pos is None:
			raise AttributeError('Error in creating Orbit object, no position data assigned')
		if self.vel is None:
			raise AttributeError('Error in creating Orbit object, no velocity data assigned')
		if self.name is None:
			raise AttributeError('Error in creating Orbit object, no name assigned')
		if self.gen_type is None:
			raise AttributeError('Error in creating Orbit object, no generator type assigned')

		if calc_astrobodies:
			logger.info('Creating ephemeris for Sun using timespan')
			# Timescale for sun position calculation should use TDB, not UTC
			# The resultant difference is likely very small
			ephem_sun = Ephem.from_body(Sun, astropyTime(self.timespan.asAstropy(scale='tdb')), attractor=Earth)
			sun_pos = ephem_sun.rv()[0]
			self.sun_pos = np.asarray(sun_pos.to(astropy_units.km))

			logger.info('Creating ephemeris for Moon using timespan')
			ephem_moon = Ephem.from_body(Moon, astropyTime(self.timespan.asAstropy(scale='tdb')), attractor=Earth)
			moon_pos = ephem_moon.rv()[0]
			self.moon_pos = np.asarray(moon_pos.to(astropy_units.km))
			self.eclipse = self._calcEclipse(self.pos, self.sun_pos)

	#pos_list
	@classmethod
	def fromListOfPositions(cls, timespan:TimeSpan, positions:np.ndarray[tuple[int, int], np.dtype[np.float_]], astrobodies:bool=True):
		"""Create an orbit by explicitly specifying the position of the
		satellite at each point in time. Useful for simplified test cases; but
		may lead to unphysical orbits.

		Parameters
		----------
		timespan : TimeSpan
			Timespan over which orbit is to be simulated

		positions: (N,3) np.array
			Position of the satellite in GCRS at each point in time.

		Returns
		-------
		satplot.Orbit
		"""
		attr_dct = createEmptyOrbitAttrDict()
		if len(positions) != len(timespan):
			raise ValueError(f'Number of supplied positions does not match timespan length: {len(positions)} =/= {len(timespan)}')

		attr_dct['timespan'] = timespan
		attr_dct['pos'] = positions
		# Assume linear motion between each position at each timestep;
		# Then assume it stops at the last timestep.
		vel = attr_dct['pos'][1:] - attr_dct['pos'][:-1]
		attr_dct['vel'] = np.concatenate((vel, np.array([[0, 0, 0]])))
		attr_dct['period'] = None
		attr_dct['period_steps'] = 1
		attr_dct['name'] = 'Sat from position list'
		attr_dct['gen_type'] = 'position list'

		return cls(attr_dct, calc_astrobodies=astrobodies)

	#TLE
	@classmethod
	def fromTLE(cls, timespan:TimeSpan, tle_path:pathlib.Path, astrobodies:bool=True, unsafe=False):
		"""Create an orbit from an existing TLE or a list of historical TLEs

		Parameters
		----------
		timespan : TimeSpan over which orbit is to be simulated
		tle_path : path to file containing TLE(s)

		Returns
		-------
		satplot.Orbit
		"""
		skyfld_earth_sats = load.tle_file(f'{tle_path}')
		epochs = np.asarray([a.epoch.utc_datetime() for a in skyfld_earth_sats])
		unq_epoch, unq_idxs = np.unique(epochs, return_index=True)
		unq_skyfld_earth_sats = list(np.asarray(skyfld_earth_sats)[unq_idxs])

		# return cls(timespan, unq_skyfld_earth_sats, type='TLE', astrobodies=astrobodies)

		tle_dates = [sat.epoch.utc_datetime() for sat in unq_skyfld_earth_sats]

		# check timespan is valid for provided TLEs
		if timespan.start < tle_dates[0] - dt.timedelta(days=14):
			logger.error("Timespan begins before provided TLEs (+14 days)")
			if not unsafe:
				raise exceptions.OutOfRange("Timespan begins before provided TLEs (+14 days)")
		elif timespan.start > tle_dates[-1] + dt.timedelta(days=14):
			logger.error("Timespan begins after provided TLEs (+14 days)")
			if not unsafe:
				raise exceptions.OutOfRange("Timespan begins after provided TLEs (+14 days)")
		elif timespan.end > tle_dates[-1] + dt.timedelta(days=14):
			logger.error("Timespan ends after provided TLEs (+14 days)")
			if not unsafe:
				raise exceptions.OutOfRange("Timespan ends after provided TLEs (+14 days)")

		attr_dct = createEmptyOrbitAttrDict()

		closest_tle_epochs = _findClosestEpochIndices(np.asarray(tle_dates), timespan[:])

		d = np.hstack((1, np.diff(closest_tle_epochs)))
		timespan_epoch_trans_idxs = np.hstack((np.where(d!=0)[0],len(timespan))) # timestep idx where TLE epoch changes, append len(timespan) to not require separate loop iteration to handle last element
		timespan_epoch_trans_tle_idxs = closest_tle_epochs[np.where(d!=0)[0]] # tle_dates idx where TLE epoch changes in timespan

		tle_epoch_idxs = []
		sub_timespans = []
		for ii, trans_start_idx in enumerate(timespan_epoch_trans_idxs[:-1]):
			tle_epoch_idxs.append(timespan_epoch_trans_tle_idxs[ii])
			sub_timespans.append(timespan[trans_start_idx,timespan_epoch_trans_idxs[ii+1]])

		skyfld_ts = load.timescale(builtin=True)

		# Calculate data for first run
		sub_timespan = sub_timespans[0]
		sub_skyfld_earthsat = unq_skyfld_earth_sats[tle_epoch_idxs[0]]
		tle_epoch = tle_dates[tle_epoch_idxs[0]]
		sat_rec = sub_skyfld_earthsat.at(skyfld_ts.utc(sub_timespan))
		pos = sat_rec.position.km.T
		vel = sat_rec.velocity.km_per_s.T * 1000
		ecef = sat_rec.frame_xyz_and_velocity(itrs)
		pos_ecef = ecef[0].km.T
		vel_ecef = ecef[1].km_per_s.T * 1000

		l, l2 = wgs84.latlon_of(sat_rec)
		lat = l.degrees
		lon = l2.degrees
		ecc = np.tile(sub_skyfld_earthsat.model.ecco,len(sub_timespan))
		inc = np.tile(sub_skyfld_earthsat.model.inclo,len(sub_timespan))
		semi_major = np.tile(sub_skyfld_earthsat.model.a * consts.R_EARTH,len(sub_timespan))
		raan = np.tile(sub_skyfld_earthsat.model.nodeo,len(sub_timespan))
		argp = np.tile(sub_skyfld_earthsat.model.argpo,len(sub_timespan))
		TLE_epochs = np.tile(tle_epoch,len(sub_timespan))


		# fill in data for any other sub timespans
		for ii in range(1,len(sub_timespans)):
			sub_timespan = sub_timespans[ii]
			sub_skyfld_earthsat = unq_skyfld_earth_sats[tle_epoch_idxs[ii]]
			tle_epoch = tle_dates[tle_epoch_idxs[ii]]
			sat_rec = sub_skyfld_earthsat.at(skyfld_ts.utc(sub_timespan))
			pos = np.vstack((pos, sat_rec.position.km.T))
			vel = np.vstack((vel, sat_rec.velocity.km_per_s.T * 1000))
			ecef = sat_rec.frame_xyz_and_velocity(itrs)
			pos_ecef = np.vstack((pos_ecef, ecef[0].km.T))
			vel_ecef = np.vstack((vel_ecef, ecef[1].km_per_s.T * 1000))

			l, l2 = wgs84.latlon_of(sat_rec)
			lat = np.concatenate((lat, l.degrees))
			lon = np.concatenate((lon, l2.degrees))
			ecc = np.concatenate((ecc, np.tile(sub_skyfld_earthsat.model.ecco, len(sub_timespan))))
			inc = np.concatenate(((inc, np.tile(sub_skyfld_earthsat.model.inclo, len(sub_timespan)))))
			semi_major = np.concatenate((semi_major, np.tile(sub_skyfld_earthsat.model.a * consts.R_EARTH, len(sub_timespan))))
			raan = np.concatenate((raan, np.tile(sub_skyfld_earthsat.model.nodeo, len(sub_timespan))))
			argp = np.concatenate((argp, np.tile(sub_skyfld_earthsat.model.argpo, len(sub_timespan))))
			TLE_epochs = np.concatenate((TLE_epochs, np.tile(tle_epoch, len(sub_timespan))))

		attr_dct['timespan'] = timespan
		attr_dct['gen_type'] = 'propagated from TLE'
		attr_dct['central_body'] = 'Earth'
		attr_dct['pos'] = pos
		attr_dct['alt'] = np.linalg.norm(pos, axis=1) - consts.R_EARTH
		attr_dct['pos_ecef'] = pos_ecef
		attr_dct['vel_ecef'] = vel_ecef
		attr_dct['vel'] = vel
		attr_dct['lat'] = lat
		attr_dct['lon'] = lon
		attr_dct['ecc'] = ecc
		attr_dct['inc'] = inc
		attr_dct['semi_major'] = semi_major
		attr_dct['raan'] = raan
		attr_dct['argp'] = argp
		attr_dct['TLE_epochs'] = TLE_epochs
		attr_dct['period'] = 2 * np.pi / sub_skyfld_earthsat.model.no_kozai * 60
		if timespan.time_step is not None:
			attr_dct['period_steps'] = int(attr_dct['period'] / timespan.time_step.total_seconds())
		else:
			attr_dct['period_steps'] = None
		attr_dct['name'] = sub_skyfld_earthsat.name

		return cls(attr_dct, calc_astrobodies=astrobodies)

	#fake TLE
	@classmethod
	def fromPropagatedOrbitalParam(cls, timespan:TimeSpan, a:float=6978, ecc:float=0, inc:float=0, raan:float=0, argp:float=0, mean_nu:float=0, name:str='Fake TLE', astrobodies:bool=True, unsafe:bool=False):
		"""Create an orbit from orbital parameters, propagated using sgp4.

		Orbits created using this class method will respect gravity corrections such as J4,
		allowing for semi-analytical sun-synchronous orbits.

		Parameters
		----------
		timespan : TimeSpan
			Timespan over which orbit is to be simulated
		a : float, optional
			semi-major axis of the orbit in km (the default is 6978 ~ 600km
			above the earth.)
		ecc : float, optional
			dimensionless number 0 < ecc < 1 (the default is 0, which is a
			circular orbit)
		inc : float, optional
			inclination of orbit in degrees (the default is 0, which represents
			an orbit around the Earth's equator)
		raan : float, optional
			right-ascension of the ascending node (the default is 0)
		argp : float, optional
			argument of the perigee in degrees (the default is 0, which
			represents an orbit with its semimajor axis in the plane of the
			Earth's equator)
		mean_nu : float, optional
			mean anomaly in degrees (the default is 0, which represents an orbit
			that is beginning at periapsis)

		Returns
		-------
		satplot.Orbit
		"""

		if a < consts.R_EARTH + consts.EARTH_MIN_ALT:
			logger.error("Semimajor axis, {}, is too close to Earth".format(a))
			if not unsafe:
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

		attr_dct = createEmptyOrbitAttrDict()
		attr_dct['name'] = name

		skyfld_ts = load.timescale(builtin=True)

		t0_epoch = epoch_u.datetime2sgp4epoch(timespan.start)
		satrec = Satrec()
		mean_motion = orbit_u.calcMeanMotion(a * 1e3)

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

		pos = np.empty((len(timespan), 3))
		vel = np.empty((len(timespan), 3))

		for ii, timestep in enumerate(timespan[:]):
			pos[ii, :] = skyfld_earthsat.at(skyfld_ts.utc(timestep)).position.km
			vel[ii, :] = skyfld_earthsat.at(skyfld_ts.utc(timestep)).velocity.km_per_s * 1000

		ecef = satrec.frame_xyz_and_velocity(itrs)
		pos_ecef = ecef[0].km.T
		vel_ecef = ecef[1].km_per_s.T * 1000

		period = 2 * np.pi / skyfld_earthsat.model.no_kozai * 60
		period_steps = int(period / timespan.time_step.total_seconds())

		lat, lon = orbit_u.convertCartesianToEllipsoidGeodetic(pos_ecef)

		attr_dct['timespan'] = timespan
		attr_dct['gen_type'] = 'propagated from orbital param'
		attr_dct['central_body'] = 'Earth'
		attr_dct['pos'] = pos
		attr_dct['alt'] = np.linalg.norm(pos, axis=1) - consts.R_EARTH
		attr_dct['vel'] = vel
		attr_dct['pos_ecef'] = pos_ecef
		attr_dct['vel_ecef'] = vel_ecef
		attr_dct['lat'] = lat
		attr_dct['lon'] = lon
		attr_dct['ecc'] = ecc * np.ones(len(timespan))
		attr_dct['inc'] = inc * np.ones(len(timespan))
		attr_dct['semi_major'] = a * np.ones(len(timespan))
		attr_dct['raan'] = raan * np.ones(len(timespan))
		attr_dct['argp'] = argp * np.ones(len(timespan))
		attr_dct['TLE_epochs'] = t0_epoch * np.ones(len(timespan))
		attr_dct['period'] = period
		attr_dct['period_steps'] = period_steps

		return cls(attr_dct, calc_astrobodies=astrobodies)

	#analytical
	@classmethod
	def fromAnalyticalOrbitalParam(cls, timespan:TimeSpan, body:str='Earth', a:float=6978, ecc:float=0, inc:float=0, raan:float=0, argp:float=0, mean_nu:float=0, name:str='Analytical', astrobodies:bool=True, unsafe:bool=False):
		"""Create an orbit from orbital parameters

		Parameters
		----------
		timespan : TimeSpan
			Timespan over which orbit is to be simulated
		body : str, optional
			Central body for the satellite: ['Earth', 'Moon', 'Mars', 'Sun'] (the default is 'Earth')
		a : float, optional
			semi-major axis of the orbit in km (the default is 6978 ~ 600km above the earth.)
		ecc : float, optional
			dimensionless number 0 < ecc < 1 (the default is 0, which is a circular orbit)
		inc : float, optional
			inclination of orbit in degrees (the default is 0, which represents
			an orbit around the Earth's equator)
		raan : float, optional
			right-ascension of the ascending node (the default is 0)
		argp : float, optional
			argument of the perigee in degrees (the default is 0, which
			represents an orbit with its semimajor axis in the plane of the
			Earth's equator)
		mean_nu : float, optional
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
			if not unsafe:
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

		attr_dct = createEmptyOrbitAttrDict()
		attr_dct['name'] = name

		logger.info("Creating analytical orbit")
		orb = hapsiraOrbit.from_classical(central_body,
											a * astropy_units.one * astropy_units.km,
											ecc * astropy_units.one,
											inc * astropy_units.one * astropy_units.deg,
											raan * astropy_units.one * astropy_units.deg,
											argp * astropy_units.one * astropy_units.deg,
											mean_nu * astropy_units.one * astropy_units.deg)


		logger.info("Creating ephemeris for orbit, using timespan")
		ephem = Ephem.from_orbit(orb, timespan.asAstropy())

		pos = np.asarray(ephem.rv()[0], dtype=np.float_)
		vel = np.asarray(ephem.rv()[1], dtype=np.float_) * 1000

		period = orb.period.unit.in_units('s') * orb.period.value
		period_steps = int(period / timespan.time_step.total_seconds())

		attr_dct['timespan'] = timespan
		attr_dct['gen_type'] = 'analytical orbit'
		attr_dct['central_body'] = body
		attr_dct['pos'] = pos
		# TODO: altitude should be sourced from central body
		attr_dct['alt'] = np.linalg.norm(pos, axis=1) - consts.R_EARTH
		attr_dct['vel'] = vel
		attr_dct['ecc'] = ecc * np.ones(len(timespan))
		attr_dct['inc'] = inc * np.ones(len(timespan))
		attr_dct['semi_major'] = a * np.ones(len(timespan))
		attr_dct['raan'] = raan * np.ones(len(timespan))
		attr_dct['argp'] = argp * np.ones(len(timespan))
		attr_dct['TLE_epochs'] = None
		attr_dct['period'] = period
		attr_dct['period_steps'] = period_steps

		return cls(attr_dct, calc_astrobodies=astrobodies)

	@classmethod
	def fromDummyConstantPosition(cls, timespan:TimeSpan,
							pos:tuple[float,float,float]|np.ndarray[tuple[int],np.dtype[np.float_]],
							sun_pos:tuple[float,float,float]|np.ndarray[tuple[int],np.dtype[np.float_]]|None=None,
							moon_pos:tuple[float,float,float]|np.ndarray[tuple[int],np.dtype[np.float_]]|None=None):
		'''
		creates an empty orbit with the motionless at the origin and the sun at +z
		useful for testing where an orbit is needed independent of any real ephemeris data
		'''
		attr_dct = createEmptyOrbitAttrDict()
		attr_dct['timespan'] = timespan
		attr_dct['name'] = 'dummy_const_pos'
		attr_dct['pos'] = np.full((len(timespan), 3), pos)
		attr_dct['vel'] = np.full((len(timespan), 3), 0)
		attr_dct['gen_type'] = 'constant position for testing'
		obj = cls(attr_dct, calc_astrobodies=False)
		if sun_pos is not None:
			obj.sun_pos = np.full((len(timespan), 3), sun_pos)
		if moon_pos is not None:
			obj.moon_pos = np.full((len(timespan), 3), moon_pos)
		return obj

	@classmethod
	def load(cls, file):
		with open(file, 'rb') as fp:
			logger.info(f'Loading orbit file from {file}')
			return pickle.load(fp)

	def getPosition(self, search_time):
		'''Return the position at the specified or closest time

		Args:
			search_time: the timestamp to search for

		Returns:
			position (km)

		Raises:
			ValueError: Raised 'search_time' is not a valid type
		'''
		return self._getAttributeClosestTime(self.pos, search_time)

	def getVelocity(self, search_time:dt.datetime|astropyTime) -> np.ndarray[tuple[int], np.dtype[np.float_]]:
		'''Return the velocity at the specified or closest time

		Args:
			search_time: the timestamp to search for

		Returns:
			velocity (m/s)

		Raises:
			ValueError: Raised 'search_time' is not a valid type
		'''
		return self._getAttributeClosestTime(self.vel, search_time)

	def _getAttributeClosestTime(self, attr:Any, search_time:dt.datetime|astropyTime) -> np.ndarray[tuple[int], np.dtype[np.float_]]:

		if isinstance(search_time, dt.datetime):
			_, closest_idx = self.timespan.getClosest(search_time)
		elif isinstance(search_time, astropyTime):
			search_dt = search_time.to_datetime()
			_, closest_idx = self.timespan.getClosest(search_dt)
		else:
			logger.error("%s cannot be used to index %s", search_time, attr)
			raise ValueError(f"{search_time} cannot be used to index {attr}")

		return attr[closest_idx, :]

	def _calcEclipse(self, pos, sun):
		earth_ang_size = np.arctan(consts.R_EARTH/consts.AU)
		neg_S_norm = np.linalg.norm(-sun,axis=1)
		neg_S_unit = -sun/neg_S_norm[:,None]
		p_neg_S = pos-sun
		p_neg_S_norm = np.linalg.norm(p_neg_S,axis=1)
		p_neg_S_unit = p_neg_S/p_neg_S_norm[:,None]
		earth_sat_ang_sep = np.arccos(np.sum(neg_S_unit*p_neg_S_unit, axis=1))
		sunlit_angsep_truth = earth_sat_ang_sep > earth_ang_size
		sunlit_dist_truth = p_neg_S_norm < neg_S_norm
		sunlit = np.logical_or(sunlit_angsep_truth, sunlit_dist_truth)

		return np.logical_not(sunlit)

	def serialise(self):
		raise NotImplementedError


def _findClosestEpochIndices(target, values):
	right_idx = np.searchsorted(target, values)
	left_idx = right_idx - 1
	# replace any idx greater than len of target with last idx in target
	right_idx[np.where(right_idx==len(target))] = len(target) - 1
	stacked_idx = np.hstack((left_idx.reshape(-1,1),right_idx.reshape(-1,1)))
	target_vals_at_idxs = target[stacked_idx]
	closest_idx_columns = np.argmin(np.abs(target_vals_at_idxs - values.reshape(-1,1)),axis=1)
	return stacked_idx[range(len(stacked_idx)),closest_idx_columns]