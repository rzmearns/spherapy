"""Class for orbital data.

This module provides:
- OrbitAttrDict: A TypedDict typing the instance attributes of an Orbit
- Orbit: A class of timestamped orbital data
"""
import datetime as dt
import logging
import pathlib

import types
import typing
from typing import Any, TypedDict, cast

from astropy import units as astropy_units
from astropy.time import Time as astropyTime
from hapsira.bodies import Earth, Mars, Moon, Sun
from hapsira.ephem import Ephem
from hapsira.twobody import Orbit as hapsiraOrbit
import numpy as np
from sgp4.api import WGS72, Satrec

# Don't need to import all of skyfield just for EarthSatellite loading.
# TODO: figure out the lightest import to use
from skyfield.api import EarthSatellite, load, wgs84
from skyfield.framelib import itrs

from spherapy.timespan import TimeSpan
from spherapy.util import epoch_u, exceptions
import spherapy.util.constants as consts
import spherapy.util.orbital_u as orbit_u

logger = logging.getLogger(__name__)


class OrbitAttrDict(TypedDict):
	"""A TypedDict providing type annotations for the Orbit class.

	Attributes:
		name:
		satcat_id:
		gen_type:
		timespan:
		TLE_epochs:
		pos:
		pos_ecef:
		vel_ecef:
		vel:
		lat:
		lon:
		sun_pos:
		moon_pos:
		alt:
		eclipse:
		central_body:
		period:
		period_steps:
		semi_major:
		ecc:
		inc:
		raan:
		argp:
	"""
	name: None|str
	satcat_id: None|int
	gen_type: None|str
	timespan: None|TimeSpan
	TLE_epochs: None|np.ndarray[tuple[int], np.dtype[np.datetime64]]
	pos: None|np.ndarray[tuple[int,int], np.dtype[np.float64]]
	pos_ecef: None|np.ndarray[tuple[int,int], np.dtype[np.float64]]
	vel_ecef: None|np.ndarray[tuple[int,int], np.dtype[np.float64]]
	vel: None|np.ndarray[tuple[int,int], np.dtype[np.float64]]
	lat: None|np.ndarray[tuple[int], np.dtype[np.float64]]
	lon: None|np.ndarray[tuple[int], np.dtype[np.float64]]
	sun_pos: None|np.ndarray[tuple[int,int], np.dtype[np.float64]]
	moon_pos: None|np.ndarray[tuple[int,int], np.dtype[np.float64]]
	alt: None|np.ndarray[tuple[int], np.dtype[np.float64]]
	eclipse: None|np.ndarray[tuple[int], np.dtype[np.bool_]]
	central_body: None|str
	period: None|float
	period_steps: None|int
	semi_major: None|np.ndarray[tuple[int], np.dtype[np.float64]]
	ecc: None|np.ndarray[tuple[int], np.dtype[np.float64]]
	inc: None|np.ndarray[tuple[int], np.dtype[np.float64]]
	raan: None|np.ndarray[tuple[int], np.dtype[np.float64]]
	argp: None|np.ndarray[tuple[int], np.dtype[np.float64]]

def _createEmptyOrbitAttrDict() -> OrbitAttrDict:
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

def _validateOrbitAttrDict(obj:OrbitAttrDict, typ: type) -> bool:
	"""Validates all fields of a TypedDict have been instantiated, and are of correct type.

	Raises:
		KeyError: If a field of the TypedDict is missing
		TypeError: If a field of the TYpedDict is incorrect
	"""
	for attr_name, spec_attr_type in typ.__annotations__.items():
		try:
			attr_value = obj.get(attr_name, None)
		except KeyError as e:
			# Check for missing keys
			raise KeyError('Typed Dict is missing key: %s', attr_name) from e
		if type(spec_attr_type) is types.UnionType:
			# check union of types
			good_type = False
			for sub_type in typing.get_args(spec_attr_type):
				# choose validator
				if sub_type in (None, types.NoneType):
					validator = _NoneValidator
				else:
					validator = _genericValidator
				if validator(attr_value, sub_type):
					good_type = True
					break
			if not good_type:
				raise TypeError(f'{attr_name} is type {type(attr_value)}, '
								f'should be {spec_attr_type}')
		elif not _genericValidator(attr_value, spec_attr_type):
			# check all other types
			raise TypeError(f'{attr_name} is type {type(attr_value)}, should be {spec_attr_type}')
	return True

def _NoneValidator(obj:object, check_type:type) -> bool:  #noqa: ARG001
	return type(obj) is types.NoneType

def _genericValidator(obj:object, check_type:type|types.GenericAlias) -> bool:
	if type(check_type) is types.GenericAlias:
		if type(obj) is check_type.__origin__:
			return True
	elif type(obj) is check_type:
			return True
	return False

class Orbit:
	"""Timestamped orbital data for a satellite, coordinate system depending on the central body.

	Attributes:
		timespan: Timespan over which orbit is to be simulated
		name: Name of the satellite
		satcat_id: NORAD satellite catelogue ID (if generated from a TLE)
		gen_type: How the orbit was generated
		TLE_epochs: Nx1 numpy array of TLE epoch used for propagation at each timestamp
				units: TLE epoch
		pos: Nx3 numpy array of cartesian coordinates of the position of the satellite
			at each timestamp
				units: km
				frame: ECI
		vel: Nx3 numpy array of cartesian velocities of the satellite at each timestamp
				units: m/s
				frame: ECI
		pos_ecef: Nx3 numpy array of cartesian coordinates of the position of the satellite
				at each timestamp
					units: km
					frame: ECEF
		vel_ecef: Nx3 numpy array of cartesian velocities of the satellite at each timestamp
					units: m/s
					frame: ECEF
		lat: Nx1 numpy array of central body latitudes of the satellite at each timestamp
				units: degrees
		lon: Nx1 numpy array of central body longitudes of the satellite at each timestamp
				units: degrees
		sun_pos: Nx3 numpy array of cartesian coordinates of the position of the Sun
				at each timestamp
					units: km
					frame: ECI
		moon_pos: Nx3 numpy array of cartesian coordinates of the position of the Moon
				at each timestamp
					units: km
					frame: ECI
		alt: Nx1 numpy array of altitudes above central body at each timestamp
				units: km
		eclipse: Nx1 numpy array of flag indicating if satellite is eclipsed at each timestamp
					units: km
		central_body: body the satellite is orbiting
		period: orbital period in secs
		period_steps: number of TimeSpan timestamps required to complete an orbit
		semi_major: Nx1 numpy array of orbit semi-major axis calculated at that timestep
						units: km
						will be constant if no orbital maneauvers
		ecc: Nx1 numpy array of orbit eccentricity calculated at that timestep
				units: unitless
				will be constant if no orbital maneauvers
		inc: Nx1 numpy array of orbit inclination calculated at that timestep
				units: degree
				will be constant if no orbital maneauvers
		raan: Nx1 numpy array of orbit RAAN calculated at that timestep
				units: degree
				will be constant if no orbital maneauvers
		argp: Nx1 numpy array of orbit Arg Perigee calculated at that timestep
				units: degree
				will be constant if no orbital maneauvers
	"""

	def __init__(self, data:OrbitAttrDict, calc_astrobodies:bool=False):
		"""The constructor should never be called directly.

		Use one of:
			Orbit.fromTLE()
			Orbit.fromListOfPositions()
			Orbit.fromPropagatedOrbitalParam()
			Orbit.fromAnalyticalOrbitalParam()
		"""
		# Should always be called from a class method, spit error if not.
		if not _validateOrbitAttrDict(data, OrbitAttrDict):
			logger.error("Orbit() should not be called directly, "
						"use one of the fromXX constructors.")
			raise ValueError("Orbit() should not be called directly, "
							"use one of the fromXX constructors")

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
			ephem_sun = Ephem.from_body(Sun,
										astropyTime(self.timespan.asAstropy(scale='tdb')),
										attractor=Earth)
			sun_pos = ephem_sun.rv()[0]
			self.sun_pos = np.asarray(sun_pos.to(astropy_units.km))

			logger.info('Creating ephemeris for Moon using timespan')
			ephem_moon = Ephem.from_body(Moon,
											astropyTime(self.timespan.asAstropy(scale='tdb')),
											attractor=Earth)
			moon_pos = ephem_moon.rv()[0]
			self.moon_pos = np.asarray(moon_pos.to(astropy_units.km))
			self.eclipse = self._calcEclipse(self.pos, self.sun_pos)

	#pos_list
	@classmethod
	def fromListOfPositions(cls, timespan:TimeSpan,
									positions:np.ndarray[tuple[int, int], np.dtype[np.float64]],
									astrobodies:bool=False) -> 'Orbit':
		"""Create an orbit from a list of positions.

		Creat an obit by explicitly specifying the position of the
		satellite at each point in time. Useful for simplified test cases; but
		may lead to unphysical orbits.

		Args:
			timespan: Timespan over which orbit is to be simulated
			positions: Nx3 numpy array of cartesian coordinates of the position of the satellite
						at each timestamp
							units: km
							frame: ECI
			astrobodies: [Optional] Flag to calculate Sun and Moon positions at timestamps
							Default is False

		Returns:
			satplot.Orbit
		"""
		attr_dct = _createEmptyOrbitAttrDict()
		if len(positions) != len(timespan):
			raise ValueError(f"Number of supplied positions does not match timespan length: "
							f"{len(positions)} =/= {len(timespan)}")

		attr_dct['timespan'] = timespan
		attr_dct['pos'] = positions
		# Assume linear motion between each position at each timestep;
		# Then assume it stops at the last timestep.
		vel = positions[1:] - positions[:-1]
		attr_dct['vel'] = np.concatenate((vel, np.array([[0, 0, 0]])))
		attr_dct['name'] = 'Sat from position list'
		attr_dct['gen_type'] = 'position list'

		return cls(attr_dct, calc_astrobodies=astrobodies)

	#TLE
	@classmethod
	def fromTLE(cls, timespan:TimeSpan,
						tle_path:pathlib.Path,
						astrobodies:bool=True,
						unsafe:bool=False) -> 'Orbit':
		"""Create an orbit from an existing TLE or a list of historical TLEs.

		Args:
			timespan : TimeSpan over which orbit is to be simulated
			tle_path : path to file containing TLEs for a satellite
			astrobodies: [Optional] Flag to calculate Sun and Moon positions at timestamps
							Default is False
			unsafe: [Optional] Flag to ignore TLE usage more than 14 days either side of timestamps
						Optional
						Default is False

		Returns:
		-------
			satplot.Orbit
		"""
		# load all TLEs and create skyfield earth sats
		skyfld_earth_sats = load.tle_file(f'{tle_path}')

		# generate list of epochs for each TLE
		epochs = np.asarray([a.epoch.utc_datetime() for a in skyfld_earth_sats])

		# ensure no duplicate skyfld sats
		_, unq_idxs = np.unique(epochs, return_index=True)
		unq_skyfld_earth_sats = list(np.asarray(skyfld_earth_sats)[unq_idxs])

		skyfld_earth_sats = unq_skyfld_earth_sats
		tle_epoch_dates = [sat.epoch.utc_datetime() for sat in skyfld_earth_sats]

		# check timespan is valid for provided TLEs
		if timespan.start < tle_epoch_dates[0] - dt.timedelta(days=14):
			logger.error("Timespan begins before provided TLEs (+14 days)")
			if not unsafe:
				raise exceptions.OutOfRangeError("Timespan begins before provided TLEs (+14 days)")
		elif timespan.start > tle_epoch_dates[-1] + dt.timedelta(days=14):
			logger.error("Timespan begins after provided TLEs (+14 days)")
			if not unsafe:
				raise exceptions.OutOfRangeError("Timespan begins after provided TLEs (+14 days)")
		elif timespan.end > tle_epoch_dates[-1] + dt.timedelta(days=14):
			logger.error("Timespan ends after provided TLEs (+14 days)")
			if not unsafe:
				raise exceptions.OutOfRangeError("Timespan ends after provided TLEs (+14 days)")

		attr_dct = _createEmptyOrbitAttrDict()

		_timespan_arr = timespan.asDatetime()
		_timespan_arr = cast("np.ndarray[tuple[int], np.dtype[np.datetime64]]", _timespan_arr)
		closest_tle_epochs = epoch_u.findClosestDatetimeIndices(_timespan_arr,
																	np.asarray(tle_epoch_dates))

		d = np.hstack((1, np.diff(closest_tle_epochs)))

		# timestep idx where TLE epoch changes, append len(timespan)
		# to not require separate loop iteration to handle last element
		timespan_epoch_trans_idxs = np.hstack((np.where(d!=0)[0],len(timespan)))
		# tle_dates idx where TLE epoch changes in timespan
		timespan_epoch_trans_tle_idxs = closest_tle_epochs[np.where(d!=0)[0]]

		tle_epoch_idxs = []
		sub_timespans = []
		for ii, trans_start_idx in enumerate(timespan_epoch_trans_idxs[:-1]):
			tle_epoch_idxs.append(timespan_epoch_trans_tle_idxs[ii])
			sub_timespans.append(timespan[trans_start_idx,timespan_epoch_trans_idxs[ii+1]])

		skyfld_ts = load.timescale(builtin=True)

		# Calculate data for first run
		sub_timespan = sub_timespans[0]
		# help out static type analysis (won't be None), no runtime effect
		sub_timespan = cast("np.ndarray[tuple[int], np.dtype[np.datetime64]]",sub_timespan)
		sub_skyfld_earthsat = skyfld_earth_sats[tle_epoch_idxs[0]]
		tle_epoch = tle_epoch_dates[tle_epoch_idxs[0]]
		sat_rec = sub_skyfld_earthsat.at(skyfld_ts.utc(sub_timespan))
		pos = sat_rec.position.km.T
		vel = sat_rec.velocity.km_per_s.T * 1000
		ecef = sat_rec.frame_xyz_and_velocity(itrs)
		pos_ecef = ecef[0].km.T
		vel_ecef = ecef[1].km_per_s.T * 1000

		skyfld_lat, skyfld_lon = wgs84.latlon_of(sat_rec)
		lat = skyfld_lat.degrees
		lon = skyfld_lon.degrees
		ecc = np.tile(sub_skyfld_earthsat.model.ecco,len(sub_timespan))
		inc = np.tile(sub_skyfld_earthsat.model.inclo,len(sub_timespan))
		semi_major = np.tile(sub_skyfld_earthsat.model.a * consts.R_EARTH,len(sub_timespan))
		raan = np.tile(sub_skyfld_earthsat.model.nodeo,len(sub_timespan))
		argp = np.tile(sub_skyfld_earthsat.model.argpo,len(sub_timespan))
		TLE_epochs = np.tile(tle_epoch,len(sub_timespan)) 		#noqa: N806


		# fill in data for any other sub timespans
		for ii in range(1,len(sub_timespans)):
			sub_timespan = sub_timespans[ii]
			# help out static type analysis (won't be None), no runtime effect
			sub_timespan = cast("np.ndarray[tuple[int], np.dtype[np.datetime64]]",sub_timespan)
			sub_skyfld_earthsat = skyfld_earth_sats[tle_epoch_idxs[ii]]
			tle_epoch = tle_epoch_dates[tle_epoch_idxs[ii]]
			sat_rec = sub_skyfld_earthsat.at(skyfld_ts.utc(sub_timespan))
			pos = np.vstack((pos, sat_rec.position.km.T))
			vel = np.vstack((vel, sat_rec.velocity.km_per_s.T * 1000))
			ecef = sat_rec.frame_xyz_and_velocity(itrs)
			pos_ecef = np.vstack((pos_ecef, ecef[0].km.T))
			vel_ecef = np.vstack((vel_ecef, ecef[1].km_per_s.T * 1000))

			skyfld_lat, skyfld_lon = wgs84.latlon_of(sat_rec)
			lat = np.concatenate((lat, skyfld_lat.degrees))
			lon = np.concatenate((lon, skyfld_lon.degrees))
			ecc = np.concatenate((ecc, np.tile(sub_skyfld_earthsat.model.ecco, len(sub_timespan))))
			inc = np.concatenate((inc,
									np.tile(sub_skyfld_earthsat.model.inclo, len(sub_timespan))))
			semi_major = np.concatenate((semi_major,
											np.tile(sub_skyfld_earthsat.model.a * consts.R_EARTH, len(sub_timespan)))) #noqa:E501
			raan = np.concatenate((raan,
									np.tile(sub_skyfld_earthsat.model.nodeo, len(sub_timespan))))
			argp = np.concatenate((argp,
									np.tile(sub_skyfld_earthsat.model.argpo, len(sub_timespan))))
			TLE_epochs = np.concatenate((TLE_epochs, np.tile(tle_epoch, len(sub_timespan)))) 		#noqa: N806

		attr_dct['timespan'] = timespan
		attr_dct['satcat_id'] = int(skyfld_earth_sats[0].target_name.split('#')[-1].split(' ')[0])
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

		period = float(2 * np.pi / sub_skyfld_earthsat.model.no_kozai * 60)
		attr_dct['period'] = period

		if timespan.time_step is not None:
			attr_dct['period_steps'] = int(period / timespan.time_step.total_seconds())
		else:
			attr_dct['period_steps'] = None
		attr_dct['name'] = sub_skyfld_earthsat.name

		return cls(attr_dct, calc_astrobodies=astrobodies)

	#fake TLE
	@classmethod
	def fromPropagatedOrbitalParam(cls, timespan:TimeSpan, 			#noqa: PLR0913
										a:float=6978,
										ecc:float=0,
										inc:float=0,
										raan:float=0,
										argp:float=0,
										mean_nu:float=0,
										name:str='Fake TLE',
										astrobodies:bool=True,
										unsafe:bool=False) -> 'Orbit':
		"""Create an orbit from orbital parameters, propagated using sgp4.

		Orbits created using this class method will respect gravity corrections such as J4,
		allowing for semi-analytical sun-synchronous orbits.

		Args:
			timespan: Timespan over which orbit is to be simulated
			a: [Optional] semi-major axis of the orbit in km
					Default is 6978 ~ 600km	above the earth.
			ecc: [Optional] eccentricty, dimensionless number 0 < ecc < 1
					Default is 0, which is a circular orbit
			inc: [Optional] inclination of orbit in degrees
					Default is 0, which represents an orbit around the Earth's equator
			raan: [Optional] right-ascension of the ascending node
					Default is 0
			argp: [Optional] argument of the perigee in degrees
					Default is 0, which	represents an orbit with its semimajor axis in
					the plane of the Earth's equator
			mean_nu: [Optional] mean anomaly in degrees
					Default is 0, which represents an orbit that is beginning at periapsis
			name: [Optional] string giving the name of the orbit
					Default is 'Fake TLE'
			astrobodies: [Optional] Flag to calculate Sun and Moon positions at timestamps
					Default is False
			unsafe: [Optional] Flag to ignore semi-major axis inside Earth's radius
					Default is False

		Returns:
			satplot.Orbit
		"""
		if a < consts.R_EARTH + consts.EARTH_MIN_ALT:
			logger.error("Semimajor axis, %s, is too close to Earth", a)
			if not unsafe:
				raise exceptions.OutOfRangeError(f"Semimajor axis, {a}, is too close to Earth")

		if ecc > 1 or ecc < 0:
			logger.error("Eccentricity, %s, is non circular or eliptical", ecc)
			raise exceptions.OutOfRangeError(f"Eccentricity, {ecc}, is non circular or eliptical")

		if inc > 180 or inc < -180: 			#noqa: PLR2004
			logger.error("Inclination, %s, is out of range, should be -180 < inc < 180", inc)
			raise exceptions.OutOfRangeError(f"Inclination, {inc}, is out of range, "
										f"should be -180 < inc < 180")

		if raan > 360 or raan < 0: 				#noqa: PLR2004
			logger.error("RAAN, %s, is out of range, should be 0 < inc < 360", raan)
			raise exceptions.OutOfRangeError(f"RAAN, {raan}, is out of range, "
												f"should be 0 < inc < 360")

		if argp > 360 or argp < 0: 				#noqa: PLR2004
			logger.error("Argument of periapsis, %s, is out of range, "
						"should be 0 < argp < 360", inc)
			raise exceptions.OutOfRangeError(f"Argument of periapsis, {inc}, is out of range, "
										f"should be 0 < argp < 360")

		if mean_nu > 360 or mean_nu < 0: 		#noqa: PLR2004
			logger.error("Mean anomaly, %s, is out of range, should be 0 < mean_nu < 360", mean_nu)
			raise exceptions.OutOfRangeError(f"Mean anomaly, {mean_nu}, is out of range, "
										f"should be 0 < mean_nu < 360")

		attr_dct = _createEmptyOrbitAttrDict()

		skyfld_ts = load.timescale(builtin=True)

		t0_epoch = epoch_u.datetime2sgp4epoch(timespan.start)
		satrec = Satrec()
		mean_motion = orbit_u.calcMeanMotion(a * 1e3)

		satrec.sgp4init(WGS72,		  # gravity model
						'i',  # mode
						1,  # satnum
						t0_epoch,  # mode
						0,  # bstar [/earth radii]
						0,  # ndot [revs/day]
						0,  # nddot [revs/day^3]
						ecc,  # ecc
						argp,  # arg perigee [radians]
						inc,  # inclination [radians]
						mean_nu,  # mean anomaly [radians]
						mean_motion * 60,  # mean motion [rad/min]
						raan)  # raan [radians]

		skyfld_earthsat = EarthSatellite.from_satrec(satrec, skyfld_ts)

		pos = np.empty((len(timespan), 3))
		vel = np.empty((len(timespan), 3))

		_timespan_arr = timespan.asDatetime()
		_timespan_arr = cast("np.ndarray[tuple[int], np.dtype[np.datetime64]]", _timespan_arr)
		for ii, timestep in enumerate(_timespan_arr):
			pos[ii, :] = skyfld_earthsat.at(skyfld_ts.utc(timestep)).position.km
			vel[ii, :] = skyfld_earthsat.at(skyfld_ts.utc(timestep)).velocity.km_per_s * 1000

		# TODO: calculate ecef and lat,lon values.
		# TODO: satrec doesn't have a frame_xyz_and_velocity in this case

		period = float(2 * np.pi / skyfld_earthsat.model.no_kozai * 60)
		if timespan.time_step is not None:
			period_steps = int(period / timespan.time_step.total_seconds())
		else:
			period_steps = None

		attr_dct['name'] = name
		attr_dct['timespan'] = timespan
		attr_dct['gen_type'] = 'propagated from orbital param'
		attr_dct['central_body'] = 'Earth'
		attr_dct['pos'] = pos
		attr_dct['alt'] = np.linalg.norm(pos, axis=1) - consts.R_EARTH
		attr_dct['vel'] = vel
		attr_dct['ecc'] = np.full(len(timespan), ecc)
		attr_dct['inc'] = np.full(len(timespan), inc)
		attr_dct['semi_major'] = np.full(len(timespan), a)
		attr_dct['raan'] = np.full(len(timespan), raan)
		attr_dct['argp'] = np.full(len(timespan), argp)
		attr_dct['TLE_epochs'] = np.full(len(timespan),timespan.start)
		attr_dct['period'] = period
		attr_dct['period_steps'] = period_steps

		return cls(attr_dct, calc_astrobodies=astrobodies)

	#analytical
	@classmethod
	def fromAnalyticalOrbitalParam(cls, timespan:TimeSpan, 					#noqa: C901, PLR0912, PLR0913
											body:str='Earth',
											a:float=6978,
											ecc:float=0,
											inc:float=0,
											raan:float=0,
											argp:float=0,
											mean_nu:float=0,
											name:str='Analytical',
											astrobodies:bool=True,
											unsafe:bool=False) -> 'Orbit':
		"""Create an analytical orbit defined by orbital parameters.

		Orbits created using this class method will NOT respect gravity corrections such as J4,
		and as such sun-synchronous orbit are not possible.

		Args:
			timespan: Timespan over which orbit is to be simulated
			body: [Optional] string indicating around what body the satellite orbits,
					Default is 'Earth'
					Options are ['Earth','Sun','Mars','Moon']
			a: [Optional] semi-major axis of the orbit in km
					Default is 6978 ~ 600km	above the earth.
			ecc: [Optional] eccentricty, dimensionless number 0 < ecc < 1
					Default is 0, which is a circular orbit
			inc: [Optional] inclination of orbit in degrees
					Default is 0, which represents an orbit around the Earth's equator
			raan: [Optional] right-ascension of the ascending node
					Default is 0
			argp: [Optional] argument of the perigee in degrees
					Default is 0, which	represents an orbit with its semimajor axis in
					the plane of the Earth's equator
			mean_nu: [Optional] mean anomaly in degrees
					Default is 0, which represents an orbit that is beginning at periapsis
			name: [Optional] string giving the name of the orbit
					Default is 'Analytical'
			astrobodies: [Optional] Flag to calculate Sun and Moon positions at timestamps
					Default is False
			unsafe: [Optional] Flag to ignore semi-major axis inside Earth's radius
					Default is False

		Returns:
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
			logger.error("Invalid central body %s. Valid options are Earth, Sun, Mars, Moon", body)
			raise ValueError(f"Invalid central body {body}.")

		if a < min_a:
			logger.error("Semimajor axis, %s, is too close to the central body, {body.upper()}", a)
			if not unsafe:
				raise exceptions.OutOfRangeError(f"Semimajor axis, {a}, is too close "
											f"to the central body, {body.upper()}")

		if ecc > 1 or ecc < 0:
			logger.error("Eccentricity, %s, is non circular or eliptical", ecc)
			raise exceptions.OutOfRangeError(f"Eccentricity, {ecc}, is non circular or eliptical")

		if inc > 180 or inc < -180: 			#noqa: PLR2004
			logger.error("Inclination, %s, is out of range, should be -180 < inc < 180", inc)
			raise exceptions.OutOfRangeError(f"Inclination, {inc}, is out of range, "
										f"should be -180 < inc < 180")

		if raan > 360 or raan < 0: 				#noqa: PLR2004
			logger.error("RAAN, %s, is out of range, should be 0 < inc < 360", inc)
			raise exceptions.OutOfRangeError(f"RAAN, {inc}, is out of range, "
												f"should be 0 < inc < 360")

		if argp > 360 or argp < 0: 				#noqa: PLR2004
			logger.error("Argument of periapsis, %s, is out of range, "
							"should be 0 < argp < 360", raan)
			raise exceptions.OutOfRangeError(f"Argument of periapsis, {raan}, is out of "
										f"range, should be 0 < argp < 360")

		if mean_nu > 360 or mean_nu < 0: 		#noqa: PLR2004
			logger.error("Mean anomaly, %s, is out of range, should be 0 < mean_nu < 360", mean_nu)
			raise exceptions.OutOfRangeError(f"Mean anomaly, {mean_nu}, is out of range, "
										f"should be 0 < mean_nu < 360")

		attr_dct = _createEmptyOrbitAttrDict()

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

		pos = np.asarray(ephem.rv()[0], dtype=np.float64)
		vel = np.asarray(ephem.rv()[1], dtype=np.float64) * 1000

		period = float(orb.period.unit.in_units('s') * orb.period.value)
		if timespan.time_step is not None:
			period_steps = int(period / timespan.time_step.total_seconds())
		else:
			period_steps = None

		attr_dct['name'] = name
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
		attr_dct['period'] = period
		attr_dct['period_steps'] = period_steps

		return cls(attr_dct, calc_astrobodies=astrobodies)

	@classmethod
	def fromDummyConstantPosition(cls, timespan:TimeSpan,
							pos:tuple[float,float,float]|np.ndarray[tuple[int],np.dtype[np.float64]],
							sun_pos:tuple[float,float,float]|np.ndarray[tuple[int],np.dtype[np.float64]]|None=None,
							moon_pos:tuple[float,float,float]|np.ndarray[tuple[int],np.dtype[np.float64]]|None=None) -> 'Orbit': 	#noqa: E501
		"""Creates an static orbit for testing.

		Satellite position is defined by pos, while sun and moon positions are optional,
		but can also be specified.

		Args:
			timespan: Timespan over which orbit is to be simulated
			pos: Nx3 numpy array of cartesian coordinates of the position of the satellite
				at each timestamp
					units: km
					frame: ECI
			sun_pos: [Optional] Nx3 numpy array of cartesian coordinates of the position of the Sun
					at each timestamp
						units: km
						frame: ECI
			moon_pos: [Optional] Nx3 numpy array of cartesian coordinates of the position of the
						Moon at each timestamp
							units: km
							frame: ECI

		Returns:
			satplot.Orbit
		"""
		attr_dct = _createEmptyOrbitAttrDict()
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

	def getPosition(self, search_time:dt.datetime|astropyTime) \
							-> np.ndarray[tuple[int],np.dtype[np.float64]]:
		"""Return the position at the specified or closest time.

		Args:
			search_time: the timestamp to search for

		Returns:
			position

		Raises:
			ValueError: orbit has no pos data
		"""
		if self.pos is None:
			raise ValueError("Orbit has no pos data")
		return self._getAttributeClosestTime(self.pos, search_time)

	def getVelocity(self, search_time:dt.datetime|astropyTime) \
							-> np.ndarray[tuple[int], np.dtype[np.float64]]:
		"""Return the velocity at the specified or closest time.

		Args:
			search_time: the timestamp to search for

		Returns:
			velocity

		Raises:
			ValueError: orbit has no vel data
		"""
		if self.vel is None:
			raise ValueError("Orbit has no vel data")
		return self._getAttributeClosestTime(self.vel, search_time)

	def _getAttributeClosestTime(self, attr:np.ndarray[Any, np.dtype[np.float64]],
										search_time:dt.datetime|astropyTime)\
										-> np.ndarray[tuple[int], np.dtype[np.float64]]:

		if self.timespan is None:
			raise ValueError("Orbit has no timespan")

		if isinstance(search_time, dt.datetime):
			_, closest_idx = self.timespan.getClosest(search_time)
		elif isinstance(search_time, astropyTime):
			search_dt = search_time.to_datetime()
			_, closest_idx = self.timespan.getClosest(search_dt)
		else:
			logger.error("%s cannot be used to index %s", search_time, attr)
			raise TypeError(f"{search_time} cannot be used to index {attr}")

		return attr[closest_idx, :]

	def _calcEclipse(self, pos:np.ndarray[tuple[int,int],np.dtype[np.float64]],
							sun:np.ndarray[tuple[int,int],np.dtype[np.float64]])\
							-> np.ndarray[tuple[int],np.dtype[np.bool_]]:
		earth_ang_size = np.arctan(consts.R_EARTH/consts.AU)
		neg_sun_norm = np.linalg.norm(-sun,axis=1)
		neg_sun_unit = -sun/neg_sun_norm[:,None]
		pos_neg_sun = pos-sun
		pos_neg_sun_norm = np.linalg.norm(pos_neg_sun,axis=1)
		pos_neg_sun_unit = pos_neg_sun/pos_neg_sun_norm[:,None]
		earth_sat_ang_sep = np.arccos(np.sum(neg_sun_unit*pos_neg_sun_unit, axis=1))
		sunlit_angsep_truth = earth_sat_ang_sep > earth_ang_size
		sunlit_dist_truth = pos_neg_sun_norm < neg_sun_norm
		sunlit = np.logical_or(sunlit_angsep_truth, sunlit_dist_truth)

		return np.logical_not(sunlit)



# def _findClosestEpochIndices(target:list[float], values:float) -> list[int]:
# 	right_idx = np.searchsorted(target, values)
# 	left_idx = right_idx - 1
# 	# replace any idx greater than len of target with last idx in target
# 	right_idx[np.where(right_idx==len(target))] = len(target) - 1
# 	stacked_idx = np.hstack((left_idx.reshape(-1,1),right_idx.reshape(-1,1)))
# 	target_vals_at_idxs = target[stacked_idx]
# 	closest_idx_columns = np.argmin(np.abs(target_vals_at_idxs - values.reshape(-1,1)),axis=1)
# 	return stacked_idx[range(len(stacked_idx)),closest_idx_columns]
