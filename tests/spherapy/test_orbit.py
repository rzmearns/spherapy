import pytest
from pytest_check import check
import pathlib
import numpy as np
from spherapy import orbit
from spherapy import timespan
from spherapy import updater
from spherapy.util import orbital_u
import spherapy.util.exceptions as exceptions
import datetime as dt

# Disable logging during testing
import logging
logging.disable(logging.CRITICAL)

class TestOrbit:
	'''
	Unit testing for the methods of the Orbit class
	'''

	# TODO: test getClosestAttribute
	# TODO: check validity of TLE Gen
	# TODO: check validity of propagated parameters Gen
	# TODO: check validity of Analytical Gen
	# TODO: test non earth orbits
	# TODO: test datetime split across multiple TLE epochs
		# TODO: check listed TLE epochs change at correct point of array
		# TODO: check less than 14 day propagation forward or back is all good

	@classmethod
	def setup_class(cls):
		"""init each test"""
		# Keep number of samples low, no need for lots of samples if only really testing attributes.
		cls.t0 = dt.datetime(2004, 12, 29, 0, 0, 1)
		cls.t = timespan.TimeSpan(cls.t0, '1S', '2M')
		cls.t02 = dt.datetime(2021, 11, 23, 0, 0, 1)
		cls.t1 = timespan.TimeSpan(cls.t02, '1S', '2M')
		cls.pos = np.tile(np.linspace(7200,8000,len(cls.t),dtype=np.float64),(3,1)).T

		# From multiple TLEs
		ISS_satcat_id = 25544 #noqa: N806
		cls.ISS_tle_path = updater.getTLEFilePaths([ISS_satcat_id])[0]
		cls.o_tle = orbit.Orbit.fromTLE(cls.t, pathlib.Path(cls.ISS_tle_path))
		cls.o_tle_astro = orbit.Orbit.fromTLE(cls.t1, pathlib.Path(cls.ISS_tle_path), astrobodies=True)

		# From fake TLE
		cls.o_ftle = orbit.Orbit.fromPropagatedOrbitalParam(cls.t, a=(6378 + 600), ecc=0, inc=45, raan=0, argp=0, mean_nu=0)
		cls.o_ftle_astro = orbit.Orbit.fromPropagatedOrbitalParam(cls.t, a=(6378 + 600), ecc=0, inc=45, raan=0, argp=0, mean_nu=0, astrobodies=True)

		# Analytic from orbital param
		cls.o_analytical = orbit.Orbit.fromAnalyticalOrbitalParam(cls.t, a=(6378 + 600), ecc=0, inc=45, raan=0, argp=0, mean_nu=0)
		cls.o_analytical_astro = orbit.Orbit.fromAnalyticalOrbitalParam(cls.t, a=(6378 + 600), ecc=0, inc=45, raan=0, argp=0, mean_nu=0, astrobodies=True)

		# From list of position
		cls.o_poslist = orbit.Orbit.fromListOfPositions(cls.t, cls.pos)
		cls.o_poslist_astro = orbit.Orbit.fromListOfPositions(cls.t, cls.pos, astrobodies=True)

		# Dummy const pos
		cls.o_static = orbit.Orbit.fromDummyConstantPosition(cls.t, np.zeros((len(cls.t),3)))
		cls.o_static_astro = orbit.Orbit.fromDummyConstantPosition(cls.t, np.zeros((len(cls.t),3)), sun_pos=np.ones((len(cls.t),3)), moon_pos=-1*np.ones((len(cls.t),3)))

	def test_TLEGenAttributes(self):
		'''
		Test for existance of required attributes in the orbit class
		'''

		self.o_tle.name
		self.o_tle.satcat_id
		self.o_tle.gen_type
		self.o_tle.timespan
		self.o_tle.TLE_epochs
		self.o_tle.pos
		self.o_tle.pos_ecef
		self.o_tle.vel_ecef
		self.o_tle.vel
		self.o_tle.lat
		self.o_tle.lon
		self.o_tle.sun_pos
		self.o_tle.moon_pos
		self.o_tle.alt
		self.o_tle.eclipse
		self.o_tle.central_body
		self.o_tle.period
		self.o_tle.period_steps
		self.o_tle.semi_major
		self.o_tle.ecc
		self.o_tle.inc
		self.o_tle.raan
		self.o_tle.argp
		assert self.o_tle.timespan is not None
		assert self.o_tle.satcat_id is not None
		assert self.o_tle.name is not None
		assert self.o_tle.gen_type is not None
		assert self.o_tle.central_body is not None
		assert self.o_tle.pos is not None
		assert self.o_tle.vel is not None
		assert self.o_tle.alt is not None
		assert self.o_tle.pos_ecef is not None
		assert self.o_tle.vel_ecef is not None
		assert self.o_tle.lat is not None
		assert self.o_tle.lon is not None
		assert self.o_tle.ecc is not None
		assert self.o_tle.inc is not None
		assert self.o_tle.semi_major is not None
		assert self.o_tle.semi_major is not None
		assert self.o_tle.raan is not None
		assert self.o_tle.argp is not None
		assert self.o_tle.TLE_epochs is not None
		assert self.o_tle.period is not None
		assert self.o_tle.period_steps is not None

		self.o_tle_astro.name
		self.o_tle_astro.satcat_id
		self.o_tle_astro.gen_type
		self.o_tle_astro.timespan
		self.o_tle_astro.TLE_epochs
		self.o_tle_astro.pos
		self.o_tle_astro.pos_ecef
		self.o_tle_astro.vel_ecef
		self.o_tle_astro.vel
		self.o_tle_astro.lat
		self.o_tle_astro.lon
		self.o_tle_astro.sun_pos
		self.o_tle_astro.moon_pos
		self.o_tle_astro.alt
		self.o_tle_astro.eclipse
		self.o_tle_astro.central_body
		self.o_tle_astro.period
		self.o_tle_astro.period_steps
		self.o_tle_astro.semi_major
		self.o_tle_astro.ecc
		self.o_tle_astro.inc
		self.o_tle_astro.raan
		self.o_tle_astro.argp
		assert self.o_tle_astro.timespan is not None
		assert self.o_tle_astro.satcat_id is not None
		assert self.o_tle_astro.name is not None
		assert self.o_tle_astro.gen_type is not None
		assert self.o_tle_astro.central_body is not None
		assert self.o_tle_astro.pos is not None
		assert self.o_tle_astro.vel is not None
		assert self.o_tle_astro.alt is not None
		assert self.o_tle_astro.pos_ecef is not None
		assert self.o_tle_astro.vel_ecef is not None
		assert self.o_tle_astro.lat is not None
		assert self.o_tle_astro.lon is not None
		assert self.o_tle_astro.ecc is not None
		assert self.o_tle_astro.inc is not None
		assert self.o_tle_astro.semi_major is not None
		assert self.o_tle_astro.semi_major is not None
		assert self.o_tle_astro.raan is not None
		assert self.o_tle_astro.argp is not None
		assert self.o_tle_astro.TLE_epochs is not None
		assert self.o_tle_astro.period is not None
		assert self.o_tle_astro.period_steps is not None
		assert self.o_tle_astro.sun_pos is not None
		assert self.o_tle_astro.moon_pos is not None
		assert self.o_tle_astro.eclipse is not None

	def test_FakeTLEGenAttributes(self):
		'''
		Test for existance of required attributes in the orbit class
		'''

		self.o_ftle.name
		self.o_ftle.satcat_id
		self.o_ftle.gen_type
		self.o_ftle.timespan
		self.o_ftle.TLE_epochs
		self.o_ftle.pos
		self.o_ftle.pos_ecef
		self.o_ftle.vel_ecef
		self.o_ftle.vel
		self.o_ftle.lat
		self.o_ftle.lon
		self.o_ftle.sun_pos
		self.o_ftle.moon_pos
		self.o_ftle.alt
		self.o_ftle.eclipse
		self.o_ftle.central_body
		self.o_ftle.period
		self.o_ftle.period_steps
		self.o_ftle.semi_major
		self.o_ftle.ecc
		self.o_ftle.inc
		self.o_ftle.raan
		self.o_ftle.argp
		assert self.o_ftle.name is not None
		assert self.o_ftle.timespan is not None
		assert self.o_ftle.gen_type is not None
		assert self.o_ftle.central_body is not None
		assert self.o_ftle.pos is not None
		assert self.o_ftle.vel is not None
		assert self.o_ftle.alt is not None
		# TODO: include pos_ecef, vel_ecef, lat, lon once these are included in the generation
		assert self.o_ftle.ecc is not None
		assert self.o_ftle.inc is not None
		assert self.o_ftle.semi_major is not None
		assert self.o_ftle.raan is not None
		assert self.o_ftle.argp is not None
		assert self.o_ftle.TLE_epochs is not None
		assert self.o_ftle.period is not None
		assert self.o_ftle.period_steps is not None


		self.o_ftle_astro.name
		self.o_ftle_astro.satcat_id
		self.o_ftle_astro.gen_type
		self.o_ftle_astro.timespan
		self.o_ftle_astro.TLE_epochs
		self.o_ftle_astro.pos
		self.o_ftle_astro.pos_ecef
		self.o_ftle_astro.vel_ecef
		self.o_ftle_astro.vel
		self.o_ftle_astro.lat
		self.o_ftle_astro.lon
		self.o_ftle_astro.sun_pos
		self.o_ftle_astro.moon_pos
		self.o_ftle_astro.alt
		self.o_ftle_astro.eclipse
		self.o_ftle_astro.central_body
		self.o_ftle_astro.period
		self.o_ftle_astro.period_steps
		self.o_ftle_astro.semi_major
		self.o_ftle_astro.ecc
		self.o_ftle_astro.inc
		self.o_ftle_astro.raan
		self.o_ftle_astro.argp
		assert self.o_ftle_astro.name is not None
		assert self.o_ftle_astro.timespan is not None
		assert self.o_ftle_astro.gen_type is not None
		assert self.o_ftle_astro.central_body is not None
		assert self.o_ftle_astro.pos is not None
		assert self.o_ftle_astro.vel is not None
		assert self.o_ftle_astro.alt is not None
		# TODO: include pos_ecef, vel_ecef, lat, lon once these are included in the generation
		assert self.o_ftle_astro.ecc is not None
		assert self.o_ftle_astro.inc is not None
		assert self.o_ftle_astro.semi_major is not None
		assert self.o_ftle_astro.raan is not None
		assert self.o_ftle_astro.argp is not None
		assert self.o_ftle_astro.TLE_epochs is not None
		assert self.o_ftle_astro.period is not None
		assert self.o_ftle_astro.period_steps is not None
		assert self.o_ftle_astro.sun_pos is not None
		assert self.o_ftle_astro.moon_pos is not None
		assert self.o_ftle_astro.eclipse is not None

	def test_AnalyticalGenAttributes(self):
		'''
		Test for existance of required attributes in the orbit class
		'''

		self.o_analytical.name
		self.o_analytical.satcat_id
		self.o_analytical.gen_type
		self.o_analytical.timespan
		self.o_analytical.TLE_epochs
		self.o_analytical.pos
		self.o_analytical.pos_ecef
		self.o_analytical.vel_ecef
		self.o_analytical.vel
		self.o_analytical.lat
		self.o_analytical.lon
		self.o_analytical.sun_pos
		self.o_analytical.moon_pos
		self.o_analytical.alt
		self.o_analytical.eclipse
		self.o_analytical.central_body
		self.o_analytical.period
		self.o_analytical.period_steps
		self.o_analytical.semi_major
		self.o_analytical.ecc
		self.o_analytical.inc
		self.o_analytical.raan
		self.o_analytical.argp
		assert self.o_analytical.name is not None
		assert self.o_analytical.timespan is not None
		assert self.o_analytical.gen_type is not None
		assert self.o_analytical.central_body is not None
		assert self.o_analytical.pos is not None
		assert self.o_analytical.vel is not None
		assert self.o_analytical.alt is not None
		assert self.o_analytical.ecc is not None
		assert self.o_analytical.inc is not None
		assert self.o_analytical.semi_major is not None
		assert self.o_analytical.raan is not None
		assert self.o_analytical.argp is not None
		assert self.o_analytical.period is not None
		assert self.o_analytical.period_steps is not None

		self.o_analytical_astro.name
		self.o_analytical_astro.satcat_id
		self.o_analytical_astro.gen_type
		self.o_analytical_astro.timespan
		self.o_analytical_astro.TLE_epochs
		self.o_analytical_astro.pos
		self.o_analytical_astro.pos_ecef
		self.o_analytical_astro.vel_ecef
		self.o_analytical_astro.vel
		self.o_analytical_astro.lat
		self.o_analytical_astro.lon
		self.o_analytical_astro.sun_pos
		self.o_analytical_astro.moon_pos
		self.o_analytical_astro.alt
		self.o_analytical_astro.eclipse
		self.o_analytical_astro.central_body
		self.o_analytical_astro.period
		self.o_analytical_astro.period_steps
		self.o_analytical_astro.semi_major
		self.o_analytical_astro.ecc
		self.o_analytical_astro.inc
		self.o_analytical_astro.raan
		self.o_analytical_astro.argp
		assert self.o_analytical_astro.name is not None
		assert self.o_analytical_astro.timespan is not None
		assert self.o_analytical_astro.gen_type is not None
		assert self.o_analytical_astro.central_body is not None
		assert self.o_analytical_astro.pos is not None
		assert self.o_analytical_astro.vel is not None
		assert self.o_analytical_astro.alt is not None
		assert self.o_analytical_astro.ecc is not None
		assert self.o_analytical_astro.inc is not None
		assert self.o_analytical_astro.semi_major is not None
		assert self.o_analytical_astro.raan is not None
		assert self.o_analytical_astro.argp is not None
		assert self.o_analytical_astro.period is not None
		assert self.o_analytical_astro.period_steps is not None
		assert self.o_analytical_astro.sun_pos is not None
		assert self.o_analytical_astro.moon_pos is not None
		assert self.o_analytical_astro.eclipse is not None

	def test_PosListGenAttributes(self):
		'''
		Test for existance of required attributes in the orbit class
		'''

		self.o_poslist.name
		self.o_poslist.satcat_id
		self.o_poslist.gen_type
		self.o_poslist.timespan
		self.o_poslist.TLE_epochs
		self.o_poslist.pos
		self.o_poslist.pos_ecef
		self.o_poslist.vel_ecef
		self.o_poslist.vel
		self.o_poslist.lat
		self.o_poslist.lon
		self.o_poslist.sun_pos
		self.o_poslist.moon_pos
		self.o_poslist.alt
		self.o_poslist.eclipse
		self.o_poslist.central_body
		self.o_poslist.period
		self.o_poslist.period_steps
		self.o_poslist.semi_major
		self.o_poslist.ecc
		self.o_poslist.inc
		self.o_poslist.raan
		self.o_poslist.argp
		assert self.o_poslist.timespan is not None
		assert self.o_poslist.name is not None
		assert self.o_poslist.gen_type is not None
		assert self.o_poslist.pos is not None
		assert self.o_poslist.vel is not None

		self.o_poslist_astro.name
		self.o_poslist_astro.satcat_id
		self.o_poslist_astro.gen_type
		self.o_poslist_astro.timespan
		self.o_poslist_astro.TLE_epochs
		self.o_poslist_astro.pos
		self.o_poslist_astro.pos_ecef
		self.o_poslist_astro.vel_ecef
		self.o_poslist_astro.vel
		self.o_poslist_astro.lat
		self.o_poslist_astro.lon
		self.o_poslist_astro.sun_pos
		self.o_poslist_astro.moon_pos
		self.o_poslist_astro.alt
		self.o_poslist_astro.eclipse
		self.o_poslist_astro.central_body
		self.o_poslist_astro.period
		self.o_poslist_astro.period_steps
		self.o_poslist_astro.semi_major
		self.o_poslist_astro.ecc
		self.o_poslist_astro.inc
		self.o_poslist_astro.raan
		self.o_poslist_astro.argp
		assert self.o_poslist_astro.timespan is not None
		assert self.o_poslist_astro.name is not None
		assert self.o_poslist_astro.gen_type is not None
		assert self.o_poslist_astro.pos is not None
		assert self.o_poslist_astro.vel is not None
		assert self.o_poslist_astro.sun_pos is not None
		assert self.o_poslist_astro.moon_pos is not None
		assert self.o_poslist_astro.eclipse is not None

	def test_StaticGenAttributes(self):
		'''
		Test for existance of required attributes in the orbit class
		'''

		self.o_static.name
		self.o_static.satcat_id
		self.o_static.gen_type
		self.o_static.timespan
		self.o_static.TLE_epochs
		self.o_static.pos
		self.o_static.pos_ecef
		self.o_static.vel_ecef
		self.o_static.vel
		self.o_static.lat
		self.o_static.lon
		self.o_static.sun_pos
		self.o_static.moon_pos
		self.o_static.alt
		self.o_static.eclipse
		self.o_static.central_body
		self.o_static.period
		self.o_static.period_steps
		self.o_static.semi_major
		self.o_static.ecc
		self.o_static.inc
		self.o_static.raan
		self.o_static.argp
		assert self.o_static.timespan is not None
		assert self.o_static.name is not None
		assert self.o_static.gen_type is not None
		assert self.o_static.pos is not None
		assert self.o_static.vel is not None

		self.o_static_astro.name
		self.o_static_astro.satcat_id
		self.o_static_astro.gen_type
		self.o_static_astro.timespan
		self.o_static_astro.TLE_epochs
		self.o_static_astro.pos
		self.o_static_astro.pos_ecef
		self.o_static_astro.vel_ecef
		self.o_static_astro.vel
		self.o_static_astro.lat
		self.o_static_astro.lon
		self.o_static_astro.sun_pos
		self.o_static_astro.moon_pos
		self.o_static_astro.alt
		self.o_static_astro.eclipse
		self.o_static_astro.central_body
		self.o_static_astro.period
		self.o_static_astro.period_steps
		self.o_static_astro.semi_major
		self.o_static_astro.ecc
		self.o_static_astro.inc
		self.o_static_astro.raan
		self.o_static_astro.argp
		assert self.o_static_astro.timespan is not None
		assert self.o_static_astro.name is not None
		assert self.o_static_astro.gen_type is not None
		assert self.o_static_astro.pos is not None
		assert self.o_static_astro.vel is not None
		assert self.o_static_astro.sun_pos is not None
		assert self.o_static_astro.moon_pos is not None

	def test_dimensions(self):
	# 	'''
	# 	Test pos and vel attributes are same dimensions as timespan
	# 	'''

		# From TLE
		check.equal(self.o_tle.pos.shape[0], len(self.t))
		check.equal(self.o_tle.vel.shape[0], len(self.t))
		check.equal(self.o_tle.pos_ecef.shape[0], len(self.t))
		check.equal(self.o_tle.vel_ecef.shape[0], len(self.t))
		check.equal(self.o_tle.lat.shape[0], len(self.t))
		check.equal(self.o_tle.lon.shape[0], len(self.t))
		check.equal(self.o_tle.alt.shape[0], len(self.t))
		check.equal(self.o_tle.eclipse.shape[0], len(self.t))
		check.equal(self.o_tle.TLE_epochs.shape[0], len(self.t))
		check.equal(self.o_tle_astro.sun_pos.shape[0], len(self.t))
		check.equal(self.o_tle_astro.moon_pos.shape[0], len(self.t))
		check.equal(self.o_tle.semi_major.shape[0], len(self.t))
		check.equal(self.o_tle.ecc.shape[0], len(self.t))
		check.equal(self.o_tle.inc.shape[0], len(self.t))
		check.equal(self.o_tle.raan.shape[0], len(self.t))
		check.equal(self.o_tle.argp.shape[0], len(self.t))

		# From fake TLE
		check.equal(self.o_ftle.pos.shape[0], len(self.t))
		check.equal(self.o_ftle.vel.shape[0], len(self.t))
		# TODO: include pos_ecef, vel_ecef, lat, lon once these are included in the generation
		check.equal(self.o_ftle.alt.shape[0], len(self.t))
		check.equal(self.o_ftle.eclipse.shape[0], len(self.t))
		check.equal(self.o_ftle.TLE_epochs.shape[0], len(self.t))
		check.equal(self.o_ftle_astro.sun_pos.shape[0], len(self.t))
		check.equal(self.o_ftle_astro.moon_pos.shape[0], len(self.t))
		check.equal(self.o_ftle.semi_major.shape[0], len(self.t))
		check.equal(self.o_ftle.ecc.shape[0], len(self.t))
		check.equal(self.o_ftle.inc.shape[0], len(self.t))
		check.equal(self.o_ftle.raan.shape[0], len(self.t))
		check.equal(self.o_ftle.argp.shape[0], len(self.t))

		# From orbital param
		check.equal(self.o_analytical.pos.shape[0], len(self.t))
		check.equal(self.o_analytical.vel.shape[0], len(self.t))
		check.equal(self.o_analytical.alt.shape[0], len(self.t))
		check.equal(self.o_analytical.eclipse.shape[0], len(self.t))
		check.equal(self.o_analytical_astro.sun_pos.shape[0], len(self.t))
		check.equal(self.o_analytical_astro.moon_pos.shape[0], len(self.t))
		check.equal(self.o_analytical.semi_major.shape[0], len(self.t))
		check.equal(self.o_analytical.ecc.shape[0], len(self.t))
		check.equal(self.o_analytical.inc.shape[0], len(self.t))
		check.equal(self.o_analytical.raan.shape[0], len(self.t))
		check.equal(self.o_analytical.argp.shape[0], len(self.t))

		# From position list
		check.equal(self.o_poslist.pos.shape[0], len(self.t))
		check.equal(self.o_poslist.vel.shape[0], len(self.t))

		# From const pos
		check.equal(self.o_static.pos.shape[0], len(self.t))
		check.equal(self.o_static.vel.shape[0], len(self.t))
		check.equal(self.o_static_astro.sun_pos.shape[0], len(self.t))
		check.equal(self.o_static_astro.moon_pos.shape[0], len(self.t))

	@pytest.mark.filterwarnings("ignore::erfa.ErfaWarning")
	def test_unsafe(self):
		# TLE prior to earliest zarya epoch
		prior_start_t = timespan.TimeSpan(dt.datetime(1997,1,1),'1S','1M')
		# timespan begins more than 14 days after last TLE epoch
		post_start_t = timespan.TimeSpan(dt.datetime.now() + dt.timedelta(days=30),'1S','1M')
		# timespan begins within given TLE epochs, but more than 14 days after last TLE epoch
		# needs to be earlier than 2100 for external library compatibility
		post_end_t = timespan.TimeSpan(dt.datetime(2005,6,1),'365d','34310d')

		with check.raises(exceptions.OutOfRangeError):
			orbit.Orbit.fromTLE(prior_start_t, pathlib.Path(self.ISS_tle_path))
			orbit.Orbit.fromTLE(post_start_t, pathlib.Path(self.ISS_tle_path))
			orbit.Orbit.fromTLE(post_end_t, pathlib.Path(self.ISS_tle_path))

		orbit.Orbit.fromTLE(prior_start_t, pathlib.Path(self.ISS_tle_path), unsafe=True)
		orbit.Orbit.fromTLE(post_start_t, pathlib.Path(self.ISS_tle_path), unsafe=True)
		orbit.Orbit.fromTLE(post_end_t, pathlib.Path(self.ISS_tle_path), unsafe=True)

		with check.raises(exceptions.OutOfRangeError):
			orbit.Orbit.fromPropagatedOrbitalParam(self.t, a=(6378 - 600), ecc=0, inc=45, raan=0, argp=0, mean_nu=0)
			orbit.Orbit.fromPropagatedOrbitalParam(self.t, a=(6378 - 600), ecc=0, inc=45, raan=0, argp=0, mean_nu=0)

		orbit.Orbit.fromPropagatedOrbitalParam(self.t, a=(6378 - 600), ecc=0, inc=45, raan=0, argp=0, mean_nu=0, unsafe=True)
		orbit.Orbit.fromPropagatedOrbitalParam(self.t, a=(6378 - 600), ecc=0, inc=45, raan=0, argp=0, mean_nu=0, unsafe=True)

	def test_analyticalValidity(self):
		'''
		Tests to check pos and vel values are sensible for an orbit
		'''

		t0 = dt.datetime(2021, 1, 1, 0, 0, 1)
		t = timespan.TimeSpan(t0, '1S', '180M')
		a = 6378+600
		o = orbit.Orbit.fromAnalyticalOrbitalParam(t, a=a, ecc=0, inc=45, raan=0, argp=0, mean_nu=0)

		# Circular orbit should have the same semi-major for its duration
		check.is_true(np.all(np.isclose(np.linalg.norm(o.pos, axis=1), a)))
		# Circular orbit should have the same orbital speed for its duration
		speed = orbital_u.calcOrbitalVel(a * 1e3, np.array((a * 1e3, 0, 0)))
		check.is_true(np.all(np.isclose(np.linalg.norm(o.vel, axis=1), speed)))
