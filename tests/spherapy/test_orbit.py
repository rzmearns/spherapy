import unittest
import numpy as np
import thermpy.solver.mission.orbit as orbit
import thermpy.solver.mission.timespan as timespan
import thermpy.util.orbital_u as orbital_u
import datetime as dt

# Disable logging during testing
import logging
logging.disable(logging.CRITICAL)

class TestOrbit(unittest.TestCase):
	'''
	Unit testing for the methods of the Orbit class
	'''

	@classmethod
	def setUpClass(cls):
		"""init each test"""
		super(TestOrbit, cls).setUpClass()
		# ##Fixtures##
		cls.t0 = dt.datetime(2004, 12, 29, 0, 0, 1)
		cls.t = timespan.TimeSpan(cls.t0, '1S', '30M')
		cls.t02 = dt.datetime(2021, 11, 23, 0, 0, 1)
		cls.t1 = timespan.TimeSpan(cls.t02, '1S', '30M')

		# From multiple TLEs
		cls.o_tle = orbit.Orbit.from_tle(cls.t, 'data/tle/zarya_20041227-20041231.tle')
		# From single TLE
		cls.o_tle2 = orbit.Orbit.from_tle(cls.t1, 'data/tle/zarya_20211124.tle')
		# From fake TLE
		cls.o_ftle = orbit.Orbit.from_tle_orbital_param(cls.t, a=(6378 + 600), ecc=0, inc=45, raan=0, argp=0, mean_nu=0)
		# From orbital param
		cls.o_param = orbit.Orbit.from_orbital_param(cls.t, a=(6378 + 600), ecc=0, inc=45, raan=0, argp=0, mean_nu=0)
		# 	# TODO: implement posvel initialisation

	def test_attributes(self):
		'''
		Test for existance of required attributes in the orbit class		
		'''
		self.o_tle.pos
		self.o_tle.vel
		self.o_tle.period
		self.o_tle.period_steps
		self.o_tle.sun
		# From fake TLE
		self.o_ftle.pos
		self.o_ftle.vel
		self.o_ftle.period
		self.o_ftle.period_steps
		self.o_ftle.sun
		# From orbital param
		self.o_param.pos
		self.o_param.vel
		self.o_param.period
		self.o_param.period_steps
		self.o_param.sun
		# TODO: implement posvel init attribute checks

	def test_dimensions(self):
		'''
		Test pos and vel attributes are same dimensions as timespan
		'''

		# From TLE
		self.assertEqual(self.o_tle.pos.shape[0], len(self.t))
		self.assertEqual(self.o_tle.vel.shape[0], len(self.t))
		# From fake TLE
		self.assertEqual(self.o_ftle.pos.shape[0], len(self.t))
		self.assertEqual(self.o_ftle.vel.shape[0], len(self.t))
		# From orbital param
		self.assertEqual(self.o_param.pos.shape[0], len(self.t))
		self.assertEqual(self.o_param.vel.shape[0], len(self.t))

	def test_get_closest(self):
		'''
		Test get_position and get_velocity methods return an indexed position
		'''

	def test_orbit_values(self):
		'''
		Tests to check pos and vel values are sensible for an orbit
		'''

		t0 = dt.datetime(2021, 1, 1, 0, 0, 1)
		t = timespan.TimeSpan(t0, '1S', '180M')

		a = 6378 + 600

		o = orbit.Orbit.from_orbital_param(t, a=a, ecc=0, inc=45, raan=0, argp=0, mean_nu=0)

		# Circular orbit should have the same semi-major for its duration
		self.assertTrue(np.all(np.isclose(np.linalg.norm(o.pos, axis=1), a)))
		# Circular orbit should have the same orbital speed for its duration
		speed = orbital_u.calc_orbital_vel(a * 1e3, np.array((a * 1e3, 0, 0)))
		self.assertTrue(np.all(np.isclose(np.linalg.norm(o.vel, axis=1), speed)))

		# Should check elliptical orbits

	def test_nonearth_orbits(self):
		'''
		Test orbits can be constructed around other celestial bodies
		'''

	def test_predef_orbit(self):
		'''
		Test orbit values are as expected for a pre-defined orbit
		'''

		# self.cls.t0 = dt.datetime(2021, 1, 1, 0, 0, 1)
		# self.cls.t = timespan.TimeSpan(self.cls.t0, '1S', '90M')

		# o_tle = orbit.Orbit.from_tle(self.cls.t, 'data/tle/zarya_.tle')

	def test_timespan_tle_mismatch(self):
		'''
		Test errors are thrown if timespan is outside the epochs of the provided TLEs
		'''

	# def tearDown(self):
	# 	"""finish any test"""
	# 	p = Stats (self.pr)
	# 	p.strip_dirs()
	# 	p.sort_stats ('cumtime')
	# 	p.print_stats ()
	# 	print("\n--->>>")
