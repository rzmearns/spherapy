import unittest
import numpy as np
import spherapy.timespan as TS
import datetime as dt
import astropy.time as astrotime
import astropy.units as u

# Disable logging during testing
import logging
logging.disable(logging.CRITICAL)


class TestTimeSpan(unittest.TestCase):
	'''
	Unit testing for the methods of the TimeSpan class
	'''
	def test_init(self):
		'''
		Tests to see that the TimeSpan can be initialised from
		an initial time, a timeperiod, and a time step.
		Makes sure the output is sensible
		'''
		birthday = dt.datetime(1996, 9, 10, tzinfo=dt.timezone.utc)
		firstweek = TS.TimeSpan(birthday, '1d', '1W')
		
		dt_week = firstweek.asDatetime()
		astropy_week = firstweek.asAstropy()
		
		self.assertEqual(len(dt_week), 8)
		self.assertEqual(len(dt_week), len(astropy_week))
		# Assert they are of the correct types (dt/astropy.Time)
		self.assertIsInstance(dt_week[0], type(birthday))
		self.assertIsInstance(astropy_week[0], type(astrotime.Time.now()))
		# Assert the time deltas of each are one day
		self.assertEqual(dt_week[1] - dt_week[0], dt.timedelta(days=1))
		self.assertTrue(astrotime.TimeDelta(1 * u.day, scale='tai').isclose(astropy_week[1] - astropy_week[0]))
		# Assert the timedelta is consistent throughout the week
		self.assertEqual(dt_week[6] - dt_week[5], dt_week[1] - dt_week[0])
		self.assertTrue(astrotime.TimeDelta(1 * u.day, scale='tai').isclose(astropy_week[6] - astropy_week[5]))
		# Assert the first day for each is Sept 10, 1996
		self.assertEqual(dt_week[0], dt.datetime(1996, 9, 10, tzinfo=dt.timezone.utc))
		self.assertTrue(astrotime.Time('1996-09-10').isclose(astropy_week[0]))
		# Assert the last day for each is Sept 17, 1996
		self.assertEqual(dt_week[-1], dt.datetime(1996, 9, 17, tzinfo=dt.timezone.utc))
		self.assertTrue(astrotime.Time('1996-09-17').isclose(astropy_week[-1]))
	
	def test_attributes(self):
		'''
		Test for existance of required attributes in the timespan class		
		'''
		birthday = dt.datetime(1996, 9, 10)
		firstweek = TS.TimeSpan(birthday, '1d', '1W')

		# Test attributes exist
		firstweek.start
		firstweek.end
		firstweek.timezone
		firstweek.timezone_str
		firstweek.time_step
		firstweek.time_period
		firstweek.num_steps
		# Test len call works
		len(firstweek)

	def test_parsing(self):
		'''
		Tests to see that this module correctly parses tricky input with 
		a variety of units
		'''
		t0 = dt.datetime(2016, 12, 30, 0, 0, 1, tzinfo=dt.timezone.utc)

		# Test all desired step and period units
		TS.TimeSpan(t0, '1uS', '10uS')
		TS.TimeSpan(t0, '1mS', '10mS')
		TS.TimeSpan(t0, '1S', '10S')
		TS.TimeSpan(t0, '1M', '10M')
		TS.TimeSpan(t0, '1H', '10H')
		TS.TimeSpan(t0, '1d', '10d')
		TS.TimeSpan(t0, '1d', '10W')
		TS.TimeSpan(t0, '1d', '10m')
		TS.TimeSpan(t0, '1d', '10y')

		with self.assertRaises(ValueError):
			# Test period smaller than timestep
			TS.TimeSpan(t0, '1d', '1M')
			# Test incorrect step unit
			TS.TimeSpan(t0, '1k', '1W')
			# Test incorrect period unit
			TS.TimeSpan(t0, '1d', '1j')
		
		# Test non integer value of timestep
		self.assertEqual(41, len(TS.TimeSpan(t0, '1.5M', '1H')))
		
		# Test non integer value of period
		self.assertEqual(91, len(TS.TimeSpan(t0, '1M', '1.5H')))

		# Test integer number of steps
		self.assertEqual(35, len(TS.TimeSpan(t0, '1.75M', '1H')))
		self.assertEqual(3570, TS.TimeSpan(t0, '1.75M', '1H').time_period.total_seconds())
		self.assertEqual(t0 + dt.timedelta(seconds=3570), TS.TimeSpan(t0, '1.75M', '1H').end)

		# Check leap second handling
		# Previous leap second Dec 31 2016 23:59:60
		astro_delta = astrotime.Time('2017-01-01T0:0:1') - astrotime.Time('2016-12-30T0:0:1')
		astro_delta.format = 'sec'
		self.assertNotEqual(len(TS.TimeSpan(t0, '1S', '2d'))-1, astro_delta.value)

		# Check leap year handling
		self.assertEqual(367, len(TS.TimeSpan(dt.datetime(2016, 1, 1), '1d', '1y')))

	def test_closest_to(self):
		'''
		Tests to see that we can find the closest date/time in the TimeSpan
		to another date
		'''
		t0 = dt.datetime(1998, 6, 5, 3, 4, 25, tzinfo=dt.timezone.utc)
		timespan = TS.TimeSpan(t0, '1S', '1W')

		inner_target = dt.datetime(1998, 6, 8, 1, 5, 0, tzinfo=dt.timezone.utc)
		pre_target = dt.datetime(1996, 1, 1, tzinfo=dt.timezone.utc)
		post_target = dt.datetime(2010, 1, 1, tzinfo=dt.timezone.utc)

		half_second = dt.timedelta(microseconds=500000)

		# Test where the target date/time is in the TimeSpan
		res_time, res_index = timespan.getClosest(inner_target)
		self.assertEqual(res_index, (inner_target - t0).total_seconds())
		# Test where the target date/time is not in the TimeSpan
		res_time, res_index = timespan.getClosest(pre_target)
		self.assertEqual(res_index, 0)
		res_time, res_index = timespan.getClosest(post_target)
		self.assertEqual(res_index, len(timespan) - 1)
		# Test where there's a tie: ensure it picks the earlier date
		res_time, res_index = timespan.getClosest(t0 + half_second)
		self.assertEqual(res_index, 0)
		
	def test_time_zones(self):
		'''
		Tests to see that the initialisation process handles time zones
		and weird conventions correctly
		'''
		
		t0 = dt.datetime(1998, 6, 5, 3, 4, 25, tzinfo=dt.timezone.utc)
		# check decimal timezone parses correctly
		timespan_tz = TS.TimeSpan(t0, '1S', '1H', timezone='+6.5')
		timespan = TS.TimeSpan(t0, '1S', '1H')

		# check difference between timespan and timezone aware timespan is always equal to the timezone
		np.testing.assert_allclose(np.vectorize(lambda x: x.total_seconds())(timespan.asDatetime() - timespan_tz.asDatetime()),
													np.vectorize(lambda x: x.total_seconds())(dt.timedelta(hours=0) * np.ones(len(timespan))))

		# check ISO timezone parses correctly
		timespan_tz = TS.TimeSpan(t0, '1S', '1H', timezone='+6:30')

		# check unsigned timezone parses correctly
		timespan_tz = TS.TimeSpan(t0, '1S', '1H', timezone='6:30')

		# check negative timezone parses correctly
		timespan_tz = TS.TimeSpan(t0, '1S', '1H', timezone='-4:45')
		self.assertIsNone(np.testing.assert_allclose(np.vectorize(lambda x: x.total_seconds())(timespan.asDatetime() - timespan_tz.asDatetime()),
													np.vectorize(lambda x: x.total_seconds())(dt.timedelta(hours=0) * np.ones(len(timespan)))))
	
	def test_seconds_since_start(self):
		'''
		Tests that the `seconds_since_start` method works as expected
		'''
		birthday = dt.datetime(1996, 9, 10, tzinfo=dt.timezone.utc)
		timespan = TS.TimeSpan(birthday, '1H', '1d')
		seconds_elapsed = timespan.secondsSinceStart()
		self.assertEqual(len(seconds_elapsed), 25)
		self.assertEqual(seconds_elapsed[0], 0)
		np.testing.assert_almost_equal(seconds_elapsed[1], 3600)
		np.testing.assert_almost_equal(seconds_elapsed[-2], 82800)
