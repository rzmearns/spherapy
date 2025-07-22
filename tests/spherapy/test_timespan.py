from pytest_check import check
import numpy as np
from spherapy.timespan import TimeSpan
import datetime as dt
import astropy.time as astrotime
import astropy.units as u

# Disable logging during testing
import logging
logging.disable(logging.CRITICAL)

## Custom checker
def checkHasattr(a:object, name:str) -> bool:
	__tracebackhide__ = True
	if hasattr(a, name):
		return True
	check.fail(f"check {hasattr({a},{name})} == 4")
	return False

class TestTimeSpan:
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
		t0 = dt.datetime(2000,1,1,12,4,5)

		first_week = TimeSpan(birthday, '1d', '7d')
		first_min = TimeSpan(birthday, '1S', '1M')
		fractional_min = TimeSpan(t0,'1S','1.5M')
		first_year = TimeSpan(birthday, '1d', '365d')
		first_decade = TimeSpan(birthday, '365d','3650d')

		minimum_ts = TimeSpan(t0, '1S', '1S')

		# check lengths are correct
		check.equal(len(first_week), 8)
		check.equal(len(first_min), 61)
		check.equal(len(fractional_min), 91)
		check.equal(len(first_year), 366)
		check.equal(len(first_decade), 11)
		check.equal(len(minimum_ts), 2)

		# check start values are correct
		check.equal(first_week[0], birthday)
		check.equal(first_min[0], birthday)
		check.equal(fractional_min[0], t0.replace(tzinfo=dt.timezone.utc))
		check.equal(first_year[0], birthday)
		check.equal(first_decade[0], birthday)
		check.equal(minimum_ts[0], t0.replace(tzinfo=dt.timezone.utc))

		# check end values are correct
		check.equal(first_week[-1], birthday + dt.timedelta(days=7))
		check.equal(first_min[-1], birthday + dt.timedelta(minutes=1))
		check.equal(fractional_min[-1], t0.replace(tzinfo=dt.timezone.utc) + dt.timedelta(minutes=1.5))
		check.equal(first_year[-1], dt.datetime(1997,9,10, tzinfo=dt.timezone.utc))
		check.equal(first_decade[-1], dt.datetime(2006,9,8, tzinfo=dt.timezone.utc))
		check.equal(minimum_ts[-1], t0.replace(tzinfo=dt.timezone.utc) + dt.timedelta(seconds=1))

	def test_len(self):
		t0 = dt.datetime(2000,1,1,12,4,5)

		first_min = TimeSpan(t0, '1S', '1M')
		non_int_sec_span = TimeSpan(t0, '1.5S', '5S')

		check.equal(len(first_min), len(first_min._timearr))
		check.equal(len(non_int_sec_span), len(non_int_sec_span._timearr))

	def test_attributes(self):
		'''
		Test for existance of required attributes in the timespan class
		'''
		birthday = dt.datetime(1996, 9, 10)
		first_week = TimeSpan(birthday, '1d', '7d')

		# Test attributes exist
		checkHasattr(first_week,'start')
		check.equal(first_week.start, birthday.replace(tzinfo=dt.timezone.utc))
		checkHasattr(first_week,'end')
		check.equal(first_week.end, dt.datetime(1996,9,17,tzinfo=dt.timezone.utc))
		checkHasattr(first_week,'time_step')
		check.equal(first_week.time_step, dt.timedelta(days=1))
		checkHasattr(first_week,'time_period')
		check.equal(first_week.time_period, dt.timedelta(days=7))

	def test_nonIntegerSpans(self):
		birthday = dt.datetime(1996, 9, 10)
		t1 = dt.datetime(2016, 12, 30, 0, 0, 1)

		non_int_day_span = TimeSpan(birthday, '1.5d', '7d')
		non_int_sec_span = TimeSpan(birthday, '1.5S', '5S')

		# check correct number of elements
		check.equal(len(non_int_day_span), 5)
		check.equal(len(non_int_sec_span), 4)

		# Test integer number of steps
		check.equal(len(TimeSpan(t1, '1.75M', '1H')), 35)
		check.equal(TimeSpan(t1, '1.75M', '1H').time_period.total_seconds(), 3570)
		check.equal(TimeSpan(t1, '1.75M', '1H').end, t1.replace(tzinfo=dt.timezone.utc) + dt.timedelta(seconds=3570))

	def test_parsing(self):
		'''
		Tests to see that this module correctly parses tricky input with
		a variety of units
		'''
		t0 = dt.datetime(2016, 12, 30, 0, 0, 1, tzinfo=dt.timezone.utc)

		# Test all desired step and period units
		TimeSpan(t0, '1uS', '10uS')
		TimeSpan(t0, '1mS', '10mS')
		TimeSpan(t0, '1S', '10S')
		TimeSpan(t0, '1M', '10M')
		TimeSpan(t0, '1H', '10H')
		TimeSpan(t0, '1d', '10d')
		TimeSpan(t0, '1d', '70d')
		TimeSpan(t0, '1d', '300d')
		TimeSpan(t0, '1d', '3650d')

		with check.raises(ValueError):
			# Test period smaller than timestep
			TimeSpan(t0, '1d', '1M')
			# Test incorrect step unit
			TimeSpan(t0, '1k', '1d')
			# Test incorrect period unit
			TimeSpan(t0, '1d', '1j')

		# Test non integer value of timestep
		check.equal(len(TimeSpan(t0, '1.5M', '1H')), 41)

		# Test non integer value of period
		check.equal(len(TimeSpan(t0, '1M', '1.5H')), 91)

	def test_leaps(self):

		t0 = dt.datetime(2016, 12, 30, 0, 0, 1, tzinfo=dt.timezone.utc)

		# Check leap second handling
		# Previous leap second Dec 31 2016 23:59:60
		astro_delta = astrotime.Time('2017-01-01T0:0:1') - astrotime.Time('2016-12-30T0:0:1')
		astro_delta.format = 'sec'
		check.not_equal(len(TimeSpan(t0, '1S', '2d'))-1, astro_delta.value)

		# Check leap year handling
		ts = TimeSpan(dt.datetime(2016, 1, 1), '1d', '366d')
		check.equal(len(ts), 367)
		check.equal(ts[-1], dt.datetime(2017,1,1, tzinfo=dt.timezone.utc))

	def test_conversions(self):

		birthday = dt.datetime(1996, 9, 10, tzinfo=dt.timezone.utc)

		first_week = TimeSpan(birthday, '1d', '7d')

		dt_week = first_week.asDatetime()
		astropy_week = first_week.asAstropy()

		check.equal(len(dt_week), len(astropy_week))
		# Assert they are of the correct types (dt/astropy.Time)
		check.is_instance(dt_week[0], type(birthday))
		check.is_instance(astropy_week[0], type(astrotime.Time.now()))
		# Assert the time deltas of each are one day
		check.equal(dt_week[1] - dt_week[0], dt.timedelta(days=1))
		check.is_true(astrotime.TimeDelta(1 * u.day, scale='tai').isclose(astropy_week[1] - astropy_week[0]))
		# Assert the timedelta is consistent throughout the week
		check.equal(dt_week[6] - dt_week[5], dt_week[1] - dt_week[0])
		check.is_true(astrotime.TimeDelta(1 * u.day, scale='tai').isclose(astropy_week[6] - astropy_week[5]))
		# Assert the first day for each is Sept 10, 1996
		check.equal(dt_week[0], dt.datetime(1996, 9, 10, tzinfo=dt.timezone.utc))
		check.is_true(astrotime.Time('1996-09-10').isclose(astropy_week[0]))
		# Assert the last day for each is Sept 17, 1996
		check.equal(dt_week[-1], dt.datetime(1996, 9, 17, tzinfo=dt.timezone.utc))
		check.is_true(astrotime.Time('1996-09-17').isclose(astropy_week[-1]))

	def test_closestTo(self):
		'''
		Tests to see that we can find the closest date/time in the TimeSpan
		to another date
		'''
		t0 = dt.datetime(1998, 6, 5, 3, 4, 25)
		timespan = TimeSpan(t0, '1S', '7d')

		inner_target = dt.datetime(1998, 6, 8, 1, 5, 0)
		pre_target = dt.datetime(1996, 1, 1)
		post_target = dt.datetime(2010, 1, 1)

		half_second = dt.timedelta(microseconds=500000)

		# Test where the target date/time is in the TimeSpan
		res_time, res_index = timespan.getClosest(inner_target)
		check.equal(res_index, (inner_target - t0).total_seconds())
		# Test where the target date/time is not in the TimeSpan
		res_time, res_index = timespan.getClosest(pre_target)
		check.equal(res_index, 0)
		res_time, res_index = timespan.getClosest(post_target)
		check.equal(res_index, len(timespan) - 1)
		# Test where there's a tie: ensure it picks the earlier date
		res_time, res_index = timespan.getClosest(t0 + half_second)
		check.equal(res_index, 0)

	def test_timeZones(self):
		"""
		Tests to see that the initialisation process handles time zonescorrectly
		"""

		# all same time
		# t == t0 == t1 == t2 (all are the same instant in time)
		t = dt.datetime(2008, 6, 5, 2, 4, 25)
		t0 = dt.datetime(2008, 6, 5, 2, 4, 25, tzinfo=dt.timezone.utc)

		tz_offset_0545 = dt.timedelta(hours=5, minutes=45) # Kathmandu
		# time @ 1998-06-05 02:04:25 UTC but expressed in timezone +05:45
		t1 = dt.datetime(2008,6,5,7,49,25, tzinfo=dt.timezone(tz_offset_0545))

		tz_offset_2130 = dt.timedelta(hours=-2, minutes=-30) # Labrador (NDT)
		# time @ 1998-06-05 02:04:25 UTC but expressed in timezone -02:30
		t2 = dt.datetime(2008,6,4,23,34,25, tzinfo=dt.timezone(tz_offset_2130))

		timespan = TimeSpan(t, '1S', '1H')
		timespan_utc = TimeSpan(t0, '1S', '1H')
		timespan_0545 = TimeSpan(t1, '1S', '1H')
		timespan_2130 = TimeSpan(t2, '1S', '1H')

		# check timespan and timespan initialised with UTC time zone have all equal elements
		check.equal(len(timespan_utc), len(timespan))
		for ii in range(len(timespan)):
			check.equal(timespan_utc[ii], timespan[ii])

		# check time offset is applied correctly
		check.equal(timespan_0545[0], t.replace(tzinfo=dt.timezone.utc))
		check.equal(timespan_2130[0], t.replace(tzinfo=dt.timezone.utc))

		# check difference between UTC derived timespan and other timezone derived timespan is always 0
		np.testing.assert_allclose(np.vectorize(lambda x: x.total_seconds())(timespan.asDatetime() - timespan_0545.asDatetime()),
													np.zeros(len(timespan)))
		np.testing.assert_allclose(np.vectorize(lambda x: x.total_seconds())(timespan.asDatetime() - timespan_2130.asDatetime()),
													np.zeros(len(timespan)))


		# all same wall time
		t0 = dt.datetime(2008, 6, 5, 2, 4, 25, tzinfo=dt.timezone.utc)
		tz_offset_0545 = dt.timedelta(hours=5, minutes=45) # Kathmandu +05:45
		t1 = dt.datetime(2008, 6, 5, 2, 4, 25, tzinfo=dt.timezone(tz_offset_0545))
		tz_offset_2130 = dt.timedelta(hours=-2, minutes=-30) # Labrador (NDT) -02:30
		t2 = dt.datetime(2008, 6, 5, 2, 4, 25, tzinfo=dt.timezone(tz_offset_2130))

		timespan_utc = TimeSpan(t0, '1S', '1H')
		timespan_0545 = TimeSpan(t1, '1S', '1H')
		timespan_2130 = TimeSpan(t2, '1S', '1H')

		# check time offset is applied correctly
		check.equal(timespan_0545[0], dt.datetime(2008, 6, 4, 20, 19, 25, tzinfo=dt.timezone.utc))
		check.equal(timespan_2130[0], dt.datetime(2008, 6, 5, 4, 34, 25, tzinfo=dt.timezone.utc))

		# check difference between timespan and timezone aware timespan is always equal to the timezone
		np.testing.assert_allclose(np.vectorize(lambda x: x.total_seconds())(timespan.asDatetime() - timespan_0545.asDatetime()),
													tz_offset_0545.total_seconds()*np.ones(len(timespan)))
		np.testing.assert_allclose(np.vectorize(lambda x: x.total_seconds())(timespan.asDatetime() - timespan_2130.asDatetime()),
													tz_offset_2130.total_seconds()*np.ones(len(timespan)))

	def test_UTC(self):
		t0 = dt.datetime(2008, 6, 5, 2, 4, 25, tzinfo=dt.timezone.utc)

		tz_offset_0545 = dt.timedelta(hours=5, minutes=45) # Kathmandu
		# time @ 1998-06-05 02:04:25 UTC but expressed in timezone +05:45
		t1 = dt.datetime(2008,6,5,7,49,25, tzinfo=dt.timezone(tz_offset_0545))

		timespan_utc = TimeSpan(t0, '1S', '1H')
		timespan_0545 = TimeSpan(t1, '1S', '1H')

		check.equal(timespan_utc[0].tzinfo, dt.timezone.utc)
		check.equal(timespan_utc.start.tzinfo, dt.timezone.utc)
		check.equal(timespan_utc.end.tzinfo, dt.timezone.utc)
		for ii in range(len(timespan_utc)):
			check.equal(timespan_utc[ii].tzinfo,timespan_0545[ii].tzinfo)


	def test_secondsSinceStart(self):
		'''
		Tests that the `seconds_since_start` method works as expected
		'''
		birthday = dt.datetime(1996, 9, 10)
		timespan = TimeSpan(birthday, '1H', '1d')
		seconds_elapsed = timespan.secondsSinceStart()
		check.equal(len(seconds_elapsed), 25)
		check.equal(seconds_elapsed[0], 0)
		np.testing.assert_almost_equal(seconds_elapsed[1], 3600)
		np.testing.assert_almost_equal(seconds_elapsed[-2], 82800)
