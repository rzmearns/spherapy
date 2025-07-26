"""Command line script to test basic functionality."""
import datetime

import spherapy.orbit
import spherapy.timespan
import spherapy.updater

if __name__ == '__main__':
	print('Updating ISS TLEs') 										# noqa: T201
	updated_TLEs = spherapy.updater.updateTLEs([25544]) 	#ISS 	# noqa: N816
	TLE_paths = spherapy.updater.getTLEFilePaths([25544]) 	#ISS
	print('Creating Timespan')			 							# noqa: T201
	t = spherapy.timespan.TimeSpan(datetime.datetime(2024,10,15,0,0,1),'1S','90M')
	print('Creating Orbit') 										# noqa: T201
	o = spherapy.orbit.Orbit.fromTLE(t, TLE_paths[0])
