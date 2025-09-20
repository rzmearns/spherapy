## About
Spherapy is a convenience wrapper around the fantastic libraries [skyfield](https://pypi.org/project/skyfield/), [spacetrack](https://pypi.org/project/spacetrack/), and [hapsira](https://pypi.org/project/hapsira/) (a maintained poliastro fork).
It provides a consistent and straightforward method to create/propagate orbits without worrying about the implementation details of each library:

- historical TLEs for known satellites (skyfield)
- propagated orbital parameters hypothetical satellites (skyfield)
- analytical orbital parameters (hapsira)
- list of positions

Addtionaly, it provides a straightforward interface for updating used TLEs via [spacetrack](https://www.space-track.org)

Created orbits contain commonly used state variables:

- satellite positions in various frames
- velocities
- etc.

## Installation
1. install the package
```bash
pip install spherapy
```
2. configure. see [Configuration](#configuration)

## Usage
- ensure spherapy has been configured, see [Configuration](#configuration)
- import the spherapy package
```python
import spherapy.updater
import spherapy.timespan
import spherapy.orbit
```
- if using real world satellites; update the desired TLEs (you can also use the TLEs supplied with the package, but these will be out of date).
	- TLE id's can be found using [NORAD's CelesTrak catalogue search](https://celestrak.org/satcat/search.php)
```python
updated_TLEs = spherapy.updater.updateTLEs([25544]) 	#ISS
TLE_paths = spherapy.updater.getTLEFilePaths([25544], use_packaged=True) 	#ISS
```
- set up a timespan
```python
t = spherapy.timespan.TimeSpan(datetime.datetime(2024,10,15,0,0,1),'1S','90M')
```
- construct an orbit
	- from a TLE (good practice to update the TLE with the most recent)
```python
o = spherapy.orbit.Orbit.fromTLE(t, TLE_paths[0])
```  
-	from orbital parameters
```python
o = spherapy.orbit.Orbit.fromAnalyticalOrbitalParam(timespan, body='Earth',
					 a=6978,
					 ecc=0,
					 inc=0,
					 raan=0,
					 argp=0,
					 mean_nu=0,
					 name='My Analytical Orbit',
					 astrobodies=True)
```  

#### Full TLE example (for copy paste)
```
import datetime
import spherapy.updater
import spherapy.timespan
import spherapy.orbit
TLE_paths = spherapy.updater.getTLEFilePaths([25544], use_packaged=True) 	#ISS
t = spherapy.timespan.TimeSpan(datetime.datetime(2024,10,15,0,0,1),'1S','90M')
o = spherapy.orbit.Orbit.fromTLE(t, TLE_paths[0])
```

#### Full Analytical exmaple (for copy paste)
```
import datetime
import spherapy.timespan
import spherapy.orbit
t = spherapy.timespan.TimeSpan(datetime.datetime(2024,10,15,0,0,1),'1S','90M')
o = spherapy.orbit.Orbit.fromAnalyticalOrbitalParam(timespan, body='Earth',
					 a=6978,
					 ecc=0,
					 inc=0,
					 raan=0,
					 argp=0,
					 mean_nu=0,
					 name='My Analytical Orbit',
					 astrobodies=True)
```

## SpaceTrack Credentials
In order to calculate the position of a historical satellite at any given time, spherapy requires [TLE information](https://en.wikipedia.org/wiki/Two-line_element_set) for each satellite which is accurate for the given time period (or 'epoch').  
TLE data can be obtained from either [Celestrak](https://celestrak.org/) or [Spacetrack](https://www.space-track.org/). 
Celestrak holds only the most recent TLE data for each satellite, while Spacetrack will provide historical TLE data. spherapy will fall back to using Celestrak if it cannot authenticate access to Spacetrack.  
In order to use Spacetrack, you must provide your [Spacetrack credentials](https://www.space-track.org/auth/createAccount) to spherapy.

## Configuration
Configuration for spherapy can either use the default settings, or settings specified in a `spherapy.conf` file.

### Default settings
If the default settings are used, spherapy will use the system user data directory and expect spacetrack credentials to be stored in the system keyring.  
The default data directories are listed in [Directories](#directories)

### Custom settings
Custom settings can be specified in a `spherapy.conf` file.  
This file can be located either in the system user's config file (described in [Directories](#directories)) or at a location specified by the environment variable `SPHERAPY_CONFIG_DIR`.

The fields of `spherapy.conf` are described below
```ini
[credentials]
SpacetrackUser = None 		# the spacetrack user, if left as "None" will source from system keyring
SpacetrackPasswd = None 	# the spacetrack password, if left as "None" will source from system keyring

[paths]
# relative paths are assumed to be relative to the location of this file
# do not quote path
TLE_path = ./spherapy/data/TLEs 	# the location where TLE files will be saved, if left empty (""), will default to system user's directory.
```

### Directories
If `SPHERAPY_CONFIG_DIR` is not set or `TLE_path` is left empty in `spherapy.conf`, the spherapy default directories will be used:

#### Linux
- TLE data dir:
``` bash
/home/{username}/.local/share/spherapy/TLEs
```
- config dir:
``` bash
/home/{username}/.config/spherapy/spherapy/{spherapy_version}/
```

#### OSX
- TLE data dir:
``` bash
/Users/{username}/Library/Application Support/spherapy/TLEs
```
- config dir:
``` bash
/Users/{username}/Library/Application Support/spherapy/{spherapy_version}/
```

#### Windows
- TLE data dir:
``` bash
'C:\\Users\\{username}\\AppData\\Local\\spherapy\\TLEs'
```
- config dir:
``` bash
'C:\\Users\\{username}\\AppData\\Local\\spherapy\\{spherapy_version}'
```

### Passwords
By default spherapy will use the system keyring to store the Spacetrack username and credentials (see [SpaceTrack Credentials](#spacetrack-credentials))  
If a config file is supplied, the credentials in the file will be used.

Credentials can be added to the keyring by running the command (in the terminal, not the python shell)
```
spherapy-create-credentials
```

### Configration File Format
The configuration file `spherapy.conf` should have the following format:
```ini
[credentials]
SpacetrackUser = None
SpacetrackPasswd = None

[paths]
# relative paths are assumed to be relative to the location of this file
# do not quote path
TLE_path =
```

## Data Storage
- TLEs
	- TLEs will be stored in the data directory, with a single file for each satellite ID.
	 `{sat_id}.tle`, containing all historical TLEs for that satellite.
	- If celestrak is used instead, the file will be saved as a temporary file `{sat_id}.temptle`, which will be overwritten on each fetch from celestrak.