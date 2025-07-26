## About
Spherapy is a convenience wrapper around the fantastic libraries [skyfield](https://pypi.org/project/skyfield/), [spacetrack](https://pypi.org/project/spacetrack/), and [hapsira](https://pypi.org/project/hapsira/) (a maintained poliastro fork). It provides a consistent and straightforward method to fetch historical TLEs for known satellites, update these TLEs as required, propagate the orbits, and turn these into commonly used variables: satellite positions in various frames, velocities, etc. without worrying about the implementation details of each library.  
Additionally, rather than use a historical TLE, orbits can be constructed from analytical orbital parameters, propagated orbital parameters, or a list of positions.

## Installation
- install to use the package:
  1. install from git
```bash
pip install git+ssh://git@gitlab.unimelb.edu.au/msl/libraries/spherapy.git
```
  2. configure. see [Configuration](#configuration)
- install to develop the package:
  1. clone the repo
  2. install as editable
```bash
pip install --editable .[dev]
```
  3. contribute, see [Contributing](https://gitlab.unimelb.edu.au/msl/libraries/spherapy/-/blob/main/CONTRIBUTING.md)
 
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
In order to calculate the position of a satellite at any given time, spherapy requires [TLE information](https://en.wikipedia.org/wiki/Two-line_element_set) for each satellite which is accurate for the given time period (or 'epoch').  
TLE data can be obtained from either [Celestrak](https://celestrak.org/) or [Spacetrack](https://www.space-track.org/). 
Celestrak holds only the most recent TLE data for each satellite, while Spacetrack will provide historical TLE data. spherapy will fall back to using Celestrak if it cannot authenticate access to Spacetrack.  
In order to use Spacetrack, you must provide your [Spacetrack credentials](https://www.space-track.org/auth/createAccount) to spherapy.

## Configuration
Configuration for spherapy can either use a supplied configuration file, or the system's user data and config directories (default).  
To use a config file, set the environment variable 'SPHERAPY_CONFIG_DIR' to the location of the `spherapy.conf` file.

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
The configuration file should have the following format:
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

## Timespan Object
The timespan object is the base time class of spherapy.

## Orbit Object
Please see [the orbit documentation](docs/orbit.md)

## Timespan Object
Please see [the timespan documentation](docs/timespan.md)
