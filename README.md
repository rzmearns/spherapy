## Installation
This package can be installed via pip directly from the git repository via:
```bash
pip install git+ssh://git@gitlab.unimelb.edu.au/msl/libraries/spherapy.git
```

## Usage
### As a submodule
- To be investigated

### As a stand alone import
- ensure the [configuration file](spheraphy.conf) has been updated with any desired paths
- ensure that spacetrack credentials are saved in the relevant credential file
	- otherwise celestrak will be used as a fallback, but historical TLEs will be unavailable
- import the spherapy package
```python
import spherapy.updater
import spherapy.timespan
import spherapy.orbit
```
- if using real satellites; update the desired TLEs
	- TLE id's can be found using [NORAD's CelesTrak catalogue search](https://celestrak.org/satcat/search.php)
```python
updated_TLEs = spherapy.updater.updateTLEs([58468])
TLE_paths = spherapy.updater.getTLEFilePaths([58468])
```
- set up a timespan
```python
t = spherapy.timespan.TimeSpan(dt.datetime(2024,10,15,0,0,1),'1S','90M')
```
- construct an orbit
	- from a TLE (good practice to update the TLE with the most recent)
```python
o = spherapy.orbit.Orbit.fromTLE(t, TLE_paths[0])
```  
-	- from orbital parameters
```python
o = spherapy.orbit.Orbit.fromAnalyticalOrbitalParam(timespan, body='Earth', a=6978, ecc=0, inc=0, raan=0, argp=0, mean_nu=0, name='Analytical', astrobodies=True)
```  

## SpaceTrack Credentials
In order to calculate the position of a satellite at any given time, Satplot requires [TLE information](https://en.wikipedia.org/wiki/Two-line_element_set) for each satellite which is accurate for the given time period.  
TLE data can be obtained from either [Celestrak](https://celestrak.org/) or [Spacetrack](https://www.space-track.org/). 
Celestrak holds only the most recent TLE data for each satellite, while Spacetrack will provide historical TLE data. Satplot will fall back to using Celestrak if it cannot authenticate access to Spacetrack.  
In order to use Spacetrack, you must provide your [Spacetrack credentials](https://www.space-track.org/auth/createAccount) to spherapy.

- use `spherapy.updater.createCredentials()` to create credentials in the configured location, this will overwrite any existing credentials

## Configuration
The configuration for spheraphy is set in an `.ini` style configuration file; `spherapy.conf`.
- All paths should be relative to the spherapy root directory.
- DO NOT EDIT the DEFAULT section of the file

## Data Storage
- TLEs
	- TLEs will be stored in the data directory (specified in the [configuration file](spheraphy.conf)), with a single file for each satellite ID.
	 `{sat_id}.tle`, containing all historical TLEs for that satellite.
	- If celestrak is used instead, the file will be saved as a temporary file `{sat_id}.temptle`, which will be overwritten on each fetch from celestrak.
- Orbits
	- The orbit object is currently not serialisable, but this is a high priority area of work to prevent needing to propagate large data sets each time it is run (investigate diskcache package?)

## Timespan Object
The timespan object is the base time class of spheraphy.


## Orbit Object
Please see [the orbit documentation](docs/orbit.md)

## Timespan Object
Please see [the timespan documentation](docs/timespan.md)
