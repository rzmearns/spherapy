## Orbit Objects

```python
class Orbit(object)
```

Contains timestepped array of orbital position and velocity data for a satellite, coordinate system depending on the central body.
Contains timestepped array of sun position.
Timing data taken from a TimeSpan object.

**Attributes**  
timespan: {satplot.TimeSpan}  
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
*Recommended: 10 < period_steps < 100*



### \_\_init\_\_

```python
def __init__(*args, **kwargs)
```

The constructor should never be called directly.
Use one of:
- Orbit.fromListOfPositions()
- Orbit.fromTLE()
- Orbit.multiFromTLE()
- Orbit.FromTLEOrbitalParam()
- Orbit.fromOrbitalParam()
- Orbit.load()

### fromListOfPositions
```python
@classmethod
def fromListOfPositions(cls, timespan, positions, astrobodies=True)
```

Create an orbit by explicitly specifying the position of the
satellite at each point in time. Useful for simplified test cases; but
may lead to unphysical orbits.

**Parameters**  
timespan : {satplot.TimeSpan}  
Timespan over which orbit is to be simulated

positions: (N,3) np.array  
Position of the satellite in GCRS at each point in time.

**Returns**  
satplot.Orbit



### fromTLE

```python
@classmethod
def fromTLE(cls, timespan, tle_path, fp=sys.stdout, astrobodies=True)
```

Create an orbit from an existing TLE or a list of historical TLEs

**Parameters**  
timespan : {satplot.TimeSpan}  
	Timespan over which orbit is to be simulated

tle_path : {path to file containing TLE(s)}  
	A plain text file containing the TLE or list of TLE(s) with each TLE line on a new line in the file.

**Returns**  
satplot.Orbit

### multiFromTLE

```python
@classmethod
def multiFromTLE(cls, timespan, tle_path, fp=sys.stdout, astrobodies=True)
```

Create an orbit from an existing TLE or a list of historical TLEs

**Parameters**  
timespan : {satplot.TimeSpan}  
Timespan over which orbit is to be simulated

tle_path : {path to file containing TLE(s)}  
A plain text file containing the TLE or list of TLE(s) with each TLE line on a new line in the file.

**Returns**  
satplot.Orbit


### fromTLEOrbitalParam

```python
@classmethod
def fromTLEOrbitalParam(cls,
                        timespan,
                        a=6978,
                        ecc=0,
                        inc=0,
                        raan=0,
                        argp=0,
                        mean_nu=0,
                        name='Fake TLE',
                        astrobodies=True)
```

Create an orbit from orbital parameters, propagated using sgp4.

Orbits created using this class method will respect gravity corrections such as J4, 
allowing for semi-analytical sun-synchronous orbits.

**Parameters**  
timespan : {satplot.TimeSpan}  
	Timespan over which orbit is to be simulated  

a : {float}, optional  
	semi-major axis of the orbit in km (the default is 6978 ~ 600km 
	above the earth.)  

ecc : {float}, optional  
	dimensionless number 0 < ecc < 1 (the default is 0, which is a 	circular orbit)

inc : {float}, optional  
	inclination of orbit in degrees (the default is 0, which represents 
	an orbit around the Earth's equator)

raan : {float}, optional  
	right-ascension of the ascending node (the default is 0)

argp : {float}, optional  
	argument of the perigee in degrees (the default is 0, which 
	represents an orbit with its semimajor axis in the plane of the 
	Earth's equator)

mean_nu : {float}, optional  
	mean anomaly in degrees (the default is 0, which represents an orbit
	that is beginning at periapsis)

**Returns**  
satplot.Orbit



### fromOrbitalParam

```python
@classmethod
def fromOrbitalParam(cls,
                     timespan,
                     body='Earth',
                     a=6978,
                     ecc=0,
                     inc=0,
                     raan=0,
                     argp=0,
                     mean_nu=0,
                     name='Analytical',
                     astrobodies=True)
```

Create an orbit from orbital parameters

**Parameters**  
timespan : {satplot.TimeSpan}  
	Timespan over which orbit is to be simulated

body : {str}, optional  
	Central body for the satellite: ['Earth', 'Moon', 'Mars', 'Sun'] (the default is 'Earth')

a : {float}, optional  
	semi-major axis of the orbit in km (the default is 6978 ~ 600km above the earth.)

ecc : {float}, optional  
	dimensionless number 0 < ecc < 1 (the default is 0, which is a circular orbit)

inc : {float}, optional  
	inclination of orbit in degrees (the default is 0, which represents 
	an orbit around the Earth's equator)

raan : {float}, optional  
	right-ascension of the ascending node (the default is 0)

argp : {float}, optional  
	argument of the perigee in degrees (the default is 0, which 
	represents an orbit with its semimajor axis in the plane of the 
	Earth's equator)

mean_nu : {float}, optional  
	mean anomaly in degrees (the default is 0, which represents an orbit
	that is beginning at periapsis)

**Returns**  
satplot.Orbit



### load

```python
@classmethod
def load(cls, file)
```



### getPosition

```python
def getPosition(time)
```

Return the position at the specified index or (closest) time

**Parameters**  
time : {int|datetime.datetime|astropy.Time}  
	The index of the timespan, or a particular time to fetch the position  
	If a particular time is specified, the nearest timestep within the timespan will be returned

**Returns**  
(3,) ndarray  
	Orbital position (km)

**Raises**  
ValueError  
	Raises a ValueError if 'time' is not an integer or a time type



### getVelocity

```python
def getVelocity(time)
```

Return the position at the specified index or (closest) time

**Parameters**  
time : {int|datetime.datetime|astropy.Time}  
	The index of the timespan, or a particular time to fetch the velocity  
	If a particular time is specified, the nearest timestep within the timespan will be returned

**Returns**  
(3,) ndarray  
	Orbital velocity (m/s)

**Raises**  
ValueError  
	Raises a ValueError if 'time' is not an integer or a time type

