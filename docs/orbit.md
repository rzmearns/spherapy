<!-- markdownlint-disable -->

<a href="../spherapy/orbit.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `orbit`
Class for orbital data. 

This module provides: 
- OrbitAttrDict: A TypedDict typing the instance attributes of an Orbit 
- Orbit: A class of timestamped orbital data 

**Global Variables**
---------------
- **WGS72**


---

<a href="../spherapy/orbit.py#L36"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `OrbitAttrDict`
A TypedDict providing type annotations for the Orbit class. 



**Attributes:**
  name:  satcat_id:  gen_type:  timespan:  TLE_epochs:  pos:  pos_ecef:  vel_ecef:  vel:  lat:  lon:  sun_pos:  moon_pos:  alt:  eclipse:  central_body:  period:  period_steps:  semi_major:  ecc:  inc:  raan:  argp: 





---

<a href="../spherapy/orbit.py#L159"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Orbit`
Timestamped orbital data for a satellite, coordinate system depending on the central body. 



**Attributes:**
 
 - <b>`timespan`</b>:  Timespan over which orbit is to be simulated 
 - <b>`name`</b>:  Name of the satellite 
 - <b>`satcat_id`</b>:  NORAD satellite catelogue ID (if generated from a TLE) 
 - <b>`gen_type`</b>:  How the orbit was generated 
 - <b>`TLE_epochs`</b>:  Nx1 numpy array of TLE epoch used for propagation at each timestamp 
   - units:  TLE epoch 
 - <b>`pos`</b>:  Nx3 numpy array of cartesian coordinates of the position of the satellite  at each timestamp 
   - units:  km 
   - frame:  ECI 
 - <b>`vel`</b>:  Nx3 numpy array of cartesian velocities of the satellite at each timestamp 
   - units:  m/s 
   - frame:  ECI 
 - <b>`pos_ecef`</b>:  Nx3 numpy array of cartesian coordinates of the position of the satellite  at each timestamp 
   - units:  km 
   - frame:  ECEF 
 - <b>`vel_ecef`</b>:  Nx3 numpy array of cartesian velocities of the satellite at each timestamp 
   - units:  m/s 
   - frame:  ECEF 
 - <b>`lat`</b>:  Nx1 numpy array of central body latitudes of the satellite at each timestamp 
   - units:  degrees 
 - <b>`lon`</b>:  Nx1 numpy array of central body longitudes of the satellite at each timestamp 
   - units:  degrees 
 - <b>`sun_pos`</b>:  Nx3 numpy array of cartesian coordinates of the position of the Sun  at each timestamp 
   - units:  km 
   - frame:  ECI 
 - <b>`moon_pos`</b>:  Nx3 numpy array of cartesian coordinates of the position of the Moon  at each timestamp 
   - units:  km 
   - frame:  ECI 
 - <b>`alt`</b>:  Nx1 numpy array of altitudes above central body at each timestamp 
   - units:  km 
 - <b>`eclipse`</b>:  Nx1 numpy array of flag indicating if satellite is eclipsed at each timestamp 
 - <b>`central_body`</b>:  body the satellite is orbiting 
 - <b>`period`</b>:  orbital period in secs 
 - <b>`period_steps`</b>:  number of TimeSpan timestamps required to complete an orbit 
 - <b>`semi_major`</b>:  Nx1 numpy array of orbit semi-major axis calculated at that timestep 
   - units:  km will be constant if no orbital maneauvers 
 - <b>`ecc`</b>:  Nx1 numpy array of orbit eccentricity calculated at that timestep 
   - units:  unitless will be constant if no orbital maneauvers 
 - <b>`inc`</b>:  Nx1 numpy array of orbit inclination calculated at that timestep 
   - units:  degree will be constant if no orbital maneauvers 
 - <b>`raan`</b>:  Nx1 numpy array of orbit RAAN calculated at that timestep 
   - units:  degree will be constant if no orbital maneauvers 
 - <b>`argp`</b>:  Nx1 numpy array of orbit Arg Perigee calculated at that timestep 
   - units:  degree will be constant if no orbital maneauvers 

<a href="../spherapy/orbit.py#L219"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(data: OrbitAttrDict, calc_astrobodies: bool = False)
```

The constructor should never be called directly. 

Use one of:  
 - Orbit.fromTLE()
 - Orbit.fromListOfPositions()
 - Orbit.fromPropagatedOrbitalParam()
 - Orbit.fromAnalyticalOrbitalParam() 




---

<a href="../spherapy/orbit.py#L621"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>classmethod</kbd> `fromAnalyticalOrbitalParam`

```python
fromAnalyticalOrbitalParam(
    timespan: TimeSpan,
    body: str = 'Earth',
    a: float = 6978,
    ecc: float = 0,
    inc: float = 0,
    raan: float = 0,
    argp: float = 0,
    mean_nu: float = 0,
    name: str = 'Analytical',
    astrobodies: bool = True,
    unsafe: bool = False
) → Orbit
```

Create an analytical orbit defined by orbital parameters. 

Orbits created using this class method will NOT respect gravity corrections such as J4, and as such sun-synchronous orbit are not possible. 



**Args:**
 
 - <b>`timespan`</b>:  Timespan over which orbit is to be simulated 
 - <b>`body`</b>:  [Optional] string indicating around what body the satellite orbits,  Default is 'Earth'  Options are ['Earth','Sun','Mars','Moon'] 
 - <b>`a`</b>:  [Optional] semi-major axis of the orbit in km  Default is 6978 ~ 600km above the earth. 
 - <b>`ecc`</b>:  [Optional] eccentricty, dimensionless number 0 < ecc < 1  Default is 0, which is a circular orbit 
 - <b>`inc`</b>:  [Optional] inclination of orbit in degrees  Default is 0, which represents an orbit around the Earth's equator 
 - <b>`raan`</b>:  [Optional] right-ascension of the ascending node  Default is 0 
 - <b>`argp`</b>:  [Optional] argument of the perigee in degrees  Default is 0, which     represents an orbit with its semimajor axis in  the plane of the Earth's equator 
 - <b>`mean_nu`</b>:  [Optional] mean anomaly in degrees  Default is 0, which represents an orbit that is beginning at periapsis 
 - <b>`name`</b>:  [Optional] string giving the name of the orbit  Default is 'Analytical' 
 - <b>`astrobodies`</b>:  [Optional] Flag to calculate Sun and Moon positions at timestamps  Default is False 
 - <b>`unsafe`</b>:  [Optional] Flag to ignore semi-major axis inside Earth's radius  Default is False 



**Returns:**
 satplot.Orbit 

---

<a href="../spherapy/orbit.py#L755"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>classmethod</kbd> `fromDummyConstantPosition`

```python
fromDummyConstantPosition(
    timespan: TimeSpan,
    pos: tuple[float, float, float] | ndarray[tuple[int], dtype[float64]],
    sun_pos: tuple[float, float, float] | ndarray[tuple[int], dtype[float64]] | None = None,
    moon_pos: tuple[float, float, float] | ndarray[tuple[int], dtype[float64]] | None = None
) → Orbit
```

Creates a static orbit for testing. 

Satellite position is defined by pos, while sun and moon positions are optional, but can also be specified. 



**Args:**
 
 - <b>`timespan`</b>:  Timespan over which orbit is to be simulated 
 - <b>`pos`</b>:  Nx3 numpy array of cartesian coordinates of the position of the satellite  at each timestamp 
   - units:  km 
   - frame:  ECI 
 - <b>`sun_pos`</b>:  [Optional] Nx3 numpy array of cartesian coordinates of the position of the Sun  at each timestamp 
   - units:  km 
   - frame:  ECI 
 - <b>`moon_pos`</b>:  [Optional] Nx3 numpy array of cartesian coordinates of the position of the  Moon at each timestamp 
   - units:  km 
   - frame:  ECI 



**Returns:**
 satplot.Orbit 

---

<a href="../spherapy/orbit.py#L298"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>classmethod</kbd> `fromListOfPositions`

```python
fromListOfPositions(
    timespan: TimeSpan,
    positions: ndarray[tuple[int, int], dtype[float64]],
    astrobodies: bool = False
) → Orbit
```

Create an orbit from a list of positions. 

Creat an obit by explicitly specifying the position of the satellite at each point in time. Useful for simplified test cases; but may lead to unphysical orbits. 



**Args:**
 
 - <b>`timespan`</b>:  Timespan over which orbit is to be simulated 
 - <b>`positions`</b>:  Nx3 numpy array of cartesian coordinates of the position of the satellite  at each timestamp 
   - units:  km 
   - frame:  ECI 
 - <b>`astrobodies`</b>:  [Optional] Flag to calculate Sun and Moon positions at timestamps  Default is False 



**Returns:**
 satplot.Orbit 

---

<a href="../spherapy/orbit.py#L489"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>classmethod</kbd> `fromPropagatedOrbitalParam`

```python
fromPropagatedOrbitalParam(
    timespan: TimeSpan,
    a: float = 6978,
    ecc: float = 0,
    inc: float = 0,
    raan: float = 0,
    argp: float = 0,
    mean_nu: float = 0,
    name: str = 'Fake TLE',
    astrobodies: bool = True,
    unsafe: bool = False
) → Orbit
```

Create an orbit from orbital parameters, propagated using sgp4. 

Orbits created using this class method will respect gravity corrections such as J4, allowing for semi-analytical sun-synchronous orbits. 



**Args:**
 
 - <b>`timespan`</b>:  Timespan over which orbit is to be simulated 
 - <b>`a`</b>:  [Optional] semi-major axis of the orbit in km  Default is 6978 ~ 600km above the earth. 
 - <b>`ecc`</b>:  [Optional] eccentricty, dimensionless number 0 < ecc < 1  Default is 0, which is a circular orbit 
 - <b>`inc`</b>:  [Optional] inclination of orbit in degrees  Default is 0, which represents an orbit around the Earth's equator 
 - <b>`raan`</b>:  [Optional] right-ascension of the ascending node  Default is 0 
 - <b>`argp`</b>:  [Optional] argument of the perigee in degrees  Default is 0, which     represents an orbit with its semimajor axis in  the plane of the Earth's equator 
 - <b>`mean_nu`</b>:  [Optional] mean anomaly in degrees  Default is 0, which represents an orbit that is beginning at periapsis 
 - <b>`name`</b>:  [Optional] string giving the name of the orbit  Default is 'Fake TLE' 
 - <b>`astrobodies`</b>:  [Optional] Flag to calculate Sun and Moon positions at timestamps  Default is False 
 - <b>`unsafe`</b>:  [Optional] Flag to ignore semi-major axis inside Earth's radius  Default is False 



**Returns:**
 satplot.Orbit 

---

<a href="../spherapy/orbit.py#L337"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>classmethod</kbd> `fromTLE`

```python
fromTLE(
    timespan: TimeSpan,
    tle_path: Path,
    astrobodies: bool = True,
    unsafe: bool = False
) → Orbit
```

Create an orbit from an existing TLE or a list of historical TLEs. 



**Args:**
 
 - <b>`timespan `</b>:  TimeSpan over which orbit is to be simulated 
 - <b>`tle_path `</b>:  path to file containing TLEs for a satellite 
 - <b>`astrobodies`</b>:  [Optional] Flag to calculate Sun and Moon positions at timestamps  Default is False 
 - <b>`unsafe`</b>:  [Optional] Flag to ignore TLE usage more than 14 days either side of timestamps  Optional  Default is False 



**Returns:**
 satplot.Orbit 

---

<a href="../spherapy/orbit.py#L796"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getPosition`

```python
getPosition(search_time: datetime | Time) → ndarray[tuple[int], dtype[float64]]
```

Return the position at the specified or closest time. 



**Args:**
 
 - <b>`search_time`</b>:  the timestamp to search for 



**Returns:**
 position 



**Raises:**
 
 - <b>`ValueError`</b>:  orbit has no pos data 

---

<a href="../spherapy/orbit.py#L813"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getVelocity`

```python
getVelocity(search_time: datetime | Time) → ndarray[tuple[int], dtype[float64]]
```

Return the velocity at the specified or closest time. 



**Args:**
 
 - <b>`search_time`</b>:  the timestamp to search for 



**Returns:**
 velocity 



**Raises:**
 
 - <b>`ValueError`</b>:  orbit has no vel data 




---
