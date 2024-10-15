
## TimeSpan Objects

```python
class TimeSpan(object)
```

### \_\_init\_\_

```python
def __init__(t0, timestep='1S', timeperiod='10S', timezone='00:00')
```

Create a TimeSpan object, times are assumed to be in UTC

An array of dates as datetimes is accessed with TimeSpan.as_datetime()  
An array of dates as astropy.time is accessed with TimeSpan.as_astropy()

Does not handle leap seconds

**Parameters**  
t0 : {datetime.datetime}  
	datetime object defining the start of the TimeSpan

num_steps: int  
	Number of timesteps

timestep : {str}, optional  
	String describing the time step of the time span. The string is constructed 
	as an integer or float, followed by a time unit: (d)ays,
	(H)ours, (M)inutes, (S)econds, (mS) milliseconds, (uS) microseconds
	(the default is '1S', which is one second.)

timeperiod : {str}, optional  
	String describing the time period of the time span. The string is constructed 
	as an integer or float, followed by a time unit: (y)ears, (m)onths, (W)eeks, (d)ays,
	(H)ours, (M)inutes, (S)econds, (mS) milliseconds, (uS) microseconds
	(the default is '1d', which is one day.)

timezone : {datetime.timezone}, optional  
	Timezone to be implicitly assumed for the timespan. Methods can check against 
	TimeSpan.timezone for instances.

**Raises**  
ValueError

### asAstropy

```python
def asAstropy(*args, scale='utc')
```

Return ndarray of TimeSpan as astropy.time objects

**Returns**  
ndarray


### asDatetime

```python
def asDatetime(*args)
```

Return ndarray of TimeSpan as datetime objects	

**Returns**  
ndarray



### asSkyfield

```python
def asSkyfield(*args)
```

Return TimeSpan element as Skyfield Time object

**Returns**  
Skyfield Time



### asText

```python
def asText(*args)
```



### secondsSinceStart

```python
def secondsSinceStart()
```

Return ndarray with the seconds of all timesteps since the beginning.



### getClosest

```python
def getClosest(t_search)
```

Find the closest time in a TimeSpan

**Parameters**  
t_search : {datetime}
	time to find

**Returns**  
datetime, int  
	Closest datetime in TimeSpan, index of closest date in TimeSpan

