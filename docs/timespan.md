<!-- markdownlint-disable -->

<a href="../spherapy/timespan.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `timespan`
Class for series of timestamps. 

This module provides: 
- TimeSpan: a series of timestamps to be used by an orbit object 



---

<a href="../spherapy/timespan.py#L21"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `TimeSpan`
A series of timestamps. 



**Attributes:**
 
 - <b>`start`</b>:  The first timestamp 
 - <b>`end`</b>:  The last timestamp 
 - <b>`time_step`</b>:  The difference between timestamps in seconds.
     - Can be None if irregular steps 
 - <b>`time_period`</b>:  The difference between end and start in seconds 

<a href="../spherapy/timespan.py#L32"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(t0: datetime, timestep: str = '1S', timeperiod: str = '10S')
```

Creates a series of timestamps in UTC. 

Difference between each timestamp = timestep Total duration = greatest integer number of timesteps less than timeperiod  
If timeperiod is an integer multiple of timestep,  then TimeSpan[-1] - TimeSpan[0] = timeperiod  
If timeperiod is NOT an integer multiple of timestep,  then TimeSpan[-1] = int(timeperiod/timestep) (Note: int cast, not round) 

Does not account for Leap seconds, similar to all Posix compliant UTC based time  representations. see: https://numpy.org/doc/stable/reference/arrays.datetime.html#datetime64-shortcomings  for equivalent shortcomings.  
Always contains at least two timestamps 

**Args**
 - t0: datetime defining the start of the TimeSpan.
     - If timezone naive, assumed to be in UTC
     - If timezone aware, will be converted to UTC
 - timestep: String describing the time step of the time span.
     - The string is constructed as an integer or float, followed by a time unit:  (d)ays, (H)ours, (M)inutes, (S)econds, (mS) milliseconds, (uS) microseconds  (the default is '1S') 
 - timeperiod: String describing the time period of the time span.
     - The string is constructed as an integer or float, followed by a time unit:  (d)ays, (H)ours, (M)inutes, (S)econds, (mS) milliseconds, (uS) microseconds  (the default is '1d') 



**Raises:**
 
ValueError 




---

<a href="../spherapy/timespan.py#L153"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `asAstropy`

```python
asAstropy(idx: None | int = None, scale: str = 'utc') → Time
```

Return ndarray of TimeSpan as astropy.time objects. 



**Args:**
 
 - <b>`idx`</b>:  timestamp index to return inside an astropyTime object  if no index supplied, returns all timestamps 
 - <b>`scale`</b>:  astropy time scale, can be one of  ('tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1', 'utc') 



**Returns:**
 
ndarray 

---

<a href="../spherapy/timespan.py#L170"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `asDatetime`

```python
asDatetime(
    idx: None | int = None
) → datetime | ndarray[tuple[int], dtype[datetime64]]
```

Return ndarray of TimeSpan as datetime objects. 



**Args:**
 
 - <b>`idx`</b>:  timestamp index to return as a datetime  If no index supplied, returns whole array 



**Returns:**
 
ndarray 

---

<a href="../spherapy/timespan.py#L186"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `asSkyfield`

```python
asSkyfield(idx: int) → Time
```

Return TimeSpan element as Skyfield Time object. 



**Args:**
 
 - <b>`idx`</b>:  timestamp index to return as skyfield time object 



**Returns:**
 
Skyfield Time 

---

<a href="../spherapy/timespan.py#L200"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `asText`

```python
asText(idx: int) → str
```

Returns a text representation of a particular timestamp. 

Timestamp will be formatted as YYYY-mm-dd HH:MM:SS 



**Args:**
 
 - <b>`idx`</b>:  timestamp index to format 



**Returns:**

 str: 

---

<a href="../spherapy/timespan.py#L321"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `cherryPickFromIndices`

```python
cherryPickFromIndices(idxs: int | tuple | slice)
```

Adjust TimeSpan to only contain the indices specified by idxs. 



**Args:**
 
 - <b>`idxs`</b>:  numpy style indexing of TimeSpan 

---

<a href="../spherapy/timespan.py#L333"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>classmethod</kbd> `fromDatetime`

```python
fromDatetime(
    dt_arr: ndarray[tuple[int], dtype[datetime64]],
    timezone: timezone = datetime.timezone.utc
) → TimeSpan
```

Create a TimeSpan from an array of datetime objects. 



**Args:**
 
 - <b>`dt_arr`</b>:  1D array of datetime objects 
 - <b>`timezone`</b>:  timezone to apply to each element of datetime array. 
 - <b>`Default`</b>:  dt.timezone.utc 

---

<a href="../spherapy/timespan.py#L225"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getClosest`

```python
getClosest(t_search: datetime) → tuple[datetime, int]
```

Find the closest time in a TimeSpan. 

Parameters 
---------- t_search : datetime to search for in TimeSpan  If timezone naive, assumed to be in UTC  If timezone aware, will be converted to UTCtime to find 



**Returns:**
 
------- datetime, int  Closest datetime in TimeSpan, index of closest date in TimeSpan 

---

<a href="../spherapy/timespan.py#L225"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `areTimesWithin`

```python
areTimesWithin(t_search:dt.datetime|np.ndarray[tuple[int], np.dtype[np.datetime64]]) → np.ndarray[tuple[int],np.dtype[np.bool_]]
```

Find if the provided times are within the timespan.

**Args**: 
 - <b>`t_search`</b>: times to check if within timespan
                        If timezone naive, assumed to be in UTC
                        If timezone aware, will be converted to UTCtime to find


**Returns:**

 ndarray of bools, True if within timespan


---

<a href="../spherapy/timespan.py#L225"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getFractionalIndices`

```python
getFractionalIndices(t_search:dt.datetime|np.ndarray[tuple[int],np.dtype[np.datetime64]]) → np.ndarray[tuple[int],np.dtype[np.float64]]
```

Find the fractional indices of timespan at wich t_search should be inserted.

Find the indices in the original timespan at which each value of t_search should be
        inserted to maintain the sorted order.
        The integer part of the index indicates the value immediately prior to the value in
        t_search, while the fractional part represents the point between the two adjacent indices
        in the timespan at which the t_search value falls.
        For example (using integers rather than datetime objects:
            timespan = [0, 1, 5, 6, 10]
            t_search = [0.5, 2, 3, 8]
            timespan.getFractionalIndices(t_search) = [0.5, 1.25, 1.5, 3.5]
        All values of t_search must be within the timespan, otherwise the output is undefined.


**Args:**
- <b>`t_search`</b>: times to locate within timespan
                If timezone naive, assumed to be in UTC
                If timezone aware, will be converted to UTCtime to find



**Returns:**
 
 ndarray of fractional indices

---

<a href="../spherapy/timespan.py#L213"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `secondsSinceStart`

```python
secondsSinceStart() → ndarray[tuple[int], dtype[float64]]
```

Return ndarray with the seconds of all timesteps since the beginning. 



**Returns:**
  array of seconds since start for each timestamp 




---

