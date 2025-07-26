<!-- markdownlint-disable -->

<a href="../spherapy/util/celestrak.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `util.celestrak`
Functions to fetch TLEs from Celestrak. 



**Attributes:**
 
 - <b>`MAX_RETRIES`</b>:  number of times to try and reach Celestrak. 
 - <b>`TIMEOUT`</b>:  timeout for connection and request servicing by Celestrak. 

**Global Variables**
---------------
- **MAX_RETRIES**
- **TIMEOUT**

---

<a href="../spherapy/util/celestrak.py#L22"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `updateTLEs`

```python
updateTLEs(sat_id_list: list[int]) → list[int]
```

Fetch most recent TLE for satcat IDs from celestrak. 

Fetch most recent TLE for provided list of satcat IDs, and store in file.  
Will try MAX_RETRIES before raising a ValueError 



**Args:**
 
 - <b>`sat_id_list`</b>:  list of satcat ids to fetch 



**Returns:**
 
 - <b>`list`</b>:  list of satcat ids successfully fetched 



**Raises:**
 TimeoutError 


---

<a href="../spherapy/util/celestrak.py#L63"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getTLEFilePath`

```python
getTLEFilePath(sat_id: int) → Path
```

Gives path to file where celestrak TLE is stored. 

Celestrak TLEs are stored in {satcadID}.temptle 



**Args:**
 
 - <b>`sat_id`</b>:  satcat ID 



**Returns:**
 path to file 


---

<a href="../spherapy/util/celestrak.py#L76"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getStoredEpochs`

```python
getStoredEpochs(sat_id: int) → None | tuple[datetime, datetime | None]
```

Return the start and end epoch for {sat_id}.temptle . 



**Args:**
 
 - <b>`sat_id`</b>:  satcat id to check 



**Returns:**
 (first epoch datetime, last epoch datetime) None if no spacetrack tle stored for sat_id 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
