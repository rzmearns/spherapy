<!-- markdownlint-disable -->

<a href="../spherapy/util/spacetrack.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `util.spacetrack`
Functions to fetch TLEs from Spacetrack. 



**Attributes:**
 
 - <b>`MAX_RETRIES`</b>:  number of times to try and reach Spacetrack. 

**Global Variables**
---------------
- **MAX_RETRIES**

---

<a href="../spherapy/util/spacetrack.py#L208"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `updateTLEs`

```python
updateTLEs(
    sat_id_list: list[int],
    user: None | str = None,
    passwd: None | str = None
) → list[int]
```

Fetch most recent TLE for satcat IDs from spacetrack. 

Fetch most recent TLEs for provided list of satcat IDs, and append to file.  
If no TLE file yet exists, will download all historical TLEs Will try MAX_RETRIES before raising a TimeoutError 



**Args:**
 
 - <b>`sat_id_list`</b>:  list of satcat ids to fetch 
 - <b>`user`</b>:  [Optional] overriding spacetrack username to use 
 - <b>`passwd`</b>:  [Optional] overriding spacetrack password to use 



**Returns:**
 
 - <b>`list`</b>: list of satcat ids successfully fetched 



**Raises:**
 TimeoutError 


---

<a href="../spherapy/util/spacetrack.py#L231"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getTLEFilePath`

```python
getTLEFilePath(sat_id: int) → Path
```

Gives path to file where spacetrack TLE is stored. 

Spacetrack TLEs are stored in {satcadID}.tle 



**Args:**
 
 - <b>`sat_id`</b>:  satcat ID 



**Returns:**
 path to file 


---

<a href="../spherapy/util/spacetrack.py#L244"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getStoredEpochs`

```python
getStoredEpochs(sat_id: int) → None | tuple[datetime, datetime | None]
```

Return the start and end epoch for {sat_id}.tle . 



**Args:**
 
 - <b>`sat_id`</b>:  satcat id to check 



**Returns:**
 (first epoch datetime, last epoch datetime) None if no spacetrack tle stored for sat_id 


---

<a href="../spherapy/util/spacetrack.py#L259"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `doCredentialsExist`

```python
doCredentialsExist() → bool
```

Checks if spacetrack credentials have been loaded into Spherapy. 



**Returns:**
  bool: 


---

<a href="../spherapy/util/spacetrack.py#L195"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `InvalidCredentialsError`
Error indicating invalid credentials used to access spacetrack. 

<a href="../spherapy/util/spacetrack.py#L197"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(message: str)
```











---


