<!-- markdownlint-disable -->

<a href="../spherapy/util/epoch_u.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `util.epoch_u`
Utility functions for converting different epochs. 



**Attributes:**
 
 - <b>`GMST_epoch`</b>:  [description] 


---

<a href="../spherapy/util/epoch_u.py#L17"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `epoch2datetime`

```python
epoch2datetime(epoch_str: str) → datetime
```

Converts a fractional epoch string to a datetime object. 



**Args:**
 epoch_str:      fractional year epoch 



**Returns:**
  equivalent datetime 


---

<a href="../spherapy/util/epoch_u.py#L40"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `epochEarlierThan`

```python
epochEarlierThan(epoch_a: str, epoch_b: str) → bool
```

Check if epoch A is earlier than epoch B. 



**Args:**
 
 - <b>`epoch_a`</b>:  TLE epoch string A 
 - <b>`epoch_b`</b>:  TLE epoch string B 



**Returns:**
 True if epoch A is earlier than epoch B 


---

<a href="../spherapy/util/epoch_u.py#L54"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `epochLaterThan`

```python
epochLaterThan(epoch_a: str, epoch_b: str) → bool
```

Check if epoch A is later than epoch B. 



**Args:**
 
 - <b>`epoch_a`</b>:  TLE epoch string A 
 - <b>`epoch_b`</b>:  TLE epoch string B 



**Returns:**
 True if epoch A is later than epoch B 


---

<a href="../spherapy/util/epoch_u.py#L68"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `datetime2TLEepoch`

```python
datetime2TLEepoch(date: datetime) → str
```

Converts a datetime to a TLE epoch string. 



**Args:**
 
 - <b>`date`</b>:  Datetime object 



**Returns:**
 TLE epoch string, with fractional seconds 


---

<a href="../spherapy/util/epoch_u.py#L88"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `datetime2sgp4epoch`

```python
datetime2sgp4epoch(date: datetime) → float
```

Converts a datetime to an sgp4 epoch. 



**Args:**
 
 - <b>`date`</b>:  Datetime object 



**Returns:**
 SGP4 epoch, with fractional seconds 


---

<a href="../spherapy/util/epoch_u.py#L103"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `findClosestDatetimeIndices`

```python
findClosestDatetimeIndices(
    test_arr: ndarray[tuple[int], dtype[datetime64]],
    source_arr: ndarray[tuple[int], dtype[datetime64]]
) → ndarray[tuple[int], dtype[int64]]
```

Find the index of the closest datetime in source arr for each datetime in test_arr. 

Both search_arr and source_arr must be sorted. 



**Args:**
 
 - <b>`test_arr`</b>:  Mx1 array of datetimes, will find closest time for each element in this array 
 - <b>`source_arr`</b>:  Nx1: array of datetimes to compare to 



**Returns:**
 Mx1 array of indices of source_arr 


---

<a href="../spherapy/util/epoch_u.py#L127"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getStoredEpochs`

```python
getStoredEpochs(tle_path: Path) → None | tuple[datetime, datetime | None]
```

Return the start and end epoch for tle_path. 



**Args:**
 
 - <b>`tle_path`</b>:  tle file 



**Returns:**
 (first epoch datetime, last epoch datetime) None if no spacetrack tle stored for sat_id 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
