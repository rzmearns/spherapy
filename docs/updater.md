<!-- markdownlint-disable -->

<a href="../spherapy/updater.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `updater`
Functions to update TLE of a satcat ID, and fetch where the file is stored. 


---

<a href="../spherapy/updater.py#L10"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `updateTLEs`

```python
updateTLEs(sat_id_list: list[int]) → list[int]
```

Fetch most recent TLE for provided list of satcat IDs. 

Fetch most recent TLE for provided list of satcat IDs,  
If credentials are stored, use Spacetrack as TLE source, this will fetch all historical TLEs for that satcat ID (from most recent stored TLE up to current)  
If no credentials are stored, use Celestrak as TLE source, this will fetch only the most  recent TLE 



**Args:**
 
 - <b>`sat_id_list`</b>:  list of satcat ids to fetch 



**Returns:**
 list of satcat ids successfully fetched 


---

<a href="../spherapy/updater.py#L32"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getTLEFilePaths`

```python
getTLEFilePaths(sat_id_list: list[int], use_packaged: bool = False) → list[Path]
```

Fetch list of paths to TLE files. 

Fetch paths to the file storing the list of TLEs for each provided satcat ID  
If credentials are stored, return path containing all historical TLEs (fetched from Spacetrack)  
If no credentials are stored, return path to .temptle containing only most recent TLE  (from Celestrak)  
You can force the use of TLE data packaged with spherapy by setting use_packaged=True  This data will be out of date, and should only be used in examples. 



**Args:**
 
 - <b>`sat_id_list`</b>:  list of satcat ids to find paths 
 - <b>`use_packaged`</b>:  use the default TLEs packaged with spherapy, these will be out of date. 



**Returns:**
 list of paths 


---

<a href="../spherapy/updater.py#L66"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getStoredEpochLimits`

```python
getStoredEpochLimits(
    sat_id_list: list[int],
    use_packaged: bool = False
) → dict[int, None | tuple[datetime, datetime | None]]
```

Returns limiting epochs for the stored TLEs for each sat in sat_id_list. 

[description] 



**Args:**
 
 - <b>`sat_id_list`</b> (list[int]):  [description] 
 - <b>`use_packaged`</b> (bool):  [description] (default: `False`) 



**Returns:**
 
 - <b>`list[tuple[dt.datetime, dt.datetime]]`</b>:  [description] 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
