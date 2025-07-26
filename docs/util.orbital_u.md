<!-- markdownlint-disable -->

<a href="../spherapy/util/orbital_u.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `util.orbital_u`
Utility functions for orbital calculations. 


---

<a href="../spherapy/util/orbital_u.py#L7"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `ssoInc`

```python
ssoInc(alt: float, e: float = 0) → float
```

Generates required inclination for given altitude [km] to maintain Sun Syncrhonous orbit. 



**Args:**
 
 - <b>`alt`</b>:  altitude of orbit in km 
 - <b>`e`</b>:  [Optional] eccentricity of orbit  Default is circular (0) 



**Returns:**
 Inclination angle in degrees 


---

<a href="../spherapy/util/orbital_u.py#L38"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calcPeriod`

```python
calcPeriod(a: float) → float
```

Returns the period of an elliptical or circular orbit. 



**Args:**
 
 - <b>`a`</b>:  semi-major axis in m 



**Returns:**
 Orbital period in s 


---

<a href="../spherapy/util/orbital_u.py#L49"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calcOrbitalVel`

```python
calcOrbitalVel(a: float, pos: ndarray[tuple[int], dtype[float64]]) → float
```

Return the instantaneous velocity magnitude for an elliptical orbit at position. 



**Args:**
 
 - <b>`a`</b>:  semi-major axis in m 
 - <b>`pos`</b>:  cartesian position, assuming the origin is at the central body. 



**Returns:**
 instantaneous velocity magnitude. 


---

<a href="../spherapy/util/orbital_u.py#L63"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calcMeanMotion`

```python
calcMeanMotion(a: float) → float
```

Returns mean motion [radians/s] for an elliptical or circular orbit with semi-major axis a. 



**Args:**
 
 - <b>`a`</b>:  semi-major axis in m 



**Returns:**
 Orbital period in s 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
