"""Utility functions for orbital calculations."""
import numpy as np

import spherapy.util.constants as consts


def ssoInc(alt:float, e:float=0) -> float:
	# TODO: update to calculate for different central bodies
	"""Generates required inclination for given altitude [km] to maintain Sun Syncrhonous orbit.

	Args:
		alt: altitude of orbit in km
		e: [Optional] eccentricity of orbit
			Default is circular (0)

	Returns:
		Inclination angle in degrees
	"""
	a = consts.R_EARTH + alt * 1e3
	# print(a)
	p = a * (1 - e**2)
	# print(p)
	period = 2 * np.pi * np.sqrt(a**3 / consts.GM_EARTH)
	# print(period)
	req_prec_rate = (2 * np.pi / 365.25) * (1 / 86400.)
	# print(req_prec_rate)
	req_prec_orb = req_prec_rate * period
	# print(req_prec_orb)

	cosi = (req_prec_orb * p**2) / (-3 * np.pi * consts.J2 * consts.R_EARTH**2)
	# print(cosi)

	inc = np.arccos(cosi)

	return np.rad2deg(inc)


def calcPeriod(a:float) -> float:
	"""Returns the period of an elliptical or circular orbit.

	Args:
		a: semi-major axis in m

	Returns:
		Orbital period in s
	"""
	return 2 * np.pi * np.sqrt(a**3 / consts.GM_EARTH)

def calcOrbitalVel(a:float, pos:np.ndarray[tuple[int],np.dtype[np.float64]]) -> float:
	"""Return the instantaneous velocity magnitude for an elliptical orbit at position.

	Args:
		a: semi-major axis in m
		pos: cartesian position, assuming the origin is at the central body.

	Returns:
		instantaneous velocity magnitude.
	"""
	r = np.linalg.norm(pos)

	return np.sqrt(consts.GM_EARTH * (2 / r - 1 / a))

def calcMeanMotion(a:float) -> float:
	"""Returns mean motion [radians/s] for an elliptical or circular orbit with semi-major axis a.

	Args:
		a: semi-major axis in m

	Returns:
		Orbital period in s
	"""
	return np.sqrt(consts.GM_EARTH / a**3)
