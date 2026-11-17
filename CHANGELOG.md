The change log should include
- theme of the upversion
- changes in API
- changes in model behaviour

# Unreleased
- **Breaking:** `Orbit.fromAnalyticalOrbitalParam()` now accepts `mean_nu` (mean anomaly at epoch, in degrees) instead of `true_nu`. Existing `true_nu=` calls raise `TypeError`; positional anomaly arguments now mean mean anomaly. To preserve an eccentric orbit's initial phase, convert its true anomaly to mean anomaly before calling. `epoch` and the zero-anomaly default are unchanged.
- Replaced Hapsira with SpiceyPy's vectorized SPICE conic propagator for analytical orbits, removing slow Numba compilation on first use.
- Calculates Sun/Moon positions directly with Astropy.
- Rejects non-finite analytical orbital parameters explicitly.

# v0.3.0
 - renamed fromAnalyticalOrbitalParams() argument mean_nu to true_nu to reflect the underlying api call to hapsira.
 - added an epoch argument for the orbital parameters of fromAnalyticalOrbitalParams(). Note the default epoch for fromAnalyticalOrbitalParams() is now the first timestep of the TimeSpan rather than astropy's J2000. This is a breaking change and thus existing calls to fromAnalyticalOrbitalParams() will almost certainly now produce different results. The previous functionality can be obtained by choosing datetime.datetime.fromisoformat('2000-01-01T12:00:00+00:00') as the epoch.
 - adding typing overloads to the TimeSpan.__getitem__() function which specify a stronger typing constraints dependent on the type of the argument.
 - fixed TimeSpan.__getitem__() adding extra layer of list nesting when index was of type list.

# v0.2.1
- minor bug fixes to the API for pulling TLEs from Celestrak

# v0.2.0
- initial release to pypi
- timespan object
- orbit object
- orbit object creation from: TLE, propagated orbital parameters, analytical orbital parameters, list of positions
- TLE updater, pulls requested TLEs from spacetrack, stores all TLEs for that satellite in a saved file
- TLE updater falls back to celestrak if no credentials provided
