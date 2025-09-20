The change log should include
- theme of the upversion
- changes in API
- changes in model behaviour

# v0.2.0
- initial release to pypi
- timespan object
- orbit object
- orbit object creation from: TLE, propagated orbital parameters, analytical orbital parameters, list of positions
- TLE updater, pulls requested TLEs from spacetrack, stores all TLEs for that satellite in a saved file
- TLE updater falls back to celestrak if no credentials provided