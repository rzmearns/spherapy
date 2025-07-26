<!-- markdownlint-disable -->

# API Overview

## Public API

- [`orbit.Orbit`](./orbit.md#class-orbit): Timestamped orbital data for a satellite
	- [`Orbit.fromAnalyticalOrbitalParam`](./orbit.md#classmethod-fromanalyticalorbitalparam)
	- [`Orbit.fromListOfPositions`](./orbit.md#classmethod-fromlistofpositions)
	- [`Orbit.fromPropagatedOrbitalParam`](./orbit.md#classmethod-frompropagatedorbitalparam)
	- [`Orbit.fromTLE`](./orbit.md#classmethod-fromtle)
- [`timespan.TimeSpan`](./timespan.md#class-timespan): A series of timestamps.
- [`updater.getStoredEpochLimits`](./updater.md#function-getstoredepochlimits): Returns limiting epochs for the stored TLEs for each sat in sat_id_list.
- [`updater.getTLEFilePaths`](./updater.md#function-gettlefilepaths): Fetch list of paths to TLE files.
- [`updater.updateTLEs`](./updater.md#function-updatetles): Fetch most recent TLE for provided list of satcat IDs.

## All Modules

- [`orbit`](./orbit.md#module-orbit): Class for orbital data.
- [`timespan`](./timespan.md#module-timespan): Class for series of timestamps.
- [`updater`](./updater.md#module-updater): Functions to update TLE of a satcat ID, and fetch where the file is stored.
- [`util`](./util.md#module-util)
- [`util.celestrak`](./util.celestrak.md#module-utilcelestrak): Functions to fetch TLEs from Celestrak.
- [`util.constants`](./util.constants.md#module-utilconstants): Namespace of useful orbital constants.
- [`util.credentials`](./util.credentials.md#module-utilcredentials): Functions to fetch and store spacetrack credentials on the system.
- [`util.elements_u`](./util.elements_u.md#module-utilelements_u): Utility functions for dealing with propagation element set strings.
- [`util.epoch_u`](./util.epoch_u.md#module-utilepoch_u): Utility functions for converting different epochs.
- [`util.exceptions`](./util.exceptions.md#module-utilexceptions): Generic exceptions that apply throughout satplot.
- [`util.orbital_u`](./util.orbital_u.md#module-utilorbital_u): Utility functions for orbital calculations.
- [`util.spacetrack`](./util.spacetrack.md#module-utilspacetrack): Functions to fetch TLEs from Spacetrack.

## All Classes

- [`orbit.Orbit`](./orbit.md#class-orbit): Timestamped orbital data for a satellite, coordinate system depending on the central body.
- [`orbit.OrbitAttrDict`](./orbit.md#class-orbitattrdict): A TypedDict providing type annotations for the Orbit class.
- [`timespan.TimeSpan`](./timespan.md#class-timespan): A series of timestamps.
- [`elements_u.ElementsLineDict`](./util.elements_u.md#class-elementslinedict): TypedDict container for lines of an element set.
- [`exceptions.DimensionError`](./util.exceptions.md#class-dimensionerror): Dimension of an array is invalid.
- [`exceptions.OutOfRangeError`](./util.exceptions.md#class-outofrangeerror): The value is out of the acceptable range.
- [`spacetrack.InvalidCredentialsError`](./util.spacetrack.md#class-invalidcredentialserror): Error indicating invalid credentials used to access spacetrack.

## All Functions

- [`updater.getStoredEpochLimits`](./updater.md#function-getstoredepochlimits): Returns limiting epochs for the stored TLEs for each sat in sat_id_list.
- [`updater.getTLEFilePaths`](./updater.md#function-gettlefilepaths): Fetch list of paths to TLE files.
- [`updater.updateTLEs`](./updater.md#function-updatetles): Fetch most recent TLE for provided list of satcat IDs.
- [`celestrak.getStoredEpochs`](./util.celestrak.md#function-getstoredepochs): Return the start and end epoch for {sat_id}.temptle .
- [`celestrak.getTLEFilePath`](./util.celestrak.md#function-gettlefilepath): Gives path to file where celestrak TLE is stored.
- [`celestrak.updateTLEs`](./util.celestrak.md#function-updatetles): Fetch most recent TLE for satcat IDs from celestrak.
- [`credentials.createCredentials`](./util.credentials.md#function-createcredentials): Script helper function to create and store credentials.
- [`credentials.fetchConfigCredentials`](./util.credentials.md#function-fetchconfigcredentials): Fetch spacetrack credentials from config file.
- [`credentials.fetchKeyringCredentials`](./util.credentials.md#function-fetchkeyringcredentials): Fetch spacetrack credentials from system keyring.
- [`credentials.storeCredentials`](./util.credentials.md#function-storecredentials): Store spacetrack credentials in system keyring.
- [`elements_u.dictify3LEs`](./util.elements_u.md#function-dictify3les): Turn list of strings into list of dicts storing TLE info.
- [`elements_u.split3LELineIntoFields`](./util.elements_u.md#function-split3lelineintofields): Create an ElementsLineDict from an element set line.
- [`elements_u.stringify3LEDict`](./util.elements_u.md#function-stringify3ledict): Turn an element set dict 3LE dict back into a \n delimited string.
- [`epoch_u.datetime2TLEepoch`](./util.epoch_u.md#function-datetime2tleepoch): Converts a datetime to a TLE epoch string.
- [`epoch_u.datetime2sgp4epoch`](./util.epoch_u.md#function-datetime2sgp4epoch): Converts a datetime to an sgp4 epoch.
- [`epoch_u.epoch2datetime`](./util.epoch_u.md#function-epoch2datetime): Converts a fractional epoch string to a datetime object.
- [`epoch_u.epochEarlierThan`](./util.epoch_u.md#function-epochearlierthan): Check if epoch A is earlier than epoch B.
- [`epoch_u.epochLaterThan`](./util.epoch_u.md#function-epochlaterthan): Check if epoch A is later than epoch B.
- [`epoch_u.findClosestDatetimeIndices`](./util.epoch_u.md#function-findclosestdatetimeindices): Find the index of the closest datetime in source arr for each datetime in test_arr.
- [`epoch_u.getStoredEpochs`](./util.epoch_u.md#function-getstoredepochs): Return the start and end epoch for tle_path.
- [`orbital_u.calcMeanMotion`](./util.orbital_u.md#function-calcmeanmotion): Returns mean motion [radians/s] for an elliptical or circular orbit with semi-major axis a.
- [`orbital_u.calcOrbitalVel`](./util.orbital_u.md#function-calcorbitalvel): Return the instantaneous velocity magnitude for an elliptical orbit at position.
- [`orbital_u.calcPeriod`](./util.orbital_u.md#function-calcperiod): Returns the period of an elliptical or circular orbit.
- [`orbital_u.ssoInc`](./util.orbital_u.md#function-ssoinc): Generates required inclination for given altitude [km] to maintain Sun Syncrhonous orbit.
- [`spacetrack.doCredentialsExist`](./util.spacetrack.md#function-docredentialsexist): Checks if spacetrack credentials have been loaded into Spherapy.
- [`spacetrack.getStoredEpochs`](./util.spacetrack.md#function-getstoredepochs): Return the start and end epoch for {sat_id}.tle .
- [`spacetrack.getTLEFilePath`](./util.spacetrack.md#function-gettlefilepath): Gives path to file where spacetrack TLE is stored.
- [`spacetrack.updateTLEs`](./util.spacetrack.md#function-updatetles): Fetch most recent TLE for satcat IDs from spacetrack.


---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
