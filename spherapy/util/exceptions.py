"""Generic exceptions that apply throughout satplot."""

class OutOfRangeError(Exception):
	"""The value is out of the acceptable range."""

class DimensionError(Exception):
	"""Dimension of an array is invalid."""
