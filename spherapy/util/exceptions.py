"""Generic exceptions that apply throughout satplot."""

class OutOfRange(Exception):
	"""The value is out of the acceptable range."""
	pass

class DimensionError(Exception):
	"""Dimension of an array is invalid."""
	pass