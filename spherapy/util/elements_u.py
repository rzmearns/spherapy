"""Utility functions for dealing with propagation element set strings."""

from typing import TypedDict


class ElementsLineDict(TypedDict):
	"""TypedDict container for lines of an element set.

	Attributes:
		fields: list of each field of an element set line as strings
		line_str: original whole line string of element set (no newline)
	"""
	fields: list[str]
	line_str: str

def split3LELineIntoFields(line:str) -> ElementsLineDict:
	"""Create an ElementsLineDict from an element set line.

	Args:
		line: string

	Returns:
		ElementsLineDict
	"""
	fields = line.split()
	# insert fields into and store original line in tle dict
	line_dict:ElementsLineDict = {'fields':fields,
							'line_str':line}

	return line_dict

def stringify3LEDict(tle_dict:dict[int, ElementsLineDict]) -> str:
	r"""Turn an element set dict 3LE dict back into a \n delimited string."""
	lines = [line_dict['line_str'] for line_dict in tle_dict.values()]
	return '\n'.join(lines)

def dictify3LEs(lines:list[str]) -> list[dict[int, ElementsLineDict]]:
	"""Turn list of strings into list of dicts storing TLE info.

	dict:
		key: line number (0-2)
		value: dict
			key: 'fields' | 'line_str'
			value: 'list of fields as strings' | original TLE line str
	"""
	if len(lines) == 0:
		raise ValueError("No data")
	if len(lines)%3 != 0:
		# If not obvious what lines relate to eachother, can't make assumption -> abort.
		raise ValueError('Incomplete TLEs present, aborting')

	list_3les = []
	tle:dict[int,ElementsLineDict] = {0:{'fields':[], 'line_str':''},
			1:{'fields':[], 'line_str':''},
			2:{'fields':[], 'line_str':''}}
	for line in lines:
		line_dict = split3LELineIntoFields(line)
		# insert fields into and store original line in tle dict
		tle_line_num = int(line_dict['fields'][0])
		tle[tle_line_num] = line_dict

		if tle_line_num == 2: 		# noqa: PLR2004
			# store parsed tle
			list_3les.append(tle.copy())
			# empty dict
			tle[0] = {'fields':[], 'line_str':''}
			tle[1] = {'fields':[], 'line_str':''}
			tle[2] = {'fields':[], 'line_str':''}

	return list_3les
