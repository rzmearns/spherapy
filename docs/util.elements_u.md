<!-- markdownlint-disable -->

<a href="../spherapy/util/elements_u.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `util.elements_u`
Utility functions for dealing with propagation element set strings. 


---

<a href="../spherapy/util/elements_u.py#L16"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `split3LELineIntoFields`

```python
split3LELineIntoFields(line: str) → ElementsLineDict
```

Create an ElementsLineDict from an element set line. 



**Args:**
 
 - <b>`line`</b>:  string 



**Returns:**
 ElementsLineDict 


---

<a href="../spherapy/util/elements_u.py#L32"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `stringify3LEDict`

```python
stringify3LEDict(tle_dict: dict[int, ElementsLineDict]) → str
```

Turn an element set dict 3LE dict back into a \n delimited string. 


---

<a href="../spherapy/util/elements_u.py#L37"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `dictify3LEs`

```python
dictify3LEs(lines: list[str]) → list[dict[int, ElementsLineDict]]
```

Turn list of strings into list of dicts storing TLE info. 

- dict:  
	- key: line number (0-2)
	- value: dict  
			- key: 'fields' | 'line_str'
			- value: 'list of fields as strings' | original TLE line str 


---

<a href="../spherapy/util/elements_u.py#L6"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ElementsLineDict`
TypedDict container for lines of an element set. 



**Attributes:**
 
 - <b>`fields`</b>:  list of each field of an element set line as strings 
 - <b>`line_str`</b>:  original whole line string of element set (no newline) 







---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
