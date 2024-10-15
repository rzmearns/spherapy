import pathlib

def isValidPath(path_str:str) -> bool:
	if path_str[-1] == '/':
		return True
	else:
		return False
	
def sanitisePath(path_str:str) -> str:

	return pathlib.Path(path_str)