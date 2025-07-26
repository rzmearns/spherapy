<!-- markdownlint-disable -->

<a href="../spherapy/util/credentials.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `util.credentials`
Functions to fetch and store spacetrack credentials on the system. 


---

<a href="../spherapy/util/credentials.py#L14"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `fetchConfigCredentials`

```python
fetchConfigCredentials(config: ConfigParser) → dict[str, str | None]
```

Fetch spacetrack credentials from config file. 



**Args:**
 
 - <b>`config`</b>:  ConfigParser object containing config fields 



**Returns:**
 dict containing username and password 


---

<a href="../spherapy/util/credentials.py#L33"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `fetchKeyringCredentials`

```python
fetchKeyringCredentials() → dict[str, str | None]
```

Fetch spacetrack credentials from system keyring. 



**Returns:**
  dict containing username and password 



**Raises:**
 
 - <b>`ValueError`</b>:  If no system keyring exists 


---

<a href="../spherapy/util/credentials.py#L66"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `storeCredentials`

```python
storeCredentials(user: None | str = None, passwd: None | str = None) → bool
```

Store spacetrack credentials in system keyring. 

Can set either username or password, but can't set a password without a username 



**Args:**
 
 - <b>`user`</b>:  [optional] username string 
 - <b>`passwd`</b>:  [optional] password string 



**Returns:**
 True if successfully stored 



**Raises:**
 
 - <b>`KeyError`</b>:  raised if trying to store a password without a username 


---

<a href="../spherapy/util/credentials.py#L96"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `createCredentials`

```python
createCredentials()
```

Script helper function to create and store credentials. 

Called by command line script. Requires user input. 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
