# Contribution Guide

## How to Contribute
Almost all communication with the spherapy maintainers should be done through the main spherapy Github repository: https://github.com/rzmearns/spherapy

- Bug reports and feature requests can be submitted through the [Issues](https://github.com/rzmearns/spherapy/issues) page on the repository. It is critical that sufficient information is provided in the bug report to recreate the problem. Please provide:
  - The spherapy version
  - The python version
  - A description of the behaviour.
  - The TLE files which resulted in the bug
  - The timespan settings

- Any changes to actual code, including documentation, should be submitted as a pull request on Github.
  - Please make sure to submit pull requests using a new branch following the [Branch Naming Convention](#branch-naming-convention). Don’t be afraid to create a pull request as soon as the branch is created. It can be always be updated later. Creating them early gives maintainers a chance to provide an early review of your code if that’s something you’re looking for. See below for more information on writing documentation and checking your changes.
  - Make sure you understand the overall architecture of spherapy before adding features.

- No matter how you contribute, spherapy maintainers will try their best to read, research, and respond to your query as soon as possible. For code changes, automated checks and tests will run on Github to provide an initial “review” of your changes.

## Submitting Bugs
Please submit problems or bugs as issues in github, they will be tagged by a maintainer as a `bug`

https://github.com/rzmearns/spherapy/issues

## Submitting Feature Ideas
Please submit feature ideas as an issue in Github, they will be tagged by a maintainer as an `enhancement`

https://github.com/rzmearns/spherapy/issues

## Asking Questions
Please submit questions as an issue in Github, they will be tagged by a maintainer as `support`

https://github.com/rzmearns/spherapy/issues

## Adding Features
If you have a contribution to make to spherapy, please clone the repo, and submit a pull request.  
In order to be considered for merging, the branch must pass the [Linting](#linting), [Testing](#testing) and [Static Typing](#typing) stages of the CI pipeline.  
For any branch you make, please follow the following branch naming conventions:

### Branch Naming Convention
- If a branch is adding a feature to spherapy: `feature/{descriptive_name}`
  - the feature must be supported by unit tests included in the branch
  - will results in a minor upversion
- If a branch is fixing a spherapy bug: `bugfix/{descriptive_name}` or `bugfix/{ISSUE#}`
  - the bugfix must be supported by regrestion tests included in the branch
  - will result in a patch upversion
- If a branch is modifying or fixing a bug in the infrastructure (CI, etc.): `infrastucture/{descriptive_name}` or `infrastructure/{ISSUE#}`
- If a branch is updating the documentation: `docs/{descriptive_name}`
- The maintainers may use `hotfix/{descriptive_name}` or `hotfix/{ISSUE#}` if an urgent fix is required to the main branch
  - will result in a patch upversion


## Development Environment
It is recommended that you setup a virtual environment running python 3.12 (spherapy is compatible with 3.10->3.12), and install the packages required within that virtual environment.

### Virtual Environment
```
pip install virtualenv
virtualenv -p=`which python3.12` .venv312
source .venv312/bin/activate
```

### Clone and install spherapy
```
git clone git@github.com:rzmearns/spherapy.git
pip install --editable .[dev]
```

### Linting
Linting is acheived using Ruff. The necessary packages will be installed as part of the \[dev\] depenedencies.  
Spherapy has an extensive list of linter warnings enabled by default, the full list can be found in the pyproject.toml under [tool.ruff.lint], while the explanation to each class of error, and specific errors can be found on Ruff's website, [Rules](https://docs.astral.sh/ruff/rules/)

If there is a good reason, individual linting errors can be excluded, but not entire error classes.

For coding style guidelines, see [Coding Style](#coding-style)

linting can be run locally using
```
ruff check
```
Tedious changes, for example sorting the import blocks, can be handled automatically by Ruff:
```
ruff check --fix
```
Linting will be automatically run whenever a commit is made to a pull request.

### Testing
Testing is achieved with pytest and pytest-check. The necessary packages will be installed as part of the \[dev\] depenedencies.  
Tests should be located in the `tests/` dir, following a similar structure to the package structure.  
If you add functionality to spherapy, there should be dedicated unit tests for each feature. If you contribute by fixing a bug, tests should be added to ensure that the bug does not get reintroduced in future merges. These regression tests should have a comment note detailing what they are protecting against.  

spherapy often uses pytest-check to run all tests in a function in parallel, so that any one failure does not stop other tests in that function from running.  

tests can be run locally using 
```
pytest
```

Tests will be automatically run whenever a commit is made to a pull request. (The tests will run under python3.10, python3.11 and python3.12)

### Typing
spherapy employs static type analysis using mypy (hopefully this will be changed to ty in the future). The necessary packages will be installed as part of the \[dev\] depenedencies.

Static typing can be run locally using
```
mypy spherapy --disable-error-code=import-untyped
```
The `--disable-error-code=import-untyped` flag is necessary to avoid errors with untyped libraries like astropy.

Static typing will be automatically run whenever a commit is made to a pull request.

### pre-commit hooks
Although not necessary, it is recommended to use [pre-commit](https://pre-commit.com/). The necessary packages will be installed as part of the \[dev\] depenedencies.  
This runs configured checks on the codebase prior to every commit. The configured checks can be seen in [.pre-commit-config.yaml](https://github.com/rzmearns/spherapy/blob/main/.pre-commit-config.yaml)  
The configured checks for spherapy are:
- `ruff check`
- `gitleaks`

To configure pre-commit, after `pip install .[dev]`, run `pre-commit install`

## Coding Style
In general, spherapy follows a modified PEP 8 style guidelines, that are enforced by the linter. The notable differences are:
- all indentations are TABS (NOT SPACES.)
- variables are snake_case (this_is_a_variable) (like PEP8)
- functions are mixedCase (thisIsAFunction)
  - if there is an acronym in the function name, it should be in UPPERCASE.