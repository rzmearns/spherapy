# Contribution Guide

## How to Contribute
Almost all communication with the spherapy maintainers should be done through the main spherapy Gitlab repository: https://gitlab.unimelb.edu.au/msl/libraries/spherapy

- Bug reports and feature requests can be submitted through the “Issues” page on the repository. It is critical that sufficient information is provided in the bug report to recreate the problem. Please provide:
  - A description of the behaviour.
  - The TLE files which resulted in the bug
  - The timespan settings

- Any changes to actual code, including documentation, should be submitted as a merge request on Gitlab.
  - Please make sure to submit merge requests using a new branch following the [Branch Naming Convention](#branch-naming-convention). Don’t be afraid to create a merge request as soon as the branch is created. It can be always be updated later. Creating them early gives maintainers a chance to provide an early review of your code if that’s something you’re looking for. See below for more information on writing documentation and checking your changes.
  - Make sure you understand the overall architecture of spherapy before adding features.

- No matter how you contribute, spherapy maintainers will try their best to read, research, and respond to your query as soon as possible. For code changes, automated checks and tests will run on Gitlab to provide an initial “review” of your changes.

## Submitting Bugs
Please submit problems or bus as issues in gitlab,

https://gitlab.unimelb.edu.au/msl/libraries/spherapy/-/issues

## Submitting Feature Ideas
Please submit feature ideas as an issue in gitlab, 

https://gitlab.unimelb.edu.au/msl/libraries/spherapy/-/issues

## Asking Questions
Please submit questions as an issue in gitlab,

https://gitlab.unimelb.edu.au/msl/libraries/spherapy/-/issues

## Adding Features
### Branch Naming Convention

## Development Environment
It is recommended that you setup a virtual environment running python 3.10 or later, and install the packages required within that virtual environment using
```
pip install -r requirements.txt
```

## Coding Style
In general, spherapy follows a modified PEP 8 style guidelines:

### Writing Tests

### Running Tests