import pathlib
import setuptools

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setuptools.setup(
    name="spherapy",
    version="0.0.1",
    description="An orbital propagator wrapper and TLE fetcher",
    long_description=README,
    long_description_content_type="text/x-rst",
    author="Robert Mearns",
    author_email="robert.mearns@unimelb.edu.au",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python"
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.10",
	install_requires=[
		'astropy ==6.1.4',
		'scipy >=1.8',
		'numpy',
		'skyfield ==1.46',
		'progressbar2 ==4.4.1',
		'hapsira ==0.18.0',
		'spacetrack ==1.2.0'
	]
)
