#!/usr/bin/env python3
"""
Copyright 2020-2024 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
https://github.com/edgardomortiz/Captus

Captus' installation script

This file is part of Captus. Captus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Captus is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Captus. If
not, see <http://www.gnu.org/licenses/>.
"""

# Make sure this is being run with Python 3.6 or later.
import sys
if sys.version_info.major != 3 or sys.version_info.minor < 6:
    sys.exit('Error: you must execute setup.py using Python 3.6 or later')

from setuptools import setup, find_packages

# Get the program version from another file.
__version__ = "0.0.0"
exec(open('captus/version.py').read())

setup(
    name = "captus",
    version = __version__,
    url = "https://github.com/edgardomortiz/Captus",
    author = "Edgardo M. Ortiz",
    author_email = "e.ortiz.v@gmail.com",
    description = "Tools for hybridization capture-based target enrichment: "
                  "Probe Design, Data Assembly, Marker Extraction and Alignment",
    long_description = open("README.md").read(),
    long_description_content_type = "text/markdown",
    packages = find_packages(),
    package_data = {
        "data": ["*"],
        "dependencies": ["scipio-1.4/*", "blat/*"]
        },
    include_package_data=True,
    entry_points = {
        "console_scripts": ["captus_assembly = captus.captus_assembly:main",
                            "captus = captus.captus_assembly:main",
                            "captus_design = captus.captus_design:main",
                            "captusd = captus.captus_design:main"]
        },
    license = "GPL",
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ]
)
