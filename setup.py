#!/usr/bin/env python3
"""
Copyright 2020-2026 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
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

from setuptools import find_packages, setup

# Get the program version from another file safely without using exec().
__version__ = "0.0.0"
with open("captus/version.py", encoding="utf-8") as f:
    for line in f:
        if line.startswith("__version__"):
            __version__ = line.split("=")[1].strip().strip('"').strip("'")
            break

with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="captus",
    version=__version__,
    url="https://github.com/edgardomortiz/Captus",
    author="Edgardo M. Ortiz",
    author_email="e.ortiz.v@gmail.com",
    description="Tools for hybridization capture-based target enrichment: "
    "Probe Design, Data Assembly, Marker Extraction and Alignment",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    # Enforces Python 3.6+ natively at the package installer level
    python_requires=">=3.6",
    package_data={
        "captus": ["data/*", "dependencies/scipio-1.4/*", "dependencies/blat/*"]
    },
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "captus_assembly = captus.captus_assembly:main",
            "captus = captus.captus_assembly:main",
            "captus_design = captus.captus_design:main",
            "captusd = captus.captus_design:main",
            "concatenate_alignments = captus.concatenate_alignments:main",
            "most_common_target_per_locus = captus.most_common_target_per_locus:main",
            "new_targets_from_alignments = captus.new_targets_from_alignments:main",
            "phylo_commands = captus.phylo_commands:main",
        ]
    },
    license="GPL",
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
    ],
)
