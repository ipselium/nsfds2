#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2016-2018 Cyril Desjouy <cyril.desjouy@univ-lemans.fr>
#
# This file is part of nsfds2
#
# nsfds2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# nsfds2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with nsfds2. If not, see <http://www.gnu.org/licenses/>.
#
#
# Creation Date : mar. 10 avril 2018 17:52:42 CEST
# Last Modified : ven. 11 mai 2018 16:13:55 CEST
"""
-----------

setup file for nsfds2

-----------
"""

from setuptools import setup, find_packages

setup(
    name='nsfds2',
    description="Finite difference solver for Navier-Stokes equations",
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    version="0.9.12",
    license="GPL",
    url='https://github.com/ipselium/nsfds2',
    author="Cyril Desjouy",
    author_email="cyril.desjouy@univ-lemans.fr",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["numpy", "scipy", "matplotlib", "ofdlib2>=0.9.5",
                      "progressbar33", "mplutils>=0.3.0", "h5py", "fdgrid>=0.8.7"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],
    entry_points={
        'console_scripts': [
            'nsfds2 = nsfds2.solver:main',
        ],
    }
)
