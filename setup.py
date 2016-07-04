#!/usr/bin/env python3
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
Installer of Dr. Disco

[License: GNU General Public License v3 (GPLv3)]
 
 This file is part of Dr. Disco.
 
 FuMa is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 FuMa is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

import drdisco

from setuptools import setup, find_packages
setup(name = "dr-disco",
        scripts=['bin/dr-disco'],
        packages = ["drdisco"],
        test_suite="tests",
        version = drdisco.__version__,
        description = "Makes discordant RNA-Seq alignments healthy, and tries to interpret intronic break points",
        author = drdisco.__author__,
        url = drdisco.__homepage__,
        keywords = ["rna-seq", "intronic", "break point"],
        classifiers = [
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
     )
