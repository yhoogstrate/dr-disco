#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""[License: GNU General Public License v3 (GPLv3)]

    Dr. Disco: fusion gene detection in random hexamer RNA-seq data
    Copyright (C) 2017  Youri Hoogstrate

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    You can contact me via the github repository at the following url:
    <https://github.com/yhoogstrate/dr-disco>

    You can e-mail me via 'y.hoogstrate' at the following webmail domain:
    gmail dot com
"""

import logging
import sys

__version_info__ = ('0', '15', '2')
__version__ = '.'.join(__version_info__) if (len(__version_info__) == 3) else '.'.join(__version_info__[0:3]) + "-" + __version_info__[3]
__author__ = 'Youri Hoogstrate'
__homepage__ = 'https://github.com/yhoogstrate/dr-disco'
__license__ = 'GNU General Public License v3 (GPLv3)'
__license_notice__ = 'License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.\nThis is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.'

__log_format__ = "[%(filename)s:%(lineno)s - %(funcName)s()] %(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(level=logging.DEBUG, format=__log_format__, stream=sys.stderr)  # bioconda seems to crash on stdout here..
log = logging.getLogger(__name__)


# parameters
MAX_ACCEPTABLE_INSERT_SIZE = 450
MAX_ACCEPTABLE_ALIGNMENT_ERROR = 3  # acceptable splice junction error or other general alignment mistakes
MAX_GENOME_DISTANCE = 999999999
MAX_SIZE_CIRCULAR_RNA = 100000

# filter settings
MIN_DISCO_PER_SUBNET_PER_NODE = 1  # minimum nodes is 2 per subnet, hence mininal 2 discordant reads are necessairy
