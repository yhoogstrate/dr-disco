#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
[License: GNU General Public License v3 (GPLv3)]

 This file is part of Dr. Disco.

 FuMa is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Dr. Disco is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

import unittest
from subprocess import Popen, PIPE
from utils import main

# Nosetests doesn't use main()


class TestIsBlacklisted(unittest.TestCase):
    def test_is_blacklisted_1pos_unstranded(self):
        test_pos = 'chrM:12345'
        regions_file = 'share/blacklist-regions.hg38.bed'
        # junctions_file = 'share/blacklist-junctions.hg38.txt'

        command = ["bin/dr-disco",
                   "is-blacklisted",
                   '--blacklist-regions', regions_file,
                   test_pos]

        p = Popen(command, stdin=PIPE, stdout=PIPE)
        stdout = p.communicate()[0].decode("utf-8")

        self.assertEqual(p.returncode, 0)
        self.assertEqual(stdout, '1\tr-0-chrM\n2\tr-1-chrM\n')

    def test_is_blacklisted_1pos_pos(self):
        test_pos = 'chrM:12345(+)'
        regions_file = 'share/blacklist-regions.hg38.bed'
        # junctions_file = 'share/blacklist-junctions.hg38.txt'

        command = ["bin/dr-disco",
                   "is-blacklisted",
                   '--blacklist-regions', regions_file,
                   test_pos]

        p = Popen(command, stdin=PIPE, stdout=PIPE)
        stdout = p.communicate()[0].decode("utf-8")

        self.assertEqual(p.returncode, 0)
        self.assertEqual(stdout, '1\tr-0-chrM\n')

    def test_is_blacklisted_1pos_neg(self):
        test_pos = 'chrM:12345(-)'
        regions_file = 'share/blacklist-regions.hg38.bed'
        # junctions_file = 'share/blacklist-junctions.hg38.txt'

        command = ["bin/dr-disco",
                   "is-blacklisted",
                   '--blacklist-regions', regions_file,
                   test_pos]

        p = Popen(command, stdin=PIPE, stdout=PIPE)
        stdout = p.communicate()[0].decode("utf-8")

        self.assertEqual(p.returncode, 0)
        self.assertEqual(stdout, '1\tr-1-chrM\n')

    def test_is_blacklisted_2pos1(self):
        test_pos1 = 'chr10:10919999'
        test_pos2 = 'chr10:11165485'
        # regions_file = 'share/blacklist-regions.hg38.bed'
        junctions_file = 'share/blacklist-junctions.hg38.txt'

        command = ["bin/dr-disco",
                   "is-blacklisted",
                   '--blacklist-junctions', junctions_file,
                   test_pos1, test_pos2]

        p = Popen(command, stdin=PIPE, stdout=PIPE)
        stdout = p.communicate()[0].decode("utf-8")

        self.assertEqual(p.returncode, 0)
        self.assertEqual(stdout, '1\tj-0\n')

    def test_is_blacklisted_2pos2(self):
        test_pos1 = 'chr10:10919999(-)'
        test_pos2 = 'chr10:11165485(-)'
        # regions_file = 'share/blacklist-regions.hg38.bed'
        junctions_file = 'share/blacklist-junctions.hg38.txt'

        command = ["bin/dr-disco",
                   "is-blacklisted",
                   '--blacklist-junctions', junctions_file,
                   test_pos1, test_pos2]

        p = Popen(command, stdin=PIPE, stdout=PIPE)
        stdout = p.communicate()[0].decode("utf-8")

        self.assertEqual(p.returncode, 0)
        self.assertEqual(stdout, '')

    def test_is_blacklisted_2pos3(self):
        test_pos1 = 'chr10:10919999(+)'
        test_pos2 = 'chr10:11165485(+)'
        # regions_file = 'share/blacklist-regions.hg38.bed'
        junctions_file = 'share/blacklist-junctions.hg38.txt'

        command = ["bin/dr-disco",
                   "is-blacklisted",
                   '--blacklist-junctions', junctions_file,
                   test_pos1, test_pos2]

        p = Popen(command, stdin=PIPE, stdout=PIPE)
        stdout = p.communicate()[0].decode("utf-8")

        self.assertEqual(p.returncode, 0)
        self.assertEqual(stdout, '')

    def test_is_blacklisted_2pos4(self):
        test_pos1 = 'chr10:10919999(+)'
        test_pos2 = 'chr10:11165485(-)'
        # regions_file = 'share/blacklist-regions.hg38.bed'
        junctions_file = 'share/blacklist-junctions.hg38.txt'

        command = ["bin/dr-disco",
                   "is-blacklisted",
                   '--blacklist-junctions', junctions_file,
                   test_pos1, test_pos2]

        p = Popen(command, stdin=PIPE, stdout=PIPE)
        stdout = p.communicate()[0].decode("utf-8")

        self.assertEqual(p.returncode, 0)
        self.assertEqual(stdout, '')

    def test_is_blacklisted_2pos5(self):
        test_pos1 = 'chr10:10919999(-)'
        test_pos2 = 'chr10:11165485(+)'
        # regions_file = 'share/blacklist-regions.hg38.bed'
        junctions_file = 'share/blacklist-junctions.hg38.txt'

        command = ["bin/dr-disco",
                   "is-blacklisted",
                   '--blacklist-junctions', junctions_file,
                   test_pos1, test_pos2]

        p = Popen(command, stdin=PIPE, stdout=PIPE)
        stdout = p.communicate()[0].decode("utf-8")

        self.assertEqual(p.returncode, 0)
        self.assertEqual(stdout, '1\tj-0\n')

    def test_is_blacklisted_2pos6(self):
        test_pos1 = 'chrM:12345(-)'
        test_pos2 = 'chr10:11165485(+)'
        regions_file = 'share/blacklist-regions.hg38.bed'
        junctions_file = 'share/blacklist-junctions.hg38.txt'

        command = ["bin/dr-disco",
                   "is-blacklisted",
                   '--blacklist-regions', regions_file,
                   '--blacklist-junctions', junctions_file,
                   test_pos1, test_pos2]

        p = Popen(command, stdin=PIPE, stdout=PIPE)
        stdout = p.communicate()[0].decode("utf-8")

        self.assertEqual(p.returncode, 0)
        self.assertEqual(stdout, '1\tr-1-chrM\n')

    def test_is_blacklisted_2pos7(self):
        test_pos1 = 'chrM:12345(-)'
        test_pos2 = 'chr10:11165485(+)'
        regions_file = 'share/blacklist-regions.hg38.bed'
        junctions_file = 'share/blacklist-junctions.hg38.txt'

        command = ["bin/dr-disco",
                   "is-blacklisted",
                   '--blacklist-regions', regions_file,
                   '--blacklist-junctions', junctions_file,
                   test_pos2, test_pos1]  # swap

        p = Popen(command, stdin=PIPE, stdout=PIPE)
        stdout = p.communicate()[0].decode("utf-8")

        self.assertEqual(p.returncode, 0)
        self.assertEqual(stdout, '1\tr-1-chrM\n')


if __name__ == '__main__':
    main()
