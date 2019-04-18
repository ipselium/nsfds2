#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2016-2019 Cyril Desjouy <cyril.desjouy@univ-lemans.fr>
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
# Creation Date : 2019-04-18 - 15:23:05
"""
-----------
DOCSTRING

@author: Cyril Desjouy
"""

import re
import itertools

def bc_combinations():
    """ List of possible combinations for bc"""

    corners = list(itertools.product(['X', 'A'], repeat=4))
    corners = [''.join(i) for i in corners if i.count('A') == 2]
    corners = [i for i in corners if any(re.match(regex, i) for regex in [r'A..A', r'AA..', r'..AA', r'.AA.'])]

    pml = list(itertools.product(['R', 'P', 'X', 'A'], repeat=4))
    pml = [''.join(i) for i in pml if i.count('A') == 1]
    pml = [i for i in pml if 'A' in i and i.count('P') < 2 and i.count('R') < 3 and i.count('X') > 0]
    pml = [i for i in pml if any(re.match(regex, i) for regex in [r'.A.X', r'.X.A', r'X.A.', r'A.X.'])]


    pml_p = [i for i in pml if 'P' in i]
    pml_r = [i for i in pml if 'R' in i and 'P' not in i]
    pml_x = [i for i in pml if 'R' not in i and 'P' not in i]

    print('\n* {} Corners : {}'.format(len(corners), ', '.join(corners)))
    print('\n* {} with only X : {}'.format(len(pml_x), ', '.join(pml_x)))
    print('\n* {} with P : {}'.format(len(pml_p), ', '.join(pml_p)))
    print('\n* {} with only R & X : {}'.format(len(pml_r), ', '.join(pml_r)))
