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

Utils: Misc.

-----------
"""

import re
import itertools
import datetime as _datetime


class colors:
    """ ANSI codes """

    PURPLE = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    END = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def bc_combinations():
    """ List of possible combinations for bc"""

    corners = list(itertools.product(['X', 'A'], repeat=4))
    corners = [''.join(i) for i in corners if i.count('A') == 2]
    corners = [i for i in corners if any(re.match(regex, i) for regex in [r'A..A', r'AA..', r'..AA', r'.AA.'])]

    pml = list(itertools.product(['W', 'P', 'X', 'A'], repeat=4))
    pml = [''.join(i) for i in pml if i.count('A') == 1]
    pml = [i for i in pml if 'A' in i and i.count('P') < 2 and i.count('W') < 3 and i.count('X') > 0]
    pml = [i for i in pml if any(re.match(regex, i) for regex in [r'.A.X', r'.X.A', r'X.A.', r'A.X.'])]


    pml_p = [i for i in pml if 'P' in i]
    pml_r = [i for i in pml if 'W' in i and 'P' not in i]
    pml_x = [i for i in pml if 'W' not in i and 'P' not in i]

    print('\n* {} Corners : {}'.format(len(corners), ', '.join(corners)))
    print('\n* {} with only X : {}'.format(len(pml_x), ', '.join(pml_x)))
    print('\n* {} with P : {}'.format(len(pml_p), ', '.join(pml_p)))
    print('\n* {} with only W & X : {}'.format(len(pml_r), ', '.join(pml_r)))


def secs_to_dhms(secs):
    """ Convert seconds to years, months, days, hh:mm:ss."""

    dhms = _datetime.datetime(1, 1, 1) + _datetime.timedelta(seconds=secs)

    year, years = f'{dhms.year-1} year, ', f'{dhms.year-1} years, '
    month, months = f'{dhms.month-1} month, ', f'{dhms.month-1} months, '
    day, days = f'{dhms.day-1} day, ', f'{dhms.day-1} days, '
    h = f'{dhms.hour}:'
    m = f'{dhms.minute:02}:'
    s = f'{dhms.second:02}'

    return (year if dhms.year == 2 else years if dhms.year > 2 else '') + \
           (month if dhms.month == 2 else months if dhms.month > 2 else '') + \
           (day if dhms.day == 2 else days if dhms.day > 2 else '') + \
           (h if dhms.hour > 0 else '') + m + s


def nearest_index(n, ns, nt):
    """ Returns nearest possible index `n`

    Parameters
    ----------
    n: int
        Index to look for
    ns: int
        Subdivision of `nt`
    nt: int
        Total number of iterations
    """

    if n > nt:
        return nt

    if n%ns == n:
        return ns

    if n%ns > ns/2:
        return (n//ns + 1)*ns

    if n%ns <= ns/2:
        return n//ns*ns

    return ns
