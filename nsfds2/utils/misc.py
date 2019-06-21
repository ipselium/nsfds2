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
import os as _os
import datetime as _datetime
import numpy as _np
import scipy.signal as _sps
import scipy.io.wavfile as _wf


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


def resample(file, target_rate, pad=None, write=False, force_mono=True):
    """ Resample target wave file with target_rate. """

    target_rate = int(target_rate)

    if not 1 < target_rate < 4.3e6:
        raise ValueError('Sampling rate must be 1 < rate < 4.3e6')

    path = _os.path.dirname(file)
    filename = _os.path.basename(file).split('.')
    rate, data = _wf.read(file)
    dtype = data.dtype
    duration = data.shape[0]/rate
    N = int(target_rate*duration)

    if len(data.shape) == 2 and force_mono:   # stereo to mono
        data = (data[:, 0] + data[:, 1])/2

    print(colors.BLUE + f'Resampling {file} at {target_rate} kHz ({N} points)...')
    print(colors.RED + f'Set nt > {N} to play the whole sound' + colors.END)
    if len(data.shape) == 1:   # mono
        data_r = _sps.resample(data, N).astype(dtype)
        if pad:
            data_r = get_padded(data_r, pad)


    if len(data.shape) == 2:   # stereo
        tmp_l = _sps.resample(data[:, 0], N).astype(dtype)
        tmp_r = _sps.resample(data[:, 1], N).astype(dtype)
        if pad:
            tmp_l = get_padded(tmp_l, pad)
            tmp_r = get_padded(tmp_r, pad)

        data_r = _np.vstack([tmp_l, tmp_r]).T


    if write:
        print(f'Writing {N} samples at {target_rate} kHz rate...')
        _wf.write(path + '{}_r.{}'.format(*filename), rate=target_rate, data=data_r)

    return data_r/abs(data_r).max()


def get_padded(s, N, value=0):
    """ Pad signal with value. """
    if N > s.shape[0]:
        return _np.concatenate([s, value*_np.ones(N - s.shape[0])])

    return s


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
