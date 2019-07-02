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
# Creation Date : 2019-06-21 - 12:30:57
"""
-----------

Utils: Sounds

-----------
"""

import os as _os
import numpy as _np
import scipy.signal as _sps
import scipy.io.wavfile as _wf
from nsfds2.utils.misc import colors as _colors
from nsfds2.utils import DataExtractor as _DataExtractor

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

    print(_colors.BLUE + f'Resampling {file} at {target_rate} kHz ({N} points)...')
    print(_colors.RED + f'Set nt > {N} to play the whole sound' + _colors.END)
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


def normalize(s, dtype=_np.int16):
    """ Normalize wav data. """
    smax = _np.iinfo(dtype).max
    return (s/abs(s).max()*smax).astype(dtype)


def probes_to_wave(datafile, dtype=_np.int16, path=None):
    """ Make wav files from probes.

    Parameters
    ----------

    datafile: the hdf5 file
    dtype: type of the output wave file (np.int8, np.int16)
    path: path to save wavefiles

    Returns
    -------

    out: list of probe signals saved in the wave files
    """

    probes = []
    home = _os.path.expanduser('~')
    datafile = datafile.replace('~', home)

    with _DataExtractor(datafile) as data:
        probes_values = data.get_dataset('probes_value')
        p0 = data.get_attr('p0')
        dt = data.get_attr('dt')

    if not path:
        path = _os.path.dirname(datafile)


    if list(probes_values):

        filename = _os.path.basename(datafile).split('.')[0]

        for i, p in enumerate(probes_values):
            name = f'{path}/{filename}_probe_{i}.wav'
            tmp = normalize(p - p0, dtype=dtype)
            print(f'Writing {name}...')
            _wf.write(name, rate=int(1/dt), data=tmp)
            probes.append(tmp)
        print('Done!')

    else:
        print('No probes found!')

    return probes
