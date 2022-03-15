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
# Creation Date : 2019-03-21 - 23:43:11
"""
-----------

Utils : Files

-----------
"""

import os
import sys
import numpy as _np
from nsfds2.utils.graphics import get_data as _get_data
from fdgrid import templates as _tplt


def get_obstacle(cfg):
    """ Get obstacle from custom file or fdgrid templates. """

    if cfg.geofile != 'None':
        try:
            sys.path.append(os.path.dirname(cfg.geofile))
            custom = __import__(os.path.basename(cfg.geofile).split('.')[0])
            obstacle = getattr(custom, cfg.geoname)(cfg.nx, cfg.nz)
            cfg.geoflag = True
        except (AttributeError, ImportError) as e:
            cfg.geoflag = False
            obstacle = []
    else:
        try:
            obstacle = getattr(_tplt, cfg.geoname)(cfg.nx, cfg.nz)
            cfg.geoflag = True
        except AttributeError as e:
            cfg.geoflag = False
            obstacle = []
    return obstacle


def get_curvilinear(cfg):
    """ Get curvilinear fonction from custom file or fdgrid templates. """

    if cfg.geofile != 'None':
        try:
            sys.path.append(os.path.dirname(cfg.geofile))
            custom = __import__(os.path.basename(cfg.geofile).split('.')[0])
            fcurv = getattr(custom, cfg.curvname)
            cfg.curvflag = True
        except (AttributeError, ImportError) as e:
            cfg.curvflag = False
            fcurv = None
    else:
        try:
            fcurv = getattr(_tplt, cfg.curvname)
            cfg.curvflag = True
        except AttributeError as e:
            cfg.curvflag = False
            fcurv = None
    return fcurv


def get_wall_function(cfg, name):
    """ Get function to apply to wall source."""

    try:
        sys.path.append(os.path.dirname(cfg.geofile))
        custom = __import__(os.path.basename(cfg.geofile).split('.')[0])
        func = getattr(custom, name)
    except (AttributeError, ImportError):
        func = None
    return func


def save_probes(filename):
    """Save probes in npz file and return file content."""
    data = _get_data(filename)
    nt, dt = data.attrs['nt'], data.attrs['dt']
    x, z = data['x'][...], data['z'][...]
    probe_values = data['probe_values'][...] - data.attrs['p0']
    probe_locs = data['probe_locations'][...]

    t = _np.arange(0, nt*dt, dt)
    _np.savez_compressed(filename.split('.')[0] + '_probes',
                         values=probe_values, locs=probe_locs, t=t, x=x, z=z)
    return _np.load(filename.split('.')[0] + '_probes.npz')
