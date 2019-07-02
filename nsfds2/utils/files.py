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
from fdgrid import templates as _tplt


def get_obstacle(cfg):
    """ Get obstacle from custom file or fdgrid templates. """

    if cfg.geofile != 'None':
        try:
            sys.path.append(os.path.dirname(cfg.geofile))
            custom = __import__(os.path.basename(cfg.geofile).split('.')[0])
            obstacle = getattr(custom, cfg.geoname)(cfg.nx, cfg.nz)
            cfg.geoflag = True
        except (AttributeError, ImportError):
            cfg.geoflag = False
            obstacle = []
    else:
        try:
            obstacle = getattr(_tplt, cfg.geoname)(cfg.nx, cfg.nz)
            cfg.geoflag = True
        except AttributeError:
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
        except (AttributeError, ImportError):
            cfg.curvflag = False
            fcurv = None
    else:
        try:
            fcurv = getattr(_tplt, cfg.curvname)
            cfg.curvflag = True
        except AttributeError:
            cfg.curvflag = False
            fcurv = None
    return fcurv
