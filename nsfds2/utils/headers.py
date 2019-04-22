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
# Creation Date : 2019-03-08 - 10:47:50
#
# pylint: disable=redefined-builtin
"""
-----------

Headers for nsfds2

@author: Cyril Desjouy
"""

import os
import sys
import platform
import textwrap
import numpy
import ofdlib2
import fdgrid
import nsfds2


def _columns():
    _, col = os.popen('stty size', 'r').read().split()
    return int(col) if int(col) < 81 else 80


def copyright():
    """ Show copyright. """
    col = 76 if _columns() > 76 else _columns
    cp = [col*'#',
          "nsfds2 v{} -- Copyright (C) 2016-2019 -- Cyril Desjouy".format(nsfds2.__version__),
          " "]
    lc = "This program comes with ABSOLUTELY NO WARRANTY. " + \
         "This is free software, and you are welcome to redistribute it " + \
         "under certain conditions; See GNU GPL v3 for more details."
    cp += textwrap.wrap(lc, int(col))
    cp += [col*"#"]
    tmp = ''
    for txt in cp:
        tmp += "# {:^{col}} #\n".format(txt, col=col)
    print(tmp)


def version():
    """ Display versions of modules in use. """

    print('System : {}\n'.format(platform.platform()))
    print('\tpython : {:<10}\tnumpy  : {:<10}'.format(sys.version.split(' ')[0], numpy.__version__))
    print('\tfdgrid : {:<10}\tofdlib : {:<10}'.format(fdgrid.__version__, ofdlib2.__version__))
    print('\n' + _columns()*"#" + "\n")


def start(cfg):
    """ Prompt start. """

    check_geo(cfg)

    log = "Starting computation for geometry : '{}' ({}x{})."
    print(log.format(cfg.geoname, cfg.nx, cfg.nz))
    key = input("Hit enter to continue (prefix 'q' to abort) ! ")
    if key == 'q':
        sys.exit(0)


def check_geo(cfg):
    """ Check geofile. """

    if not cfg.geoflag:
        print('Warning: Unable to load {} from {}'.format(cfg.geoname, cfg.geofile))
        cfg.geoname = 'empty'


def parameters(cfg):
    """ Show simulation parameters. """

    s = "geometry : '{}' ({}x{})".format(cfg.geoname, cfg.nx, cfg.nz)
    s += "\n\t* {} points PML : ".format(cfg.Npml)
    s += "sigma=({}, {}) and beta={}".format(cfg.sigmax, cfg.sigmaz, cfg.beta)
    s += "\n\t* source : {} at ({}, {}).".format(cfg.typ, cfg.ixS, cfg.izS)
    s += "\n\t* bc = {} (PML: {} points)".format(cfg.bc, cfg.Npml)
    s += "\n\t* dx = {} m\tdz = {} m.".format(cfg.dx, cfg.dz)
    s += "\n\t* dt = {} s\tnt = {}.".format(cfg.dt, cfg.nt)

    print(s)
