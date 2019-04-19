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
#
# Creation Date : 2019-03-07 - 23:54:38
"""
-----------

Some tools

@author: Cyril Desjouy
"""

import re
import time
import warnings
import numpy as np


class Check:
    """ Check FDTD parameters. """

    def __init__(self, cfg, msh):

        self.cfg = cfg
        self.msh = msh

        self.source()
        self.obstacles()
        self.domain()


    def source(self):
        """ Check if source is not in an obstacle. """
        for o in self.msh.obstacles:
            if o.ix[0] < self.cfg.ixS < o.ix[1] and o.iz[0] < self.cfg.izS < o.iz[1]:
                raise ValueError('source cannot be in an obstacle')

    def obstacles(self):
        """ Check validity of obstacles boundary conditions. """

        for obs in self.msh.obstacles:
            if obs.bc != 'RRRR':
                s = "Obstacles can only be 'RRRR' for now. "
                s += "Fix bcs to 'RRRR'."
                warnings.warn(s, stacklevel=8)
                obs.bc = 'RRRR'

    def domain(self):
        """ Check validity of the bcs. """

        if not re.match(r'[PRA][PRA][PRA][PRA]', self.msh.bc):
            s = "Only 'R', 'P', and 'A' bc are implemented for now. "
            s += "Fix bcs to 'RRRR'."
            warnings.warn(s, stacklevel=8)
            self.msh.bc = 'RRRR'


def timed(name):
    """ Time method of a class containing 'bench' attribute. """
    def layer(func):
        def wrapper(*args, **kwargs):
            start = time.perf_counter()
            func(*args, **kwargs)
            args[0].bench[name].append(time.perf_counter() - start)
        return wrapper
    return layer


def disp_bench(bench, it, residu):
    """ Display time spent at each step. """

    template = "Iteration : {0:5} | Res. : {1:.8f} | Time : {2:.4f} s."
    print(template.format(it, residu, np.mean(bench['total'])))
    template = "\t {:6} : {:.4f} s."
    if 'vfluxes' not in bench.keys():
        bench['vfluxes'] = 0
    if not bench['vfluxes']:
        bench['vfluxes'] = 0
    if not bench['save']:
        bench['save'] = 0
    print(template.format('EFlux', np.mean(bench['efluxes'])))
    if 0 not in bench['vfluxes']:
        print(template.format('VFlux', np.mean(bench['vfluxes'])))
    print(template.format('Filt', np.mean(bench['sfilt'])))
    print(template.format('Capt', np.mean(bench['scapt'])))
    print(template.format('Save', np.mean(bench['save'])))
    bench = {key: [] for key in bench.keys()}

    return bench
