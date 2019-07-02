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

Timing utils

-----------
"""

import time
import numpy as np


def proceed(name):
    """ Time method of a class containing 'bench' attribute. """
    def layer(func):
        def wrapper(*args, **kwargs):
            start = time.perf_counter()
            func(*args, **kwargs)
            args[0].bench[name].append(time.perf_counter() - start)
        return wrapper
    return layer


def disp(bench, it, residu):
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
