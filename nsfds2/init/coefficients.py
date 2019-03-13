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
# Creation Date : 2019-03-07 - 11:43:44
"""
-----------

Initialization of finite difference coefficients

@author: Cyril Desjouy
"""

import time
import ofdlib.coefficients as cf


class Coefficients:
    """ Finite difference scheme coefficients. """

    def __init__(self, cfg, stencil):

        ti = time.perf_counter()
        self.stencil = stencil
        self._cfg = cfg
        self.ac, self.ad = self.a_id(stencil)()
        self.dtrk = self._cfg.dt*cf.rk4o()
        self.xnu = cf.fxnu()*2.5                   # [0.7, 0.1] / Olivier : 0.8
        self.d11c, self.d11d = cf.d11om()
        _, self.d7d15 = cf.d7o()
        self.d4d03 = cf.d4o()
        self.c_sc = cf.shock()
        msg = 'Finite difference coefficients initialized in {:.4f} s.'
        print(msg.format(time.perf_counter() - ti))

    @staticmethod
    def a_id(stencil):
        """ Return coefficients corresponding to stencil. """
        return getattr(cf, 'a{}o'.format(stencil))
