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
# Creation Date : 2018-04-14 01:46:17
"""
-----------

Compute Eulerian fluxes

@author: Cyril Desjouy
"""


import numpy as np
from ofdlib.fdtdc import time_advance, pnl
from nsfds2.utils.array import empty_like
from .cin import Cin

class EulerianFluxes:
    """ Compute Eulerian fluxes. """

    def __init__(self, msh, fld, cfg, cff):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg
        self.cff = cff
        self.p, self.rho, self.rhou, self.rhov, self.rhoe = empty_like(msh.shape, 5)
        self.ccin = Cin(msh, fld, cff)

    def rk4(self, p, rho, rhou, rhov, rhoe):
        """
          Avancement de la solution en temps à l'aide d'un algortihme de
          Runge-Kutta à 6 étapes
        """

        self.p = p.copy()
        self.rho = rho.copy()
        self.rhou = rhou.copy()
        self.rhov = rhov.copy()
        self.rhoe = rhoe.copy()

        for irk in range(1, 7):
            self.cin()

            self.rho = time_advance(np.empty_like(p), rho, self.fld.K, self.cff.rk[irk],
                                    self.cfg.dt, self.msh.nx, self.msh.nz)
            self.rhou = time_advance(np.empty_like(p), rhou, self.fld.Ku, self.cff.rk[irk],
                                     self.cfg.dt, self.msh.nx, self.msh.nz)
            self.rhov = time_advance(np.empty_like(p), rhov, self.fld.Kv, self.cff.rk[irk],
                                     self.cfg.dt, self.msh.nx, self.msh.nz)
            self.rhoe = time_advance(np.empty_like(p), rhoe, self.fld.Ke, self.cff.rk[irk],
                                     self.cfg.dt, self.msh.nx, self.msh.nz)

            self.cout()

            self.p = pnl(self.msh.nx, self.msh.nz,
                         self.p, self.rho, self.rhou, self.rhov, self.rhoe, self.cfg.gamma)



        return self.p, self.rho, self.rhou, self.rhov, self.rhoe

    def cin(self):
        """ Interior domain. """
        self.ccin.dispatch(self.p, self.rho, self.rhou, self.rhov, self.rhoe)

    def cout(self):
        """ Boundaries. """

        for s in self.msh.obstacles:

            self.rhou[s.xz[0]:s.xz[2]+1, s.xz[1]] = 0
            self.rhov[s.xz[0]:s.xz[2]+1, s.xz[1]] = 0

            self.rhou[s.xz[0]:s.xz[2]+1, s.xz[3]] = 0
            self.rhov[s.xz[0]:s.xz[2]+1, s.xz[3]] = 0

            self.rhou[s.xz[0], s.xz[1]:s.xz[3]+1] = 0
            self.rhov[s.xz[0], s.xz[1]:s.xz[3]+1] = 0

            self.rhou[s.xz[2], s.xz[1]:s.xz[3]+1] = 0
            self.rhov[s.xz[2], s.xz[1]:s.xz[3]+1] = 0

        self.rhou[[0, -1], :] = 0
        self.rhov[[0, -1], :] = 0

        self.rhou[:, [0, -1]] = 0
        self.rhov[:, [0, -1]] = 0
