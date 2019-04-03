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
from ofdlib2.fdtd import time_advance, comp_p
from nsfds2.utils.array import empty_like
from .cin import Cin
import ofdlib.coefficients as cf


class EulerianFluxes:
    """ Compute Eulerian fluxes. """

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg
        self.p, self.rho, self.rhou, self.rhov, self.rhoe = empty_like(msh.shape, 5)
        self.ccin = Cin(msh, fld)
        self.dtrk = self.cfg.dt*cf.rk4o()

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

            # Eulerian fluxes
            self.cin()

            # Intregration of eulerian fluxes : update self.rhox
            time_advance(self.rho, rho, self.fld.K, self.dtrk[irk])
            time_advance(self.rhou, rhou, self.fld.Ku, self.dtrk[irk])
            time_advance(self.rhov, rhov, self.fld.Kv, self.dtrk[irk])
            time_advance(self.rhoe, rhoe, self.fld.Ke, self.dtrk[irk])

            # Boundary conditions
            self.cout()

            # Compute p
            comp_p(self.p, self.rho, self.rhou, self.rhov, self.rhoe, self.cfg.gamma)

        return self.p, self.rho, self.rhou, self.rhov, self.rhoe

    def cin(self):
        """ Interior domain. """
        self.ccin.dispatch(self.p, self.rho, self.rhou, self.rhov, self.rhoe)

    def cout(self):
        """ Boundaries. """

        for s in self.msh.obstacles:

            self.rhou[s.sx, s.iz[0]] = 0
            self.rhov[s.sx, s.iz[0]] = 0

            self.rhou[s.sx, s.iz[1]] = 0
            self.rhov[s.sx, s.iz[1]] = 0

            self.rhou[s.ix[0], s.sz] = 0
            self.rhov[s.ix[0], s.sz] = 0

            self.rhou[s.ix[1], s.sz] = 0
            self.rhov[s.ix[1], s.sz] = 0

        if self.msh.bc[0] == 'R':
            self.rhou[0, :] = 0
            self.rhov[0, :] = 0

        if self.msh.bc[2] == 'R':
            self.rhou[-1, :] = 0
            self.rhov[-1, :] = 0

        if self.msh.bc[1] == 'R':
            self.rhou[:, 0] = 0
            self.rhov[:, 0] = 0

        if self.msh.bc[3] == 'R':
            self.rhou[:, -1] = 0
            self.rhov[:, -1] = 0


        self.rhou[:103, -1] = 0
        self.rhov[:103, -1] = 0

        self.rhou[153:, -1] = 0
        self.rhov[153:, -1] = 0
