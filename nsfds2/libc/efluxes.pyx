# -*- coding: utf-8 -*-
#cython: language_level=3
#cython: wraparound=False
#cython: boundscheck=False
#cython: nonecheck=False
#cython: cdivision=True
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

cimport cython
import numpy as np
cimport numpy as np
import ofdlib2.fdtd as fdtd
import ofdlib2.coefficients as cf
from nsfds2.libc.cin import Cin


cdef class EulerianFluxes:
    """ Compute Eulerian fluxes. """

    cdef public object msh, fld, cfg, ccin
    cdef double[:] dtrk
    cdef double[:,::1] p, rho, rhou, rhov, rhoe

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg

        self.p = np.empty(msh.shape)
        self.rho = np.empty(msh.shape)
        self.rhou = np.empty(msh.shape)
        self.rhov = np.empty(msh.shape)
        self.rhoe = np.empty(msh.shape)

        self.ccin = Cin(msh, fld)
        self.dtrk = self.cfg.dt*cf.rk4o()

    def rk4(self, double[:,::1] p, double[:,::1] rho, double[:,::1] rhou,
            double[:,::1] rhov, double[:,::1] rhoe):
        """
          Avancement de la solution en temps à l'aide d'un algortihme de
          Runge-Kutta à 6 étapes
        """

        cdef int irk

        self.p = p.copy()
        self.rho = rho.copy()
        self.rhou = rhou.copy()
        self.rhov = rhov.copy()
        self.rhoe = rhoe.copy()

        for irk in range(1, 7):

            # Eulerian fluxes
            self.cin()

            # Intregration of eulerian fluxes : update self.rhox
            fdtd.time_advance(self.rho, rho, self.fld.K, self.dtrk[irk])
            fdtd.time_advance(self.rhou, rhou, self.fld.Ku, self.dtrk[irk])
            fdtd.time_advance(self.rhov, rhov, self.fld.Kv, self.dtrk[irk])
            fdtd.time_advance(self.rhoe, rhoe, self.fld.Ke, self.dtrk[irk])

            # Boundary conditions
            self.cout()

            # Compute p
            fdtd.comp_p(self.p, self.rho, self.rhou, self.rhov, self.rhoe, self.cfg.gamma)

        return np.asarray(self.p), np.asarray(self.rho), np.asarray(self.rhou), np.asarray(self.rhov), np.asarray(self.rhoe)

    def cin(self):
        """ Interior domain. """
        self.ccin.dispatch(self.p, self.rho, self.rhou, self.rhov, self.rhoe)

    def cout(self):
        """ Boundaries. """

        cdef object s

        for s in self.msh.obstacles:

            self.rhou[s.ix[0]:s.ix[1]+1, s.iz[0]] = 0
            self.rhov[s.ix[0]:s.ix[1]+1, s.iz[0]] = 0

            self.rhou[s.ix[0]:s.ix[1]+1, s.iz[1]] = 0
            self.rhov[s.ix[0]:s.ix[1]+1, s.iz[1]] = 0

            self.rhou[s.ix[0], s.iz[0]:s.iz[1]+1] = 0
            self.rhov[s.ix[0], s.iz[0]:s.iz[1]+1] = 0

            self.rhou[s.ix[1], s.iz[0]:s.iz[1]+1] = 0
            self.rhov[s.ix[1], s.iz[0]:s.iz[1]+1] = 0

        self.rhou[0, :] = 0
        self.rhov[0, :] = 0

        self.rhou[:, 0] = 0
        self.rhov[:, 0] = 0

        self.rhou[self.msh.nx-1, :] = 0
        self.rhov[self.msh.nx-1, :] = 0

        self.rhou[:, self.msh.nz-1] = 0
        self.rhov[:, self.msh.nz-1] = 0
