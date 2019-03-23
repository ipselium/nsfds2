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
# Creation Date : 2018-04-13 - 00:48:58
"""
-----------

Compute viscous fluxes

@author: Cyril Desjouy
"""


import numpy as np
from ofdlib2.fdtd import integrate
from ofdlib2.derivation import dudx3RR, dudz3RR, dudxz3c
from ofdlib2.vflux import cErhou, cErhov, cErhoe, cFrhov, cFrhoe
from ofdlib2.utils import cmult, cdiv


class ViscousFluxes:
    """ Compute viscous flux : update rhou, rhov and rhoe only. """

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg

    def dispatch(self):
        """ Dispatch domains to integrate. """
        pass

    def integrate(self):
        """
            Viscous flux integration : interior points [optimized]
        """

        idx = [0, self.msh.nx-1]
        idz = [0, self.msh.nz-1]

        self.fld.Ku[:, :] = 0
        self.fld.Kv[:, :] = 0
        self.fld.Ke[:, :] = 0
        self.fld.Fu[:, :] = 0
        self.fld.Fv[:, :] = 0
        self.fld.Fe[:, :] = 0
        self.fld.Ee[:, :] = 0

        # dE/dx term
        self.fld.Eu = cdiv(self.fld.rhou, self.fld.rho)
        self.fld.Ev = cdiv(self.fld.rhov, self.fld.rho)

        tau11 = np.zeros_like(self.fld.p)
        tau22 = np.zeros_like(self.fld.p)
        tau12 = np.zeros_like(self.fld.p)

        # Strain tensor : WARNING : CHECK THAT DUDX3RR CORRESPONDS TO THIS NEED !
        dudx3RR(self.fld.Eu, tau11, self.msh.one_dx, *idx, *idz)
        dudz3RR(self.fld.Ev, tau22, self.msh.one_dz, *idx, *idz, False)
        dudxz3c(self.fld.Eu, self.fld.Ev, tau12, self.msh.one_dx, self.msh.dz, *idx, *idz)

        # Dynamic viscosity
        mu = cmult(self.fld.rho, self.cfg.nu)

        # dE/dx term
        self.fld.Eu = cErhou(mu, self.fld.Eu, tau11, tau22)
        self.fld.Ev = cErhov(mu, self.fld.Ev, tau12)
        self.fld.Ee = cErhoe(mu, self.fld.Ee, tau11, tau12, tau22, self.fld.rho, self.fld.rhou, self.fld.rhov)

        # dF/dz term
        self.fld.Fu = self.fld.Ev
        self.fld.Fv = cFrhov(mu, self.fld.Fv, tau11, tau22)
        self.fld.Fe = cFrhoe(mu, self.fld.Fe, tau11, tau12, tau22, self.fld.rho, self.fld.rhou, self.fld.rhov)

        # viscous flux : order 2 centered scheme
        idx = [1, self.msh.nx-1]
        idz = [1, self.msh.nz-1]

        dudx3RR(self.fld.Eu, self.fld.Ku, self.msh.one_dx, *idx, *idz)
        dudx3RR(self.fld.Ev, self.fld.Kv, self.msh.one_dx, *idx, *idz)
        dudx3RR(self.fld.Ee, self.fld.Ke, self.msh.one_dx, *idx, *idz)

        dudz3RR(self.fld.Fu, self.fld.Ku, self.msh.one_dz, *idx, *idz, True)
        dudz3RR(self.fld.Fv, self.fld.Kv, self.msh.one_dz, *idx, *idz, True)
        dudz3RR(self.fld.Fe, self.fld.Ke, self.msh.one_dz, *idx, *idz, True)

        # Integrate in time
        integrate(self.fld.rhou, self.fld.Ku, self.cfg.dt, *idx, *idz)
        integrate(self.fld.rhov, self.fld.Kv, self.cfg.dt, *idx, *idz)
        integrate(self.fld.rhoe, self.fld.Ke, self.cfg.dt, *idx, *idz)
