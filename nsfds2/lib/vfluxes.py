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
import ofdlib2.fdtd as fdtd
import ofdlib2.derivation as drv
import ofdlib2.vflux as vf
from ofdlib2.utils import cmult, cdiv


class ViscousFluxes:
    """ Compute viscous flux : update rhou, rhov and rhoe only. """

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg

        cls = getattr(drv, 'du{}'.format(self.cfg.vsc_stencil))
        self.du = cls(msh.nx, msh.nz, msh.one_dx, msh.one_dz)

        for sub in self.msh.dmdomains:

            bc = sub.bc.replace('.', '')
            name = 'dud{}_{}'.format(sub.axname, bc)
            sub.du = getattr(self.du, name)

    def integrate(self):
        """
            Viscous flux integration : interior points [optimized]
        """

        # dE/dx term
        self.fld.Eu = self.fld.ru/self.fld.r
        self.fld.Ev = self.fld.rv/self.fld.r

        # Strain tensor : WARNING : CHECK THAT DUDX3RR CORRESPONDS TO THIS NEED !
        for sub in self.msh.dxdomains:
            sub.du(self.fld.Eu, self.fld.tau11, *sub.ix, *sub.iz)
            sub.du(0.5*self.fld.Ev, self.fld.tau12, *sub.ix, *sub.iz)

        # Important : Init tau22 (because of adding previous tau22 in dudz !)
        self.fld.tau22[:, :] = 0
        for sub in self.msh.dzdomains:
            sub.du(self.fld.Ev, self.fld.tau22, *sub.ix, *sub.iz)
            sub.du(0.5*self.fld.Eu, self.fld.tau12, *sub.ix, *sub.iz)


        # Dynamic viscosity
        mu = self.fld.r*self.cfg.nu

        # dE/dx and dF/dz
        vf.dEF(self.fld.Eu, self.fld.Ev, self.fld.Ee,
               self.fld.Fu, self.fld.Fv, self.fld.Fe, mu,
               self.fld.tau11, self.fld.tau12, self.fld.tau22,
               self.fld.r, self.fld.ru, self.fld.rv)


        # viscous flux : order 2 centered scheme
        for sub in self.msh.dxdomains:
            sub.du(self.fld.Eu, self.fld.Ku, *sub.ix, *sub.iz)
            sub.du(self.fld.Ev, self.fld.Kv, *sub.ix, *sub.iz)
            sub.du(self.fld.Ee, self.fld.Ke, *sub.ix, *sub.iz)

        for sub in self.msh.dzdomains:
            sub.du(self.fld.Fu, self.fld.Ku, *sub.ix, *sub.iz)
            sub.du(self.fld.Fv, self.fld.Kv, *sub.ix, *sub.iz)
            sub.du(self.fld.Fe, self.fld.Ke, *sub.ix, *sub.iz)

        # Integrate in time
        for sub in self.msh.dsdomains:
            fdtd.integrate(self.fld.ru, self.fld.Ku, self.cfg.dt, *sub.ix, *sub.iz)
            fdtd.integrate(self.fld.rv, self.fld.Kv, self.cfg.dt, *sub.ix, *sub.iz)
            fdtd.integrate(self.fld.re, self.fld.Ke, self.cfg.dt, *sub.ix, *sub.iz)
