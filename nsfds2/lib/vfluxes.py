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

        for sub in self.msh.all_domains:

            bc = sub.bc.replace('.', '')
            name = 'dud{}{}{}'.format(sub.axname, self.cfg.vsc_stencil, bc)
            sub.tau = getattr(drv, name)
            name = 'dud{}{}{}'.format(sub.axname, self.cfg.vsc_stencil, bc)
            sub.visc = getattr(drv, name)

    def integrate(self):
        """
            Viscous flux integration : interior points [optimized]
        """

        # dE/dx term
        self.fld.Eu = self.fld.rhou/self.fld.rho
        self.fld.Ev = self.fld.rhov/self.fld.rho

        # Strain tensor : WARNING : CHECK THAT DUDX3RR CORRESPONDS TO THIS NEED !
        for sub in self.msh.xdomains:
            sub.tau(self.fld.Eu, self.fld.tau11, self.msh.one_dx, *sub.ix, *sub.iz)
            sub.tau(0.5*self.fld.Ev, self.fld.tau12, self.msh.one_dx, *sub.ix, *sub.iz)

        for sub in self.msh.zdomains:
            sub.tau(self.fld.Ev, self.fld.tau22, self.msh.one_dz, *sub.ix, *sub.iz, False)
            sub.tau(0.5*self.fld.Eu, self.fld.tau12, self.msh.one_dz, *sub.ix, *sub.iz, True)


        # Dynamic viscosity
        mu = self.fld.rho*self.cfg.nu

        # dE/dx and dF/dz
        vf.dEF(self.fld.Eu, self.fld.Ev, self.fld.Ee,
               self.fld.Fu, self.fld.Fv, self.fld.Fe, mu,
               self.fld.tau11, self.fld.tau12, self.fld.tau22,
               self.fld.rho, self.fld.rhou, self.fld.rhov)


        # viscous flux : order 2 centered scheme
        for sub in self.msh.xdomains:
            sub.visc(self.fld.Eu, self.fld.Ku, self.msh.one_dx, *sub.ix, *sub.iz)
            sub.visc(self.fld.Ev, self.fld.Kv, self.msh.one_dx, *sub.ix, *sub.iz)
            sub.visc(self.fld.Ee, self.fld.Ke, self.msh.one_dx, *sub.ix, *sub.iz)

        for sub in self.msh.zdomains:
            sub.visc(self.fld.Fu, self.fld.Ku, self.msh.one_dz, *sub.ix, *sub.iz, True)
            sub.visc(self.fld.Fv, self.fld.Kv, self.msh.one_dz, *sub.ix, *sub.iz, True)
            sub.visc(self.fld.Fe, self.fld.Ke, self.msh.one_dz, *sub.ix, *sub.iz, True)

        # Integrate in time
        for i in self.msh.xdomains:
            fdtd.integrate(self.fld.rhou, self.fld.Ku, self.cfg.dt, *sub.ix, *sub.iz)
            fdtd.integrate(self.fld.rhov, self.fld.Kv, self.cfg.dt, *sub.ix, *sub.iz)
            fdtd.integrate(self.fld.rhoe, self.fld.Ke, self.cfg.dt, *sub.ix, *sub.iz)
