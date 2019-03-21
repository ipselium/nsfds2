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
# Creation Date : 2019-03-01 - 10:29:55
"""
-----------

Compute Derivatives

@author: Cyril Desjouy
"""


import numpy as np
from ofdlib2.fdtd import cEuv, cEvu, cEe
import ofdlib2.derivation as drv

class Cin:
    """ Dispatch domains to associated computation functions. """

    def __init__(self, msh, fld):

        self.msh = msh
        self.one_dx = msh.one_dx
        self.one_dz = msh.one_dz
        self.fld = fld

        for sub in self.msh.all_domains:
            fname = self.cin_id(sub, self.msh.stencil)
            sub.cin_method = getattr(drv, fname)

    @staticmethod
    def cin_id(sub, stencil):
        """ Identify which computation function to use for subdomain. """
        return 'dud{}{}{}'.format(sub.axis, stencil, sub.bc.replace('.', ''))

    def dispatch(self, p, rho, rhou, rhov, rhoe):
        """ Dispatch the domains to the functions. """

        # dE/dz ---------------------------------------------------------------
        self.fld.E = rhou
        self.fld.Eu = cEuv(self.fld.Eu, rho, rhou, p)
        self.fld.Ev = cEvu(self.fld.Ev, rho, rhou, rhov)
        self.fld.Ee = cEe(self.fld.Ee, rho, rhou, rhoe, p)

        for sub in self.msh.xdomains:
            sub.cin_method(self.fld.E, self.fld.K, self.one_dx, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Eu, self.fld.Ku, self.one_dx, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Ev, self.fld.Kv, self.one_dx, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Ee, self.fld.Ke, self.one_dx, *sub.ix, *sub.iz)

        # dF/dz ---------------------------------------------------------------
        self.fld.F = rhov
        self.fld.Fu = self.fld.Ev
        self.fld.Fv = cEuv(self.fld.Fv, rho, rhov, p)
        self.fld.Fe = cEe(self.fld.Fe, rho, rhov, rhoe, p)

        for sub in self.msh.zdomains:
            sub.cin_method(self.fld.F, self.fld.K, self.one_dz, *sub.ix, *sub.iz, True)
            sub.cin_method(self.fld.Fu, self.fld.Ku, self.one_dz, *sub.ix, *sub.iz, True)
            sub.cin_method(self.fld.Fv, self.fld.Kv, self.one_dz, *sub.ix, *sub.iz, True)
            sub.cin_method(self.fld.Fe, self.fld.Ke, self.one_dz, *sub.ix, *sub.iz, True)
