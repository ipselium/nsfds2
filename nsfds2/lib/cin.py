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
from ofdlib2.derivation import dudx11c, dudx11p, dudx11m
from ofdlib2.derivation import dudz11c, dudz11p, dudz11m


class Cin:
    """ Dispatch domains to associated computation functions. """

    def __init__(self, msh, fld, cff):

        self.msh = msh
        self.nx, self.nz = msh.shape
        self.one_dx = msh.one_dx
        self.one_dz = msh.one_dz
        self.fld = fld
        self.ac = cff.ac
        self.ad = cff.ad

        for subdomain in self.msh.all_domains:
            fname = self.cin_id(subdomain, self.msh.stencil)
            subdomain.cin_method = getattr(self, fname)

    def dispatch(self, p, rho, rhou, rhov, rhoe):
        """ Dispatch the domains to the functions. """

        #######################################################################
        # dE/dz
        #######################################################################
        self.fld.E = rhou
        self.fld.Eu = cEuv(self.fld.Eu, rho, rhou, p)
        self.fld.Ev = cEvu(self.fld.Ev, rho, rhou, rhov)
        self.fld.Ee = cEe(self.fld.Ee, rho, rhou, rhoe, p)

        for sub in self.msh.xdomains:
            sub.cin_method(self.fld.E, self.fld.K, sub)
            sub.cin_method(self.fld.Eu, self.fld.Ku, sub)
            sub.cin_method(self.fld.Ev, self.fld.Kv, sub)
            sub.cin_method(self.fld.Ee, self.fld.Ke, sub)

        #######################################################################
        # dF/dz
        #######################################################################
        self.fld.F = rhov
        self.fld.Fu = self.fld.Ev
        self.fld.Fv = cEuv(self.fld.Fv, rho, rhov, p)
        self.fld.Fe = cEe(self.fld.Fe, rho, rhov, rhoe, p)

        for sub in self.msh.zdomains:
            sub.cin_method(self.fld.F, self.fld.K, sub)
            sub.cin_method(self.fld.Fu, self.fld.Ku, sub)
            sub.cin_method(self.fld.Fv, self.fld.Kv, sub)
            sub.cin_method(self.fld.Fe, self.fld.Ke, sub)

    @staticmethod
    def cin_id(subdomain, stencil):
        """ Identify which computation function to use for subdomain. """
        return 'dud{}{}{}'.format(subdomain.axis, stencil, subdomain.bc.replace('.', ''))

    def dudx11RR(self, u, K, subdomain):
        """ Rigid-Rigid following x with a 11 points scheme """

        idx = [subdomain.xz[0]+5, subdomain.xz[2]-4]
        idz = [subdomain.xz[1], subdomain.xz[3]+1]

        dudx11c(u, K, self.one_dx, *idx, *idz)
        dudx11p(u, K, self.one_dx, subdomain.xz[1], *idz)
        dudx11m(u, K, self.one_dx, subdomain.xz[2], *idz)

    def dudz11RR(self, u, K, subdomain):
        """ Rigid-Rigid following z with a 11 points scheme """

        idx = [subdomain.xz[0], subdomain.xz[2]+1]
        idz = [subdomain.xz[1]+5, subdomain.xz[3]-4]

        dudz11c(u, K, self.one_dz, *idx, *idz, True)
        dudz11p(u, K, self.one_dz, *idx, subdomain.xz[1], True)
        dudz11m(u, K, self.one_dz, *idx, subdomain.xz[3], True)

    def dudx11RX(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudx11XR(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudx11XX(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudx11PR(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudx11RP(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudx11PX(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudx11XP(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudz11RX(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudz11XR(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudz11XX(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudz11PR(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudz11RP(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudz11PX(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass

    def dudz11XP(self):
        """ Rigid-Rigid following z with a 11 points scheme """
        pass
