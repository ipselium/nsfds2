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
from ofdlib.fdtdc import cErhouv, cErhovu, cErhoe
from ofdlib.dschmc import dudx11c, dudx11dp, dudx11dm
from ofdlib.dschmc import dudz11c, dudz11dp, dudz11dm


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

        ###########################################################################
        # Terme dE/dz (no numba)
        ###########################################################################
        self.fld.E = rhou
        self.fld.Eu = cErhouv(self.nx, self.nz, self.fld.Eu, rho, rhou, p)
        self.fld.Ev = cErhovu(self.nx, self.nz, self.fld.Ev, rho, rhou, rhov)
        self.fld.Ee = cErhoe(self.nx, self.nz, self.fld.Ee, rho, rhou, rhoe, p)

        for sub in self.msh.xdomains:
            self.fld.K = sub.cin_method(self.fld.E, self.fld.K, sub)
            self.fld.Ku = sub.cin_method(self.fld.Eu, self.fld.Ku, sub)
            self.fld.Kv = sub.cin_method(self.fld.Ev, self.fld.Kv, sub)
            self.fld.Ke = sub.cin_method(self.fld.Ee, self.fld.Ke, sub)

        ###########################################################################
        # Terme dF/dz (no numba)
        ###########################################################################
        self.fld.F = rhov
        self.fld.Fu = self.fld.Ev
        self.fld.Fv = cErhouv(self.nx, self.nz, self.fld.Fv, rho, rhov, p)
        self.fld.Fe = cErhoe(self.nx, self.nz, self.fld.Fe, rho, rhov, rhoe, p)

        for sub in self.msh.zdomains:
            self.fld.K = sub.cin_method(self.fld.F, self.fld.K, sub)
            self.fld.Ku = sub.cin_method(self.fld.Fu, self.fld.Ku, sub)
            self.fld.Kv = sub.cin_method(self.fld.Fv, self.fld.Kv, sub)
            self.fld.Ke = sub.cin_method(self.fld.Fe, self.fld.Ke, sub)

    @staticmethod
    def cin_id(subdomain, stencil):
        """ Identify which computation function to use for subdomain. """
        return 'dud{}{}{}'.format(subdomain.axis, stencil, subdomain.bc.replace('.', ''))

    def dudx11RR(self, u, K, subdomain):
        """ Rigid-Rigid following x with a 11 points scheme """

        idx = np.array([subdomain.xz[0]+5, subdomain.xz[2]-4])
        idz = np.array([subdomain.xz[1], subdomain.xz[3]+1])

        K = dudx11c(u, self.ac, self.one_dx, idx, idz, K)
        K = dudx11dp(u, self.ad, self.one_dx, subdomain.xz[0], idz, K)
        K = dudx11dm(u, self.ad, self.one_dx, subdomain.xz[2], idz, K)

        return K

    def dudz11RR(self, u, K, subdomain):
        """ Rigid-Rigid following z with a 11 points scheme """

        idx = np.array([subdomain.xz[0], subdomain.xz[2]+1])
        idz = np.array([subdomain.xz[1]+5, subdomain.xz[3]-4])

        K = dudz11c(u, self.ac, self.one_dz, idx, idz, True, K)
        K = dudz11dp(u, self.ad, self.one_dz, idx, subdomain.xz[1], True, K)
        K = dudz11dm(u, self.ad, self.one_dz, idx, subdomain.xz[3], True, K)

        return K

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

    def linear(self, p, rho, rhou, rhov, rhoe):
        """ Cin linear test. """

        ###########################################################################
        # Terme dE/dz (no numba)
        ###########################################################################
        self.fld.E = rhou
        self.fld.Eu = cErhouv(self.nx, self.nz, self.fld.Eu, rho, rhou, p)
        self.fld.Ev = cErhovu(self.nx, self.nz, self.fld.Ev, rho, rhou, rhov)
        self.fld.Ee = cErhoe(self.nx, self.nz, self.fld.Ee, rho, rhou, rhoe, p)
        ###########################################################################
        # Calcul intérieur suivant x : schéma centré sur 11 points de iz = 5 à self.msh.nz-6
        ###########################################################################
        idx = np.array([5, self.nx-5])
        idz = np.array([0, self.nz])
        self.fld.K = dudx11c(self.fld.E, self.ac, self.one_dx, idx, idz, self.fld.K)
        self.fld.Ku = dudx11c(self.fld.Eu, self.ac, self.one_dx, idx, idz, self.fld.Ku)
        self.fld.Kv = dudx11c(self.fld.Ev, self.ac, self.one_dx, idx, idz, self.fld.Kv)
        self.fld.Ke = dudx11c(self.fld.Ee, self.ac, self.one_dx, idx, idz, self.fld.Ke)
        ###########################################################################
        # BC Bottom : schéma décentré sur 11 points
        ###########################################################################
        self.fld.K = dudx11dp(self.fld.E, self.ad, self.one_dx, 0, idz, self.fld.K)
        self.fld.Ku = dudx11dp(self.fld.Eu, self.ad, self.one_dx, 0, idz, self.fld.Ku)
        self.fld.Kv = dudx11dp(self.fld.Ev, self.ad, self.one_dx, 0, idz, self.fld.Kv)
        self.fld.Ke = dudx11dp(self.fld.Ee, self.ad, self.one_dx, 0, idz, self.fld.Ke)
        ###########################################################################
        # BC Top : schéma décentré sur 11 points
        ###########################################################################
        self.fld.K = dudx11dm(self.fld.E, self.ad, self.one_dx, self.nx-1, idz, self.fld.K)
        self.fld.Ku = dudx11dm(self.fld.Eu, self.ad, self.one_dx, self.nx-1, idz, self.fld.Ku)
        self.fld.Kv = dudx11dm(self.fld.Ev, self.ad, self.one_dx, self.nx-1, idz, self.fld.Kv)
        self.fld.Ke = dudx11dm(self.fld.Ee, self.ad, self.one_dx, self.nx-1, idz, self.fld.Ke)
        ###########################################################################
        # Terme dF/dz (no numba)
        ###########################################################################
        self.fld.F = rhov
        self.fld.Fu = self.fld.Ev
        self.fld.Fv = cErhouv(self.nx, self.nz, self.fld.Fv, rho, rhov, p)
        self.fld.Fe = cErhoe(self.nx, self.nz, self.fld.Fe, rho, rhov, rhoe, p)
        ###########################################################################
        # BC Bottom : schéma décentré sur 11 points
        ###########################################################################
        idx = np.array([0, self.nx])
        idz = np.array([5, self.nz-5])
        self.fld.K = dudz11dp(self.fld.F, self.ad, self.one_dz, idx, 0, True, self.fld.K)
        self.fld.Ku = dudz11dp(self.fld.Fu, self.ad, self.one_dz, idx, 0, True, self.fld.Ku)
        self.fld.Kv = dudz11dp(self.fld.Fv, self.ad, self.one_dz, idx, 0, True, self.fld.Kv)
        self.fld.Ke = dudz11dp(self.fld.Fe, self.ad, self.one_dz, idx, 0, True, self.fld.Ke)
        ###########################################################################
        # BC Top : schéma décentré sur 11 points
        ###########################################################################
        self.fld.K = dudz11dm(self.fld.F, self.ad, self.one_dz, idx, self.nz-1, True, self.fld.K)
        self.fld.Ku = dudz11dm(self.fld.Fu, self.ad, self.one_dz, idx, self.nz-1, True, self.fld.Ku)
        self.fld.Kv = dudz11dm(self.fld.Fv, self.ad, self.one_dz, idx, self.nz-1, True, self.fld.Kv)
        self.fld.Ke = dudz11dm(self.fld.Fe, self.ad, self.one_dz, idx, self.nz-1, True, self.fld.Ke)
        ###########################################################################
        # Calcul intérieur suivant z : schéma centré sur 11 points de iz = 5 à self.msh.nz-6
        ###########################################################################
        self.fld.K = dudz11c(self.fld.F, self.ac, self.one_dz, idx, idz, True, self.fld.K)
        self.fld.Ku = dudz11c(self.fld.Fu, self.ac, self.one_dz, idx, idz, True, self.fld.Ku)
        self.fld.Kv = dudz11c(self.fld.Fv, self.ac, self.one_dz, idx, idz, True, self.fld.Kv)
        self.fld.Ke = dudz11c(self.fld.Fe, self.ac, self.one_dz, idx, idz, True, self.fld.Ke)
