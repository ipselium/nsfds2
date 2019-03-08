#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2016-2019 Cyril Desjouy <ipselium@free.fr>
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
# Creation Date : 2018-04-13 00:19:03
"""
-----------

The SelectiveFilter class propose a filtering procedure to reduce Gibbs
oscillations.

@author: Cyril Desjouy
"""

import numpy as np
from ofdlib.fschmc import fxd, fzd, fx11c, fz11c, fx11d, fz11d
from ofdlib.fschmc import filt, xsingle_filt, zsingle_filt


class SelectiveFilter:
    """ Filter rho, rhou, rhov and rhoe. """


    def __init__(self, msh, fld, cff):

        self.msh = msh
        self.fld = fld
        self.xnu = cff.xnu
        self.d11c = cff.d11c
        self.d11d = cff.d11d
        self.d7d15 = cff.d7d15
        self.d4d03 = cff.d4d03

        for subdomain in self.msh.all_domains:
            fname = self.filt_id(subdomain, self.msh.stencil)
            subdomain.filt_method = getattr(self, fname)

    @staticmethod
    def filt_id(subdomain, stencil):
        """ Identify which filter function to use for subdomain. """
        return 'f{}{}'.format(stencil, subdomain.bc.replace('.', ''))

    def apply(self):
        """ Dispatch filtering. """

        for sub in self.msh.xdomains:
            sub.filt_method(self.fld.rho, sub)
            sub.filt_method(self.fld.rhou, sub)
            sub.filt_method(self.fld.rhov, sub)
            sub.filt_method(self.fld.rhoe, sub)

    def f11RR(self, u, sub):
        """
            Selective filter : attenuate high frequencies
            4, 7 and 11 points filters
        """

        idx = np.array([sub.xz[0], sub.xz[2]])
        idz = np.array([sub.xz[1], sub.xz[3]])

        # X FILTERING #########################################################
        K = np.zeros_like(u)
        # 1st and last points : 4 points scheme
        K = fxd(u, self.d4d03, idx[0], idz, idx[0], 4, 1, K)
        K = fxd(u, self.d4d03, idx[1], idz, idx[1], 4, -1, K)
        u = xsingle_filt(u, K, idx[0], idz, self.xnu[1])
        u = xsingle_filt(u, K, idx[1], idz, self.xnu[1])
        # Other points
        K = fx11d(u, K, self.d7d15, self.d11d, idx, idz)
        K = fx11c(u, self.d11c, np.array([idx[0]+5, idx[1]-4]), idz, K)
        u = filt(u, K, np.array([idx[0]+1, idx[1]]), idz, self.xnu[0])
        # Z FILTERING #########################################################
        K = np.zeros_like(u)
        # 1st and last points : 4 points scheme
        K = fzd(u, self.d4d03, idx, idz[0], idz[0], 4, 1, K)
        K = fzd(u, self.d4d03, idx, idz[1], idz[1], 4, -1, K)
        u = zsingle_filt(u, K, idx, idz[0], self.xnu[1])
        u = zsingle_filt(u, K, idx, idz[1], self.xnu[1])
        # Other points
        K = fz11d(u, K, self.d7d15, self.d11d, idx, idz)
        K = fz11c(u, self.d11c, idx, np.array([idz[0]+5, idz[1]-4]), K)
        u = filt(u, K, idx, np.array([idz[0]+1, idz[1]]), self.xnu[0])
