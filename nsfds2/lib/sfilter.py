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
from ofdlib2.filters import filt, fx11c, fz11c
from ofdlib2.filters import fx711d, fz711d, fx4d, fz4d
from ofdlib2.filters import fx11xx, fz11xx

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
        return 'f{}{}{}'.format(stencil, subdomain.axis, subdomain.bc.replace('.', ''))

    def apply(self):
        """ Dispatch filtering. """

        for sub in self.msh.xdomains:
            sub.filt_method(self.fld.rho, self.fld.K, sub)
            sub.filt_method(self.fld.rhou, self.fld.Ku, sub)
            sub.filt_method(self.fld.rhov, self.fld.Kv, sub)
            sub.filt_method(self.fld.rhoe, self.fld.Ke, sub)

        self.xupdate()

        for sub in self.msh.zdomains:
            sub.filt_method(self.fld.rho, self.fld.K, sub)
            sub.filt_method(self.fld.rhou, self.fld.Ku, sub)
            sub.filt_method(self.fld.rhov, self.fld.Kv, sub)
            sub.filt_method(self.fld.rhoe, self.fld.Ke, sub)

        self.zupdate()

    def xupdate(self):
        pass

    def zupdate(self):
        pass

    def f11xRR(self, u, K, sub):
        """
            Selective filter : attenuate high frequencies
            4, 7 and 11 points filters
        """

        fx4d(u, K, self.xnu[1], *sub.ix, *sub.iz) # 1st and last points : 4pts scheme
        fx711d(u, K, *sub.ix, *sub.iz)             # Other points : 7 & 11pts scheme
        fx11c(u, K, *sub.ix, *sub.iz)
        # Update all but the first and last points
        filt(u, K, self.xnu[0], *[sub.ix[0]+1, sub.ix[1]], *sub.iz)

    def f11zRR(self, u, K, sub):
        """ Selective filter. """

        fz4d(u, K, self.xnu[1], *sub.ix, *sub.iz) # 1st and last points : 4pts scheme
        fz711d(u, K, *sub.ix, *sub.iz)             # Other points : 7 & 11pts scheme
        fz11c(u, K, *sub.ix, *sub.iz)
        # Update all but the first and last points
        filt(u, K, self.xnu[0], *sub.ix, *[sub.iz[0]+1, sub.iz[1]])

    def f11xXX(self, u, K, sub):
        fx11xx(u, K, self.xnu[0], *sub.ix, sub.iz[0])

    def f11zXX(self, u, K, sub):
        fz11xx(u, K, self.xnu[0], sub.ix[0], *sub.iz)

