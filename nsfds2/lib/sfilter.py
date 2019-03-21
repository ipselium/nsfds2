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
from ofdlib2.filters import fx11rr, fz11rr
from ofdlib2.filters import fx11xx, fz11xx


class SelectiveFilter:
    """ Filter rho, rhou, rhov and rhoe. """


    def __init__(self, msh, fld, cff):

        self.msh = msh
        self.fld = fld
        self.xnu = cff.xnu

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

        for sub in self.msh.zdomains:
            sub.filt_method(self.fld.rho, self.fld.K, sub)
            sub.filt_method(self.fld.rhou, self.fld.Ku, sub)
            sub.filt_method(self.fld.rhov, self.fld.Kv, sub)
            sub.filt_method(self.fld.rhoe, self.fld.Ke, sub)

    def f11xRR(self, u, K, sub):
        """ Rigid-Rigid : 11 points selective filtering following x. """
        fx11rr(u, K, self.xnu[0], *sub.ix, *sub.iz)

    def f11zRR(self, u, K, sub):
        """ Rigid-Rigid : 11 points selective filtering following z. """
        fz11rr(u, K, self.xnu[0], *sub.ix, *sub.iz)

    def f11xXX(self, u, K, sub):
        """ Open-Open : 11 points selective filtering following x. """
        fx11xx(u, K, self.xnu[0], *sub.ix, sub.iz[0])

    def f11zXX(self, u, K, sub):
        """ Open-Open : 11 points selective filtering following z. """
        fz11xx(u, K, self.xnu[0], sub.ix[0], *sub.iz)
