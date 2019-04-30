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
import ofdlib2.filters as filters

class SelectiveFilter:
    """ Filter rho, rhou, rhov and rhoe. """


    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        cls = getattr(filters, f'f{cfg.flt_stencil}')
        self.flt = cls(msh.nx, msh.nz, cfg.xnu)


        for subdomain in self.msh.fmdomains:
            fname = self.filt_id(subdomain, self.msh.stencil)
            subdomain.filt_method = getattr(self.flt, fname)

    @staticmethod
    def filt_id(sub, stencil):
        """ Identify which filter function to use for subdomain. """
        return f"f{sub.axname}_{sub.bc.replace('.', '')}"

    def apply(self):
        """ Dispatch filtering. """

        for sub in self.msh.fxdomains:
            sub.filt_method(self.fld.r, self.fld.K, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.ru, self.fld.Ku, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.rv, self.fld.Kv, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.re, self.fld.Ke, *sub.ix, *sub.iz)

        self.update(self.msh.fxdomains)

        for sub in self.msh.fzdomains:
            sub.filt_method(self.fld.r, self.fld.K, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.ru, self.fld.Ku, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.rv, self.fld.Kv, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.re, self.fld.Ke, *sub.ix, *sub.iz)

        self.update(self.msh.fzdomains)

    def update(self, domains):
        """ Update fields. """
        for sub in domains:
            self.flt.apply(self.fld.r, self.fld.K, *sub.ix, *sub.iz)
            self.flt.apply(self.fld.ru, self.fld.Ku, *sub.ix, *sub.iz)
            self.flt.apply(self.fld.rv, self.fld.Kv, *sub.ix, *sub.iz)
            self.flt.apply(self.fld.re, self.fld.Ke, *sub.ix, *sub.iz)
