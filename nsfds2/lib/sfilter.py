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
# Creation Date : 2018-04-13 00:19:03
"""
-----------

The SelectiveFilter class provides a filtering procedure to reduce Gibbs
oscillations.

-----------
"""

import ofdlib2.filters as flt


class SelectiveFilter:
    """ Filter rho, rhou, rhov and rhoe. """

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.flt = flt.sfilter(msh.nx, msh.nz,
                               cfg.flt_xnu, cfg.flt_xnu0,
                               stencil=cfg.flt_stencil,
                               cpu=cfg.cpu)

        for sub in self.msh.fmdomains:
            bc = sub.bc.replace('.', '').replace('V', 'W')
            sub.filt_method = getattr(self.flt, f"f{sub.axname}_{bc}")

    def apply(self):
        """ Apply selective filter. """

        # Calculate filter following x
        for sub in self.msh.fxdomains:
            sub.filt_method(self.fld.r, self.fld.K, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.ru, self.fld.Ku, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.rv, self.fld.Kv, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.re, self.fld.Ke, *sub.ix, *sub.iz)

        # Apply filtering following x
        for sub in self.msh.fxdomains:
            self.flt.update(self.fld.r, self.fld.K, *sub.ix, *sub.iz)
            self.flt.update(self.fld.ru, self.fld.Ku, *sub.ix, *sub.iz)
            self.flt.update(self.fld.rv, self.fld.Kv, *sub.ix, *sub.iz)
            self.flt.update(self.fld.re, self.fld.Ke, *sub.ix, *sub.iz)

        # Calculate filter following x
        for sub in self.msh.fzdomains:
            sub.filt_method(self.fld.r, self.fld.K, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.ru, self.fld.Ku, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.rv, self.fld.Kv, *sub.ix, *sub.iz)
            sub.filt_method(self.fld.re, self.fld.Ke, *sub.ix, *sub.iz)

        # Apply filtering following z
        for sub in self.msh.fzdomains:
            self.flt.update(self.fld.r, self.fld.K, *sub.ix, *sub.iz)
            self.flt.update(self.fld.ru, self.fld.Ku, *sub.ix, *sub.iz)
            self.flt.update(self.fld.rv, self.fld.Kv, *sub.ix, *sub.iz)
            self.flt.update(self.fld.re, self.fld.Ke, *sub.ix, *sub.iz)
