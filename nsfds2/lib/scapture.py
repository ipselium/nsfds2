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
# Creation Date : 2019-03-06 - 10:31:53
"""
-----------

This module provides the ShockCapture class implementing the shock capturing
procedure proposed by Bogey & al -- JCP228 -- 2009

-----------
"""

import ofdlib2.derivation as drv
import ofdlib2.filters as flt


class ShockCapture:
    """ Shock Capturing procedure. (Bogey & al. -- JCP 228 -- 2009)"""

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg

        sg_cls = getattr(flt, 'sigma_{}'.format(cfg.cpt_meth[0]))

        self.du = drv.du(msh.x, msh.z, cfg.cpt_stencil)
        self.sg = sg_cls(msh.nx, msh.nz, cfg.rth, cfg.gamma)
        self.lpl = flt.lplf3(msh.nx, msh.nz)
        self.cpt = flt.capture(msh.nx, msh.nz)

        for sub in self.msh.fmdomains:
            bc = sub.bc.replace('.', '').replace('V', 'W')
            sub.dltn = getattr(self.du, f'dud{sub.axname}_{bc}')
            sub.lpl = getattr(self.lpl, f'lplf{sub.axname}_{bc}')
            sub.cpt = getattr(self.cpt, f'cpt{sub.axname}_{bc}')
            sub.sg = getattr(self.sg, f'sg{sub.axname}_{bc}')

    def apply(self):
        """ Run shock capture. """

        for direction in [self.msh.fxdomains, self.msh.fzdomains]:

            self.update_reference()

            for sub in direction:
                self.laplacian(sub)

            for sub in direction:
                self.filter_magnitude(sub)

            for sub in direction:
                self.filter(sub)

            for sub in direction:
                self.update(sub)

    def update_reference(self):
        """ Update pressure / dilatation. """

        self.fld.fdtools.p(self.fld.p, self.fld.r, self.fld.ru,
                           self.fld.rv, self.fld.re)

        if self.cfg.cpt_meth == 'dilatation':
            self.dilatation()

    def dilatation(self):
        """ Compute dilatation with 7 points scheme : dltn = div.v. """

        self.fld.E = self.fld.ru/self.fld.r
        self.fld.F = self.fld.rv/self.fld.r

        for sub in self.msh.fxdomains:
            sub.dltn(self.fld.E, self.fld.dltn, *sub.ix, *sub.iz)

        for sub in self.msh.fzdomains:
            sub.dltn(self.fld.F, self.fld.dltn, *sub.ix, *sub.iz)

    def laplacian(self, sub):
        """  Determine the high frequency components
        of pressure (or dilatation) with a high pass filter (laplacian filter)
        """

        if self.cfg.cpt_meth == 'pressure':
            sub.lpl(self.fld.p, self.fld.dp, *sub.ix, *sub.iz)

        elif self.cfg.cpt_meth == "dilatation":
            sub.lpl(self.fld.dltn, self.fld.dp, *sub.ix, *sub.iz)

    def filter_magnitude(self, sub):
        """ Determine the filtering magnitude (sg). """

        if self.cfg.cpt_meth == 'pressure':
            sub.sg(self.fld.p, self.fld.sg, self.fld.dp, *sub.ix, *sub.iz)

        elif self.cfg.cpt_meth == "dilatation":
            sub.sg(self.fld.p, self.fld.r, self.fld.sg, self.fld.dp, *sub.ix, *sub.iz)

    def filter(self, sub):
        """ Compute filter. """

        sub.cpt(self.fld.r, self.fld.K, self.fld.sg, *sub.ix, *sub.iz)
        sub.cpt(self.fld.ru, self.fld.Ku, self.fld.sg, *sub.ix, *sub.iz)
        sub.cpt(self.fld.rv, self.fld.Kv, self.fld.sg, *sub.ix, *sub.iz)
        sub.cpt(self.fld.re, self.fld.Ke, self.fld.sg, *sub.ix, *sub.iz)

    def update(self, sub):
        """ Apply filter to conservative variables. """

        self.cpt.update(self.fld.r, self.fld.K, *sub.ix, *sub.iz)
        self.cpt.update(self.fld.ru, self.fld.Ku, *sub.ix, *sub.iz)
        self.cpt.update(self.fld.rv, self.fld.Kv, *sub.ix, *sub.iz)
        self.cpt.update(self.fld.re, self.fld.Ke, *sub.ix, *sub.iz)
