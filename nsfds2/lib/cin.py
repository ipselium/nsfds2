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
#
# pylint: disable=too-few-public-methods
"""
-----------

Compute Derivatives

-----------
"""


import ofdlib2.derivation as drv


class Cin:
    """ Dispatch domains to associated computation functions. """

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg
        self.du = drv.du(msh.x, msh.z, msh.stencil)

        for sub in self.msh.dmdomains:
            bc = sub.bc.replace('.', '').replace('V', 'W')
            sub.cin_method = getattr(self.du, f'dud{sub.axname}_{bc}')

    def dispatch(self):
        """ Dispatch the domains to the functions. """

        # dE/dz ---------------------------------------------------------------
        if self.cfg.mesh in ['regular', 'adaptative']:
            self.fld.fdtools.Eu(self.fld.E, self.fld.Eu, self.fld.Ev, self.fld.Ee,
                                self.fld.r, self.fld.ru, self.fld.rv, self.fld.re,
                                self.fld.p)

        elif self.cfg.mesh == 'curvilinear':
            self.fld.fdtools.EuJ(self.fld.E, self.fld.Eu, self.fld.Ev, self.fld.Ee,
                                 self.fld.r, self.fld.ru, self.fld.rv, self.fld.re,
                                 self.fld.p, self.msh.dxn_dxp, self.msh.dxn_dzp)

        for sub in self.msh.dxdomains:
            sub.cin_method(self.fld.E, self.fld.K, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Eu, self.fld.Ku, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Ev, self.fld.Kv, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Ee, self.fld.Ke, *sub.ix, *sub.iz)

        # dF/dz ---------------------------------------------------------------
        if self.cfg.mesh in ['regular', 'adaptative']:
            self.fld.fdtools.Fu(self.fld.F, self.fld.Fu, self.fld.Fv, self.fld.Fe,
                                self.fld.r, self.fld.ru, self.fld.rv, self.fld.re,
                                self.fld.p)

        elif self.cfg.mesh == 'curvilinear':
            self.fld.fdtools.FuJ(self.fld.F, self.fld.Fu, self.fld.Fv, self.fld.Fe,
                                 self.fld.r, self.fld.ru, self.fld.rv, self.fld.re,
                                 self.fld.p, self.msh.dzn_dxp, self.msh.dzn_dzp)

        for sub in self.msh.dzdomains:
            sub.cin_method(self.fld.F, self.fld.K, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Fu, self.fld.Ku, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Fv, self.fld.Kv, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Fe, self.fld.Ke, *sub.ix, *sub.iz)
