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


from ofdlib2 import fdtd
import ofdlib2.derivation as drv

class Cin:
    """ Dispatch domains to associated computation functions. """

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg
        cls = getattr(drv, 'du{}'.format(msh.stencil))
        self.du = cls(msh.x, msh.z)

        for sub in self.msh.dmdomains:
            fname = self.cin_id(sub)
            sub.cin_method = getattr(self.du, fname)

    @staticmethod
    def cin_id(sub):
        """ Identify which computation function to use for subdomain. """

        return 'dud{}_{}'.format(sub.axname, sub.bc.replace('.', ''))

    def dispatch(self):
        """ Dispatch the domains to the functions. """

        # dE/dz ---------------------------------------------------------------
        if self.cfg.mesh in ['regular', 'adaptative']:
            fdtd.Eu(self.fld.E, self.fld.Eu, self.fld.Ev, self.fld.Ee,
                    self.fld.r, self.fld.ru, self.fld.rv, self.fld.re, self.fld.p)

        elif self.cfg.mesh == 'curvilinear':
            fdtd.EuJ(self.fld.E, self.fld.Eu, self.fld.Ev, self.fld.Ee,
                     self.fld.F, self.fld.Fu, self.fld.Fv, self.fld.Fe,
                     self.fld.r, self.fld.ru, self.fld.rv, self.fld.re, self.fld.p,
                     self.msh.dxn_dxp, self.msh.dxn_dzp)

        for sub in self.msh.dxdomains:
            sub.cin_method(self.fld.E, self.fld.K, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Eu, self.fld.Ku, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Ev, self.fld.Kv, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Ee, self.fld.Ke, *sub.ix, *sub.iz)

        # dF/dz ---------------------------------------------------------------
        if self.cfg.mesh in ['regular', 'adaptative']:
            fdtd.Fu(self.fld.F, self.fld.Fu, self.fld.Fv, self.fld.Fe,
                    self.fld.r, self.fld.ru, self.fld.rv, self.fld.re, self.fld.p)

        elif self.cfg.mesh == 'curvilinear':
            fdtd.FuJ(self.fld.E, self.fld.Eu, self.fld.Ev, self.fld.Ee,
                     self.fld.F, self.fld.Fu, self.fld.Fv, self.fld.Fe,
                     self.fld.r, self.fld.ru, self.fld.rv, self.fld.re, self.fld.p,
                     self.msh.dzn_dxp, self.msh.dzn_dzp)

        for sub in self.msh.dzdomains:
            sub.cin_method(self.fld.F, self.fld.K, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Fu, self.fld.Ku, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Fv, self.fld.Kv, *sub.ix, *sub.iz)
            sub.cin_method(self.fld.Fe, self.fld.Ke, *sub.ix, *sub.iz)
