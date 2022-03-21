#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2016-2020 Cyril Desjouy <cyril.desjouy@univ-lemans.fr>
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
# Creation Date : 2022-02-24 - 20:50:15
#
# pylint: disable=too-few-public-methods
"""
-----------

Compute Vorticity field

-----------
"""

import ofdlib2.derivation as drv


class Vorticity:
    """ Dispatch domains to computate vorticity. """

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg
        self.du = drv.du(msh.x, msh.z, cfg.stencil, cpu=cfg.cpu, add=False)

        for sub in self.msh.dmdomains:
            bc = sub.bc.replace('.', '').replace('V', 'W')
            sub.du = getattr(self.du, f'dud{sub.axname}_{bc}')

    def compute(self):
        """ Compute vorticity. """

        # dvz/dx
        for sub in self.msh.dxdomains:
            sub.du(self.fld.rv/self.fld.r, self.fld.E, *sub.ix, *sub.iz)

        # dvx/dz
        for sub in self.msh.dzdomains:
            sub.du(self.fld.ru/self.fld.r, self.fld.F, *sub.ix, *sub.iz)

        self.fld.vxz = self.fld.E - self.fld.F
