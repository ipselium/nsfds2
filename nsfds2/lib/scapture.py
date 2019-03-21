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

@author: Cyril Desjouy
"""

import numpy as np
from ofdlib2.capture import sigma_p, sigma_d, fo2
from ofdlib2.capture import xcapture, zcapture, update
from ofdlib2.fdtd import comp_p
import ofdlib.coefficients as cf

class ShockCapture:
    """ Shock Capturing procedure. (Bogey & al. -- JCP 228 -- 2009)"""

    def __init__(self, msh, fld, cfg, cin):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg
        self.cin = cin
        self.c_sc = cf.shock()

    def apply(self):
        """ Run shock capture. """
        sigmax, sigmaz = self.detector()
        self.new_fields(sigmax, sigmaz)

    def dilatation(self):
        """ Compute dilatation : Theta = div.v. """

        self.fld.E = self.fld.rhou/self.fld.rho
        self.fld.F = self.fld.rhov/self.fld.rho

        for sub in self.msh.xdomains:
            name = self.cin.cin_id(sub, self.msh.stencil)
            self.fld.K = getattr(self.cin, name)(self.fld.E, self.fld.K, sub)

        for sub in self.msh.zdomains:
            name = self.cin.cin_id(sub, self.msh.stencil)
            self.fld.K = getattr(self.cin, name)(self.fld.F, self.fld.K, sub)

        return self.fld.K

    def detector(self):
        """ Detector determination.

        * Determine the high frequency components of pressure (or dilatation)
        with a high pass filter (fo2)
        * Determine the filtering magnitude (sigmax, sigmaz)
        """

        if self.cfg.scapt_meth == 'pressure':
            comp_p(self.fld.p, self.fld.rho, self.fld.rhou,
                   self.fld.rhov, self.fld.rhoe, self.cfg.gamma)
            dpx, dpz = fo2(self.fld.p)
            sigmax, sigmaz = sigma_p(self.fld.p, dpx, dpz,
                                     self.cfg.rth, self.cfg.eps_machine)
        elif self.cfg.scapt_meth == 'dilatation':
            theta = self.dilatation()
            dpx, dpz = fo2(theta)
            sigmax, sigmaz = sigma_d(self.fld.p, self.fld.rho, dpx, dpz,
                                     self.cfg.gamma, self.cfg.rth,
                                     self.cfg.eps_machine)
        return sigmax, sigmaz

    def new_fields(self, sigmax, sigmaz):
        """ Apply shock capture. """

        # Capture following x
        xcapture(self.fld.rho, self.fld.K, sigmax, self.c_sc)
        xcapture(self.fld.rhou, self.fld.Ku, sigmax, self.c_sc)
        xcapture(self.fld.rhov, self.fld.Kv, sigmax, self.c_sc)
        xcapture(self.fld.rhoe, self.fld.Ke, sigmax, self.c_sc)
        self.update()

        # Capture following z
        zcapture(self.fld.rho, self.fld.K, sigmaz, self.c_sc)
        zcapture(self.fld.rhou, self.fld.Ku, sigmaz, self.c_sc)
        zcapture(self.fld.rhov, self.fld.Kv, sigmaz, self.c_sc)
        zcapture(self.fld.rhoe, self.fld.Ke, sigmaz, self.c_sc)
        self.update()

    def update(self):
        """ Apply shock capture. """

        idx = [0, self.msh.nx]
        idz = [0, self.msh.nz]
        update(self.fld.rho, self.fld.K, *idx, *idz)
        update(self.fld.rhou, self.fld.Ku, *idx, *idz)
        update(self.fld.rhov, self.fld.Kv, *idx, *idz)
        update(self.fld.rhoe, self.fld.Ke, *idx, *idz)
