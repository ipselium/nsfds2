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
# Creation Date : 2019-03-07 - 21:31:14
"""
-----------

Post-treatment tools for nsfds2

@author: Cyril Desjouy
"""

import sys
import numpy as np
import numba as nb
from ofdlib2.coefficients import a7o
from ofdlib2 import fdtd
from ofdlib2.utils import cdiv

@nb.jit
def rot(data, iref, nx, nz, one_dx, one_dz, a7):
    """ Rotational. """
    name = "{}_it{}"
    rho = data[name.format('rho', iref)][:, :]
    rhou = data[name.format('rhou', iref)][:, :]
    rhov = data[name.format('rhov', iref)][:, :]

    u = rhou/rho
    v = rhov/rho

    vort = np.zeros_like(u)

    for i in range(3, nx-3):
        for j in range(3, nz-3):
            for l in range(-3, 4):
                vort[i, j] = vort[i, j] + a7[l]*v[i+l, j]*one_dx - a7[l]*u[i, j+l]*one_dz

    for i in range(3):
        for j in range(3, nz-3):
            vort[i, j] = vort[i, j] + (v[i+1, j] - v[i, j])*one_dx
            for l in range(-3, 4):
                vort[i, j] = vort[i, j] - a7[l]*u[i, j+l]*one_dz

    for i in range(nx-3, nx):
        for j in range(3, nz-3):
            vort[i, j] = vort[i, j] + (v[i, j] - v[i-1, j])*one_dx
            for l in range(-3, 4):
                vort[i, j] = vort[i, j] - a7[l]*u[i, j+l]*one_dz

    return vort


class FrameGenerator:
    """ Frame Genrator """

    def __init__(self, data, view, iref=20):

        self.data = data
        self.view = view
        self.imin = 0
        self.iref = iref
        self.ns = data['ns'][...]
        self.icur = self.imin
        self.var = {'p': 'p',
                    'rho': 'rho',
                    'vx': 'rhou',
                    'vz': 'rhov',
                    'e': 'rhoe'}
        self.a7, _ = a7o()
        self.one_dx = 1/self.data['dx'][...]
        self.one_dz = 1/self.data['dz'][...]
        self.nx = self.data['nx'][...]
        self.nz = self.data['nz'][...]

    def reference(self):
        """ Generate the reference for min/max colormap values """

        if self.view == "p":
            ref = self.p_from_rhoX(self.iref).T
        elif self.view in ['rho', 'vx', 'vz', 'e']:
            ref = self.data["{}_it{}".format(self.var[self.view], self.iref)][:, :].T
        elif self.view in ['vort']:
            ref = rot(self.data, self.iref, self.nx, self.nz, self.one_dx, self.one_dz, self.a7).T
        else:
            print("Only 'p', 'rho', 'vx', 'vz' and 'e' available !")
            sys.exit(1)

        return ref.max(), ref.min()

    def p_from_rhoX(self, i):
        """ Compute p from rho, rhou, rhov and rhoe. """

        rho = self.data["{}_it{}".format('rho', i)][:, :]
        rhou = self.data["{}_it{}".format('rhou', i)][:, :]
        rhov = self.data["{}_it{}".format('rhov', i)][:, :]
        rhoe = self.data["{}_it{}".format('rhoe', i)][:, :]
        p = np.empty_like(rho)
        fdtd.p(p, rho, rhou, rhov, rhoe, self.data['gamma'][...])
        return p - self.data['p0'][...]

    def next_item(self):
        """ Generate next value of variable """

        if self.view == 'p':
            return self.p_from_rhoX(self.icur).T

        elif self.view == 'rho':
            return self.data["{}_it{}".format(self.var[self.view], self.icur)][:, :].T

        elif self.view in ['vx', 'vz', 'e']:
            vX = self.data["{}_it{}".format(self.var[self.view], self.icur)][:, :]
            rho = self.data["{}_it{}".format('rho', self.icur)][:, :]
            return cdiv(vX, rho).T

        elif self.view in ['vort']:
            return rot(self.data, self.icur, self.nx, self.nz, self.one_dx, self.one_dz, self.a7).T

    def __iter__(self):
        """ Iterator """

        return self

    def __next__(self):
        """ Next element of iterator : (frame_number, variable) """

        try:
            self.icur += self.ns
            self.item = self.next_item()
            return self.icur, self.item
        except KeyError:
            raise StopIteration
