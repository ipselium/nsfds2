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
# Creation Date : 2019-03-01 - 14:09:42
"""
-----------

Initialize all fields

@author: Cyril Desjouy
"""

import os
import h5py
import numpy as np


class Fields:
    """ Fields initialization. """

    def __init__(self, msh, cfg):

        self._msh = msh
        self._cfg = cfg

        self.init_fields()
        self.init_derivatives()
        if self._cfg.save:
            self.init_save()

    def init_fields(self):
        """ Setup initial fields. """

        self.p = np.zeros(self._msh.shape)
        self.rho = np.zeros_like(self.p) + self._cfg.rho0
        self.rhou = np.zeros_like(self.p)
        self.rhov = np.zeros_like(self.p)

        # Location
        ixS = self._cfg.ixS
        izS = self._cfg.izS
        Bx = 5*self._msh.dx                            # spatial bandwidth
        S0 = self._cfg.S0

        for iz in range(0, self._msh.nz):
            self.p[:, iz] = self._cfg.p0 + S0*np.exp(-np.log(2)*((self._msh.x-self._msh.x[ixS])**2 +
                                                                 (self._msh.z[iz]-self._msh.z[izS])**2)/Bx**2)
        self.rhoe = self.p/(self._cfg.gamma-1.)

    def init_derivatives(self):
        """ Init derivatives. """

        self.K = np.zeros_like(self.p)
        self.Ku = np.zeros_like(self.p)
        self.Kv = np.zeros_like(self.p)
        self.Ke = np.zeros_like(self.p)

        self.E = np.empty_like(self.p)
        self.Eu = np.empty_like(self.p)
        self.Ev = np.empty_like(self.p)
        self.Ee = np.empty_like(self.p)

        self.F = np.empty_like(self.p)
        self.Fu = np.empty_like(self.p)
        self.Fv = np.empty_like(self.p)
        self.Fe = np.empty_like(self.p)

    def init_save(self):
        """ Init save. """
        self.PATH = "results/"
        self.sfile = h5py.File(self.PATH + self._cfg.filename + '.hdf5', 'w')
        self.sfile.create_dataset('x', data=self._msh.x, compression=self._cfg.comp)
        self.sfile.create_dataset('z', data=self._msh.z, compression=self._cfg.comp)
        self.sfile.create_dataset('dx', data=self._msh.dx)
        self.sfile.create_dataset('dz', data=self._msh.dz)
        self.sfile.create_dataset('dt', data=self._cfg.dt)
        self.sfile.create_dataset('nx', data=self._cfg.nx)
        self.sfile.create_dataset('nz', data=self._cfg.nz)
        self.sfile.create_dataset('nt', data=self._cfg.nt)
        self.sfile.create_dataset('ns', data=self._cfg.ns)
        self.sfile.create_dataset('p0', data=self._cfg.p0)
        self.sfile.create_dataset('rho0', data=self._cfg.rho0)
        self.sfile.create_dataset('gamma', data=self._cfg.gamma)
        self.sfile.create_dataset('obstacles', data=self._msh.get_obstacles())

    def get(self):
        """ Get initial fields as a tuple. """
        return self.p, self.rho, self.rhou, self.rhov, self.rhoe
