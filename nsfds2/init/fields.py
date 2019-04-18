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

import h5py
import numpy as np


class Fields:
    """ Fields initialization. """

    def __init__(self, msh, cfg):

        self._msh = msh
        self._cfg = cfg

        self.init_fields()
        self.init_derivatives()
        self.init_filters()
        if self._cfg.save:
            self.init_save()
        if 'A' in self._msh.bc:
            self.init_pml()

    def init_fields(self):
        """ Setup initial fields. """

        self.p = np.zeros(self._msh.shape)
        self.r = np.zeros_like(self.p) + self._cfg.rho0
        self.ru = np.zeros_like(self.p)
        self.rv = np.zeros_like(self.p)
        self.dltn = np.zeros_like(self.p)

        # Location
        ixS = self._cfg.ixS
        izS = self._cfg.izS
        Bx = 5*self._msh.dx                            # spatial bandwidth
        S0 = self._cfg.S0

        for iz in range(0, self._msh.nz):
            self.p[:, iz] = self._cfg.p0 + \
                    S0*np.exp(-np.log(2)*((self._msh.x-self._msh.x[ixS])**2 +
                                          (self._msh.z[iz]-self._msh.z[izS])**2)/Bx**2)
        self.re = self.p/(self._cfg.gamma-1.)

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

    def init_filters(self):
        """ Init variables used for filtering. """

        self.dp = np.zeros_like(self.p)
        self.sg = np.zeros_like(self.p)
        self.tau11 = np.zeros_like(self.p)
        self.tau22 = np.zeros_like(self.p)
        self.tau12 = np.zeros_like(self.p)

    def init_pml(self):
        """ Init PMLs. """

        # sigmax
        Dx = self._msh.Npml*self._msh.dx                      # Width of the PML
        self.sx = np.zeros(self._msh.nx)
        self.sx[:self._msh.Npml] = self._cfg.sigmax*abs((self._msh.x[:self._msh.Npml] \
                                 - self._msh.x[self._msh.Npml])/Dx)**self._cfg.alpha
        self.sx[self._msh.nx - self._msh.Npml:] = self.sx[self._msh.Npml-1::-1]

        # sigmaz
        Dz = self._msh.Npml*self._msh.dz                      # Width of the PML
        self.sz = np.zeros(self._msh.nz)
        self.sz[:self._msh.Npml] = self._cfg.sigmaz*abs((self._msh.z[:self._msh.Npml] \
                                 - self._msh.z[self._msh.Npml])/Dz)**self._cfg.alpha
        self.sz[self._msh.nz - self._msh.Npml:] = self.sz[self._msh.Npml-1::-1]


        # Init Q
        self.qx = np.zeros_like(self.p)
        self.qux = np.zeros_like(self.p)
        self.qvx = np.zeros_like(self.p)
        self.qex = np.zeros_like(self.p)
        self.qz = np.zeros_like(self.p)
        self.quz = np.zeros_like(self.p)
        self.qvz = np.zeros_like(self.p)
        self.qez = np.zeros_like(self.p)

        # Init K
        self.Kx = np.zeros_like(self.p)
        self.Kux = np.zeros_like(self.p)
        self.Kvx = np.zeros_like(self.p)
        self.Kex = np.zeros_like(self.p)
        self.Kz = np.zeros_like(self.p)
        self.Kuz = np.zeros_like(self.p)
        self.Kvz = np.zeros_like(self.p)
        self.Kez = np.zeros_like(self.p)

        # Initial E & F
        self.Ei = self.ru
        self.Eui = self.ru**2/self.r + self.p
        self.Evi = self.ru*self.rv/self.r
        self.Eei = self.ru/self.r*(self.re + self.p)

        self.Fi = self.rv
        self.Fui = self.ru*self.rv/ self.r
        self.Fvi = self.rv**2/self.r + self.p
        self.Fei = self.rv/self.r*(self.re + self.p)

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
        return self.p, self.r, self.ru, self.rv, self.re
