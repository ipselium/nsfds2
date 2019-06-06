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
#
# pylint: disable=too-many-instance-attributes
"""
-----------

Initialize all fields

@author: Cyril Desjouy
"""

import h5py
import numpy as np
from ofdlib2.fdtd import fdtools


def set_sigma(cfg, msh):
    """ Set sigmax and sigmaz. """

    try:
        cfg.sigmax = float(cfg.sigmax)
    except ValueError:
        D = msh.x[msh.Npml] - msh.x[0]
        cfg.sigmax = 8*cfg.c0*(cfg.alpha+1)/D

    try:
        cfg.sigmaz = float(cfg.sigmaz)
    except ValueError:
        D = msh.z[msh.Npml] - msh.z[0]
        cfg.sigmaz = 8*cfg.c0*(cfg.alpha+1)/D


class Fields:
    """ Fields initialization. """

    def __init__(self, msh, cfg):


        self._msh = msh
        self._cfg = cfg

        # Fields
        self.init_fields()
        self.init_derivatives()
        self.init_filters()

        self.fdtools = fdtools(self._cfg.nx, self._cfg.nz,
                               self._cfg.p0, self._cfg.gamma, self._cfg.dt)

        # PMLS
        if 'A' in self._msh.bc:
            self.init_pml()

        # Obstacles to 0
        self.zero_obstacles()

        # Save file
        if self._cfg.save:
            self.init_save()

    def init_fields(self):
        """ Setup initial fields. """

        # Fields
        self.p = np.zeros(self._msh.shape) + self._cfg.p0
        self.r = np.zeros_like(self.p) + self._cfg.rho0
        self.ru = np.zeros_like(self.p) + self._cfg.U0
        self.rv = np.zeros_like(self.p) + self._cfg.V0
        self.dltn = np.zeros_like(self.p)

        # Source
        self.Bx = self._cfg.B0*self._cfg.dx
        self.src = np.empty_like(self.p)
        self.update_source = None
        self.source_select()

        # Flow
        self.flow_select()

        # Curvilinear
        if self._cfg.mesh == 'curvilinear':
            self.r = self.r/self._msh.J
            self.ru = self.ru/self._msh.J
            self.rv = self.rv/self._msh.J
            self.p = self.p/self._msh.J

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

        set_sigma(self._cfg, self._msh)

        # sigmax
        Dx = self._msh.x[self._msh.Npml] - self._msh.x[0]             # Width of the PML
        self.sx = np.zeros(self._msh.nx)
        self.sx[:self._msh.Npml] = self._cfg.sigmax*abs((self._msh.x[:self._msh.Npml] \
                                 - self._msh.x[self._msh.Npml])/Dx)**self._cfg.alpha
        self.sx[self._msh.nx - self._msh.Npml:] = self.sx[self._msh.Npml-1::-1]

        # sigmaz
        Dz = self._msh.z[self._msh.Npml] - self._msh.z[0]             # Width of the PML
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
        self.Ei = np.zeros_like(self.p)
        self.Eui = np.zeros_like(self.p)
        self.Evi = np.zeros_like(self.p)
        self.Eei = np.zeros_like(self.p)
        self.Fi = np.zeros_like(self.p)
        self.Fui = np.zeros_like(self.p)
        self.Fvi = np.zeros_like(self.p)
        self.Fei = np.zeros_like(self.p)

        if self._cfg.mesh in ['regular', 'adaptative']:
            self.fdtools.Eu(self.Ei, self.Eui, self.Evi, self.Eei,
                            self.r, self.ru, self.rv, self.re, self.p)

        elif self._cfg.mesh == 'curvilinear':
            self.fdtools.EuJ(self.Ei, self.Eui, self.Evi, self.Eei,
                             self.Fi, self.Fui, self.Fvi, self.Fei,
                             self.r, self.ru, self.rv, self.re, self.p,
                             self._msh.dxn_dxp, self._msh.dxn_dzp)

        if self._cfg.mesh in ['regular', 'adaptative']:
            self.fdtools.Fu(self.Fi, self.Fui, self.Fvi, self.Fei,
                            self.r, self.ru, self.rv, self.re, self.p)

        elif self._cfg.mesh == 'curvilinear':
            self.fdtools.FuJ(self.Ei, self.Eui, self.Evi, self.Eei,
                             self.Fi, self.Fui, self.Fvi, self.Fei,
                             self.r, self.ru, self.rv, self.re, self.p,
                             self._msh.dzn_dxp, self._msh.dzn_dzp)

    def init_save(self):
        """ Init save. """

        self.sfile = h5py.File(self._cfg.datafile, 'w')

        self.sfile.create_dataset('x', data=self._msh.x, compression=self._cfg.comp)
        self.sfile.create_dataset('z', data=self._msh.z, compression=self._cfg.comp)
        self.sfile.create_dataset('rho_init', data=self.r, compression=self._cfg.comp)
        self.sfile.create_dataset('rhou_init', data=self.ru, compression=self._cfg.comp)
        self.sfile.create_dataset('rhov_init', data=self.rv, compression=self._cfg.comp)
        self.sfile.create_dataset('rhoe_init', data=self.re, compression=self._cfg.comp)

        self.sfile.attrs['obstacles'] = self._msh.get_obstacles()
        self.sfile.attrs['dx'] = self._msh.dx
        self.sfile.attrs['dz'] = self._msh.dz
        self.sfile.attrs['dt'] = self._cfg.dt
        self.sfile.attrs['nx'] = self._cfg.nx
        self.sfile.attrs['nz'] = self._cfg.nz
        self.sfile.attrs['nt'] = self._cfg.nt
        self.sfile.attrs['ns'] = self._cfg.ns
        self.sfile.attrs['p0'] = self._cfg.p0
        self.sfile.attrs['rho0'] = self._cfg.rho0
        self.sfile.attrs['gamma'] = self._cfg.gamma
        self.sfile.attrs['Npml'] = self._cfg.Npml
        self.sfile.attrs['mesh'] = self._cfg.mesh
        self.sfile.attrs['bc'] = self._cfg.bc

        probes = np.zeros((len(self._cfg.probes), self._cfg.nt))
        self.sfile.create_dataset('probes_location', data=self._cfg.probes)
        self.sfile.create_dataset('probes_value', data=probes,
                                  compression=self._cfg.comp)

        if self._cfg.mesh == 'curvilinear':
            self.sfile.create_dataset('J', data=self._msh.J, compression=self._cfg.comp)
            self.sfile.create_dataset('xn', data=self._msh.xn, compression=self._cfg.comp)
            self.sfile.create_dataset('zn', data=self._msh.zn, compression=self._cfg.comp)
            self.sfile.create_dataset('xp', data=self._msh.xp, compression=self._cfg.comp)
            self.sfile.create_dataset('zp', data=self._msh.zp, compression=self._cfg.comp)

    def source_select(self):
        """ Source selection. """

        if self._cfg.stype in ["", "None", "none"]:
            pass

        elif self._cfg.stype == "pulse":
            self.pulse()

        elif self._cfg.stype == "harmonic":
            self.update_source = self.update_harmonic
            self.harmonic()

        elif self._cfg.stype == "white":
            self.update_source = self.update_white
            self.white()

        else:
            raise ValueError('Only pulse, harmonic and white supported for now')

    def flow_select(self):
        """ Flow selection. """

        if self._cfg.ftype in ["", "None", "none"]:
            pass
        elif self._cfg.ftype == "custom":
            pass
        elif self._cfg.ftype == "vortex":
            self.isentropic_vortex()
        elif self._cfg.ftype == "poiseuille":
            self.poiseuille()
        elif self._cfg.ftype == "kh":
            self.kelvin_helmholtz()
        else:
            raise ValueError('Only custom and vortex supported for now')

    def pulse(self):
        """ Pulse source. """

        # Location
        ixS = self._cfg.ixS
        izS = self._cfg.izS

        for iz, z in enumerate(self._msh.z):
            self.p[:, iz] += self._cfg.S0*np.exp(-np.log(2)*((self._msh.x-self._msh.x[ixS])**2 +
                                                             (z - self._msh.z[izS])**2)/self.Bx**2)

    def harmonic(self):
        """ Harmonic ponctual source. """

        for iz, z in enumerate(self._msh.z):
            self.src[:, iz] = np.exp(-np.log(2)*((self._msh.x - self._msh.x[self._cfg.ixS])**2 +
                                                 (z - self._msh.z[self._cfg.izS])**2)/self.Bx**2)
        self.src = self._cfg.S0*self.src

    def white(self):
        """ white noise. """

        for iz, z in enumerate(self._msh.z):
            tmp = (1 - np.log(2)*((self._msh.x - self._msh.x[self._cfg.ixS])**2 +
                                  (z - self._msh.z[self._cfg.izS])**2)/self.Bx**2)
            self.src[:, iz] = tmp*np.exp(-np.log(2)*((self._msh.x - self._msh.x[self._cfg.ixS])**2 +
                                                     (z - self._msh.z[self._cfg.izS])**2)/self.Bx**2)
        self.src = self._cfg.S0*self.src
        self.noise = np.random.normal(size=self._cfg.nt)

    def update_harmonic(self, it):
        """ Harmonic source time evolution. """

        return self.src*np.sin(2*np.pi*self._cfg.f0*it*self._cfg.dt)

    def update_wall(self, it, f0=None):
        """ Harmonic source time evolution for walls. """

        if not f0:
            f0 = self._cfg.f0

        return np.sin(2*np.pi*f0*it*self._cfg.dt)

    def update_white(self, it):
        """ White noise time evolution. """

        return self.src*self.noise[it-1]

    def isentropic_vortex(self):
        """ Initialize fields with an isentropic vortex [Hu, JCP, 2008] """

        Umax = 0.5*self._cfg.U0
        b = 0.2

        theta = np.zeros(self._cfg.nx)

        # initialisation des fluctuations
        for iz in range(self._cfg.nz):

            r = np.sqrt(self._msh.x**2 + self._msh.z[iz]**2)

            for ix in range(self._cfg.nx):
                if self._msh.x[ix] == 0 and self._msh.z[iz] == 0:
                    theta[ix] = 0
                elif self._msh.x[ix] >= 0:
                    theta[ix] = np.arcsin(self._msh.z[iz]/r[ix])
                elif self._msh.x[ix] < 0:
                    theta[ix] = - np.arcsin(self._msh.z[iz]/r[ix]) + np.pi

            self.r[:, iz] = (1 - 0.5*(self._cfg.gamma - 1)*Umax**2 \
                             * np.exp(1-r**2/b**2))**(1/(self._cfg.gamma - 1))
            self.ru[:, iz] -= Umax*r*np.exp(0.5*(1-r**2/b**2))/b*np.sin(theta)
            self.rv[:, iz] += Umax*r*np.exp(0.5*(1-r**2/b**2))/b*np.cos(theta)

        self.p = self.r**self._cfg.gamma/self._cfg.gamma
        self.ru *= self.r
        self.rv *= self.r

    def poiseuille(self):
        """ Poiseuille """

        # Poiseuille
        tmp_z = np.arange(self._msh.nz)*self._msh.dz
        self._cfg.U0 *= 0.5*(4*tmp_z/(tmp_z[-1]-tmp_z[0]) - 4*tmp_z**2/(tmp_z[-1]-tmp_z[0])**2)

    def kelvin_helmholtz(self):
        """
            Initialize fields with a mixing layer with roll-up vortices induced
            by the Kelvin–Helmholtz instability.
        """

        U1 = 0.8
        U2 = 0.2
        delta = 0.4
        T1 = 1.
        T2 = 0.8

        T = np.zeros((self._cfg.nx, self._cfg.nz))

        for iz in range(self._cfg.nz):
            self.ru[:, iz] = 0.5*((U1+U2) +
                                  (U1-U2)*np.tanh(2*self._msh.z[iz]/delta))
            T[:, iz] = T1*(self.ru[:, iz] - U2)/(U1 - U2) \
                     + T2*(U1-self.ru[:, iz])/(U1 - U2) \
                     + 0.5*(self._cfg.gamma-1)*(U1 - self.ru[:, iz])*(self.ru[:, iz]-U2)
            self.r[:, iz] = 1/T[:, iz]

        self.p = np.ones((self._msh.nx, self._msh.nz))/self._cfg.gamma
        self.ru = self.ru*self.r

    def zero_obstacles(self):
        """ Set velocity to 0 in obstacles. """

        for sub in self._msh.obstacles:
            self.ru[sub.sx, sub.sz] = 0
            self.rv[sub.sx, sub.sz] = 0

    def get(self):
        """ Get initial fields as a tuple. """
        return self.p, self.r, self.ru, self.rv, self.re
