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

To initialize all the fields, :py:class:`Fields` needs :

    * the mesh (from :py:class:fdgrid.mesh.Mesh)
    * the configuration of the simulation as an instance of :py:class:`nsfds2.init.config.CfgSetup`


Example
-------

::

    fld = Fields(msh, cfg)

-----------
"""

import os as _os
import sys as _sys
import h5py
import numpy as _np
from ofdlib2.fdtd import fdtools
from nsfds2.utils import headers, misc
from nsfds2.utils.sounds import resample


def set_sigma(cfg, msh):
    """ Set sigmax and sigmaz for PMLs. """

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
    """ Initialize fields.

    Parameters
    ----------
    msh : :py:class:fdgrid.mesh.Mesh
        The mesh
    cfg : :py:class:`nsfds2.init.config.CfgSetup`
        The simulation configuration

    """

    def __init__(self, msh, cfg):

        self._msh = msh
        self._cfg = cfg

        self._nx = self._cfg.nx
        self._nz = self._cfg.nz
        self._x = self._msh.x
        self._z = self._msh.z

        # Headers
        if not self._cfg.quiet:
            headers.copyright()
            headers.version()

        # Fields
        self._init_fields()
        self._init_derivatives()
        self._init_filters()
        self._init_wav_wall()

        self.fdtools = fdtools(self._nx, self._nz,
                               self._cfg.p0, self._cfg.gamma, self._cfg.dt)

        # PMLS
        if 'A' in self._msh.bc:
            self._init_pml()
            self._init_radiation()

        # Obstacles to 0
        self.zero_obstacles()

    def _init_fields(self):
        """ Setup initial fields. """

        if self._cfg.mesh == 'curvilinear':
            self.p = _np.zeros(self._msh.shape) + self._cfg.p0/self._msh.J
            self.r = _np.zeros_like(self.p) + self._cfg.rho0/self._msh.J
            self.ru = _np.zeros_like(self.p) + self._cfg.U0/self._msh.J
            self.rv = _np.zeros_like(self.p) + self._cfg.V0/self._msh.J
        else:
            self.p = _np.zeros(self._msh.shape) + self._cfg.p0
            self.r = _np.zeros_like(self.p) + self._cfg.rho0
            self.ru = _np.zeros_like(self.p) + self._cfg.U0
            self.rv = _np.zeros_like(self.p) + self._cfg.V0

        if self._cfg.cpt_meth == 'dilatation':
            self.dltn = _np.zeros_like(self.p)

        # Source
        self.source_select()

        # Flow
        self.flow_select()

        # Update rho.e
        self.re = self.p/(self._cfg.gamma-1.)

    def _init_derivatives(self):
        """ Init derivatives. """

        self.K = _np.zeros_like(self.p)
        self.Ku = _np.zeros_like(self.p)
        self.Kv = _np.zeros_like(self.p)
        self.Ke = _np.zeros_like(self.p)

        self.E = _np.empty_like(self.p)
        self.Eu = _np.empty_like(self.p)
        self.Ev = _np.empty_like(self.p)
        self.Ee = _np.empty_like(self.p)

        self.F = _np.empty_like(self.p)
        self.Fu = _np.empty_like(self.p)
        self.Fv = _np.empty_like(self.p)
        self.Fe = _np.empty_like(self.p)

    def _init_filters(self):
        """ Init variables used for filtering. """

        self.dp = _np.zeros_like(self.p)
        self.sg = _np.zeros_like(self.p)
        self.tau11 = _np.zeros_like(self.p)
        self.tau22 = _np.zeros_like(self.p)
        self.tau12 = _np.zeros_like(self.p)

    def _init_pml(self):
        """ Init PMLs. """

        set_sigma(self._cfg, self._msh)

        # sigmax
        Dx = self._x[self._msh.Npml] - self._x[0]             # Width of the PML
        self.sx = _np.zeros(self._nx)
        self.sx[:self._msh.Npml] = self._cfg.sigmax*abs((self._x[:self._msh.Npml] \
                                 - self._x[self._msh.Npml])/Dx)**self._cfg.alpha
        self.sx[self._nx - self._msh.Npml:] = self.sx[self._msh.Npml-1::-1]

        # sigmaz
        Dz = self._z[self._msh.Npml] - self._z[0]             # Width of the PML
        self.sz = _np.zeros(self._nz)
        self.sz[:self._msh.Npml] = self._cfg.sigmaz*abs((self._z[:self._msh.Npml] \
                                 - self._z[self._msh.Npml])/Dz)**self._cfg.alpha
        self.sz[self._nz - self._msh.Npml:] = self.sz[self._msh.Npml-1::-1]

        # Init Q
        self.qx = _np.zeros_like(self.p)
        self.qux = _np.zeros_like(self.p)
        self.qvx = _np.zeros_like(self.p)
        self.qex = _np.zeros_like(self.p)
        self.qz = _np.zeros_like(self.p)
        self.quz = _np.zeros_like(self.p)
        self.qvz = _np.zeros_like(self.p)
        self.qez = _np.zeros_like(self.p)

        # Init K
        self.Kx = _np.zeros_like(self.p)
        self.Kux = _np.zeros_like(self.p)
        self.Kvx = _np.zeros_like(self.p)
        self.Kex = _np.zeros_like(self.p)
        self.Kz = _np.zeros_like(self.p)
        self.Kuz = _np.zeros_like(self.p)
        self.Kvz = _np.zeros_like(self.p)
        self.Kez = _np.zeros_like(self.p)

        # Initial E & F
        self.Ei = _np.zeros_like(self.p)
        self.Eui = _np.zeros_like(self.p)
        self.Evi = _np.zeros_like(self.p)
        self.Eei = _np.zeros_like(self.p)
        self.Fi = _np.zeros_like(self.p)
        self.Fui = _np.zeros_like(self.p)
        self.Fvi = _np.zeros_like(self.p)
        self.Fei = _np.zeros_like(self.p)

        if self._cfg.mesh in ['regular', 'adaptative']:
            self.fdtools.Eu(self.Ei, self.Eui, self.Evi, self.Eei,
                            self.r, self.ru, self.rv,
                            self._cfg.p0*_np.ones_like(self.p)/(self._cfg.gamma-1.),
                            self._cfg.p0*_np.ones_like(self.p))

        elif self._cfg.mesh == 'curvilinear':
            self.fdtools.EuJ(self.Ei, self.Eui, self.Evi, self.Eei,
                             self.r, self.ru, self.rv, self.re, self.p,
                             self._msh.dxn_dxp, self._msh.dxn_dzp)

        if self._cfg.mesh in ['regular', 'adaptative']:
            self.fdtools.Fu(self.Fi, self.Fui, self.Fvi, self.Fei,
                            self.r, self.ru, self.rv,
                            self._cfg.p0*_np.ones_like(self.p)/(self._cfg.gamma-1.),
                            self._cfg.p0*_np.ones_like(self.p))

        elif self._cfg.mesh == 'curvilinear':
            self.fdtools.FuJ(self.Fi, self.Fui, self.Fvi, self.Fei,
                             self.r, self.ru, self.rv, self.re, self.p,
                             self._msh.dzn_dxp, self._msh.dzn_dzp)

    def _init_radiation(self):
        """ Init radiation conditions. """

        # Init K
        self.Kx = _np.zeros_like(self.p)
        self.Kux = _np.zeros_like(self.p)
        self.Kvx = _np.zeros_like(self.p)
        self.Kex = _np.zeros_like(self.p)
        self.Kz = _np.zeros_like(self.p)
        self.Kuz = _np.zeros_like(self.p)
        self.Kvz = _np.zeros_like(self.p)
        self.Kez = _np.zeros_like(self.p)

        xCL = self._x[self._cfg.ixS]
        zCL = self._z[self._cfg.izS]

        self.ixCL = _np.zeros(11)
        self.izCL = _np.zeros(11)

        self.ixCL[-5] = self._nx - 5
        self.ixCL[-4] = self._nx - 4
        self.ixCL[-3] = self._nx - 3
        self.ixCL[-2] = self._nx - 2
        self.ixCL[-1] = self._nx - 1
        self.ixCL[0] = 1
        self.ixCL[1] = 2
        self.ixCL[2] = 3
        self.ixCL[3] = 4
        self.ixCL[4] = 5

        self.izCL[-5] = self._nz - 5
        self.izCL[-4] = self._nz - 4
        self.izCL[-3] = self._nz - 3
        self.izCL[-2] = self._nz - 2
        self.izCL[-1] = self._nz - 1
        self.izCL[0] = 0
        self.izCL[1] = 1
        self.izCL[2] = 2
        self.izCL[3] = 3
        self.izCL[4] = 4

        self.r_v = _np.zeros((10, self._nz))
        self.cosv = _np.zeros((10, self._nz))
        self.sinv = _np.zeros((10, self._nz))
        self.r_h = _np.zeros((self._nx, 10))
        self.cosh = _np.zeros((self._nx, 10))
        self.sinh = _np.zeros((self._nx, 10))


        # prendre (xCL,yCL) pour origine pour le calcul de one_r, cosphi et sinphi
        # calcul pour les bandes verticales
        for ix in range(-5, 5):
            self.r_v[ix, :] = 0.5/_np.sqrt((self._x[ix] - xCL)**2 + (self._z - zCL)**2)
            self.cosv[ix, :] = (self._x[ix] - xCL)*2*self.r_v[ix, :]
            self.sinv[ix, :] = (self._z - zCL)*2*self.r_v[ix, :]

        # calcul pour les bandes horizontales
        for iz in range(-5, 5):
            self.r_h[:, iz] = 0.5/_np.sqrt((self._x - xCL)**2 + (self._z[iz] - zCL)**2)
            self.cosh[:, iz] = (self._x - xCL)*2.*self.r_h[:, iz]
            self.sinh[:, iz] = (self._z[iz] - zCL)*2.*self.r_h[:, iz]

    def _init_wav_wall(self):
        """ Initialize wall with a time evolution built from a wavfile. """

        for obs in self._msh.obstacles:
            for bc in obs.edges:
                if isinstance(bc.f0_n, str):
                    filename = bc.f0_n.replace('~', self._cfg.home)
                    if filename in self.wav_lst:
                        bc.wav = self.wav_lst[filename]
                    else:
                        bc.wav = resample(filename, 1/self._cfg.dt, pad=self._cfg.nt+1)
                        self.wav_lst[filename] = bc.wav

    def init_save(self):
        """ Init save. """

        if _os.path.isfile(self._cfg.datafile):
            msg = f'{self._cfg.datafile} already exists. Overwrite ? [yes]/no '
            overwrite = input(misc.colors.RED + msg + misc.colors.END)
            if overwrite.lower() in ['n', 'no']:
                _sys.exit(1)

        self.sfile = h5py.File(self._cfg.datafile, 'w')

        self.sfile.create_dataset('x', data=self._x, compression=self._cfg.comp)
        self.sfile.create_dataset('z', data=self._z, compression=self._cfg.comp)
        self.sfile.create_dataset('rho_init', data=self.r, compression=self._cfg.comp)
        self.sfile.create_dataset('rhou_init', data=self.ru, compression=self._cfg.comp)
        self.sfile.create_dataset('rhov_init', data=self.rv, compression=self._cfg.comp)
        self.sfile.create_dataset('rhoe_init', data=self.re, compression=self._cfg.comp)

        self.sfile.attrs['obstacles'] = self._msh.get_obstacles()
        self.sfile.attrs['dx'] = self._msh.dx
        self.sfile.attrs['dz'] = self._msh.dz
        self.sfile.attrs['dt'] = self._cfg.dt
        self.sfile.attrs['nx'] = self._nx
        self.sfile.attrs['nz'] = self._nz
        self.sfile.attrs['nt'] = self._cfg.nt
        self.sfile.attrs['ns'] = self._cfg.ns
        self.sfile.attrs['p0'] = self._cfg.p0
        self.sfile.attrs['rho0'] = self._cfg.rho0
        self.sfile.attrs['gamma'] = self._cfg.gamma
        self.sfile.attrs['Npml'] = self._cfg.Npml
        self.sfile.attrs['mesh'] = self._cfg.mesh
        self.sfile.attrs['bc'] = self._cfg.bc

        probes = _np.zeros((len(self._cfg.probes), self._cfg.nt))
        self.sfile.create_dataset('probes_location', data=self._cfg.probes)
        self.sfile.create_dataset('probes_value', data=probes,
                                  compression=self._cfg.comp)

        if self._cfg.mesh == 'curvilinear':
            self.sfile.create_dataset('J', data=self._msh.J, compression=self._cfg.comp)
            self.sfile.create_dataset('xn', data=self._xn, compression=self._cfg.comp)
            self.sfile.create_dataset('zn', data=self._zn, compression=self._cfg.comp)
            self.sfile.create_dataset('xp', data=self._xp, compression=self._cfg.comp)
            self.sfile.create_dataset('zp', data=self._zp, compression=self._cfg.comp)

    def source_select(self):
        """ Source selection. """

        self.Bx = self._cfg.B0*self._cfg.dx
        self.src = _np.empty_like(self.p)
        self.update_source = None
        self.wav_lst = dict()

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

        elif self._cfg.stype == "wav" and self._cfg.wavfile:
            self.wav = resample(self._cfg.wavfile, 1/self._cfg.dt, pad=self._cfg.nt+1)
            self.wav_lst[self._cfg.wavfile] = self.wav
            self.update_source = self.update_wav
            self.harmonic()

        else:
            raise ValueError('Only pulse, harmonic, white, wave supported for now')

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

        if self._cfg.mesh == 'curvilinear':
            for iz in range(self._nz):
                for ix in range(self._nx):
                    self.p[ix, iz] = self._cfg.p0/self._msh.J[ix, iz] \
                             + self._cfg.S0*_np.exp(-_np.log(2) \
                             * ((self._msh.xp[ix, iz] - self._msh.xp[ixS, izS])**2 +
                                (self._msh.zp[ix, iz] - self._msh.zp[ixS, izS])**2)/self.Bx**2)
        else:
            for iz, z in enumerate(self._z):
                self.p[:, iz] += self._cfg.S0*_np.exp(-_np.log(2) \
                            * ((self._x-self._x[ixS])**2 +
                               (z - self._z[izS])**2)/self.Bx**2)

    def harmonic(self):
        """ Harmonic ponctual source. """

        for iz, z in enumerate(self._z):
            self.src[:, iz] = _np.exp(-_np.log(2)*((self._x - self._x[self._cfg.ixS])**2 +
                                                   (z - self._z[self._cfg.izS])**2)/self.Bx**2)
        self.src = self._cfg.S0*self.src

    def white(self):
        """ white noise. """

        for iz, z in enumerate(self._z):
            tmp = (1 - _np.log(2)*((self._x - self._x[self._cfg.ixS])**2 +
                                   (z - self._z[self._cfg.izS])**2)/self.Bx**2)
            self.src[:, iz] = tmp*_np.exp(-_np.log(2)*((self._x - self._x[self._cfg.ixS])**2 +
                                                       (z - self._z[self._cfg.izS])**2)/self.Bx**2)
        self.src = self._cfg.S0*self.src
        _np.random.seed(self._cfg.seed)
        self.noise = _np.random.normal(size=self._cfg.nt)

    def update_wav(self, it):
        """ Wav source time evolution. """

        return self.src*self.wav[it]

    def update_harmonic(self, it):
        """ Harmonic source time evolution. """

        if it <= self._cfg.off:
            return self.src*_np.sin(2*_np.pi*self._cfg.f0*it*self._cfg.dt)

        return 0

    def update_wall(self, it, f=None, phi=0):
        """ Harmonic source time evolution for walls. """

        if not f:
            f = self._cfg.f0

        return _np.sin(2*_np.pi*f*it*self._cfg.dt + phi)

    def update_white(self, it):
        """ White noise time evolution. """

        if it <= self._cfg.off:
            return self.src*self.noise[it-1]

        return 0

    def isentropic_vortex(self):
        """ Initialize fields with an isentropic vortex [Hu, JCP, 2008] """

        Umax = 0.5*self._cfg.U0
        b = 0.2

        theta = _np.zeros(self._nx)

        # initialisation des fluctuations
        for iz in range(self._nz):

            r = _np.sqrt(self._x**2 + self._z[iz]**2)

            for ix in range(self._nx):
                if self._x[ix] == 0 and self._z[iz] == 0:
                    theta[ix] = 0
                elif self._x[ix] >= 0:
                    theta[ix] = _np.arcsin(self._z[iz]/r[ix])
                elif self._x[ix] < 0:
                    theta[ix] = - _np.arcsin(self._z[iz]/r[ix]) + _np.pi

            self.r[:, iz] = (1 - 0.5*(self._cfg.gamma - 1)*Umax**2 \
                             * _np.exp(1-r**2/b**2))**(1/(self._cfg.gamma - 1))
            self.ru[:, iz] -= Umax*r*_np.exp(0.5*(1-r**2/b**2))/b*_np.sin(theta)
            self.rv[:, iz] += Umax*r*_np.exp(0.5*(1-r**2/b**2))/b*_np.cos(theta)

        self.p = self.r**self._cfg.gamma/self._cfg.gamma
        self.ru *= self.r
        self.rv *= self.r

    def poiseuille(self):
        """ Poiseuille """

        # Poiseuille
        tmp_z = _np.arange(self._nz)*self._msh.dz
        self._cfg.U0 *= (4*tmp_z/(tmp_z[-1]-tmp_z[0]) -
                         4*tmp_z**2/(tmp_z[-1]-tmp_z[0])**2)

        for ix in range(self._nx):
            self.ru[ix, :] = self._cfg.U0

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

        T = _np.zeros((self._nx, self._nz))

        for iz in range(self._nz):
            self.ru[:, iz] = 0.5*((U1+U2) +
                                  (U1-U2)*_np.tanh(2*self._z[iz]/delta))
            T[:, iz] = T1*(self.ru[:, iz] - U2)/(U1 - U2) \
                     + T2*(U1-self.ru[:, iz])/(U1 - U2) \
                     + 0.5*(self._cfg.gamma-1)*(U1 - self.ru[:, iz])*(self.ru[:, iz]-U2)
            self.r[:, iz] = 1/T[:, iz]

        self.p = _np.ones((self._nx, self._nz))/self._cfg.gamma
        self.ru = self.ru*self.r

    def zero_obstacles(self):
        """ Set velocity to 0 in obstacles. """

        for sub in self._msh.obstacles:
            self.ru[sub.sx, sub.sz] = 0
            self.rv[sub.sx, sub.sz] = 0

    def get(self):
        """ Get initial fields as a tuple. """
        return self.p, self.r, self.ru, self.rv, self.re
