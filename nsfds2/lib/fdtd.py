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
# Creation Date : 2019-03-01 - 14:43:51
#
# pylint: disable=too-many-instance-attributes
"""
-----------

Finite Difference Time Domain class

@author: Cyril Desjouy
"""

import os
import time
import numpy as np
from ofdlib2 import fdtd
from nsfds2.utils import headers
from .utils import disp_bench, timed, Check
from .efluxes import EulerianFluxes
from .vfluxes import ViscousFluxes
from .sfilter import SelectiveFilter
from .scapture import ShockCapture


class FDTD:
    """ FDTD. """

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg

        # Headers
        headers.copyright()
        headers.version()

        # Check some of the simulation parameters
        Check(self.cfg, self.msh)

        # Init libs
        self.efluxes = EulerianFluxes(self.msh, self.fld, self.cfg)
        self.vfluxes = ViscousFluxes(self.msh, self.fld, self.cfg)
        self.sfilter = SelectiveFilter(self.msh, self.fld, self.cfg)
        self.scapture = ShockCapture(self.msh, self.fld, self.cfg)

        # Init probes
        if cfg.probes and cfg.probes_loc:
            self.probes = np.zeros((len(cfg.probes_loc), cfg.ns))
        else:
            cfg.probes = False

        # Prompt user for start
        headers.start(cfg)

        # Init timings
        self.it = 0
        timing = ['total', 'efluxes', 'vfluxes', 'sfilt',
                  'scapt', 'save', 'pressure', 'probe']
        self.bench = {i: [] for i in timing}
        self.tloopi = time.perf_counter()

    @property
    def columns(self):
        """ Return max terminal width. """
        _, col = os.popen('stty size', 'r').read().split()
        return int(col) if int(col) < 81 else 80

    def run(self):
        """ Main loop. """

        print('-'*int(self.columns))
        print('# Start main loop ')

        for self.it in range(self.cfg.nt+1):

            tt = time.perf_counter()

            # FDTD
            self.eulerian_fluxes()
            self.viscous_flux()
            self.selective_filter()
            self.shock_capture()
            self.update_pressure()

            # Break when computation diverges
            res = fdtd.residual(self.fld.p, self.cfg.p0)
            if (abs(res) > 100*self.cfg.S0) or np.any(np.isnan(self.fld.p)):
                print('Stop simulation at iteration ', self.it)
                if np.any(np.isnan(self.fld.p)):
                    print('Nan : {}'.format(np.argwhere(np.isnan(self.fld.p))))
                if self.cfg.save:
                    self.fld.sfile.close()
                break

            # Display log
            if self.it % self.cfg.ns == 0:
                self.save()
                self.bench['total'].append(time.perf_counter() - tt)
                self.bench = disp_bench(self.bench, self.it, res)
            else:
                self.bench['total'].append(time.perf_counter() - tt)

            if self.cfg.probes:
                self.update_probes()

        if self.cfg.save:
            self.fld.sfile.close()

        print('-'*int(self.columns))
        msg = '# Simulation completed in {:.2f} s.\n'
        msg += '# End at t = {:.4f} sec.'
        print(msg.format(time.perf_counter() - self.tloopi, self.cfg.dt*self.it))
        print('-'*int(self.columns))

    @timed('efluxes')
    def eulerian_fluxes(self):
        """ Compute Eulerian fluxes. """

        self.efluxes.rk4()

    @timed('vfluxes')
    def viscous_flux(self):
        """ Viscous flux """
        if self.cfg.vsc:
            self.vfluxes.integrate()

    @timed('sfilt')
    def selective_filter(self):
        """ Selective filter """
        if self.cfg.flt:
            self.sfilter.apply()

    @timed('scapt')
    def shock_capture(self):
        """ Shock Capture """
        if self.cfg.cpt:
            self.scapture.apply()

    @timed('save')
    def save(self):
        """ Save data """

        if self.cfg.save and self.cfg.probes:
            self.fld.sfile['probes'][:, self.it-self.cfg.ns:self.it] = self.probes

        if self.cfg.save and self.cfg.onlyp:
            self.fld.sfile.create_dataset('p_it' + str(self.it),
                                          data=self.fld.p, compression=self.cfg.comp)

        elif self.cfg.save:
            self.fld.sfile.create_dataset('rho_it' + str(self.it),
                                          data=self.fld.r, compression=self.cfg.comp)
            self.fld.sfile.create_dataset('rhou_it' + str(self.it),
                                          data=self.fld.ru, compression=self.cfg.comp)
            self.fld.sfile.create_dataset('rhov_it' + str(self.it),
                                          data=self.fld.rv, compression=self.cfg.comp)
            self.fld.sfile.create_dataset('rhoe_it' + str(self.it),
                                          data=self.fld.re, compression=self.cfg.comp)

    @timed('pressure')
    def update_pressure(self):
        """ Update pressure field """

        fdtd.p(self.fld.p, self.fld.r, self.fld.ru,
               self.fld.rv, self.fld.re, self.cfg.gamma)

    @timed('probe')
    def update_probes(self):
        """ Update probes """
        for n, c in enumerate(self.cfg.probes_loc):
            self.probes[n, self.it%self.cfg.ns] = self.fld.p[c[0], c[1]]

    def __str__(self):

        rep = '-'*int(self.columns)
        rep += '\n'
        rep += 'Mesh shape    : {}\n'.format(self.msh.shape)
        rep += 'Spacial steps : dx={}, dz={}\n'.format(self.msh.dx, self.msh.dz)
        rep += 'CFL           : {}\n'.format(self.cfg.CFL)
        rep += 'Time step     : {}\n'.format(self.cfg.dt)
        rep += '-'*int(self.columns)
        return rep

    def __repr__(self):
        return self.__str__()
