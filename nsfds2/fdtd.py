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
import sys
import time
import numpy as np
from progressbar import ProgressBar, Bar, ReverseBar, ETA
from nsfds2.utils import headers, check, timing
from nsfds2.lib.efluxes import EulerianFluxes
from nsfds2.lib.vfluxes import ViscousFluxes
from nsfds2.lib.sfilter import SelectiveFilter
from nsfds2.lib.scapture import ShockCapture


class FDTD:
    """ FDTD. """

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg

        # Headers
        if not self.cfg.quiet:
            headers.copyright()
            headers.version()

        # Check some of the simulation parameters
        check.Check(self.cfg, self.msh)

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
        if not self.cfg.quiet:
            headers.start(cfg)

        # Init timings and residuals
        self.res = 0
        self.it = 0
        self.tt = 0
        timed_methods = ['total', 'efluxes', 'vfluxes', 'sfilt',
                         'scapt', 'save', 'pressure', 'probe']
        self.bench = {i: [] for i in timed_methods}
        self.tloopi = time.perf_counter()

    def run(self):
        """ Main loop. """

        if not self.cfg.quiet:
            print('-'*int(self.columns))
            print('# Start main loop ')

        if not self.cfg.quiet and not self.cfg.timings:
            widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
            pbar = ProgressBar(widgets=widgets, maxval=self.cfg.nt,
                               term_width=self.columns).start()

        for self.it in range(self.cfg.nt+1):

            self.tt = time.perf_counter()

            # FDTD
            self.eulerian_fluxes()
            self.viscous_flux()

            if self.cfg.mesh == 'curvilinear':
                self.num2phys()

            self.selective_filter()
            self.shock_capture()

            if self.cfg.mesh == 'curvilinear':
                self.phys2num()

            self.update_pressure()

            # Break when computation diverges
            self.check_results()

            # Display log
            self.log()

            # Update probes
            if self.cfg.probes:
                self.update_probes()

            # Progress bar
            if not self.cfg.timings and not self.cfg.quiet:
                pbar.update(self.it)


        if self.cfg.save:
            self.fld.sfile.close()

        if not self.cfg.quiet and not self.cfg.timings:
            pbar.finish()

        if not self.cfg.quiet:
            print('-'*int(self.columns))
            msg = '# Simulation completed in {:.2f} s.\n'
            msg += '# End at t = {:.4f} sec.'
            print(msg.format(time.perf_counter() - self.tloopi, self.cfg.dt*self.it))
            print('-'*int(self.columns))

    @timing.proceed('efluxes')
    def eulerian_fluxes(self):
        """ Compute Eulerian fluxes. """

        self.efluxes.rk4()

    @timing.proceed('vfluxes')
    def viscous_flux(self):
        """ Viscous flux """
        if self.cfg.vsc:
            self.vfluxes.integrate()

    @timing.proceed('sfilt')
    def selective_filter(self):
        """ Selective filter """
        if self.cfg.flt:
            self.sfilter.apply()

    @timing.proceed('scapt')
    def shock_capture(self):
        """ Shock Capture """
        if self.cfg.cpt:
            self.scapture.apply()

    @timing.proceed('save')
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

    @timing.proceed('pressure')
    def update_pressure(self):
        """ Update pressure field """

        self.fld.fdtools.p(self.fld.p, self.fld.r, self.fld.ru,
                           self.fld.rv, self.fld.re)

    @timing.proceed('probe')
    def update_probes(self):
        """ Update probes. """

        for n, c in enumerate(self.cfg.probes_loc):
            self.probes[n, self.it%self.cfg.ns] = self.fld.p[c[0], c[1]]

    def check_results(self):
        """ Check if computation diverges. """

        if self.cfg.mesh == 'curvilinear':
            self.res = self.fld.fdtools.residual(self.fld.p*self.msh.J)
        else:
            self.res = self.fld.fdtools.residual(self.fld.p)

        if (abs(self.res) > 100*self.cfg.S0) or np.any(np.isnan(self.fld.p)):
            print('Stop simulation at iteration ', self.it)
            if np.any(np.isnan(self.fld.p)):
                print('Nan : {}'.format(np.argwhere(np.isnan(self.fld.p))))
            if self.cfg.save:
                self.fld.sfile.close()
            sys.exit(1)

    def phys2num(self):
        """ Convert curvilinear coordinates : from physical to numeric. """

        self.fld.r = self.fld.r/self.msh.J
        self.fld.ru = self.fld.ru/self.msh.J
        self.fld.rv = self.fld.rv/self.msh.J
        self.fld.re = self.fld.re/self.msh.J

    def num2phys(self):
        """ Convert curvilinear coordinates : from numeric to physical. """

        self.fld.r = self.fld.r*self.msh.J
        self.fld.ru = self.fld.ru*self.msh.J
        self.fld.rv = self.fld.rv*self.msh.J
        self.fld.re = self.fld.re*self.msh.J

    def log(self):
        """ Display timings. """

        if self.it % self.cfg.ns == 0:
            self.save()
            if self.cfg.timings and not self.cfg.quiet:
                self.bench['total'].append(time.perf_counter() - self.tt)
                self.bench = timing.disp(self.bench, self.it, self.res)
        elif self.cfg.timings and not self.cfg.quiet:
            self.bench['total'].append(time.perf_counter() - self.tt)

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

    @property
    def columns(self):
        """ Return max terminal width. """
        try:
            _, col = os.popen('stty size', 'r').read().split()
            return int(col) if int(col) < 81 else 80
        except ValueError:
            return 80
