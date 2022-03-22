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

The :py:class:`FDTD` class is the heart of `nsfds2`. It coordinates all the
calculations that have to be led to perform the Finite Differences Time Domain
(FDTD) simulation. The :py:class:`FDTD` works with :

    * a configuration file (parsed by :py:class:`nsfds2.init.config.CfgSetup`)
    * initial fields (created by :py:class:`nsfds2.init.fields.Fields`)
    * a mesh (created with fdgrid)


Example
-------

::

    from nsfds2.init import CfgSetup, Fields
    from nsfds2.fdtd import FDTD
    from fdgrid import mesh

    # Read simulation parameter in config file ~/nsfds2/nsfds2.conf
    cfg = CfgSetup()

    # Define the mesh
    msh = mesh.build('regular', (cfg.nx, cfg.nz), (cfg.dx, cfg.dz),
                     origin=(cfg.ix0, cfg.iz0), obstacles=[])

    # Create simulation
    fdtd = FDTD(msh, cfg)
    fdtd.run()



-----------
"""

import os
import sys
import time
import numpy as np
from progressbar import ProgressBar, Bar, ReverseBar, ETA
from nsfds2.utils import headers, check, timing, misc
from nsfds2.init.fields import Fields
from nsfds2.lib.efluxes import EulerianFluxes
from nsfds2.lib.vfluxes import ViscousFluxes
from nsfds2.lib.sfilter import SelectiveFilter
from nsfds2.lib.scapture import ShockCapture
from nsfds2.lib.vorticity import Vorticity


class FDTD:
    """

    The Finite Differences Time Domain Solver

    Parameters
    ----------
        msh : :py:class:`nsfds2.fdgrid.mesh.Mesh`
            Mesh
        cfg : :py:class:`nsfds2.init.config.CfgSetup`
            Configuration

    """

    def __init__(self, msh, cfg):

        self.msh = msh
        self.cfg = cfg
        self.fld = Fields(msh, cfg)

        # Check some of the simulation parameters
        check.Check(self.cfg, self.msh)

        # Init timings and residuals
        self.res = 0
        self.tt = 0
        timed_methods = ['total', 'efluxes', 'save', 'pressure']

        # Init libs
        self.efluxes = EulerianFluxes(self.msh, self.fld, self.cfg)
        if self.cfg.vsc:
            self.vfluxes = ViscousFluxes(self.msh, self.fld, self.cfg)
            timed_methods += ['vfluxes']
        if self.cfg.flt:
            self.sfilter = SelectiveFilter(self.msh, self.fld, self.cfg)
            timed_methods += ['sfilt']
        if self.cfg.cpt:
            self.scapture = ShockCapture(self.msh, self.fld, self.cfg)
            timed_methods += ['scapt']
        if self.cfg.save_vortis:
            self.vorticity = Vorticity(self.msh, self.fld, self.cfg)
            timed_methods += ['vorticity']

        # Init probes
        if self.cfg.probes:
            self.probes = np.zeros((len(cfg.probes), cfg.ns))
            timed_methods += ['probes']

        # Prompt user for start
        if not self.cfg.quiet:
            headers.start(cfg)

        # Max pressure
        v = (self.fld.ru.max()/self.cfg.rho0)**2 + (self.fld.rv.max()/self.cfg.rho0)**2
        self.pmax = - 0.5*self.cfg.rho0*v \
                    *(self.cfg.gamma - 2)/(self.cfg.gamma - 1) + self.cfg.S0

        # Start timings
        self.bench = {i: [] for i in timed_methods}
        self.tloopi = time.perf_counter()

    def _pre_run(self):

        # Save file
        if self.cfg.resume:
            self.fld.init_resume()
        else:
            self.fld.init_save()

        if not self.cfg.quiet:
            print('-'*int(self.columns))
            print('# Start main loop ')

    def run(self):
        """ Run the main loop. """

        self._pre_run()

        if not self.cfg.quiet and not self.cfg.timings:
            widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
            pbar = ProgressBar(widgets=widgets, maxval=self.cfg.nt,
                               term_width=self.columns).start()

        for self.cfg.it in range(self.cfg.it, self.cfg.nt+1):

            self.tt = time.perf_counter()

            # FDTD
            self.eulerian_fluxes()
            self.viscous_flux()

            self.num2phys()
            self.selective_filter()
            self.shock_capture()
            self.phys2num()

            self.update_pressure()
            self.update_vorticity()

            # Break when computation diverges
            self.check_results()

            # Display log
            self.log()

            # Update probes
            self.update_probes()

            # Progress bar
            if not self.cfg.timings and not self.cfg.quiet:
                pbar.update(self.cfg.it)

        self.fld.sfile.close()

        if not self.cfg.quiet and not self.cfg.timings:
            pbar.finish()

        if not self.cfg.quiet:
            print('-'*int(self.columns))
            msg = '# Simulation completed in {}.\n'
            msg += '# End at t = {:.4f} sec.'
            print(msg.format(misc.secs_to_dhms(time.perf_counter() - self.tloopi),
                             self.cfg.dt*self.cfg.it))
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
            self.efluxes.cout()

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

        self.fld.sfile.attrs['itmax'] = self.cfg.it

        if self.cfg.save_fields:
            self.fld.sfile.create_dataset('r_it' + str(self.cfg.it),
                                          data=self.fld.r,
                                          compression=self.cfg.comp)
            self.fld.sfile.create_dataset('ru_it' + str(self.cfg.it),
                                          data=self.fld.ru,
                                          compression=self.cfg.comp)
            self.fld.sfile.create_dataset('rv_it' + str(self.cfg.it),
                                          data=self.fld.rv,
                                          compression=self.cfg.comp)
            self.fld.sfile.create_dataset('re_it' + str(self.cfg.it),
                                          data=self.fld.re,
                                          compression=self.cfg.comp)

        if self.cfg.save_vortis:
            self.fld.sfile.create_dataset('vxz_it' + str(self.cfg.it),
                                          data=self.fld.vxz,
                                          compression=self.cfg.comp)

        if self.cfg.probes:
            self.fld.sfile['probe_values'][:, self.cfg.it-self.cfg.ns:self.cfg.it] = self.probes

    @timing.proceed('pressure')
    def update_pressure(self):
        """ Update pressure field """

        self.fld.fdt.p(self.fld.p, self.fld.r, self.fld.ru,
                       self.fld.rv, self.fld.re)

    @timing.proceed('vorticity')
    def update_vorticity(self):
        """ Update vorticity field """

        if self.cfg.save_vortis:
            self.vorticity.compute()

    @timing.proceed('probes')
    def update_probes(self):
        """ Update probes. """

        if self.cfg.probes:
            for n, c in enumerate(self.cfg.probes):
                self.probes[n, self.cfg.it % self.cfg.ns] = \
                        self.fld.p[c[0], c[1]]

    def check_results(self):
        """ Check if computation diverges. """

        if self.cfg.mesh == 'curvilinear':
            self.res = self.fld.fdt.residual(self.fld.p*self.msh.J)
        else:
            self.res = self.fld.fdt.residual(self.fld.p)

#        if (abs(self.res) > 100*self.pmax) or np.any(np.isnan(self.fld.p)):
        if np.any(np.isnan(self.fld.p)):

            print('Stop simulation at iteration ', self.cfg.it)
            if np.any(np.isnan(self.fld.p)):
                print(f'Nan : {np.argwhere(np.isnan(self.fld.p))}')

            self.fld.sfile.close()
            sys.exit(1)

    def phys2num(self):
        """ Convert curvilinear coordinates : from physical to numeric. """

        if self.cfg.mesh == 'curvilinear':
            self.fld.r = self.fld.r/self.msh.J
            self.fld.ru = self.fld.ru/self.msh.J
            self.fld.rv = self.fld.rv/self.msh.J
            self.fld.re = self.fld.re/self.msh.J

    def num2phys(self):
        """ Convert curvilinear coordinates : from numeric to physical. """

        if self.cfg.mesh == 'curvilinear':
            self.fld.r = self.fld.r*self.msh.J
            self.fld.ru = self.fld.ru*self.msh.J
            self.fld.rv = self.fld.rv*self.msh.J
            self.fld.re = self.fld.re*self.msh.J

    def log(self):
        """ Display timings. """

        if self.cfg.it % self.cfg.ns == 0:
            self.save()
            if self.cfg.timings and not self.cfg.quiet:
                self.bench['total'].append(time.perf_counter() - self.tt)
                self.bench = timing.disp(self.bench, self.cfg.it, self.res)
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
