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
# Creation Date : 2016-11-29 - 23:18:27
"""
-----------

Parse config file and set all simulation parameters

@author: Cyril Desjouy
"""


import sys
import os
import time
try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser


class CfgSetup:
    """ Handle configuration file. """

    # pylint: disable=too-many-instance-attributes

    def __init__(self):

        ti = time.perf_counter()
        self.cfg = ConfigParser.RawConfigParser()
        self.home = os.path.expanduser("~")
        self.path = self.home + '/.nsfds2/'

        # Check if config dir exists. If not create it.
        self.check_dir(self.path)

        # Check if config file exist. If not create it
        self.init_cfg()

        # Read config file
        self.cfg.read(self.path + 'nsfds2.conf')

        # Parse arguments
        self.run()

        # Log
        msg = 'Config parsed in {:.4f} s.'
        print(msg.format(time.perf_counter() - ti))

    @staticmethod
    def check_dir(directory):
        """ Check if dir exists. If not, create it."""

        if not os.path.isdir(directory):
            os.makedirs(directory)
            print("Create directory :", directory)
            time.sleep(0.5)

    def init_cfg(self):
        """ Check if nsfds2.conf exists. If not create it. """

        if not os.path.exists(self.path + 'nsfds2.conf'):
            open(self.path + 'nsfds2.conf', 'a').close()
            print("Create configuration file : {}nsfds2.conf".format(self.path))
            time.sleep(0.5)
            self.write_default()

    def write_default(self):
        """ Write default configuration file. """

        self.cfg.add_section('simulation')
        self.cfg.set('simulation', 'viscosity', 'True')
        self.cfg.set('simulation', 'probes', 'True')
        self.cfg.set('simulation', 'nt', '150')
        self.cfg.set('simulation', 'ns', '10')
        self.cfg.set('simulation', 'nx', '256')
        self.cfg.set('simulation', 'nz', '256')
        self.cfg.set('simulation', 'ix0', '0')
        self.cfg.set('simulation', 'iz0', '0')
        self.cfg.set('simulation', 'dx', '1')
        self.cfg.set('simulation', 'dz', '1')
        self.cfg.set('simulation', 'CFL', '0.5')
        self.cfg.set('simulation', 'bc', 'RRRR')
        self.cfg.set('simulation', 'Npml', '15')
        self.cfg.set('simulation', 'Stencil', '11')
        self.cfg.set('simulation', 'viscosity', 'True')

        self.cfg.add_section('filtering')
        self.cfg.set('filtering', 'filter', 'True')

        self.cfg.add_section('shock capture')
        self.cfg.set('shock capture', 'shock capture', 'True')
        self.cfg.set('shock capture', 'method', 'pressure')

        self.cfg.add_section('source')
        self.cfg.set('source', 'type', 'pulse')
        self.cfg.set('source', 'ixS', '64')
        self.cfg.set('source', 'izS', '64')
        self.cfg.set('source', 'S0', '1e3')

        self.cfg.add_section('probes')
        self.cfg.set('probes', 'probes', 'False')

        self.cfg.add_section('save')
        self.cfg.set('save', 'save', 'True')
        self.cfg.set('save', 'filename', 'tmp')
        self.cfg.set('save', 'compression', 'lzf')
        self.cfg.set('save', 'only p', 'False')
        self.cfg.set('save', 'view', 'p')


        with open(self.path + 'nsfds2.conf', 'w') as configfile:
            self.cfg.write(configfile)

    def run(self):
        """ Run configuration. """

        # Thermophysic parameters
        self.gamma = 1.4
        self.nu = 1.5e-5
        self.c0 = 340.
        self.rho0 = 1.22
        self.p0 = self.rho0*self.c0**2/self.gamma

        try:
            SIM = self.cfg['simulation']
            self.viscosity = SIM.getboolean('viscosity', True)
            self.nt = SIM.getint('nt', 150)
            self.ns = SIM.getint('ns', 10)
            self.nx = SIM.getint('nx', 256)
            self.nz = SIM.getint('nz', 256)
            self.ix0 = SIM.getint('ix0', 0)
            self.iz0 = SIM.getint('iz0', 0)
            self.dx = SIM.getfloat('dx', 1)
            self.dz = SIM.getfloat('dz', 1)
            self.CFL = SIM.getfloat('CFL', 0.5)
            self.bc = SIM.get('bc', 'RRRR')
            self.stencil = SIM.getint('stencil', 11)
            self.Npml = SIM.getint('Npml', 15)
            self.dt = min(self.dx, self.dz)*self.CFL/self.c0

            FILT = self.cfg['filtering']
            self.filt = FILT.getboolean('filter', True)

            SCAPT = self.cfg['shock capture']
            self.scapt = SCAPT.getboolean('shock capture', True)
            self.scapt_meth = SCAPT.get('method', 'pressure')
            self.rth = 1e-6
            self.eps_machine = 6e-8                 # epsilon single precision

            SRC = self.cfg['source']
            self.typ = SRC.get('type', 'pulse')
            self.ixS = SRC.getint('ixS', 64)
            self.izS = SRC.getint('izS', 64)
            self.S0 = SRC.getfloat('S0', 1e3)

            SAVE = self.cfg['save']
            self.save = SAVE.getboolean('save', True)
            self.filename = SAVE.get('filename', 'tmp')
            self.comp = SAVE.get('compression', 'lzf')
            self.view = SAVE.get('view', 'p')
            self.onlyp = SAVE.getboolean('only p', False)
            if self.comp == 'None':
                self.comp = None

            PRBS = self.cfg['probes']
            self.probes = PRBS.getboolean('probes', False)


        except ConfigParser.Error as err:
            print('Bad cfg file : ', err)
            sys.exit(1)
