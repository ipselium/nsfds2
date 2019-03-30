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
try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser


class CfgSetup:
    """ Handle configuration file. """

    # pylint: disable=too-many-instance-attributes

    def __init__(self):

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
        self.cfg.set('simulation', 'nt', '150')
        self.cfg.set('simulation', 'ns', '10')
        self.cfg.set('simulation', 'CFL', '0.5')
        self.cfg.set('simulation', 'Npml', '15')
        self.cfg.set('simulation', 'Stencil', '11')

        self.cfg.add_section('geometry')
        self.cfg.set('geometry', 'file', 'None')
        self.cfg.set('geometry', 'name', 'square')
        self.cfg.set('geometry', 'bc', 'RRRR')
        self.cfg.set('geometry', 'nx', '256')
        self.cfg.set('geometry', 'nz', '256')
        self.cfg.set('geometry', 'ix0', '0')
        self.cfg.set('geometry', 'iz0', '0')
        self.cfg.set('geometry', 'dx', '1')
        self.cfg.set('geometry', 'dz', '1')

        self.cfg.add_section('source')
        self.cfg.set('source', 'type', 'pulse')
        self.cfg.set('source', 'ixS', '32')
        self.cfg.set('source', 'izS', '32')
        self.cfg.set('source', 'S0', '1e3')

        self.cfg.add_section('filtering')
        self.cfg.set('filtering', 'filter', 'True')
        self.cfg.set('filtering', 'stength', '0.75')

        self.cfg.add_section('viscous fluxes')
        self.cfg.set('viscous fluxes', 'viscosity', 'True')

        self.cfg.add_section('shock capture')
        self.cfg.set('shock capture', 'shock capture', 'True')
        self.cfg.set('shock capture', 'method', 'pressure')

        self.cfg.add_section('probes')
        self.cfg.set('probes', 'probes', 'False')

        self.cfg.add_section('figures')
        self.cfg.set('figures', 'figures', 'True')

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
            self.nt = SIM.getint('nt', 150)
            self.ns = SIM.getint('ns', 10)
            self.CFL = SIM.getfloat('CFL', 0.5)
            self.stencil = SIM.getint('stencil', 11)
            self.Npml = SIM.getint('Npml', 15)

            GEO = self.cfg['geometry']
            self.geofile = GEO.get('file', 'None')
            self.geoname = GEO.get('name', 'square')
            self.bc = GEO.get('bc', 'RRRR')
            self.nx = GEO.getint('nx', 256)
            self.nz = GEO.getint('nz', 256)
            self.ix0 = GEO.getint('ix0', 0)
            self.iz0 = GEO.getint('iz0', 0)
            self.dx = GEO.getfloat('dx', 1)
            self.dz = GEO.getfloat('dz', 1)
            self.dt = min(self.dx, self.dz)*self.CFL/self.c0

            SRC = self.cfg['source']
            self.typ = SRC.get('type', 'pulse')
            self.ixS = SRC.getint('ixS', 32)
            self.izS = SRC.getint('izS', 32)
            self.S0 = SRC.getfloat('S0', 1e3)

            FILT = self.cfg['filtering']
            self.filt = FILT.getboolean('filter', True)
            self.xnu = FILT.getfloat('strength', 0.75)

            VISC = self.cfg['viscous fluxes']
            self.vsc = VISC.getboolean('viscosity', True)

            SCAPT = self.cfg['shock capture']
            self.scapt = SCAPT.getboolean('shock capture', True)
            self.scapt_meth = SCAPT.get('method', 'pressure')
            self.rth = 1e-6

            SAVE = self.cfg['save']
            self.save = SAVE.getboolean('save', True)
            self.filename = SAVE.get('filename', 'tmp')
            self.comp = SAVE.get('compression', 'lzf')
            self.view = SAVE.get('view', 'p')
            self.onlyp = SAVE.getboolean('only p', False)
            if self.comp == 'None':
                self.comp = None

            FIGS = self.cfg['figures']
            self.figures = FIGS.getboolean('figures', True)

            PRBS = self.cfg['probes']
            self.probes = PRBS.getboolean('probes', False)

        except ConfigParser.Error as err:
            print('Bad cfg file : ', err)
            sys.exit(1)
