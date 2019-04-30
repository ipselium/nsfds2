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
#
# pylint: disable=too-many-statements
"""
-----------

Parse config file and set all simulation parameters

@author: Cyril Desjouy
"""

import json
import time
import sys
import os
import configparser


class CfgSetup:
    """ Handle configuration file. """

    # pylint: disable=too-many-instance-attributes

    def __init__(self, args=None):

        # Command line arguments
        self.args = args

        # Create config parser
        self.cfg = configparser.RawConfigParser()
        self.home = os.path.expanduser("~")
        self.path = self.home + '/.nsfds2/'

        # Check if config dir exists. If not create it.
        self.check_dir(self.path)

        # Check if config file exist. If not create it
        self.init_cfg()

        # Read config file (can be overridden by command line)
        cf = getattr(self.args, 'cfgfile', None)
        if cf:
            self.cfg.read(cf)
        else:
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

        self.cfg.add_section('configuration')
        self.cfg.set('configuration', 'timings', 'True')
        self.cfg.set('configuration', 'quiet', 'False')

        self.cfg.add_section('simulation')
        self.cfg.set('simulation', 'nt', '500')
        self.cfg.set('simulation', 'ns', '10')
        self.cfg.set('simulation', 'CFL', '0.5')

        self.cfg.add_section('geometry')
        self.cfg.set('geometry', 'mesh', 'regular')
        self.cfg.set('geometry', 'file', 'None')
        self.cfg.set('geometry', 'geoname', 'helmholtz_double')
        self.cfg.set('geometry', 'curvname', 'None')
        self.cfg.set('geometry', 'bc', 'PPPP')
        self.cfg.set('geometry', 'nx', '256')
        self.cfg.set('geometry', 'nz', '256')
        self.cfg.set('geometry', 'ix0', '0')
        self.cfg.set('geometry', 'iz0', '0')
        self.cfg.set('geometry', 'dx', '1')
        self.cfg.set('geometry', 'dz', '1')

        self.cfg.add_section('PML')
        self.cfg.set('PML', 'beta', 0.)
        self.cfg.set('PML', 'alpha', 4.)
        self.cfg.set('PML', 'sigmax', 20.)
        self.cfg.set('PML', 'sigmaz', 20.)
        self.cfg.set('PML', 'Npml', 15)

        self.cfg.add_section('source')
        self.cfg.set('source', 'type', 'pulse')
        self.cfg.set('source', 'ixS', '32')
        self.cfg.set('source', 'izS', '128')
        self.cfg.set('source', 'S0', '1e6')

        self.cfg.add_section('eulerian fluxes')
        self.cfg.set('eulerian fluxes', 'stencil', '11')

        self.cfg.add_section('filtering')
        self.cfg.set('filtering', 'filter', 'True')
        self.cfg.set('filtering', 'stencil', '11')
        self.cfg.set('filtering', 'stength', '0.75')

        self.cfg.add_section('viscous fluxes')
        self.cfg.set('viscous fluxes', 'viscosity', 'True')
        self.cfg.set('viscous fluxes', 'stencil', '7')

        self.cfg.add_section('shock capture')
        self.cfg.set('shock capture', 'shock capture', 'True')
        self.cfg.set('shock capture', 'stencil', '7')
        self.cfg.set('shock capture', 'method', 'pressure')

        self.cfg.add_section('figures')
        self.cfg.set('figures', 'figures', 'True')

        self.cfg.add_section('save')
        self.cfg.set('save', 'save', 'True')
        self.cfg.set('save', 'path', 'results/')
        self.cfg.set('save', 'filename', 'tmp')
        self.cfg.set('save', 'compression', 'lzf')
        self.cfg.set('save', 'only p', 'False')
        self.cfg.set('save', 'probes', 'False')
        self.cfg.set('save', 'probes_locations', '[]')

        with open(self.path + 'nsfds2.conf', 'w') as cf:
            self.cfg.write(cf)

    def run(self):
        """ Run configuration. """

        # Thermophysic parameters
        self.gamma = 1.4
        self.nu = 1.5e-5
        self.c0 = 340.
        self.rho0 = 1.22
        self.p0 = self.rho0*self.c0**2/self.gamma

        try:
            CFG = self.cfg['configuration']

            self.timings = getattr(self.args, 'timings', None)
            if not isinstance(self.timings, bool):
                self.timings = CFG.getboolean('timings', True)

            self.quiet = getattr(self.args, 'quiet', None)
            if not isinstance(self.quiet, bool):
                self.quiet = CFG.getboolean('quiet', False)

            SIM = self.cfg['simulation']
            self.nt = getattr(self.args, 'nt', None)
            if self.nt == None:
                self.nt = SIM.getint('nt', 500)
            self.ns = SIM.getint('ns', 10)
            self.CFL = SIM.getfloat('CFL', 0.5)

            GEO = self.cfg['geometry']
            self.mesh = GEO.get('mesh', 'regular')
            self.geofile = getattr(self.args, 'geofile', None)
            if  self.geofile:
                self.geofile, self.geoname = self.args.geofile
            else:
                self.geofile = GEO.get('file', 'None')
                self.geoname = GEO.get('geoname', 'square')
            self.curvname = GEO.get('curvname', 'None')
            self.bc = GEO.get('bc', 'RRRR')
            self.nx = GEO.getint('nx', 256)
            self.nz = GEO.getint('nz', 256)
            self.ix0 = GEO.getint('ix0', 0)
            self.iz0 = GEO.getint('iz0', 0)
            self.dx = GEO.getfloat('dx', 1)
            self.dz = GEO.getfloat('dz', 1)
            self.dt = min(self.dx, self.dz)*self.CFL/self.c0

            PML = self.cfg['PML']
            self.beta = PML.getfloat('beta', 0.)
            self.alpha = PML.getfloat('alpha', 4.)
            self.sigmax = PML.getfloat('sigmax', 20.)
            self.sigmaz = PML.getfloat('sigmaz', 20.)
            self.Npml = PML.getint('Npml', 15)

            SRC = self.cfg['source']
            self.typ = SRC.get('type', 'pulse')
            self.ixS = SRC.getint('ixS', 32)
            self.izS = SRC.getint('izS', 32)
            self.S0 = SRC.getfloat('S0', 1e3)

            EUL = self.cfg['eulerian fluxes']
            self.stencil = EUL.getint('stencil', 11)

            FLT = self.cfg['filtering']
            self.flt = FLT.getboolean('filter', True)
            self.flt_stencil = FLT.getint('stencil', 11)
            self.xnu = FLT.getfloat('strength', 0.75)

            VSC = self.cfg['viscous fluxes']
            self.vsc = VSC.getboolean('viscosity', True)
            self.vsc_stencil = VSC.getint('stencil', 3)

            CPT = self.cfg['shock capture']
            self.cpt = CPT.getboolean('shock capture', True)
            self.cpt_stencil = CPT.getint('stencil', 7)
            self.cpt_meth = CPT.get('method', 'pressure')
            self.rth = 1e-6

            SAVE = self.cfg['save']
            self.save = SAVE.getboolean('save', True)
            self.savepath = SAVE.get('path', 'results/')
            self.savefile = SAVE.get('filename', 'tmp') + '.hdf5'
            self.comp = SAVE.get('compression', 'lzf')
            self.onlyp = SAVE.getboolean('only p', False)
            self.probes = SAVE.getboolean('probes', False)
            self.probes_loc = json.loads(SAVE.get('probes_locations', '[]'))

            # datapath and datafile
            if self.comp == 'None':
                self.comp = None

            if self.savepath and not self.savepath.endswith('/'):
                self.savepath += '/'

            # if self.savepath does not exist, create it
            self.check_dir(self.savepath)

            self.datafile = getattr(self.args, 'datafile', None)
            if not self.datafile:
                self.datafile = self.savepath + self.savefile

            FIGS = self.cfg['figures']
            self.figures = FIGS.getboolean('figures', True)

        except configparser.Error as err:
            print('Bad cfg file : ', err)
            sys.exit(1)
