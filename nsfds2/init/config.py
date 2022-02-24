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
# pylint: disable=attribute-defined-outside-init
# pylint: disable=too-many-instance-attributes
"""
-----------

This module contains the :py:class:`CfgSetup` that read the configuration file
and set all simulation parameters.

Example
-------

::

    from nsfds2.init import CfgSetup

    cfg = CfgSetup()

-----------
"""

import json
import time
import sys
import os
import shutil
import datetime
import pathlib
import configparser
from pkg_resources import parse_version
import nsfds2
from nsfds2.utils import files


class CfgSetup:
    """ Handle configuration file. """

    def __init__(self, args=None):

        # Minimal version of the config file
        self.base_version = '0.9.15'

        # Command line arguments + home
        self.args = args
        self.home = pathlib.Path.home()
        self.path_default = self.home / '.nsfds2'
        self.cfgfile_default = self.path_default / 'nsfds2.conf'

        # Create config parser
        self.cfg = configparser.ConfigParser(allow_no_value=True)

        # Load config file
        if isinstance(self.args, str):
            self.cfgfile = pathlib.Path(self.args)
        else:
            self.cfgfile = getattr(self.args, 'cfgfile', None)

        # Check cfg file
        if not self.cfgfile:
            self.path = self.path_default
            self.cfgfile = self.cfgfile_default
            self.check_dir(self.path) # Check if cfg dir exists. If not create it.
            self.init_cfg()           # Check if cfg file exist. If not create it
        else:
            self.cfgfile = pathlib.Path(self.cfgfile)
            self.path = self.cfgfile.absolute().parent

        # Check if config file version is ok
        self.check_config_file()

        # Read config file (can be overridden by command line)
        self.cfg.read(self.cfgfile)

        # Parse arguments
        self.run()

    @staticmethod
    def check_dir(directory):
        """ Check if dir exists. If not, create it."""

        if not directory.is_dir():
            directory.mkdir()
            print("Create directory :", directory)
            time.sleep(0.5)

    def check_config_file(self):
        """ Check version of the config file. Overwrite it if too old. """

        cfg = configparser.ConfigParser(allow_no_value=True)
        cfg.read(self.cfgfile)

        try:
            CFG = cfg['configuration']
            version = CFG.get('version')
            is_default_cfg = self.cfgfile == self.cfgfile_default
            is_version_ok = parse_version(version) >= parse_version(self.base_version)

            if not is_version_ok and is_default_cfg:
                self._overwrite_config_file()
            elif not is_version_ok:
                print(f'Config file version must be >= {self.base_version}')
                sys.exit(1)

        except (KeyError, TypeError):
            print(f'Config file version must be >= {self.base_version}')
            if is_default_cfg:
                self._overwrite_config_file()
            else:
                sys.exit(1)

    def _overwrite_config_file(self):

        # Backup old config file
        now = datetime.datetime.now()
        name = f'{now.year}{now.month}{now.day}{now.hour}{now.minute}{now.second}'
        shutil.move(self.path / 'nsfds2.conf', self.path / f'nsfds2_{name}.conf')

        print(f'Current configfile backup : nsfds2_{name}.conf')
        time.sleep(1)

        # Create new config file
        self.init_cfg()

    def init_cfg(self):
        """ Check if config file exists. If not create it. """

        if not (self.path / 'nsfds2.conf').is_file():
            open(self.path / 'nsfds2.conf', 'a').close()
            print("Create configuration file : {}/nsfds2.conf".format(self.path))
            time.sleep(0.5)
            self.write_default()

    def write_default(self):
        """ Write default configuration file. """

        self.cfg.add_section('configuration')
        self.cfg.set('configuration', 'version', str(nsfds2.__version__))
        self.cfg.set('configuration', 'timings', 'False')
        self.cfg.set('configuration', 'quiet', 'False')
        self.cfg.set('configuration', 'cpu', '1')

        self.cfg.add_section('simulation')
        self.cfg.set('simulation', 'nt', '500')
        self.cfg.set('simulation', 'ns', '10')
        self.cfg.set('simulation', 'CFL', '0.5')

        self.cfg.add_section('thermophysic')
        self.cfg.set('thermophysic', 'P0', '101325.0')
        self.cfg.set('thermophysic', 'T0', '20.0')
        self.cfg.set('thermophysic', 'gamma', '1.4')
        self.cfg.set('thermophysic', 'prandtl', '0.7')

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
        self.cfg.set('PML', 'beta', '0.')
        self.cfg.set('PML', 'alpha', '4.')
        self.cfg.set('PML', 'sigmax', 'auto')
        self.cfg.set('PML', 'sigmaz', 'auto')
        self.cfg.set('PML', 'Npml', '15')

        self.cfg.add_section('source')
        self.cfg.set('source', 'type', 'pulse')
        self.cfg.set('source', 'ixS', '32')
        self.cfg.set('source', 'izS', '128')
        self.cfg.set('source', 'S0', '1e6')
        self.cfg.set('source', 'B0', '5')
        self.cfg.set('source', 'f0', '20000')
        self.cfg.set('source', 'wavfile', 'None')

        self.cfg.add_section('flow')
        self.cfg.set('flow', 'type', 'None')
        self.cfg.set('flow', 'U0', '5')
        self.cfg.set('flow', 'V0', '5')

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
        self.cfg.set('figures', 'probes', 'True')
        self.cfg.set('figures', 'pml', 'True')
        self.cfg.set('figures', 'bc_profiles', 'True')
        self.cfg.set('figures', 'fps', '24')

        self.cfg.add_section('save')
        self.cfg.set('save', 'path', 'results/')
        self.cfg.set('save', 'filename', 'tmp')
        self.cfg.set('save', 'compression', 'lzf')
        self.cfg.set('save', 'fields', 'True')
        self.cfg.set('save', 'probes', '[]')

        with open(self.path / 'nsfds2.conf', 'w') as cf:
            self.cfg.write(cf)

    def run(self):
        """ Run configuration. """

        self.none = ['none', 'false', '']

        try:
            self._cfg()
            self._sim()
            self._thp()
            self._geo()
            self._pml()
            self._src()
            self._flw()
            self._eul()
            self._flt()
            self._vsc()
            self._cpt()
            self._save()
            self._figs()

        except configparser.Error as err:
            print('Bad cfg file : ', err)
            sys.exit(1)

        self.dt = min(self.dx, self.dz)*self.CFL/(self.c0 +
                                                  max(abs(self.U0), abs(self.V0)))

    def _cfg(self):

        CFG = self.cfg['configuration']
        self.timings = getattr(self.args, 'timings', None)
        self.quiet = getattr(self.args, 'quiet', None)
        self.cpu = CFG.getint('cpu', 1)

        if not isinstance(self.timings, bool):
            self.timings = CFG.getboolean('timings', False)

        if not isinstance(self.quiet, bool):
            self.quiet = CFG.getboolean('quiet', False)

    def _sim(self):

        SIM = self.cfg['simulation']
        self.nt = getattr(self.args, 'nt', None)
        self.ns = SIM.getint('ns', 10)
        self.CFL = SIM.getfloat('CFL', 0.5)

        if self.nt is None:
            self.nt = SIM.getint('nt', 500)

        if self.nt % self.ns:
            self.nt -= self.nt % self.ns

        self.it = 0

    def _thp(self):

        THP = self.cfg['thermophysic']

        self.Ssu = 111.0  # Sutherland constant
        self.T0 = 273.0
        self.T = self.T0 + THP.getfloat('T0', 20.0)

        self.gamma = THP.getfloat('gamma', 1.4)
        self.cv = 717.5
        self.cp = self.cv*self.gamma

        self.p0 = THP.getfloat('P0', 101325.0)
        self.rho0 = self.p0/(self.T*(self.cp - self.cv))
        self.c0 = (self.gamma*self.p0/self.rho0)**0.5
        self.mu0 = 0.00001716
        self.mu = (self.mu0*(self.T/self.T0)**(3./2.) *
                   (self.T0 + self.Ssu)/(self.T + self.Ssu))
        self.nu = self.mu/self.rho0
        self.prandtl = THP.getfloat('prandtl', 0.7)

#        self.nu = THP.getfloat('nu', 1.5e-5)
#        self.p0 = self.rho0*self.c0**2/self.gamma

        if self.c0 < 1:
            raise ValueError('c0 must be >= 1')

    def _geo(self):

        GEO = self.cfg['geometry']
        self.mesh = GEO.get('mesh', 'regular').lower()
        self.curvflag = True if self.mesh == 'curvilinear' else False
        self.curvname = GEO.get('curvname', 'None')
        self.geofile = GEO.get('file', 'None')
        self.geoname = GEO.get('geoname', 'square')
        self.geoflag = True
        self.bc = GEO.get('bc', 'WWWW').upper()
        self.nx = GEO.getint('nx', 256)
        self.nz = GEO.getint('nz', 256)
        self.ix0 = GEO.getint('ix0', 0)
        self.iz0 = GEO.getint('iz0', 0)
        self.dx = GEO.getfloat('dx', 1)
        self.dz = GEO.getfloat('dz', 1)

        if self.geofile != "None":
            self.geofile = self.path / self.geofile

        if self.mesh not in ['regular', 'adaptative', 'curvilinear']:
            raise ValueError('mesh must be regular, adaptative, or curvilinear')

        self.obstacles = files.get_obstacle(self)

    def _pml(self):

        PML = self.cfg['PML']
        self.beta = PML.getfloat('beta', 0.)
        self.alpha = PML.getfloat('alpha', 4.)
        self.sigmax = PML.get('sigmax', 'auto')
        self.sigmaz = PML.get('sigmaz', 'auto')
        self.Npml = PML.getint('Npml', 15)

    def _src(self):

        SRC = self.cfg['source']
        self.stype = SRC.get('type', 'pulse').lower()
        self.ixS = SRC.getint('ixS', 32)
        self.izS = SRC.getint('izS', 32)
        self.S0 = SRC.getfloat('S0', 1e3)
        self.B0 = SRC.getfloat('B0', 5)
        self.f0 = SRC.getfloat('f0', 20000)
        self.wavfile = SRC.get('wavfile', None)
        self.seed = SRC.get('seed', None)
        self.off = SRC.getint('off', self.nt)

        if self.wavfile:
            self.wavfile = pathlib.Path(self.wavfile).expanduser()

        if self.seed:
            try:
                self.seed = int(self.seed)
            except ValueError:
                raise ValueError('Seed must be int or None')

        if self.stype in self.none:
            self.S0 = 0

    def _flw(self):

        FLW = self.cfg['flow']
        self.ftype = FLW.get('type', 'None').lower()
        self.U0 = FLW.getfloat('U0', 5)
        self.V0 = FLW.getfloat('V0', 5)

        if self.ftype in self.none:
            self.U0, self.V0 = 0., 0.

    def _eul(self):

        EUL = self.cfg['eulerian fluxes']
        self.stencil = EUL.getint('stencil', 11)

        if self.stencil not in [3, 7, 11]:
            raise ValueError('stencil must be 3, 7 or 11')

    def _flt(self):

        FLT = self.cfg['filtering']
        self.flt = FLT.getboolean('filter', True)
        self.flt_stencil = FLT.getint('stencil', 11)
        self.xnu = FLT.getfloat('strength', 0.75)

        if self.flt_stencil not in [7, 11]:
            raise ValueError('only 7 and 11 pts filters implemented for now')

    def _vsc(self):

        VSC = self.cfg['viscous fluxes']
        self.vsc = VSC.getboolean('viscosity', True)
        self.vsc_stencil = VSC.getint('stencil', 3)

        if self.vsc_stencil not in [3, 7, 11]:
            raise ValueError('viscous fluxes only available with 3, 7 or 11 pts')

    def _cpt(self):

        CPT = self.cfg['shock capture']
        self.cpt = CPT.getboolean('shock capture', True)
        self.cpt_stencil = CPT.getint('stencil', 7)
        self.cpt_meth = CPT.get('method', 'pressure').lower()
        self.rth = 1e-6

        if self.cpt_stencil not in [3, 7, 11]:
            raise ValueError('capture only available with 3, 7 or 11 pts')

        if self.cpt_meth not in ['pressure', 'dilatation']:
            raise ValueError('capture method must be pressure or dilatation')

    def _save(self):

        SAVE = self.cfg['save']
        self.save = SAVE.getboolean('fields', True)
        if self.path == self.path_default:
            self.savepath = pathlib.Path(SAVE.get('path', 'results'))
        else:
            self.savepath = self.path / SAVE.get('path', 'results')
        self.savefile = SAVE.get('filename', 'tmp') + '.hdf5'
        self.comp = SAVE.get('compression', 'lzf')
        try:
            self.probes = json.loads(SAVE.get('probes', '[]'))
        except:
            raise ValueError('probe and probe_location have been merged. Update config file.')

        # Check probes
        if self.probes:
            for c in self.probes:
                if not 0 <= c[0] < self.nx or not 0 <= c[1] < self.nz:
                    raise ValueError('probes must be in the domain')

        # datapath and datafile
        if self.comp == 'None':
            self.comp = None

        # if self.savepath does not exist, create it
        self.check_dir(self.savepath)

        self.datafile = getattr(self.args, 'datafile', None)
        if not self.datafile:
            self.datafile = self.savepath / self.savefile
        else:
            self.datafile = pathlib.Path(self.datafile).expanduser()

    def _figs(self):

        FIGS = self.cfg['figures']
        self.figures = FIGS.getboolean('figures', True)
        self.show_probes = FIGS.getboolean('probes', True)
        self.show_pml = FIGS.getboolean('pml', True)
        self.bc_profiles = FIGS.getboolean('bc_profiles', True)
        self.fps = FIGS.getint('fps', 24)
