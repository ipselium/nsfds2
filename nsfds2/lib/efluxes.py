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
# Creation Date : 2018-04-14 01:46:17
#
# pylint: disable=too-many-instance-attributes
"""
-----------

Compute Eulerian fluxes

@author: Cyril Desjouy
"""


import re
import ofdlib2.derivation as drv
from nsfds2.utils.array import empty_like
from .cin import Cin


class EulerianFluxes:
    """ Compute Eulerian fluxes. """

    def __init__(self, msh, fld, cfg):

        self.msh = msh
        self.fld = fld
        self.cfg = cfg
        self.p, self.r, self.ru, self.rv, self.re = empty_like(msh.shape, 5)
        self.ccin = Cin(msh, fld, cfg)
        if 'A' in msh.bc:
            self.init_pml()

    def rk4(self):
        """
          Avancement de la solution en temps à l'aide d'un algortihme de
          Runge-Kutta à 6 étapes
        """

        self.p = self.fld.p.copy()
        self.r = self.fld.r.copy()
        self.ru = self.fld.ru.copy()
        self.rv = self.fld.rv.copy()
        self.re = self.fld.re.copy()

        for irk in range(1, 7):

            # Eulerian fluxes
            self.cin()
            if 'A' in self.msh.bc:
                self.pml()

            self.fld.fdtools.adtime(self.fld.r, self.r, self.fld.K, irk)
            self.fld.fdtools.adtime(self.fld.ru, self.ru, self.fld.Ku, irk)
            self.fld.fdtools.adtime(self.fld.rv, self.rv, self.fld.Kv, irk)
            self.fld.fdtools.adtime(self.fld.re, self.re, self.fld.Ke, irk)

            # Boundary conditions
            self.cout()

            # Compute p
            self.fld.fdtools.p(self.fld.p, self.fld.r, self.fld.ru,
                               self.fld.rv, self.fld.re)

            if self.cfg.stype in ['harmonic', 'white']:
                self.fld.p += self.fld.update_source(self.cfg.it)

        if 'A' in self.msh.bc:
            self.update_pml()

    def cin(self):
        """ Interior domain. """
        self.ccin.dispatch()

    def cout(self):
        """ Boundaries. """

        self.cout_obstacles()
        self.cout_bc()
        self.cout_periodic()
        self.cout_source()

    def cout_obstacles(self):
        """ Obstacle walls. """

        for s in self.msh.obstacles:

            self.fld.ru[s.sx, s.iz[0]] = 0
            self.fld.rv[s.sx, s.iz[0]] = 0

            self.fld.ru[s.sx, s.iz[1]] = 0
            self.fld.rv[s.sx, s.iz[1]] = 0

            self.fld.ru[s.ix[0], s.sz] = 0
            self.fld.rv[s.ix[0], s.sz] = 0

            self.fld.ru[s.ix[1], s.sz] = 0
            self.fld.rv[s.ix[1], s.sz] = 0

    def cout_bc(self):
        """ Boundary conditions.

        Note
        ----

        For PML, following Hu 2002 :
            * ymin, ymax, xmax : p=0
            * and xmin : p=v=0
            * or periodic bc
        """

        if self.msh.bc[0] in ['R', 'A']:
            self.fld.ru[0, :] = 0
            self.fld.rv[0, :] = 0

        if self.msh.bc[1] in ['R', 'A']:
            self.fld.ru[:, 0] = 0
            self.fld.rv[:, 0] = 0

        if self.msh.bc[2] in ['R', 'A']:
            self.fld.ru[-1, :] = 0
            self.fld.rv[-1, :] = 0

        if self.msh.bc[3] in ['R', 'A']:
            self.fld.ru[:, -1] = 0
            self.fld.rv[:, -1] = 0

    def cout_periodic(self):
        """ Additional rigid boundaries for periodic condition. """

        if re.match(r'P.P.', self.msh.bc):
            for s in self.msh.dxdomains.additional_rigid_bc:
                self.fld.ru[s] = 0
                self.fld.rv[s] = 0

        if re.match(r'.P.P', self.msh.bc):
            for s in self.msh.dzdomains.additional_rigid_bc:
                self.fld.ru[s] = 0
                self.fld.rv[s] = 0

    def cout_source(self):
        """ Wall sources. """

        for obs in self.msh.obstacles:
            for bc in obs.edges:
                if bc.type == 'U':
                    self.fld.ru[bc.sx, bc.sz] = bc.prf*self.fld.update_wall(self.cfg.it, bc.f)
                elif bc.type == 'V':
                    self.fld.rv[bc.sx, bc.sz] = bc.prf*self.fld.update_wall(self.cfg.it, bc.f)

    def init_pml(self):
        """ Initialize PMLs. """

        self.du = drv.duA(self.msh.x, self.msh.z,
                          self.fld.sx, self.fld.sz, self.cfg.beta,
                          self.cfg.dt, stencil=self.cfg.stencil)

        for sub in self.msh.adomains:
            sub.pml_method = getattr(self.du, f'du_{sub.bc}')
            sub.pml_intgrt = getattr(self.du, f'update_{sub.axname}')

    def pml(self):
        """ PMLs."""

        # Remove initial field
        self.fld.E -= self.fld.Ei
        self.fld.Eu -= self.fld.Eui
        self.fld.Ev -= self.fld.Evi
        self.fld.Ee -= self.fld.Eei

        self.fld.F -= self.fld.Fi
        self.fld.Fu -= self.fld.Fui
        self.fld.Fv -= self.fld.Fvi
        self.fld.Fe -= self.fld.Fei

        for sub in self.msh.adomains:
            sub.pml_method(self.fld.E, self.fld.F, self.fld.qx, self.fld.qz,
                           self.fld.Kx, self.fld.Kz, self.fld.K, *sub.ix, *sub.iz)
            sub.pml_method(self.fld.Eu, self.fld.Fu, self.fld.qux, self.fld.quz,
                           self.fld.Kux, self.fld.Kuz, self.fld.Ku, *sub.ix, *sub.iz)
            sub.pml_method(self.fld.Ev, self.fld.Fv, self.fld.qvx, self.fld.qvz,
                           self.fld.Kvx, self.fld.Kvz, self.fld.Kv, *sub.ix, *sub.iz)
            sub.pml_method(self.fld.Ee, self.fld.Fe, self.fld.qex, self.fld.qez,
                           self.fld.Kex, self.fld.Kez, self.fld.Ke, *sub.ix, *sub.iz)

    def update_pml(self):
        """ Integrate PMLs. """

        qx = self.fld.qx.copy()
        qux = self.fld.qux.copy()
        qvx = self.fld.qvx.copy()
        qex = self.fld.qex.copy()

        qz = self.fld.qz.copy()
        quz = self.fld.quz.copy()
        qvz = self.fld.qvz.copy()
        qez = self.fld.qez.copy()

        for irk in range(1, 7):

            for sub in self.msh.adomains:

                sub.pml_intgrt(self.fld.qx, self.fld.qz, qx, qz,
                               self.fld.Kx, self.fld.Kz, self.fld.E, irk,
                               *sub.ix, *sub.iz)
                sub.pml_intgrt(self.fld.qux, self.fld.quz, qux, quz,
                               self.fld.Kux, self.fld.Kuz, self.fld.Eu, irk,
                               *sub.ix, *sub.iz)
                sub.pml_intgrt(self.fld.qvx, self.fld.qvz, qvx, qvz,
                               self.fld.Kvx, self.fld.Kvz, self.fld.Ev, irk,
                               *sub.ix, *sub.iz)
                sub.pml_intgrt(self.fld.qex, self.fld.qez, qex, qez,
                               self.fld.Kex, self.fld.Kez, self.fld.Ee, irk,
                               *sub.ix, *sub.iz)
