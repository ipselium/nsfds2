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
# Creation Date : 2018-04-13 00:48:58
"""
-----------

Compute viscous fluxes

@author: Cyril Desjouy
"""


import numpy as np
from .init import Param
from ofdlib.fdtd import integrate
from ofdlib.dschm import dudx3cG, dudz3cG
from ofdlib.viscosity import cErhou, cErhov, cErhoe, cFrhov, cFrhoe, cMu
from ofdlib.viscosity import tau11x, tau12x, tau12z, tau22z


class ViscousFluxes(Param):
    """ Compute viscous flux : update rhou, rhov and rhoe only. """

    def __init__(self):

        super(ViscousFluxes, self).__init__()
        super(ViscousFluxes, self).init_coefficients()

    def base(self, rho, rhou, rhov, rhoe):
        """
            Viscous flux integration : interior points [optimized]
        """

        # Array initialization
        tau11 = np.zeros_like(rho)
        tau12 = np.zeros_like(rho)
        tau22 = np.zeros_like(rho)
        mu = np.zeros_like(rho)
        Krhou = np.zeros((self.nbx, self.nbz))
        Krhov = np.zeros((self.nbx, self.nbz))
        Krhoe = np.zeros((self.nbx, self.nbz))
        Frhou = np.zeros((self.nbx, self.nbz))
        Frhov = np.zeros((self.nbx, self.nbz))
        Frhoe = np.zeros((self.nbx, self.nbz))
        Erhoe = np.zeros((self.nbx, self.nbz))

        # dE/dx term
        Erhou = rhou/rho
        Erhov = rhov/rho

        # Strain tensor
        tau11 = tau11x(self.nbx, self.nbz, tau11, Erhou, Erhov, self.one_dx)
        tau12 = tau12x(self.nbx, self.nbz, tau12, Erhou, Erhov, self.one_dx)
        tau12 = tau12z(self.nbx, self.nbz, tau12, Erhou, Erhov, self.one_dz)
        tau22 = tau22z(self.nbx, self.nbz, tau22, Erhou, Erhov, self.one_dz)

        # Dynamic viscosity
        mu = cMu(self.nbx, self.nbz, rho, self.nu, mu)

        # dE/dx term
        Erhou = cErhou(self.nbx, self.nbz, mu, Erhou, tau11, tau22)
        Erhov = cErhov(self.nbx, self.nbz, mu, Erhov, tau12)
        Erhoe = cErhoe(self.nbx, self.nbz, mu, Erhoe, tau11, tau12, tau22, rho, rhou, rhov)

        # dF/dz term
        Frhou = Erhov
        Frhov = cFrhov(self.nbx, self.nbz, mu, Frhov, tau11, tau22)
        Frhoe = cFrhoe(self.nbx, self.nbz, mu, Frhoe, tau11, tau12, tau22, rho, rhou, rhov)

        # viscous flux : order 2 centered scheme
        Krhou = dudx3cG(Erhou, self.one_dx, self.nbx, self.nbz, Krhou)
        Krhov = dudx3cG(Erhov, self.one_dx, self.nbx, self.nbz, Krhov)
        Krhoe = dudx3cG(Erhoe, self.one_dx, self.nbx, self.nbz, Krhoe)

        Krhou = dudz3cG(Frhou, self.one_dz, self.nbx, self.nbz, Krhou)
        Krhov = dudz3cG(Frhov, self.one_dz, self.nbx, self.nbz, Krhov)
        Krhoe = dudz3cG(Frhoe, self.one_dz, self.nbx, self.nbz, Krhoe)

        # Integrate in time
        rhou = integrate(self.nbx, self.nbz, self.dt, rhou, Krhou)
        rhov = integrate(self.nbx, self.nbz, self.dt, rhov, Krhov)
        rhoe = integrate(self.nbx, self.nbz, self.dt, rhoe, Krhoe)
