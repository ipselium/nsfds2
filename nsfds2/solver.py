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
# Creation Date : 2019-03-01 - 12:05:08
"""
-----------

Navier Stokes Finite Differences Solver

@author: Cyril Desjouy
"""

from fdgrid import mesh, templates
from nsfds2.init import CfgSetup, Coefficients, Fields
from nsfds2.lib import FDTD
from nsfds2.utils import figures, headers


def main():
    """ Main """

    # Headers
    headers.copyright()
    headers.version()

    # Parse Config
    cfg = CfgSetup()

    # Mesh
    obstacle = templates.square(cfg.nx, cfg.nz)
    msh = mesh.Mesh((cfg.nx, cfg.nz),
                    (cfg.dx, cfg.dz),
                    origin=(cfg.ix0, cfg.iz0),
                    bc=cfg.bc, obstacles=obstacle,
                    Npml=cfg.Npml,
                    stencil=cfg.stencil)

    # Simulation parameters
    fld = Fields(msh, cfg)
    cff = Coefficients(cfg, msh.stencil)

    # Wait for user input
    _ = input('Hit enter to continue...')


    # Simulation
    fdtd = FDTD(msh, fld, cff, cfg)
    p, rho, rhou, rhov, rhoe = fdtd.run()

    # Figures
    if cfg.figures:
        figures.fields(p, rhou/rho, rhov/rho, rhoe/rho, msh, cfg)

if __name__ == "__main__":
    import os
    os.nice(20)
    main()
