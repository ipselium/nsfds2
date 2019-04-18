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

import re
import os
import warnings
from fdgrid import mesh
from nsfds2.init import CfgSetup, Fields
from nsfds2.lib import FDTD
from nsfds2.utils import figures, headers, files


def check_source(xs, zs, obstacles):
    """ Check if source is not in an obstacle. """
    for obs in obstacles:
        if obs.ix[0] < xs < obs.ix[1] and obs.iz[0] < zs < obs.iz[1]:
            raise ValueError('source cannot be in an obstacle')


def check_obstacles(obstacles):
    """ Check validity of obstacles boundary conditions. """
    flag = False

    for obs in obstacles:
        if obs.bc is not 'RRRR':
            s = "Obstacles can only be 'RRRR' for now. "
            s += "Fix bcs to 'RRRR'."
            warnings.warn(s, stacklevel=8)
            obs.bc = 'RRRR'


def check_domain(domain):
    """ Check validity of the bcs. """
    if not re.match(r'[PRA][PRA][PRA][PRA]', domain.bc):
        s = "Only 'R' and 'P' bc are implemented for now. "
        s += "Fix bcs to 'RRRR'."
        warnings.warn(s, stacklevel=8)
        domain.bc = 'RRRR'


def main():
    """ Main """

    # Headers
    headers.copyright()
    headers.version()

    # Parse Config
    cfg = CfgSetup()

    # Geometry
    obstacles = files.get_obstacle(cfg)
    check_obstacles(obstacles)

    # Check source location
    check_source(cfg.ixS, cfg.izS, obstacles)

    # Mesh
    msh = mesh.Mesh((cfg.nx, cfg.nz),
                    (cfg.dx, cfg.dz),
                    origin=(cfg.ix0, cfg.iz0),
                    bc=cfg.bc, obstacles=obstacles,
                    Npml=cfg.Npml,
                    stencil=cfg.stencil)
    check_domain(msh)

    # Simulation parameters
    fld = Fields(msh, cfg)

    # Prompt user for start
    headers.start(cfg)

    # Simulation
    fdtd = FDTD(msh, fld, cfg)
    p, rho, rhou, rhov, rhoe = fdtd.run()

    # Figures
    if cfg.figures:
        figures.fields(p, rhou/rho, rhov/rho, rhoe/rho, msh, cfg)


if __name__ == "__main__":

    os.nice(20)
    main()
