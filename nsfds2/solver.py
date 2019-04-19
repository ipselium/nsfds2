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

import os
from fdgrid import mesh
from nsfds2.init import CfgSetup, Fields
from nsfds2.lib import FDTD
from nsfds2.utils import figures, files


def main():
    """ Main """

    # Parse Config
    cfg = CfgSetup()

    # Geometry
    obstacles = files.get_obstacle(cfg)

    # Mesh
    msh = mesh.Mesh((cfg.nx, cfg.nz), (cfg.dx, cfg.dz), origin=(cfg.ix0, cfg.iz0),
                    bc=cfg.bc, obstacles=obstacles, Npml=cfg.Npml, stencil=cfg.stencil)

    # Simulation parameters
    fld = Fields(msh, cfg)

    # Simulation
    fdtd = FDTD(msh, fld, cfg)
    fdtd.run()

    # Figures
    figures.fields(cfg)
    figures.probes(cfg)
    figures.show()


if __name__ == "__main__":

    os.nice(20)
    main()
