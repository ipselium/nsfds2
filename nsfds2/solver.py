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
import argparse
from fdgrid import mesh
from nsfds2.init import CfgSetup, Fields
from nsfds2.fdtd import FDTD
from nsfds2.utils import files, headers, graphics


def parse_args():
    """ Parse arguments. """

    parser = argparse.ArgumentParser(prog='nsfds',
                                     description='A Navier-Stokes Finite Difference Solver')

    subparsers = parser.add_subparsers(dest='command')
    solve_parser = subparsers.add_parser("solve")
    movie_parser = subparsers.add_parser("movie")
    show_parser = subparsers.add_parser("show")

    # Common
    for p in [solve_parser, show_parser]:
        p.add_argument('-g', '--geo-file', metavar='GF', dest='geofile', nargs=2,
                       help='file and pattern to use')
        p.add_argument('-c', '--cfg-file', metavar='CF', dest='cfgfile',
                       help='path to config file')

    for p in [solve_parser, movie_parser]:
        p.add_argument('-q', '--quiet', action="store_true", help='Quiet mode')
        p.add_argument('-i', dest='nt', type=int,
                       help='Number of time iterations')
        p.add_argument('-d', '--dat-file', metavar='DF', dest='datafile',
                       help='path to hdf5 data file')


    # Movie parser
    movie_parser.add_argument('view', nargs='?', default='p',
                              choices=['p', 'rho', 'vx', 'vz', 'e'])
    movie_parser.add_argument('-r', dest='ref', type=int,
                              help='Reference frame for colormap')

    # Show parser
    show_parser.add_argument('view', nargs='?', default='grid',
                             choices=['grid', 'domains', 'all'])

    # Solver parser
    solve_parser.add_argument('-t', '--timings', action="store_true", default=None,
                              help='Display complete timings')

    return parser.parse_args()


def solve(cfg, msh):
    """ Solve NS equations. """

    # Simulation parameters
    fld = Fields(msh, cfg)

    # Simulation
    fdtd = FDTD(msh, fld, cfg)
    fdtd.run()

    # Figures
    graphics.fields(cfg)
    graphics.probes(cfg)
    graphics.show()


def show(cfg, msh):
    """ Show simulation parameters and grid. """

    headers.check_geo(cfg)
    headers.parameters(cfg)

    if cfg.args.view == 'grid':
        msh.plot_grid()

    elif cfg.args.view == 'domains':
        print('\n')
        msh.plot_domains(legend=True)

    elif cfg.args.view == 'all':
        print('\n')
        msh.plot_grid()
        msh.plot_domains(legend=True)

    msh.show_figures()


def movie(cfg, _):
    """ Create a movie from a dataset. """

    movie = graphics.Movie(cfg.datafile, view=cfg.args.view,
                           ref=cfg.args.ref, nt=cfg.nt, quiet=cfg.quiet)
    movie.make()


def main():
    """ Main """

    # Parse arguments
    args = parse_args()

    # Parse config file
    cfg = CfgSetup(args=args)

    # Geometry
    obstacles = files.get_obstacle(cfg)

    # Mesh
    if cfg.mesh == 'regular':
        msh = mesh.Mesh((cfg.nx, cfg.nz), (cfg.dx, cfg.dz), origin=(cfg.ix0, cfg.iz0),
                        bc=cfg.bc, obstacles=obstacles, Npml=cfg.Npml,
                        stencil=cfg.stencil)
    elif cfg.mesh == 'curvilinear':
        fcurv = files.get_curvilinear(cfg)
        msh = mesh.CurvilinearMesh((cfg.nx, cfg.nz), (cfg.dx, cfg.dz),
                                   origin=(cfg.ix0, cfg.iz0),
                                   bc=cfg.bc, obstacles=obstacles, Npml=cfg.Npml,
                                   stencil=cfg.stencil, fcurvxz=fcurv)

    if args.command:
        globals()[args.command](cfg, msh)
    else:
        headers.copyright()
        print('Must specify an action among solve/movie/show')
        print('See nsfds2 -h for help')


if __name__ == "__main__":

    os.nice(20)
    main()
