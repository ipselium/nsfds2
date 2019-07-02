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

-----------
"""

import os
import argparse
from fdgrid import mesh
from nsfds2.init import CfgSetup, Fields
from nsfds2.fdtd import FDTD
from nsfds2.utils import files, headers, graphics, sounds


def parse_args():
    """ Parse arguments. """

    # Options gathered in some parsers
    commons = argparse.ArgumentParser(add_help=False)
    commons.add_argument('-q', '--quiet', action="store_true", help='quiet mode')
    commons.add_argument('-c', '--cfg-file', metavar='CF', dest='cfgfile',
                         help='path to config file')

    view = argparse.ArgumentParser(add_help=False)
    view.add_argument('-i', dest='nt', type=int, help='number of time iterations')
    view.add_argument('-r', dest='ref', type=int, help='reference frame for colormap')
    view.add_argument('view', nargs='*', default='p',
                      choices=['p', 'rho', 'vx', 'vz', 'vort', 'e'],
                      help='variable(s) to plot')

    data = argparse.ArgumentParser(add_help=False)
    data.add_argument('-d', '--dat-file', metavar='DF', dest='datafile',
                      help='path to hdf5 data file')

    time = argparse.ArgumentParser(add_help=False)
    time.add_argument('-t', '--timings', action="store_true", default=None,
                      help='Display complete timings')

    description = 'A Navier-Stokes Finite Difference Solver'
    root = argparse.ArgumentParser(prog='nsfds2', description=description)

    # Subparsers : solve/movie/show commands
    commands = root.add_subparsers(dest='command',
                                   help='see nsfds2 `command` -h for further help')


    commands.add_parser("solve", parents=[commons, view, data, time],
                        description="Navier-Stokes equation solver",
                        help="solve NS equations with given configuration")
    shw = commands.add_parser("show",
                              description="Helper commands for parameters/results analysis",
                              help="show results and simulation configuration")
    mak = commands.add_parser("make",
                              description="Make movie/sound files",
                              help="make movie/sound files")


    # show section subsubparsers : frame/probe/
    shw_cmds = shw.add_subparsers(dest='show_command',
                                  help='see -h for further help')
    shw_cmds.add_parser('frame', parents=[commons, view, data],
                        description="Extract frame from hdf5 file and display it",
                        help="show results at a given iteration")
    shw_cmds.add_parser('probes', parents=[commons, data],
                        description="Display pressure at probe locations",
                        help="plot pressure at probes locations")
    shw_cmds.add_parser('spectrogram', parents=[commons, data],
                        description="Display spectrograms at probe locations",
                        help="plot spectrograms at probes locations")
    shw_cmds.add_parser('grid', parents=[commons],
                        description="Display numerical grid mesh",
                        help="show numerical grid mesh")
    shw_cmds.add_parser('pgrid', parents=[commons],
                        description="Display physical grid mesh",
                        help="show physical grid mesh")
    shw_cmds.add_parser('domains', parents=[commons],
                        description="Display subdomains",
                        help="show domain decomposition")
    shw_cmds.add_parser('parameters', parents=[commons],
                        description="Display some simulation parameters",
                        help="display some simulation parameters")

    # make section subsubparsers : movie/wav
    mak_cmds = mak.add_subparsers(dest='make_command',
                                  help='see -h for further help')
    mak_cmds.add_parser("movie", parents=[commons, view, data],
                        description="Make a movie from existing results",
                        help="make movie from existing results")
    mak_cmds.add_parser("sound", parents=[commons, data],
                        description="Make sound files from existing results",
                        help="make sound files from existing results")

    return root.parse_args()


def solve(cfg, msh):
    """ Solve NS equations. """

    # Simulation parameters
    fld = Fields(msh, cfg)

    # Simulation
    fdtd = FDTD(msh, fld, cfg)
    fdtd.run()

    if cfg.figures:
        plt = graphics.Plot(cfg.datafile, quiet=cfg.quiet)
        if cfg.save:
            plt.fields(view=cfg.args.view, ref=cfg.args.ref,
                       show_pml=cfg.show_pml, show_probes=cfg.show_probes)
        if cfg.probes:
            plt.probes()

        plt.show()


def show(cfg, msh):
    """ Show simulation parameters and grid. """


    if cfg.args.show_command == 'parameters':
        headers.version()
        headers.check_geo(cfg)
        headers.parameters(cfg, msh)

    elif cfg.args.show_command == 'grid':
        msh.plot_grid(axis=True, pml=cfg.show_pml, bc_profiles=cfg.bc_profiles,
                      probes=cfg.probes if cfg.show_probes else False)

    elif cfg.args.show_command == 'pgrid':
        msh.plot_physical(pml=cfg.show_pml, bc_profiles=cfg.bc_profiles,
                          probes=cfg.probes if cfg.show_probes else False)

    elif cfg.args.show_command == 'domains':
        if cfg.quiet:
            msh.plot_domains(legend=False)
        else:
            msh.plot_domains(legend=True)

    elif cfg.args.show_command == 'frame':
        plt = graphics.Plot(cfg.datafile, quiet=cfg.quiet)
        plt.fields(view=cfg.args.view, iteration=cfg.args.nt, ref=cfg.args.ref,
                   show_pml=cfg.show_pml, show_probes=cfg.show_probes)

    elif cfg.args.show_command == 'probes':
        plt = graphics.Plot(cfg.datafile, quiet=cfg.quiet)
        plt.probes()

    elif cfg.args.show_command == 'spectrogram':
        plt = graphics.Plot(cfg.datafile, quiet=cfg.quiet)
        plt.spectrogram()


    else:
        headers.copyright()
        headers.version()

    msh.show_figures()


def make(cfg, _):
    """ Create a movie from a dataset. """

    if cfg.args.make_command == 'movie':

        plt = graphics.Plot(cfg.datafile, quiet=cfg.quiet)
        plt.movie(view=cfg.args.view, nt=cfg.args.nt, ref=cfg.args.ref,
                  show_pml=cfg.show_pml, show_probes=cfg.show_probes)
        plt.show()

    elif cfg.args.make_command == 'sound':
        _ = sounds.probes_to_wave(cfg.datafile)


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
    elif cfg.mesh == 'adaptative':
        msh = mesh.AdaptativeMesh((cfg.nx, cfg.nz), (cfg.dx, cfg.dz),
                                  origin=(cfg.ix0, cfg.iz0),
                                  bc=cfg.bc, obstacles=obstacles, Npml=cfg.Npml,
                                  stencil=cfg.stencil)

    if args.command:
        globals()[args.command](cfg, msh)
    else:
        headers.copyright()
        print('Must specify an action among solve/movie/show')
        print('See nsfds2 -h for help')


if __name__ == "__main__":

    os.nice(20)
    main()
