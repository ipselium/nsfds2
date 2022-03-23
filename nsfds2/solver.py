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
import time
import argparse
import pathlib
from multiprocessing import cpu_count, Pool, current_process
from fdgrid import mesh
from nsfds2.init import CfgSetup, create_template
from nsfds2.fdtd import FDTD
from nsfds2.utils import files, headers, graphics, sounds


class VirtualArguments:

    def __init__(self, cfgfile, command='solve', quiet=True):
        """Generate a set of virtual arguments as for configparser"""
        self.cfgfile = cfgfile
        self.command = command
        self.quiet = quiet


def _wrapper(run, args):
    """Wrapper function for ploop jobs."""
    cur = current_process().name
    ti = time.perf_counter()
    print(f'{cur} - Solving {args.cfgfile}')
    run(args)
    ts = time.perf_counter() - ti
    print(f'{cur} - {args.cfgfile} solved in {ts:.2f} s')
    ti = time.perf_counter()
    cfg = CfgSetup(args=args)
    plt = graphics.Plot(cfg.datafile, quiet=args.quiet)
    plt.movie()
    tm = time.perf_counter() - ti
    print(f'{cur} - {args.cfgfile} movie made in {tm:.2f} s')


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
    view.add_argument('-l', dest='logscale', action="store_true", help='Plot in log scale')
    view.add_argument('view', nargs='*', default='p',
                      choices=['p', 'rho', 'vx', 'vz', 'vxz', 'e'],
                      help='variable(s) to plot')

    data = argparse.ArgumentParser(add_help=False)
    data.add_argument('-d', '--dat-file', metavar='DF', dest='datafile',
                      help='path to hdf5 data file')

    path = argparse.ArgumentParser(add_help=False)
    path.add_argument('-p', dest='path', type=str, required=True,
                      help='loop over this path')

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
    loop = commands.add_parser("loop", parents=[path],
                               description="Loop over .conf in a directory",
                               help="Solve multiple configurations")
    ploop = commands.add_parser("ploop", parents=[path],
                                description="Loop over .conf in a directory [parallel version]",
                                help="Solve multiple configurations")

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

    # make section subsubparsers : movie/wav/template
    mak_cmds = mak.add_subparsers(dest='make_command',
                                  help='see -h for further help')
    mak_cmds.add_parser("movie", parents=[commons, view, data],
                        description="Make a movie from existing results",
                        help="make movie from existing results")
    mak_cmds.add_parser("sound", parents=[commons, data],
                        description="Make sound files from existing results",
                        help="make sound files from existing results")
    mak_cmds.add_parser("template", parents=[commons, data],
                        description="Create basic configuration file",
                        help="Create basic configuration file")

    return root.parse_args()


def solve(cfg, msh):
    """ Solve NS equations. """

    # Simulation
    fdtd = FDTD(msh, cfg)
    fdtd.run()

    if cfg.figures:
        plt = graphics.Plot(cfg.datafile, quiet=cfg.quiet)
        if cfg.save_fields:
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
                   show_pml=cfg.show_pml, show_probes=cfg.show_probes,
                   logscale=cfg.args.logscale)

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
                  show_pml=cfg.show_pml, show_probes=cfg.show_probes,
                  fps=cfg.fps, logscale=cfg.args.logscale)
        plt.show()

    elif cfg.args.make_command == 'sound':
        _ = sounds.probes_to_wave(cfg.datafile)


def template(args):
    """Make template."""
    if not args.cfgfile:
        print('Path/filename must be specified with -c option')
    else:
        cfgfile = pathlib.Path(args.cfgfile).expanduser()
        path, filename = cfgfile.parent, cfgfile.stem + cfgfile.suffix
        create_template(path=path, filename=filename)
        print(f"{cfgfile} created")


def loop(path):
    """Loop simulations over .conf files"""
    headers.copyright()
    path = pathlib.Path(path)
    files = [f for f in pathlib.os.listdir(path.expanduser())
             if f.endswith('conf')]
    ti = time.perf_counter()
    for file in files:
        print(f'Processing {file}...')
        args = VirtualArguments(path.expanduser() / file)
        run(args=args)
        print(f'Making movie for {file}...')
        cfg = CfgSetup(args=args)
        plt = graphics.Plot(cfg.datafile, quiet=args.quiet)
        plt.movie()
    print(f'Simulations ended in {time.perf_counter() - ti:.2f}')


def ploop(path):
    """Loop simulations over .conf files [separate processes]."""
    headers.copyright()
    path = pathlib.Path(path)
    files = [f for f in pathlib.os.listdir(path.expanduser())
             if f.endswith('conf')]
    ti = time.perf_counter()
    with Pool(processes=(cpu_count() - 1)) as pool:
        for file in files:
            args = VirtualArguments(path.expanduser() / file)
            pool.apply_async(_wrapper, args=(run, args))
        pool.close()
        pool.join()
    print(f'Simulations ended in {time.perf_counter() - ti:.2f}')


def run(args):
    """Run nsfds2."""

    # Parse config file
    cfg = CfgSetup(args=args)

    # Mesh
    msh = mesh.build(cfg.mesh, (cfg.nx, cfg.nz), (cfg.dx, cfg.dz),
                     origin=(cfg.ix0, cfg.iz0),
                     bc=cfg.bc, obstacles=cfg.obstacles,
                     Npml=cfg.Npml, stencil=cfg.stencil,
                     dilatation=cfg.Rx, Nd=cfg.Nd, only_pml=cfg.only_pml,
                     fcurvxz=files.get_curvilinear(cfg))

    if args.command:
        globals()[args.command](cfg, msh)
    else:
        headers.copyright()
        print('Must specify an action among solve/make/show/loop')
        print('See nsfds2 -h for help')


def main(args=None):
    """ Main """

    # Parse arguments
    if not args:
        args = parse_args()

    # Command
    if getattr(args, 'make_command', None) == 'template':
        template(args)
    elif getattr(args, 'command', None) == 'loop':
        loop(args.path)
    elif getattr(args, 'command', None) == 'ploop':
        ploop(args.path)
    else:
        run(args)


if __name__ == "__main__":

    os.nice(20)
    main()
