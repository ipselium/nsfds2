#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2016-2019 Cyril Desjouy <cyril.desjouy@univ-lemans.fr>
#
# This file is part of {name}
#
# {name} is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# {name} is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with {name}. If not, see <http://www.gnu.org/licenses/>.
#
# Creation Date : 2019-06-04 - 11:40:41
"""
-----------
DOCSTRING

@author: Cyril Desjouy
"""

import argparse



def parse():

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
    commands.add_parser("movie", parents=[commons, view, data],
                        description="Make a movie from existing results",
                        help="make movie from existing results")
    shw = commands.add_parser("show",
                              description="Helper commands for parameters/results analysis",
                              help="show results and simulation configuration")

    # show section subsubparsers : frame/probe/
    show_cmds = shw.add_subparsers(dest='show_command',
                                   help='see -h for further help')
    show_cmds.add_parser('frame', parents=[commons, view, data],
                         description="Extract frame from hdf5 file and display it",
                         help="show results at a given iteration")
    show_cmds.add_parser('probes', parents=[commons, data],
                         description="Display pressure at probe locations",
                         help="plot pressure at probes locations")
    show_cmds.add_parser('grid', parents=[commons],
                         description="Display grid mesh",
                         help="show grid mesh")
    show_cmds.add_parser('domains', parents=[commons],
                         description="Display subdomains",
                         help="show domain decomposition")
    show_cmds.add_parser('parameters', parents=[commons],
                         description="Display some simulation parameters",
                         help="display some simulation parameters")

    return root.parse_args()


if __name__ == '__main__':
    args = parse()
    print(args)
