#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2016-2020 Cyril Desjouy <cyril.desjouy@univ-lemans.fr>
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
# Creation Date : 2022-03-15 - 22:17:54
"""
-----------
DOCSTRING

-----------
"""

import pathlib
import nsfds2.solver as ns
from nsfds2.utils import graphics
from nsfds2.init import CfgSetup

class Args:
    def __init__(self, cfgfile):
        self.cfgfile = cfgfile
        self.command = 'solve'
        self.quiet = True


if __name__ == '__main__':

    path = pathlib.Path('~/.nsfds2/')

    for file in pathlib.os.listdir(path.expanduser()):
        if file.endswith('conf'):
            print(f'Processing {file}...')
            args = Args(path.expanduser() / file)
            ns.main(args=args)
            print(f'Making movie for {file}...')
            cfg = CfgSetup(args=args)
            plt = graphics.Plot(cfg.datafile, quiet=args.quiet)
            plt.movie()
