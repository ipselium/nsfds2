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
# Creation Date : 2019-04-19 - 15:14:22
"""
-----------

Utils

-----------
"""

import re
import warnings


class Check:
    """ Check FDTD parameters. """

    def __init__(self, cfg, msh):

        self.cfg = cfg
        self.msh = msh

        self.source()
        self.obstacles()
        self.domain()

    def source(self):
        """ Check if source is not in an obstacle. """
        for o in self.msh.obstacles:
            if o.ix[0] < self.cfg.ixS < o.ix[1] and o.iz[0] < self.cfg.izS < o.iz[1]:
                raise ValueError('source cannot be in an obstacle')

    def obstacles(self):
        """ Check validity of obstacles boundary conditions. """

        for obs in self.msh.obstacles:
            if not re.match(r'^[WV]+$', obs.bc):
                s = "Obstacles can only be combination of 'WV' for now. "
                raise ValueError(s)

    def domain(self):
        """ Check validity of the bcs. """

        if not re.match(r'[PRAW][PRAW][PRAW][PRAW]', self.msh.bc):
            s = "Only 'W', 'R', 'P', and 'A' bc are implemented for now. "
            s += "Fix bcs to 'WWWW'."
            warnings.warn(s, stacklevel=8)
            self.msh.bc = 'WWWW'
