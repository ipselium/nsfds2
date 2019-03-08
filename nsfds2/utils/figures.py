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
# Creation Date : 2019-03-07 - 22:59:23
"""
-----------

Plotting library for nsfds2

@author: Cyril Desjouy
"""


import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpltools import modified_jet, MidpointNormalize


def fields(p, u, v, e, msh, cfg):
    """ Make figure """

    cm = modified_jet()
    norm = MidpointNormalize(midpoint=0)

    _, axes = plt.subplots(2, 2, figsize=(12, 9))

    im1 = axes[0, 0].pcolorfast(msh.x, msh.z, (p-cfg.p0).T, cmap=cm, norm=norm)
    im2 = axes[0, 1].pcolorfast(msh.x, msh.z, u.T, cmap=cm)
    im3 = axes[1, 0].pcolorfast(msh.x, msh.z, v.T, cmap=cm)
    im4 = axes[1, 1].pcolorfast(msh.x, msh.z, e.T, cmap=cm)
    ims = [im1, im2, im3, im4]

    for ax, im in zip(axes.ravel(), ims):
        msh.plot_obstacles(msh.x, msh.z, ax, msh.get_obstacles(), facecolor='y')
        ax.set_xlabel(r'$x$ [m]')
        ax.set_ylabel(r'$z$ [m]')
        ax.set_aspect('equal')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)

    plt.show()
