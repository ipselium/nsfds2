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
#
# pylint: disable=too-many-locals
"""
-----------

Plotting library for nsfds2

@author: Cyril Desjouy
"""


import h5py
import numpy as _np
import matplotlib.pyplot as _plt
from ofdlib2 import fdtd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from fdgrid.mesh import plot_obstacles
from mpltools import modified_jet, MidpointNormalize


def show():
    """ Show all figures. """
    _plt.show()


def fields(cfg):
    """ Make figure """

    if not cfg.figures or not cfg.save:
        return None

    with h5py.File(cfg.savepath + cfg.filename + '.hdf5', 'r') as sfile:

        nt = sfile['nt'][...]
        x = sfile['x'][...]
        z = sfile['z'][...]
        obstacles = sfile['obstacles'][...]
        if cfg.onlyp:
            p = sfile[f'p_it{nt}'][...]
        else:
            rho = sfile[f'rho_it{nt}'][...]
            rhou = sfile[f'rhou_it{nt}'][...]
            rhov = sfile[f'rhov_it{nt}'][...]
            rhoe = sfile[f'rhoe_it{nt}'][...]
            u = rhou/rho
            v = rhov/rho
            e = rhoe/rho
            p = _np.empty_like(rho)
            fdtd.p(p, rho, rhou, rhov, rhoe, cfg.gamma)

    cm = modified_jet()
    norm = MidpointNormalize(midpoint=0)

    if cfg.onlyp:
        _, ax = _plt.subplots(figsize=(12, 9))
        im = ax.pcolorfast(x, z, (p-cfg.p0).T, cmap=cm, norm=norm)
        plot_obstacles(x, z, ax, obstacles)
        ax.set_xlabel(r'$x$ [m]')
        ax.set_ylabel(r'$z$ [m]')
        ax.set_aspect('equal')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        _plt.colorbar(im, cax=cax)
        if cfg.probes:
            ax.plot(*cfg.probes_loc, 'ro')

    else:
        _, axes = _plt.subplots(2, 2, figsize=(12, 9))

        im1 = axes[0, 0].pcolorfast(x, z, (p-cfg.p0).T, cmap=cm, norm=norm)
        im2 = axes[0, 1].pcolorfast(x, z, u.T, cmap=cm)
        im3 = axes[1, 0].pcolorfast(x, z, v.T, cmap=cm)
        im4 = axes[1, 1].pcolorfast(x, z, e.T, cmap=cm)
        ims = [im1, im2, im3, im4]

        for ax, im in zip(axes.ravel(), ims):
            plot_obstacles(x, z, ax, obstacles)
            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            ax.set_aspect('equal')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            _plt.colorbar(im, cax=cax)
            if cfg.probes and cfg.probes_loc:
                ax.plot(*cfg.probes_loc, 'ro')

    return None


def probes(cfg):
    """ Plot probes. """

    if cfg.figures and cfg.probes:

        with h5py.File(cfg.savepath + cfg.filename + '.hdf5', 'r') as sfile:

            nt = sfile['nt'][...]
            dt = sfile['dt'][...]
            p0 = sfile['p0'][...]
            pressure = sfile['probes'][...] - p0
            probes_loc = sfile['probes_location'][...]

        t = _np.arange(nt)*dt

        _, ax = _plt.subplots(figsize=(9, 4))
        for i, c in enumerate(probes_loc):
            ax.plot(t, pressure[i, :], label=f'@{tuple(c)}')
        ax.set_xlim(t.min(), t.max())
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Pressure [Pa]')
        ax.legend()
        ax.grid()
