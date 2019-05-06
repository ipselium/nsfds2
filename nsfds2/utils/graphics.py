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
# Creation Date : 2019-05-02 - 14:09:56
#
# pylint: disable=too-many-locals
# pylint: disable=too-many-statements
# pylint: disable=too-many-instance-attributes

"""
-----------
DOCSTRING

@author: Cyril Desjouy
"""

import os
import sys
import getpass
import h5py
import numpy as _np
import matplotlib.pyplot as _plt
import matplotlib.animation as _ani
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ofdlib2 import fdtd as _fdtd
from progressbar import ProgressBar, Bar, ReverseBar, ETA
from mplutils import modified_jet, MidpointNormalize
from fdgrid.mesh import plot_subdomains as _plt_subdomains


__all__ = ['get_data', 'show', 'fields', 'probes', 'FrameGenerator', 'Movie']


def get_data(filename):
    """ Get data from filename. """

    try:
        data = h5py.File(filename, 'r')
    except OSError:
        print('You must provide a valid hdf5 file : moviemaker mydir/myfile.hdf5')
        sys.exit(1)
    else:
        return data


def show():
    """ Show all figures. """
    _plt.show()


def fields(cfg):
    """ Make figure """

    if not cfg.figures or not cfg.save:
        return None

    with h5py.File(cfg.datafile, 'r') as sfile:

        nt = sfile['nt'][...]
        x = sfile['x'][...]
        z = sfile['z'][...]
        obstacles = sfile['obstacles'][...]
        if cfg.mesh == 'curvilinear':
            J = sfile['J'][...]
        else:
            J = _np.ones((x.size, z.size))

        if cfg.onlyp:
            p = sfile[f'p_it{nt}'][...]
        else:
            rho = sfile[f'rho_it{nt}'][...]*J
            rhou = sfile[f'rhou_it{nt}'][...]*J
            rhov = sfile[f'rhov_it{nt}'][...]*J
            rhoe = sfile[f'rhoe_it{nt}'][...]*J
            u = rhou/rho
            v = rhov/rho
            e = rhoe/rho
            p = _np.empty_like(rho)
            _fdtd.p(p, rho, rhou, rhov, rhoe, cfg.gamma)


    cm = modified_jet()
    norm = MidpointNormalize(midpoint=0)

    if cfg.onlyp:
        _, ax = _plt.subplots(figsize=(12, 9))
        im = ax.pcolorfast(x, z, (p-cfg.p0).T, cmap=cm, norm=norm)
        _plt_subdomains(ax, x, z, obstacles)
        ax.set_xlabel(r'$x$ [m]')
        ax.set_ylabel(r'$z$ [m]')
        ax.set_aspect('equal')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        _plt.colorbar(im, cax=cax)
        if cfg.probes and cfg.probes_loc:
            ax.plot(*[[x[i], z[j]] for i, j in cfg.probes_loc], 'ro')

    else:
        _, axes = _plt.subplots(2, 2, figsize=(12, 9))

        im1 = axes[0, 0].pcolorfast(x, z, (p-cfg.p0).T, cmap=cm, norm=norm)
        im2 = axes[0, 1].pcolorfast(x, z, u.T, cmap=cm)
        im3 = axes[1, 0].pcolorfast(x, z, v.T, cmap=cm)
        im4 = axes[1, 1].pcolorfast(x, z, e.T, cmap=cm)
        ims = [im1, im2, im3, im4]

        for ax, im in zip(axes.ravel(), ims):
            _plt_subdomains(ax, x, z, obstacles)
            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            ax.set_aspect('equal')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            _plt.colorbar(im, cax=cax)
            if cfg.probes and cfg.probes_loc:
                ax.plot(*[[x[i], z[j]] for i, j in cfg.probes_loc], 'ro')

    return None


def probes(cfg):
    """ Plot probes. """

    if cfg.figures and cfg.probes:

        with h5py.File(cfg.datafile, 'r') as sfile:

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


class FrameGenerator:
    """ Frame Genrator

    Parameters
    ----------

    data : Path to hdf5 file (str) or data from hdf5 file
    view : The variable to display. string.
    ref : The reference frame for colormap. int.
    nt : The last frame to consider. int

    """

    def __init__(self, data, view='p', ref=20, nt=None):

        if isinstance(data, str):
            self.data = get_data(data)
        else:
            self.data = data
        self.view = view
        self.ref = ref
        self.nt = nt
        self.ns = data['ns'][...]
        self.icur = 0
        self.nx = self.data['nx'][...]
        self.nz = self.data['nz'][...]
        if self.data['mesh'][...] == 'curvilinear':
            self.J = self.data['J'][...]
        else:
            self.J = _np.ones((self.nx, self.nz))

    def reference(self):
        """ Generate the reference for min/max colormap values """

        if self.view == "p":
            ref = self.compute_p(self.ref)
        elif self.view in ['rho', 'vx', 'vz', 'e']:
            ref = self.data[f"{self.view}_it{self.ref}"][...]
        else:
            print("Only 'p', 'rho', 'vx', 'vz' and 'e' available !")
            sys.exit(1)

        ref = ref*self.J

        return ref[ref > 0].mean()*2, ref[ref < 0].mean()*2

    def compute_p(self, it):
        """ Compute p from rho, rhou, rhov and rhoe at iteration it. """

        rho = self.data[f"rho_it{it}"][...]*self.J
        rhou = self.data[f"rhou_it{it}"][...]*self.J
        rhov = self.data[f"rhov_it{it}"][...]*self.J
        rhoe = self.data[f"rhoe_it{it}"][...]*self.J
        p = _np.empty_like(rho)
        _fdtd.p(p, rho, rhou, rhov, rhoe, self.data['gamma'][...])

        return p - self.data['p0'][...]

    def _next_item(self):
        """ Generate next value of variable """

        if self.view == 'p':
            return self.compute_p(self.icur).T

        if self.view == 'rho':
            return self.data[f"{self.view}_it{self.icur}"][...].T

        if self.view in ['vx', 'vz', 'e']:
            vx = self.data[f"{self.view}_it{self.icur}"][...].T
            rho = self.data[f"rho_it{self.icur}"][...].T
            return vx/rho

        return None

    def __iter__(self):
        """ Iterator """

        return self

    def __next__(self):
        """ Next element of iterator : (frame_number, variable) """

        try:
            self.icur += self.ns
            if self.nt:
                if self.icur > self.nt:
                    raise StopIteration

            return self.icur, self._next_item()

        except KeyError:
            raise StopIteration


class Movie:
    """ Make movie from nsfds2 hdf5 ouput file.

    Parameters
    ----------

    filename : hdf5 file
    view : The variable to display. string.
    ref : The reference frame for colormap. int.
    nt : The last frame to consider. int.
    quiet : quiet mode. Boolean.
    """

    def __init__(self, filename, view='p', ref=None, nt=None, quiet=False):

        self.filename = filename
        self.path = os.path.dirname(filename) + '/'
        self.data = get_data(filename)
        self.view = view
        self.ref = ref
        self.nt = nt
        self.quiet = quiet
        self.curvilinear = True if self.data['mesh'][...] == 'curvilinear' else False
        self.obstacles = self.data['obstacles'][...]

        self.check()
        self._init_coordinates()
        self._init_generator()
        self._init_cmap()
        self._init_pbar()
        self._init_movie()

    def check(self):
        """ Check parameters. """

        if not self.nt:
            self.nt = self.data['nt'][...]

        if not self.ref:
            self.ref = int(self.nt/2)

    def _init_coordinates(self):
        """ Init coordinate system. """

        if self.curvilinear:
            self.x = self.data['xp'][...]
            self.z = self.data['zp'][...]
        else:
            self.x = self.data['x'][...]
            self.z = self.data['z'][...]

    def _init_generator(self):
        """ Init FrameGenerator. """

        try:
            self.frames = FrameGenerator(self.data, self.view, ref=self.ref, nt=self.nt)
            self.maxp, self.minp = self.frames.reference()
        except KeyError:
            self.frames = FrameGenerator(self.data, self.view, nt=self.nt)
            self.maxp, self.minp = self.frames.reference()

    def _init_cmap(self):
        """ Init cmap. """

        self.mycm = modified_jet()
        if self.minp < 0:
            self.norm = MidpointNormalize(vmax=self.maxp, vmin=self.minp, midpoint=0)
        else:
            self.norm = None

    def _init_pbar(self):
        """ Init progress bar. """

        if not self.quiet:
            widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
            self.pbar = ProgressBar(widgets=widgets, maxval=self.nt).start()

    def _init_movie(self):
        """ Init movie parameters. """

        title = os.path.basename(self.filename).split('.')[0]
        metadata = dict(title=title, artist=getpass.getuser(), comment='From nsfds2')
        self.writer = _ani.FFMpegWriter(fps=24, metadata=metadata, bitrate=-1, codec="libx264")
        self.movie_filename = '{}.mkv'.format(title)

    def make(self):
        """ Make movie. """

        # 1st frame
        _, p = next(self.frames)
        movie = _plt.figure(figsize=(20, 9))
        movie.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15)
        axm = movie.add_subplot(111)
        # Labels
        title = r'{} -- iteration : {}'
        axm.set_xlabel(r'$x$ [m]', fontsize=22)
        axm.set_ylabel(r'$y$ [m]', fontsize=22)
        axm.set_title(title.format(0, self.view))
        axm.set_aspect('equal')
        # plot
        if self.curvilinear:
            movie_plt = axm.pcolormesh(self.x, self.z, p.T, cmap=self.mycm, norm=self.norm)
            axm.plot(self.x[0, :], self.z[0, :], 'k', linewidth=3)
            axm.plot(self.x[-1, :], self.z[-1, :], 'k', linewidth=3)
            axm.plot(self.x[:, 0], self.z[:, 0], 'k', linewidth=3)
            axm.plot(self.x[:, -1], self.z[:, -1], 'k', linewidth=3)
            _plt_subdomains(axm, self.x, self.z, self.obstacles, edgecolor='k', curvilinear=True)
        else:
            movie_plt = axm.pcolorfast(self.x, self.z, p, cmap=self.mycm, norm=self.norm)
            _plt_subdomains(axm, self.x, self.z, self.obstacles)

        _plt.colorbar(movie_plt)

        # Start Video
        with self.writer.saving(movie, self.path + self.movie_filename, dpi=100):
            for i, var in self.frames:
                axm.set_title(r'{} -- iteration : {}'.format(self.view, i))
                if self.curvilinear:
                    # StackOv : using-set-array-with-pyplot-pcolormesh-ruins-figure
                    movie_plt.set_array(var[:-1, :-1].T.ravel())
                else:
                    movie_plt.set_data(var)
                self.writer.grab_frame()
                if not self.quiet:
                    self.pbar.update(i)
            if not self.quiet:
                self.pbar.finish()
        _plt.show()
