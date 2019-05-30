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
# pylint: disable=no-member
# pylint: disable=invalid-unary-operand-type
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
import fdgrid.graphics as _graphics


__all__ = ['get_data', 'show', 'fields', 'probes', 'DataGenerator', 'Movie']


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


def initial_fields(cfg, msh, fld):
    """ Show initial fields. """

    _, axes = _plt.subplots(2, 2, figsize=(12, 9))

    cm = modified_jet()

    im1 = axes[0, 0].pcolormesh(msh.x, msh.z, fld.p.T, cmap=cm)
    im2 = axes[0, 1].pcolormesh(msh.x, msh.z, fld.re.T, cmap=cm)
    im3 = axes[1, 0].pcolormesh(msh.x, msh.z, fld.ru.T, cmap=cm)
    im4 = axes[1, 1].pcolormesh(msh.x, msh.z, fld.rv.T, cmap=cm)
    ims = [im1, im2, im3, im4]

    axes[0, 0].set_title(r'$p_a$ [Pa]')
    axes[0, 1].set_title(r'$\rho e$ [kg.m$^2$.s$^{-2}$]')
    axes[1, 0].set_title(r'$\rho v_x$ [m/s]')
    axes[1, 1].set_title(r'$\rho v_z$ [m/s]')

    for ax, im in zip(axes.ravel(), ims):
        _graphics.plot_subdomains(ax, msh.x, msh.z, msh.obstacles)
        if cfg.show_pml:
            _graphics.plot_pml(ax, msh.x, msh.z, msh.bc, msh.Npml)

        ax.set_xlabel(r'$x$ [m]')
        ax.set_ylabel(r'$z$ [m]')
        ax.set_aspect('equal')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        _plt.colorbar(im, cax=cax)

    _plt.tight_layout()


def fields(cfg):
    """ Make figure """

    if not cfg.figures or not cfg.save:
        return None

    data = get_data(cfg.datafile)
    var = DataGenerator(data)

    obstacles = data.attrs['obstacles']
    nt = data.attrs['nt']

    Npml = data.attrs['Npml']
    bc = data.attrs['bc']
    x = data['x'][...]
    z = data['z'][...]

    p = var.get(view='p', iteration=nt)
    pmin, pmax = var.reference(view='p')
    if not cfg.onlyp:
        u = var.get(view='vx', iteration=nt)
        umin, umax = var.reference(view='vx')
        v = var.get(view='vz', iteration=nt)
        vmin, vmax = var.reference(view='vz')
        e = var.get(view='e', iteration=nt)
        emin, emax = var.reference(view='e')

    cm = modified_jet()
    unorm = MidpointNormalize(vmin=umin, vmax=umax, midpoint=0)
    vnorm = MidpointNormalize(vmin=vmin, vmax=vmax, midpoint=0)
    pnorm = MidpointNormalize(vmin=pmin, vmax=pmax, midpoint=0)
    enorm = MidpointNormalize(vmin=emin, vmax=emax, midpoint=0)

    if cfg.onlyp:
        _, ax = _plt.subplots(figsize=(12, 9))
        im = ax.pcolormesh(x, z, p, cmap=cm, norm=pnorm)
        _graphics.plot_subdomains(ax, x, z, obstacles)
        ax.set_xlabel(r'$x$ [m]')
        ax.set_ylabel(r'$z$ [m]')
        ax.set_title('Pressure')
        ax.set_aspect('equal')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        _plt.colorbar(im, cax=cax)
        if cfg.probes and cfg.probes_loc:
            ax.plot(*[[x[i], z[j]] for i, j in cfg.probes_loc], 'ro')

    else:
        _, axes = _plt.subplots(2, 2, figsize=(12, 9))

        im1 = axes[0, 0].pcolormesh(x, z, p, cmap=cm, norm=pnorm)
        im2 = axes[0, 1].pcolormesh(x, z, e, cmap=cm, norm=enorm)
        im3 = axes[1, 0].pcolormesh(x, z, u, cmap=cm, norm=unorm)
        im4 = axes[1, 1].pcolormesh(x, z, v, cmap=cm, norm=vnorm)
        ims = [im1, im2, im3, im4]

        axes[0, 0].set_title(r'$p_a$ [Pa]')
        axes[0, 1].set_title(r'$e$ [kg.m$^2$.s$^{-2}$]')
        axes[1, 0].set_title(r'$v_x$ [m/s]')
        axes[1, 1].set_title(r'$v_z$ [m/s]')

        for ax, im in zip(axes.ravel(), ims):
            _graphics.plot_subdomains(ax, x, z, obstacles)
            if cfg.show_pml:
                _graphics.plot_pml(ax, x, z, bc, Npml)

            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            ax.set_aspect('equal')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            _plt.colorbar(im, cax=cax)
            if cfg.probes and cfg.probes_loc:
                ax.plot(*[[x[i], z[j]] for i, j in cfg.probes_loc], 'ro')

        _plt.tight_layout()

    return None


def probes(cfg):
    """ Plot probes. """

    if cfg.figures and cfg.probes:

        with h5py.File(cfg.datafile, 'r') as sfile:

            nt = sfile.attrs['nt']
            dt = sfile.attrs['dt']
            p0 = sfile.attrs['p0']
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


class DataGenerator:
    """ Data Generator

    Parameters
    ----------

    data : Path to hdf5 file (str) or data from hdf5 file
    view : The variable to display. string.
    ref : The reference frame for colormap. int.
    nt : The last frame to consider. int

    """

    def __init__(self, data, view='p', ref=None, nt=None):

        if isinstance(data, str):
            self.data = get_data(data)
        else:
            self.data = data
        self.view = view
        self.var = {'p':'p', 'rho':'rho', 'vx':'rhou', 'vz':'rhov', 'e':'rhoe'}
        self.ref = ref
        self.nt = nt
        self.ns = self.data.attrs['ns']
        self.icur = -self.ns
        self.nx = self.data.attrs['nx']
        self.nz = self.data.attrs['nz']
        if self.data.attrs['mesh'] == 'curvilinear':
            self.J = self.data['J'][...]
        else:
            self.J = _np.ones((self.nx, self.nz))

    def reference(self, view=None, ref=None):
        """ Generate the reference for min/max colormap values """

        if not view:
            view = self.view

        if ref is None:
            ref = self.autoref(view=view)

        if view == "p" and isinstance(ref, int):
            var = self.compute_p(ref)*self.J
            return var.min(), var.max()

        if view == "p" and isinstance(ref, tuple):
            varmin = self.compute_p(ref[0])*self.J
            varmax = self.compute_p(ref[1])*self.J
            return varmin.min(), varmax.max()

        if view in ['rho', 'vx', 'vz', 'e'] and isinstance(ref, int):
            var = self.data[f"{self.var[view]}_it{ref}"][...]*self.J
            if view != 'rho':
                rho = self.data[f"rho_it{ref}"][...]*self.J
                return (var/rho).min(), (var/rho).max()
            if view == 'rho':
                return var.min(), var.max()

        if view in ['rho', 'vx', 'vz', 'e'] and isinstance(ref, tuple):
            varmin = self.data[f"{self.var[view]}_it{ref[0]}"][...]*self.J
            varmax = self.data[f"{self.var[view]}_it{ref[1]}"][...]*self.J
            if view != 'rho':
                rhomin = self.data[f"rho_it{ref[0]}"][...]*self.J
                rhomax = self.data[f"rho_it{ref[1]}"][...]*self.J
                return (varmin/rhomin).min(), (varmax/rhomax).max()
            if view == 'rho':
                return varmin.min(), varmax.max()

        print("Only 'p', 'rho', 'vx', 'vz' and 'e' available !")
        sys.exit(1)

    def autoref(self, view='p'):
        """ Autoset reference. """

        var = DataGenerator(self.data, view=view)

        maxs = []
        mins = []
        for _, v in var:
            maxs.append(v.max())
            mins.append(v.min())

        maxs = _np.array(maxs)
        mins = _np.array(mins)

        refmax = abs(maxs - maxs.mean()).argmin()*self.ns
        refmin = abs(mins - mins.mean()).argmin()*self.ns

        return refmin, refmax

    def close(self):
        """ Close hdf5 file """

        self.data.close()

    def get(self, view='p', iteration=20):
        """ Get data at iteration. """

        if view == 'p':
            return self.compute_p(iteration).T

        if view == 'rho':
            return (self.data[f"{view}_it{iteration}"][...]*self.J).T

        if view in ['vx', 'vz', 'e']:
            vx = (self.data[f"{self.var[view]}_it{iteration}"][...]*self.J).T
            rho = (self.data[f"rho_it{iteration}"][...]*self.J).T
            return vx/rho

        return None

    def compute_p(self, it):
        """ Compute p from rho, rhou, rhov and rhoe at iteration it. """

        rho = self.data[f"rho_it{it}"][...]*self.J
        rhou = self.data[f"rhou_it{it}"][...]*self.J
        rhov = self.data[f"rhov_it{it}"][...]*self.J
        rhoe = self.data[f"rhoe_it{it}"][...]*self.J
        p = _np.empty_like(rho)
        _fdtd.p(p, rho, rhou, rhov, rhoe, self.data.attrs['gamma'])

        return p - self.data.attrs['p0']

    def _next_item(self):
        """ Generate next value of variable """

        if self.view == 'p':
            return self.compute_p(self.icur).T

        if self.view == 'rho':
            return self.data[f"{self.view}_it{self.icur}"][...].T

        if self.view in ['vx', 'vz', 'e']:
            vx = self.data[f"{self.var[self.view]}_it{self.icur}"][...].T
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
        self.ns = self.data.attrs['ns']
        self.nt = nt
        self.quiet = quiet
        self.obstacles = self.data.attrs['obstacles']

        self.Npml = self.data.attrs['Npml']
        self.mesh = self.data.attrs['mesh']
        self.bc = self.data.attrs['bc']

        self.check()
        self._init_coordinates()
        self._init_generator()
        self._init_cmap()
        self._init_pbar()
        self._init_movie()

    def check(self):
        """ Check parameters. """

        if not self.nt:
            self.nt = self.data.attrs['nt']

    def _init_coordinates(self):
        """ Init coordinate system. """

        if self.mesh == 'curvilinear':
            self.x = self.data['xp'][...]
            self.z = self.data['zp'][...]
        else:
            self.x = self.data['x'][...]
            self.z = self.data['z'][...]

    def _init_generator(self):
        """ Init DataGenerator. """

        try:
            self.frames = DataGenerator(self.data, self.view, ref=self.ref, nt=self.nt)
            self.pmin, self.pmax = self.frames.reference(ref=self.ref)
        except KeyError:
            self.frames = DataGenerator(self.data, self.view, nt=self.nt)
            self.pmin, self.pmax = self.frames.reference()

    def _init_cmap(self):
        """ Init cmap. """

        self.mycm = modified_jet()
        self.norm = MidpointNormalize(vmax=self.pmax, vmin=self.pmin, midpoint=0)

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

    def make(self, show_pml=False):
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
        if self.mesh == 'curvilinear':
            movie_plt = axm.pcolormesh(self.x, self.z, p.T, cmap=self.mycm, norm=self.norm)
        elif self.mesh == 'adaptative':
            movie_plt = axm.pcolormesh(self.x, self.z, p, cmap=self.mycm, norm=self.norm)
        else:
            movie_plt = axm.pcolorfast(self.x, self.z, p, cmap=self.mycm, norm=self.norm)

        _graphics.plot_subdomains(axm, self.x, self.z, self.obstacles)
        if show_pml:
            _graphics.plot_pml(axm, self.x, self.z, self.bc, self.Npml)
        _plt.colorbar(movie_plt)

        # Start Video
        with self.writer.saving(movie, self.path + self.movie_filename, dpi=100):
            for i, var in self.frames:
                axm.set_title(title.format(self.view, i))
                if self.mesh == 'curvilinear':
                    # StackOv : using-set-array-with-pyplot-pcolormesh-ruins-figure
                    movie_plt.set_array(var[:-1, :-1].T.ravel())
                elif self.mesh == 'adaptative':
                    movie_plt.set_array(var[:-1, :-1].flatten())
                else:
                    movie_plt.set_data(var)
                self.writer.grab_frame()
                if not self.quiet:
                    self.pbar.update(i)
            if not self.quiet:
                self.pbar.finish()
        _plt.show()
