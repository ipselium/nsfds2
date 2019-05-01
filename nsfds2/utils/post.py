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
# Creation Date : 2019-03-07 - 21:31:14
#
# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-locals
#
"""
-----------

Post-treatment tools for nsfds2

@author: Cyril Desjouy
"""

import os
import sys
import getpass
import h5py
import numpy as np
import numba as nb
from ofdlib2.coefficients import a7o
from ofdlib2 import fdtd
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from progressbar import ProgressBar, Bar, ReverseBar, ETA
from fdgrid.mesh import plot_subdomains
from mpltools.custom_cmap import MidpointNormalize, modified_jet


@nb.jit
def rot(data, ref, nx, nz, one_dx, one_dz, a7):
    """ Rotational. """

    name = "{}_it{}"
    rho = data[name.format('rho', ref)][:, :]
    rhou = data[name.format('rhou', ref)][:, :]
    rhov = data[name.format('rhov', ref)][:, :]

    u = rhou/rho
    v = rhov/rho

    vort = np.zeros_like(u)

    for i in range(3, nx-3):
        for j in range(3, nz-3):
            for l in range(-3, 4):
                vort[i, j] = vort[i, j] + a7[l]*v[i+l, j]*one_dx - a7[l]*u[i, j+l]*one_dz

    for i in range(3):
        for j in range(3, nz-3):
            vort[i, j] = vort[i, j] + (v[i+1, j] - v[i, j])*one_dx
            for l in range(-3, 4):
                vort[i, j] = vort[i, j] - a7[l]*u[i, j+l]*one_dz

    for i in range(nx-3, nx):
        for j in range(3, nz-3):
            vort[i, j] = vort[i, j] + (v[i, j] - v[i-1, j])*one_dx
            for l in range(-3, 4):
                vort[i, j] = vort[i, j] - a7[l]*u[i, j+l]*one_dz

    return vort.T


class FrameGenerator:
    """ Frame Genrator """

    def __init__(self, data, view='p', ref=20, nt=None):

        self.data = data
        self.view = view
        self.nt = nt
        self.imin = 0
        self.ref = ref
        self.ns = data['ns'][...]
        self.icur = self.imin
        self.var = {'p': 'p',
                    'rho': 'rho',
                    'vx': 'rhou',
                    'vz': 'rhov',
                    'e': 'rhoe'}
        self.a7, _ = a7o()
        self.one_dx = 1/self.data['dx'][...]
        self.one_dz = 1/self.data['dz'][...]
        self.nx = self.data['nx'][...]
        self.nz = self.data['nz'][...]
        if self.data['mesh'][...] == 'curvilinear':
            self.J = self.data['J'][...]
        else:
            self.J = np.ones((self.nx, self.nz))

    def reference(self):
        """ Generate the reference for min/max colormap values """

        if self.view == "p":
            ref = self.p_from_rhoX(self.ref)
        elif self.view in ['rho', 'vx', 'vz', 'e']:
            ref = self.data["{}_it{}".format(self.var[self.view], self.ref)][:, :]
        elif self.view in ['vort']:
            ref = rot(self.data, self.ref, self.nx, self.nz,
                      self.one_dx, self.one_dz, self.a7)
        else:
            print("Only 'p', 'rho', 'vx', 'vz' and 'e' available !")
            sys.exit(1)

        return (ref*self.J).max(), (ref*self.J).min()

    def p_from_rhoX(self, i):
        """ Compute p from rho, rhou, rhov and rhoe. """

        rho = self.data["{}_it{}".format('rho', i)][:, :]*self.J
        rhou = self.data["{}_it{}".format('rhou', i)][:, :]*self.J
        rhov = self.data["{}_it{}".format('rhov', i)][:, :]*self.J
        rhoe = self.data["{}_it{}".format('rhoe', i)][:, :]*self.J
        p = np.empty_like(rho)
        fdtd.p(p, rho, rhou, rhov, rhoe, self.data['gamma'][...])
        return p - self.data['p0'][...]

    def next_item(self):
        """ Generate next value of variable """

        if self.view == 'p':
            return self.p_from_rhoX(self.icur).T

        if self.view == 'rho':
            return self.data["{}_it{}".format(self.var[self.view], self.icur)][:, :].T

        if self.view in ['vx', 'vz', 'e']:
            vX = self.data["{}_it{}".format(self.var[self.view], self.icur)][:, :]
            rho = self.data["{}_it{}".format('rho', self.icur)][:, :]
            return vX.T/rho.T

        if self.view in ['vort']:
            return rot(self.data, self.icur, self.nx, self.nz,
                       self.one_dx, self.one_dz, self.a7)

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

            return self.icur, self.next_item()

        except KeyError:
            raise StopIteration


def get_data(filename):
    """ Get data from filename. """

    try:
        data = h5py.File(filename, 'r')
    except OSError:
        print('You must provide a valid hdf5 file : moviemaker mydir/myfile.hdf5')
        sys.exit(1)
    else:
        return data


def make_movie(filename, view='p', ref=None, nt=None, quiet=False):
    """ Movie maker main. """

    path = os.path.dirname(filename) + '/'
    data = get_data(filename)

    mesh = data['mesh'][...]

    # Mesh and time parameters
    if mesh == 'curvilinear':
        x = data['xp'][...]
        z = data['zp'][...]
    else:
        x = data['x'][...]
        z = data['z'][...]

    obstacles = data['obstacles'][...]

    if not nt:
        nt = data['nt'][...]
    if not ref:
        ref = int(nt/2)

    # Movie Parameters
    title = os.path.basename(filename).split('.')[0]
    metadata = dict(title=title, artist=getpass.getuser(), comment='From nsfds2')
    writer = ani.FFMpegWriter(fps=24, metadata=metadata, bitrate=-1, codec="libx264")
    movie_filename = '{}.mkv'.format(title)
    try:
        frames = FrameGenerator(data, view, ref=ref, nt=nt)
        maxp, minp = frames.reference()
    except KeyError:
        frames = FrameGenerator(data, view, nt=nt)
        maxp, minp = frames.reference()

    # CMAP
    mycm = modified_jet()
    if minp < 0:
        norm = MidpointNormalize(vmax=maxp, vmin=minp, midpoint=0)
    else:
        norm = None

    # Progress bar
    if not quiet:
        widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
        pbar = ProgressBar(widgets=widgets, maxval=nt).start()

    # 1st frame
    _, p = next(frames)
    movie = plt.figure(figsize=(20, 9))
    movie.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15)
    axm = movie.add_subplot(111)
    # Labels
    title = r'{} -- iteration : {}'
    axm.set_xlabel(r'$x$ [m]', fontsize=22)
    axm.set_ylabel(r'$y$ [m]', fontsize=22)
    axm.set_title(title.format(0, view))
    axm.set_aspect('equal')
    # plot
    if mesh == 'curvilinear':
        movie_plt = axm.pcolormesh(x, z, p.T, cmap=mycm, norm=norm)
        axm.plot(x[0, :], z[0, :], 'k', linewidth=3)
        axm.plot(x[-1, :], z[-1, :], 'k', linewidth=3)
        axm.plot(x[:, 0], z[:, 0], 'k', linewidth=3)
        axm.plot(x[:, -1], z[:, -1], 'k', linewidth=3)
        plot_subdomains(axm, x, z, obstacles, edgecolor='k', curvilinear=True)
    else:
        movie_plt = axm.pcolorfast(x, z, p, cmap=mycm, norm=norm)
        plot_subdomains(axm, x, z, obstacles)

    plt.colorbar(movie_plt)

    # Start Video
    with writer.saving(movie, path + movie_filename, dpi=100):
        for i, var in frames:
            axm.set_title(r'{} -- iteration : {}'.format(view, i))
            if mesh == 'curvilinear':
                # StackOv : using-set-array-with-pyplot-pcolormesh-ruins-figure
                movie_plt.set_array(var[:-1, :-1].T.ravel())
            else:
                movie_plt.set_data(var)
            writer.grab_frame()
            if not quiet:
                pbar.update(i)
        if not quiet:
            pbar.finish()
    plt.show()
