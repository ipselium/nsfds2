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
# pylint: disable=too-many-instance-attributes
# pylint: disable=no-member
"""
--------------

Graphic utilities for nsfds2

--------------
"""

import os
import sys
import getpass
import h5py
import numpy as _np
from scipy import signal as _signal
import matplotlib.pyplot as _plt
import matplotlib.animation as _ani
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ofdlib2 import fdtd as _fdtd, derivation as _derivation
from progressbar import ProgressBar, Bar, ReverseBar, ETA
from mplutils import modified_jet, MidpointNormalize, set_figsize, get_subplot_shape
import fdgrid.graphics as _graphics
from nsfds2.utils.misc import nearest_index as _ne

__all__ = ['get_data', 'DataIterator', 'DataExtractor', 'Plot']


def get_data(filename):
    """ Get data from `filename` (hdf5 file). """

    try:
        home = os.path.expanduser("~")
        filename = filename.replace('~', home)
        data = h5py.File(filename, 'r')
    except OSError:
        print('You must provide a valid hdf5 file')
        sys.exit(1)
    else:
        return data


class DataIterator:
    """ Data Generator

    Parameters
    ----------
    data :
        DataIterator instance.
    view : tuple
        The variable to display.
    nt : int
        The last frame to consider.

    """

    def __init__(self, data_extractor, view=('p'), nt=None):

        self.data = data_extractor
        self.view = view
        self.ns = self.data.get_attr('ns')
        self.icur = -self.ns
        if nt is None:
            self.nt = self.data.get_attr('nt')
        else:
            self.nt = nt

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

            tmp = [self.icur]
            for var in self.view:
                tmp.append(self.data.get(view=var, iteration=self.icur))

            return tmp

        except KeyError:
            raise StopIteration


class DataExtractor:
    """ Extract data from hdf5 file

    Parameters
    ----------
    data : str, hdf5file
        Path to hdf5 file or data from hdf5 file.

    """

    def __init__(self, data):

        if isinstance(data, str):
            self.data = get_data(data)
        else:
            self.data = data

        self.var = {'p':'p', 'rho':'rho', 'vx':'rhou', 'vz':'rhov', 'e':'rhoe'}
        self.nt = self.get_attr('nt')
        self.ns = self.get_attr('ns')

        if self.get_attr('mesh') == 'curvilinear':
            self.J = self.get_dataset('J')
        else:
            self.J = _np.ones((self.get_attr('nx'), self.get_attr('nz')))

    def __enter__(self):
        return self

    def __exit__(self, mtype, value, traceback):
        self.close()

    def reference(self, view='p', ref=None):
        """ Generate the references  for min/max colormap values """

        if ref is None:
            ref = self.autoref(view=view)

        if isinstance(ref, int):
            var = self.get(view=view, iteration=_ne(ref, self.ns, self.nt))
            return var.min(), var.max()

        if isinstance(ref, tuple):
            varmin = self.get(view=view, iteration=_ne(ref[0], self.ns, self.nt))
            varmax = self.get(view=view, iteration=_ne(ref[1], self.ns, self.nt))
            return varmin.min(), varmax.max()

        print("Only 'p', 'rho', 'vx', 'vz', 'vort' and 'e' available !")
        sys.exit(1)

    def autoref(self, view='p'):
        """ Autoset reference. """

        var = DataIterator(self, view=(view, ))

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

    def list(self):
        """ List all datasets and attributes. """

        datasets = [i for i in self.data.keys() if '_it' not in i]
        print('datasets: ', *datasets)
        print('attrs: ', *list(self.data.attrs))

    def get(self, view='p', iteration=0):
        """ Get data at iteration. """

        if view == 'p':
            rho = self.data[f"rho_it{iteration}"][...]*self.J
            rhou = self.data[f"rhou_it{iteration}"][...]*self.J
            rhov = self.data[f"rhov_it{iteration}"][...]*self.J
            rhoe = self.data[f"rhoe_it{iteration}"][...]*self.J
            p = _np.empty_like(rho)
            _fdtd.p(p, rho, rhou, rhov, rhoe, self.data.attrs['gamma'])

            return p.T - self.data.attrs['p0']

        if view == 'rho':
            return (self.data[f"{view}_it{iteration}"][...]*self.J).T

        if view in ['vx', 'vz', 'e']:
            vx = (self.data[f"{self.var[view]}_it{iteration}"][...]*self.J)
            rho = (self.data[f"rho_it{iteration}"][...]*self.J)
            return (vx/rho).T

        if view == 'vort':
            rho = (self.data[f"rho_it{iteration}"][...]*self.J)
            vx = (self.data[f"rhou_it{iteration}"][...]*self.J)/rho
            vz = (self.data[f"rhov_it{iteration}"][...]*self.J)/rho
            return _derivation.vorticity(self.data['x'][...], self.data['x'][...],
                                         vx, vz).T

        raise ValueError("view must be 'p', 'rho', 'vx', 'vz', 'e', or 'vort'")

    def get_attr(self, attr):
        """ Get attribute from hdf5 file. attr must be string."""
        return self.data.attrs[attr]

    def get_dataset(self, dataset):
        """ Get dataset from hdf5 file. attr must be string."""

        return self.data[dataset][...]


class Plot:
    """ Helper class to plot results from nsfds2.

    Parameters
    ----------
    filename : str
        hdf5 file
    quiet : bool, optional
        Quiet mode.

    """

    def __init__(self, filename, quiet=False):

        self.filename = filename
        self.path = os.path.dirname(filename) + '/'
        self.quiet = quiet

        self.data = DataExtractor(self.filename)
        self.nt = self.data.get_attr('nt')
        self.ns = self.data.get_attr('ns')

        self._init_geo()
        self._init_fig()

    def _init_geo(self):
        """ Init coordinate system. """

        self.obstacles = self.data.get_attr('obstacles')
        self.Npml = self.data.get_attr('Npml')
        self.mesh = self.data.get_attr('mesh')
        self.bc = self.data.get_attr('bc')

        if self.mesh == 'curvilinear':
            self.x = self.data.get_dataset('xp')
            self.z = self.data.get_dataset('zp')
        else:
            self.x, self.z = _np.meshgrid(self.data.get_dataset('x'), self.data.get_dataset('z'))
            self.x = _np.ascontiguousarray(self.x.T)
            self.z = _np.ascontiguousarray(self.z.T)

    def _init_fig(self):
        """ Init figure parameters. """

        self.cm = modified_jet()
        self.title = r'{} -- iteration : {}'
        self.titles = {'p': r'$p_a$ [Pa]',
                       'e': r'$e$ [kg.m$^2$.s$^{-2}$]',
                       'vx': r'$v_x$ [m/s]',
                       'vz': r'$v_z$ [m/s]',
                       'rho': r'$\rho$ [kg.m$^3$]',
                       'vort': r'$\omega$ [m/s]'}

    def movie(self, view=('p', 'e', 'vx', 'vz'), nt=None, ref=None,
              figsize='auto', show_pml=False, show_probes=False, dpi=100):
        """ Make movie. """

        # Progress bar
        if not self.quiet:
            widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
            pbar = ProgressBar(widgets=widgets, maxval=self.nt).start()

        # Movie parameters
        title = os.path.basename(self.filename).split('.')[0]
        metadata = dict(title=title, artist=getpass.getuser(), comment='From nsfds2')
        writer = _ani.FFMpegWriter(fps=24, metadata=metadata, bitrate=-1, codec="libx264")
        movie_filename = f'{title}.mkv'

        # Nb of iterations
        if nt is None:
            nt = self.nt
        else:
            nt = _ne(nt, self.ns, self.nt)

        # Create Iterator and make 1st frame
        data = DataIterator(self.data, view=view, nt=nt)
        i, *var = next(data)
        fig, axes, ims = self.fields(view=view, iteration=i, ref=ref,
                                     show_pml=show_pml,
                                     show_probes=show_probes,
                                     figsize=figsize)

        with writer.saving(fig, self.path + movie_filename, dpi=dpi):

            writer.grab_frame()

            for i, *var in data:

                # StackOv : using-set-array-with-pyplot-pcolormesh-ruins-figure
                for ax, mesh, v, j in zip(axes.ravel(), ims, var, range(len(ims))):
                    mesh.set_array(v[:-1, :-1].T.flatten())
                    ax.set_title(self.titles[view[j]] + f' (n={i})')

                writer.grab_frame()

                if not self.quiet:
                    pbar.update(i)

            if not self.quiet:
                pbar.finish()

    def probes(self):
        """ Plot pressure at probes. """

        probes = self.data.get_dataset('probes_location').tolist()

        if not probes:
            return None

        p = self.data.get_dataset('probes_value')
        t = _np.arange(self.nt)*self.data.get_attr('dt')

        _, ax = _plt.subplots(figsize=(9, 4))
        for i, c in enumerate(probes):
            if self.data.get_attr('mesh') == 'curvilinear':
                p0 = self.data.get_attr('p0')/self.data.get_dataset('J')[c[0], c[1]]
            else:
                p0 = self.data.get_attr('p0')
            ax.plot(t, p[i, :] - p0, label=f'@{tuple(c)}')
        ax.set_xlim(t.min(), t.max())
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Pressure [Pa]')
        ax.legend()
        ax.grid()

        return None

    def spectrogram(self):
        """ Plot spectograms at probes. """

        probes = self.data.get_dataset('probes_location').tolist()

        if not probes:
            return None

        p = self.data.get_dataset('probes_value')

        M = 1024

        fig, ax = _plt.subplots(p.shape[0], figsize=(9, 4))
        for i, c in enumerate(probes):

            if self.data.get_attr('mesh') == 'curvilinear':
                p0 = self.data.get_attr('p0')/self.data.get_dataset('J')[c[0], c[1]]
            else:
                p0 = self.data.get_attr('p0')

            freqs, times, Sx = _signal.spectrogram(p[i, :] - p0,
                                                   fs=1/self.data.get_attr('dt'),
                                                   window='hanning',
                                                   nperseg=M, noverlap=M-100,
                                                   detrend=False, scaling='spectrum')
            im = ax[i].pcolormesh(times, freqs/1000, 10*_np.log10(Sx), cmap='viridis')
            ax[i].set_ylabel('Frequency [kHz]')
            if i != len(probes) - 1:
                ax[i].set_xticks([])

            fig.colorbar(im, ax=ax[i], label=f'probe {i}')

        ax[-1].set_xlabel('Time [s]')
        ax[0].set_title('Square spectrum magitude')
        _plt.tight_layout()

        return None

    def fields(self, view=('p', 'e', 'vx', 'vz'), iteration=None, ref=None,
               show_pml=False, show_probes=True, figsize='auto'):
        """ Make figure """

        if iteration is None:
            iteration = self.nt
        else:
            iteration = _ne(iteration, self.ns, self.nt)

        var = []
        norm = []
        ims = []

        for v in view:
            var.append(self.data.get(view=v, iteration=iteration).T)
            vmin, vmax = self.data.reference(view=v, ref=ref)
            norm.append(MidpointNormalize(vmin=vmin, vmax=vmax, midpoint=0))

        fig, axes = _plt.subplots(*get_subplot_shape(len(var)))

        if not isinstance(axes, _np.ndarray):   # if only 1 varible in view
            axes = _np.array(axes)

        for i, ax in enumerate(axes.ravel()):
            if i < len(var):
                ims.append(ax.pcolormesh(self.x, self.z, var[i],
                                         cmap=self.cm, norm=norm[i]))
                ax.set_title(self.titles[view[i]] + f' (n={iteration})')
                ax.set_xlabel(r'$x$ [m]')
                ax.set_ylabel(r'$z$ [m]')
                ax.set_aspect('equal')
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                _plt.colorbar(ims[i], cax=cax)

                probes = self.data.get_dataset('probes_location').tolist()
                if probes and show_probes:
                    _ = [ax.plot(self.x[i, j], self.z[i, j], 'ro') for i, j in probes]

                _graphics.plot_subdomains(ax, self.x, self.z, self.obstacles)
                if show_pml:
                    _graphics.plot_pml(ax, self.x, self.z, self.bc, self.Npml)
            else:
                ax.remove()

        fig.set_size_inches(*set_figsize(axes, figsize))
        _plt.tight_layout()

        return fig, axes, ims

    @staticmethod
    def show():
        """ Show all figures. """
        _plt.show()
