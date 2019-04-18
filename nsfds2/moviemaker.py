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
# Creation Date : 2015-02-11 - 01:05:29
"""
-----------

Movie Maker for nsfds2

To make a movie from nsfds2 results :

    python moviemaker.py Mydir/myfile.hdf5

@author: Cyril Desjouy
"""

import os
import sys
import getpass
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from progressbar import ProgressBar, Bar, ReverseBar, ETA
from fdgrid.mesh import Mesh
from mpltools.custom_cmap import MidpointNormalize, modified_jet
from nsfds2.init import CfgSetup
from nsfds2.utils import FrameGenerator


def get_arg():
    """ Get filename. """
    try:
        filename = sys.argv[1]
        if not os.path.exists(filename):
            sys.stderr.write('You must provide a valid hdf5 file : moviemaker mydir/myfile.hdf5/\n')
            sys.exit(1)
    except TypeError:
        sys.stderr.write('You must provide the working directory : moviemaker mydir/myfile.hdf5/\n')
        sys.exit(1)
    else:
        return filename


def get_data(filename):
    """ Get data from filename. """

    try:
        data = h5py.File(filename, 'r')
    except OSError:
        print('You must provide a valid hdf5 file : moviemaker mydir/myfile.hdf5')
        sys.exit(1)
    else:
        return data


def main():
    """ Movie maker main. """

    filename = get_arg()
    path = os.path.dirname(filename) + '/'
    data = get_data(filename)

    cfg = CfgSetup()
    cfg.run()

    # Mesh and time parameters
    x = data['x'][:]
    z = data['z'][:]
    nt = data['nt'][...]
    ns = data['nt'][...]
    obstacles = data['obstacles'][:]

    # Movie Parameters
    title = os.path.basename(filename).split('.')[0]
    metadata = dict(title=title, artist=getpass.getuser(), comment='From nsfds2')
    writer = ani.FFMpegWriter(fps=24, metadata=metadata, bitrate=-1, codec="libx264")
    movie_filename = '{}.mkv'.format(title)
    try:
        frames = FrameGenerator(data, cfg.view, iref=int(nt/2))
        maxp, minp = frames.reference()
    except KeyError:
        frames = FrameGenerator(data, cfg.view)
        maxp, minp = frames.reference()

    # CMAP
    mycm = modified_jet()
    if minp < 0:
        norm = MidpointNormalize(vmax=maxp, vmin=minp, midpoint=0)
    else:
        norm = None

    # Progress bar
    widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
    pbar = ProgressBar(widgets=widgets, maxval=nt).start()

    # Init figure with 1st image
    _, p = next(frames)
    movie = plt.figure(figsize=(20, 9))
    movie.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15)
    axm = movie.add_subplot(111)
    # Labels
    title = r'{} -- iteration : {}'
    axm.set_xlabel(r'$x$ [m]', fontsize=22)
    axm.set_ylabel(r'$y$ [m]', fontsize=22)
    axm.set_title(title.format(0, cfg.view))
    axm.set_aspect('equal')
    # plot
    movie_plt = axm.pcolorfast(x, z, p, cmap=mycm, norm=norm)
    cbar = plt.colorbar(movie_plt)
    Mesh.plot_obstacles(x, z, axm, obstacles)

    # Start Video
    with writer.saving(movie, path + movie_filename, dpi=100):
        for i, var in frames:
            axm.set_title(r'{} -- iteration : {}'.format(cfg.view, i))
            movie_plt.set_data(var)
            Mesh.plot_obstacles(x, z, axm, obstacles)
            writer.grab_frame()
            pbar.update(i)
        pbar.finish()
    plt.show()


if __name__ == "__main__":
    main()
