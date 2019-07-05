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
# Creation Date : 2019-02-13 - 10:58:09
"""
-----------

Examples of obstacle arangements.

@author: Cyril Desjouy
"""


import numpy as _np
from fdgrid import Domain, Obstacle


def logo(xn, zn):
    """ Curvlinear coordinates : test case 1. Physical == numerical """

    xp = xn.copy()
    zp = zn \
        + _np.linspace(0.5, 0, zn.shape[1])*(100*xp**2)
    return xp, zp


def circle(xn, zn):
    """ circle """

    dx = xn[1, 0] - xn[0, 0]
    width = xn[-1, 0] - xn[0, 0]
    R = (width + dx)/(2*_np.pi)
    xp = (zn + R)*_np.sin(xn/R)
    zp = (zn + R)*_np.cos(xn/R)

    return xp, zp


def letter_a(nx, nz, size=10):
    """ letter a.

    Parameters:
    ----------

    size : percentage of domain
    """

    size = size/100
    xc, zc = int(nx/2), int(nz/2)
    xmin, xmax = int(xc-size*nx), int(xc+size*nx)
    zmin, zmax = int(zc-2*size*nz), int(zc+2*size*nz)

    geo = [Obstacle([xmin, zc-6, xmax, zc+6], 'WWWW'),
           Obstacle([xmin, zmax-12, xmax, zmax], 'WWWW'),
           Obstacle([xmin-12, zmin, xmin, zmax], 'WWWW'),
           Obstacle([xmax, zmin, xmax+12, zmax], 'WWWW')]

    return Domain((nx, nz), data=geo)


def moving_square(nx, nz, size_percent=20):
    """ Square in the middle.

    Parameters:
    -----------

    size_percent (float): size of the square in percent of the largest
    dimension of the domain.
    """

    size = int(min(nx, nz)*size_percent/100)
    obs1 = Obstacle([int(nx/2)-size, int(nz/2)-size,
                     int(nx/2)+size, int(nz/2)+size], 'VWWWW')
    obs2 = Obstacle([nx-11, 0, nx-1, nz-1], 'VWWWW')

    obs1.set_moving_bc({'f': 70000, 'A': 1, 'func': 'sine'})
    obs2.set_moving_bc({'f': 73000, 'A': 1, 'func': 'sine'})

    return Domain((nx, nz), data=[obs1, obs2])


def ceperley(nx, nz, source_size=5):
    """ Ceperley resonator.
    """

    dx = 0.01

    source_size = int(nx*source_size/100)
    loc1 = int(nx/4)
    loc2 = int(nx/2)

    obs1 = Obstacle([0, 0, nx-1, 11], 'WWWW')
    obs2 = Obstacle([0, nz-11, loc1 , nz-1], 'WWWW')
    obs3 = Obstacle([loc1, nz-11, loc1+source_size, nz-1], 'WVWW')
    obs4 = Obstacle([loc1+source_size, nz-11, loc2, nz-1], 'WWWW')
    obs5 = Obstacle([loc2, nz-11, loc2+source_size, nz-1], 'WVWW')
    obs6 = Obstacle([loc2+source_size, nz-11, nx-1, nz-1], 'WWWW')

    obs3.set_moving_bc({'f': 340/(nx*dx), 'A': 2, 'func': 'sine', 'kwargs':{'n':2}})
    obs5.set_moving_bc({'f': 340/(nx*dx), 'A': 2, 'phi': _np.pi/4,  'func': 'sine', 'kwargs':{'n':2}})

    return Domain((nx, nz), data=[obs1, obs2, obs3, obs4, obs5, obs6])


def wav_squares(nx, nz, size_percent=5):
    """ Square in the middle.

    Parameters:
    -----------

    size_percent (float): size of the square in percent of the largest
    dimension of the domain.
    """

    size = int(min(nx, nz)*size_percent/100)

    obs1 = Obstacle([int(nx/2)-size, int(nz/4),
                     int(nx/2)+size, int(nz/4)+2*size], 'WWVW')
    obs2 = Obstacle([int(nx/2)-size, int(3*nz/4)-2*size,
                     int(nx/2)+size, int(3*nz/4)], 'WWVW')

    obs1.set_moving_bc({'f': '~/.nsfds2/spam3.wav', 'A': 0.05, 'func': 'tukey'})
    obs2.set_moving_bc({'f': '~/.nsfds2/spam3.wav', 'A': -0.05, 'func': 'tukey'})

    return Domain((nx, nz), data=[obs1, obs2])

