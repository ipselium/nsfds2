=============
Classical use
=============

Basics
======

**nsfds2** can be used from a terminal with::

   nsfds2 solve|make|show|loop

* *solve* : solves Navier-Stokes equation using default config file *~/.nsfds2/nsfds2.conf*
* *make* : makes movie or sound files from results obtained with *solve* subcommand
* *show* : set of commands for simulations parameters and grid and results inspection
* *loop* : Run simulations for a set of config files in a path
* *ploop* : Run simulations for a set of config files in a path in parallel processes

See `-h` option for further help::

   nsfds2 solve -h
   nsfds2 make -h 		# 'movie' and 'sound' subcommands
   nsfds2 show -h 		# 'parameters', 'grid', 'pgrid', 'frame', 'probes' subcommands
   nsfds2 loop -h
   nsfds2 ploop -h

The solve subcommand
====================

The `solve` subcommand launch the solver. By default, the solver uses
`~/.nsfds2/nsfds2.conf` as configuration file. To target another configuration
file, use the `-c` argument as follows::

   nsfds2 solve -c my_confif_file.conf


It is possible to display timings with::

   nsfds2 solve -t

or make the solver quiet with::

   nsfds2 solver -q


To see other options::

   nsfds2 solver -h

The make subcommand
===================

The `make` subcommand can either generate movie (mp4 file), sound (wav file)
from an hdf5 file or a config file template::

   nsfds2 make sound
   nsfds2 make movie
   nsfds2 make template

Work with sounds
----------------

Sounds are generated from the probes (pressure) saved in an hdf5 file. If no
probe has been set in the computation domain, no sound will be generated from
`nsfds2 make sound`.  If this command is called without any parameters, the
default configuration file is used. It is also possible to call the command
with a target configuration file::

   nsfds2 make sound -c my_confif_file.conf

or with target data file::

   nsfds2 make sound -d my_data_file.hdf5


Work with movies
----------------

Movies can be created from an hdf5 file if the `save` option has been selected.
Then, the following variables are allowed as argument of the movie subcommand:

+------+-----------------------------+
| var  | variable                    |
+======+=============================+
| p    | pressure                    |
+------+-----------------------------+
| rho  | density                     |
+------+-----------------------------+
| vx   | x component of the velocity |
+------+-----------------------------+
| vz   | z component of the velocity |
+------+-----------------------------+
| e    | energy                      |
+------+-----------------------------+
| vxz  | vorticity                   |
+------+-----------------------------+

For instance, to create a movie with both components of velocity::

   nsfds2 make movie vx vz

It is also possible to make movie from a target configuration file::

   nsfds2 make movie e -c my_confif_file.conf

or from a data file::

   nsfds2 make movie rho -d my_data_file.hdf5


The show subcommand
===================

The `show` subcommand provides a set of commands to show simulation parameters,
results, or grid configuration. The main `show` subcommands are:

:frame:  Display an acoustic variable at a given iteration (works like `movie`)
:probes: Display pressure field(s) as a function of time at probe location(s)
:spectrogram: Display spectrogram(s) at probe location(s)
:pgrid: Display physical grid
:parameters: Display some simulation parameters

For instance, to display the density at iteration 100::

   nsfds2 show frame rho -i 100


Config file
===========

By default, **nsfds2** create the config file `~/.nfds2/nsfds2.conf`. This
file contains simulation parameters that are used by the solver.

::

   [configuration]
   version = 0.13.0
   timings = False      # Display timing detail
   quiet = False        # Quiet mode
   cpu = 1              # Number of cpu used by the solver

   [simulation]
   nt = 500             # Number of time iterations
   ns = 10              # Save fields each ns iterations
   cfl = 0.5             # Courant–Friedrichs–Lewy number

   [thermophysic]
   norm = True|False    # Normalize p0, rho0, c0 and T0 (Override other values).
   p0 = 101325.0        # Atmospheric pressure (Pa)
   t0 = 20.0            # Ambiant temperature (°C)
   gamma = 1.4          # Heat capacity ratio
   prandtl = 0.7        # Prandtl number

   [geometry]
   mesh = regular|curvilinear|adaptative	 # Mesh type
   file = None|path             	 # Path to .py file (for geo/curvname)
   geoname = helmholtz_double 		# Python function for geometry
   curvname = curvz			# Python function for curvilinear coordinates
   only_pml = False                     # Adaptative only in PML
   Nd = 23                              # Adaptative over Nd points
   Rx = 3.                              # Dilatation rate [adaptative mesh]
   bc = PPPP                            # Boundary conditions. Must be a mix of AWP
   nx = 256                             # Number of grid points along x-axis
   nz = 256                             # Number of grid points along z-axis
   ix0 = 0                              # Origin of the grid
   iz0 = 0                              # Origin of the grid
   dx = 1                               # Spatial x-step
   dz = 1                               # Spatial z-step

   [PML]
   beta = 0.0 				# Depends on pseudo mean flow profile
   alpha = 4.0 				# Order of the spatial repartition law
   sigmax = 20|auto 			# Filter strength along x. Can be 'auto'
   sigmaz = 20|auto 			# Filter strength along z. Can be 'auto'
   npml = 15				# Number of points of the PML

   [source]
   type = None|pulse|harmonic|white|wav     # Source type
   ixs = 64                                 # Source x-location
   izs = 128                                # Source z-location
   s0 = 1e6                                 # Sources strength [Pa]
   b0 = 2                                   # Half spatial bandwidth
   f0 = 60000                               # Frequency [harmonic only] [Hz]
   seed = None                              # Seed [white noise only]. Must be integer.
   off = 100                                 # Stop source at iteration 100. nt by default.
   wavfile = None|path                       # Path to wavfile (for wav only)

   [flow]
   type = None                  # Flow type [custom/vortex]
   U0 = 5                       # Flow velocity following x [m/s] [only for custom]
   V0 = 5                       # Flow velocity following z [m/s] [only for custom]

   [eulerian fluxes]
   stencil = 3|7|11             # Number of points of stencil

   [filtering]
   filter = True|False           # Activate selective filter
   stencil = 11                 # Number of points of stencil used by filters
   strength = 0.75              # Strength of the filter
   strength_on_walls = 0.01     # Strength on the nearest point from a wall

   [viscous fluxes]
   viscosity = True|False       # Activate viscosity
   stencil = 7                  # Number of points of stencil used for viscosity

   [shock capture]
   shock capture = True|False   # Activate shock capture procedure
   stencil = 7                  # Number of points of stencil for capture
   method = pressure|dilatation # Capture based on pressure or dilatation

   [figures]
   figures = True|False          # Activate figures
   probes = True|False          # Show probes in maps
   pml = True|False             # Show PML in maps
   bc_profiles = True            # Show bc profiles
   fps = 24                     # Framerate for movies

   [save]
   resume = True|False          # Resume older simulation
   path = results/              # path to data file
   filename = tmp                # data filename
   compression = None|lzf       # Activate data compression
   fields = True                 # Save fields
   vorticity = False            # Save vorticity
   probes = []                  # Probe locations. Must be list of lists


Customize geometry
==================

To customize geometry, one can provide a set of custom obstacles to the `Mesh`
constructor. To learn more about this, see `fdgrid documentation
<http://perso.univ-lemans.fr/~cdesjouy/fdgrid>`_.

Note on Wav sources
===================

**Important:** When using wav source, pay attention to the spatial steps (*dx*,
*dz*). To resolve frequencies until 20 kHz, *dx* and *dz* must be < 0.017 m.
