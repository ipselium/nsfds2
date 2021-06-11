=============
Classical use
=============

Basics
======

**nsfds2** can be used from a terminal with::

   nsfds2 solve|make|show

* *solve* : solves Navier-Stokes equation using default config file *~/.nsfds2/nsfds2.conf*
* *make* : makes movie or sound files from results obtained with *solve* subcommand
* *show* : set of commands for simulations parameters and grid and results inspection

See `-h` option for further help::

   nsfds2 solve -h
   nsfds2 make -h 		# 'movie' and 'sound' subcommands
   nsfds2 show -h 		# 'parameters', 'grid', 'pgrid', 'frame', 'probes' subcommands

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

The `make` subcommand can either generate movie (mp4 file) or sound (wav file)
from an hdf5 file::

   nsfds2 make sound
   nsfds2 make movie

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
| vort | vorticity                   |
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
   timings = True|False 			# Display timings
   quiet = True|False 			# Quiet mode

   [simulation]
   nt = 500 				# Number of time iterations
   ns = 10 				# Save each ns iterations
   cfl = 0.5 				# Courant–Friedrichs–Lewy number

   [geometry]
   mesh = regular|curvilinear|adaptative	# Mesh type
   file = None|path 			# Path to python file (for geoname and curvname)
   geoname = helmholtz_double 		# Python function for geometry
   curvname = curvz			# Python function for curvilinear coordinates
   bc = ARRA 				# Boundary conditions. Must be a mix of ARP
   nx = 256				# Number of grid points along x-axis
   nz = 256				# Number of grid points along z-axis
   ix0 = 128 				# Origin of the grid
   iz0 = 0					# Origin of the grid
   dx = 1e-4				# Spatial x-step
   dz = 1e-4 				# Spatial z-step

   [PML]
   beta = 0.0 				# Depends on pseudo mean flow profile
   alpha = 4.0 				# Order of the spatial repartition law
   sigmax = 20|auto 			# Filter strength along x. Can be 'auto'
   sigmaz = 20|auto 			# Filter strength along z. Can be 'auto'
   npml = 15				# Number of points of the PML

   [source]
   type = None|pulse|harmonic|white|wav 	# Source type
   ixs = 64				# Source x-location
   izs = 128 				# Source z-location
   s0 = 1e6 				# Sources strength [Pa]
   B0 = 2 					# Half spatial bandwidth
   f0 = 60000 				# Frequency (for harmonic only) [Hz]
   wavfile = None|path 			# path to wavfile (for wav only)

   [flow]
   type = None 				# Flow type
   U0 = 5 					# Flow velocity following x [m/s]
   V0 = 5 					# Flow velocity following z [m/s]

   [eulerian fluxes]
   stencil = 3|7|11 			# Number of points of stencil

   [filtering]
   filter = True|False 			# Activate selective filter
   stencil = 11 				# Number of points of stencil used by filters
   stength = 0.75 				# Strength of the filter

   [viscous fluxes]
   viscosity = True|False 			# Activate viscosity
   stencil = 7 				# Number of points of stencil used for viscosity

   [shock capture]
   shock capture = True|False 		# Activate shock capture procedure
   stencil = 7 				# Number of points of stencil for capture
   method = pressure|dilatation 		# Capture based on pressure or dilatation

   [figures]
   figures = True|False 			# Activate figures
   probes = True|False 			# Show probes in maps
   pml = True|False 			# Show PML in maps

   [save]
   save = True|False 			# Activate save
   path = results 				# Path to data file
   filename = tmp 				# Data filename
   compression = None|lzf 			# Activate compression
   probes = [[128, 128], [128, 192]] 	# Probe locations. Must be a list of lists


Customize geometry
==================

To customize geometry, one can provide a set of custom obstacles to the `Mesh`
constructor. To learn more about this, see `fdgrid documentation
<http://perso.univ-lemans.fr/~cdesjouy/fdgrid>`_.

Note on Wav sources
===================

**Important:** When using wav source, pay attention to the spatial steps (*dx*,
*dz*). To resolve frequencies until 20 kHz, *dx* and *dz* must be < 0.017 m.
