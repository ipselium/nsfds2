# nsfds2 : 2D Navier-Stokes Finite Differences Solver

## Introducing ***nsfds2***

***nsfds2*** is 2D Navier-Stokes Solver that uses finite difference method. In particular, ***nsfds2*** is specialized in acoustic simulations.

***nsfds2*** is still in developpement. It is still full of bugs and comes with ***ABSOLUTELY NO WARRANTY***.


## Dependencies

* python > 3.6
* numpy
* matplotlib
* h5py
* progressbar33
* ofdlib2 >= 0.9.3
* fdgrid >= 0.8.0
* mplutils >= 0.3.0

***Important:*** To create animations using `nsfds2 make movie`, you also need to have ***ffmpeg*** installed on your system.


## Installation

```
python setup.py install
```

or

```
pip install nsfds2
```


**Note:** To compile *ofdlib2*, OS X users may recquire :

```
xcode-select --install
```

## Classical use

***nsfds2*** can be used from a terminal with :

```
nsfds2 solve|make|show
```

* *solve* : solves Navier-Stokes equation using default config file *~/.nsfds2/nsfds.conf*
* *make* : makes movie or sound files from results obtained with *solve* subcommand
* *show* : set of commands for simulations parameters and grid inspection

See `-h` option for further help :

```
nsfds2 solve -h
nsfds2 make -h 		# 'movie' and 'sound' subcommands
nsfds2 show -h 		# 'parameters', 'grid', 'pgrid', 'frame', 'probes' subcommands
```


## Custom Use

***nsfds2*** can also be used as a classical Python package. It provides
several main objects to perform numerical simulations :

* `init` package povides :

	* `CfgSetup` class : Parses all simulation parameters from
		*~/.nsfds2/nsfds2.conf* file.
	* `Fields` class : Initialize all simulation parameters and fields

* `fdtd` module provides `FDTD` class : Execute simulation

* `utils` package provides in paricular `graphics` module to make result post-treatment


 The following example gives the general philosophy to use ***nsfds2*** :


```python
import matplotlib.pyplot as plt
from nsfds2.init import CfgSetup, Fields
from nsfds2.fdtd import FDTD
from fdgrid.mesh import Mesh

# Read simulation parameter in config file ~/nsfds2/nsfds2.conf
cfg = CfgSetup()

# Define the mesh
msh = Mesh((cfg.nx, cfg.nz), (cfg.dx, cfg.dz), origin=(cfg.ix0, cfg.iz0), obstacles=[])

# Init acoustic fields
fld = Fields(msh, cfg)

# Create simulation
fdtd = FDTD(msh, fld, cfg)
fdtd.run()

# Figures
plt.figure()
plt.imshow(fld.p)
plt.show()
```

## Customize geometry

To customize geometry, one can provide a set of custom obstacles to the `Mesh`
constructor. To learn more about this, see [`fdgrid` documentation](https://github.com/ipselium/fdgrid).


## Wav sources

**Important:** When using wav source, pay attention to the spatial steps (*dx*,
*dz*). To resolve frequencies until 20 kHz, *dx* and *dz* must be < 0.017 m.

## Config file

By default, ***nsfds2*** create the config file `~/.nfds2/nsfds2.conf`. This
file contains simulation parameters that are used by the solver.


```python
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
```

## Examples

The *docs* folder gathers some configuration examples.
Copy one of these files to *.nsfds2/nsfds2.conf* and :

```
nsfds2 solve
```

## Known Bugs

* Mean flows are not yet fully supported.
* Curvilinear geometries are still experimental : PML, viscous flux, and moving boundaries are
  not properly calculated. Be careful with the validity of the results !
