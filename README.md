# nsfds2 : 2D Navier-Stokes Finite Differences Solver

## Introducing ***nsfds2***

***nsfds2*** is 2D Navier-Stokes Solver that uses finite difference method.

***nsfds2*** is still in developpement.


## Dependencies

* python > 3.6
* numpy
* matplotlib
* h5py
* progressbar33
* ofdlib2
* fdgrid
* mplutils

## Installation

`python setup.py install`

or

`pip install nsfds2`


**Note:** MAC users may recquire :

```
xcode-select --install
```

## Classical use

**nsfds2** can be used from a terminal with :

```
nsfds2 solve|movie|show
```

* *solve* : solves Navier-Stokes equation using default config file *~/.nsfds2/nsfds.conf*
* *movie* : makes movie from results obtained with *solve* subcommand
* *show* : shows simulations parameters and grid

## Custom Use


```python
import matplotlib.pyplot as plt
from nsfds2.init import CfgSetup, Fields
from nsfds2.fdtd import FDTD
from nsfds2.utils import graphics
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


## Config file

```python
[configuration]
timings = True|False 			# Display timings
quiet = True|False 			# Quiet mode

[simulation]
nt = 500 				# Number of time iterations
ns = 10 				# Save each ns iterations
cfl = 0.5 				# Courant–Friedrichs–Lewy number

[geometry]
mesh = regular|curvilinear		# Mesh type
file = None|path 			# Path to python file (geometry)
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
sigmax = 20 				# Filter strength along x
sigmaz = 20 				# Filter strength along z
npml = 15				# Number of points of the PML

[source]
type = None|pulse|harmonic|white 	# Source type
ixs = 64				# Source x-location
izs = 128 				# Source z-location
s0 = 1e6 				# Sources strength
B0 = 2 					# Half spatial bandwidth
f0 = 60000 				# Frequency (only for harmonic)

[flow]
type = None 				# Flow type
U0 = 5 					# Flow velocity following x
V0 = 5 					# Flow velocity following z

[eulerian fluxes]
stencil = 3|7|11 			# Number of points of stencil

[filtering]
filter = True|False 			# Activate selective filter
stencil = 11 				# Number of points of stencil used by filter
stength = 0.75 				# Strength of the filter

[viscous fluxes]
viscosity = True|False 			# Activate viscosity
stencil = 7 				# Number of points of stencil

[shock capture]
shock capture = True|False 		# Activate shock capture procedure
stencil = 7 				# Number of points of stencil
method = pressure|dilatation 		# Capture based on pressure or dilatation

[figures]
figures = True|False 			# Activate figures
pml = True|False 			# Show PML

[save]
save = True|False 			# Activate save
path = results 				# Path to data file
filename = tmp 				# Data filename
compression = None|lzf 			# Activate compression
only p = True|False			# Save only pressure
probes = True|False                     # Activate probes
probes_locations = [[128, 128], [128, 192]] # Probe locations. Must be list of lists
```

## Changelog

### 0.8.5

* fix: raise ValueError when using old version of config file
* fix: Some bugs in graphics chg: argparse improved in solver

* chg: pcolorfast replaced by pcolormesh
* chg: graphics module rewritten
* chg: Probe definied with only 1 entry in cfg file
* chg: only_p option removed

* fix: residual calculation now also based on flow
* fix: config parser improved
* fix: Bug in cin
* fix: CFL with flow
* fix: Graphic improvements : colormaps

### 0.8.4
* new: white noise pressure source

### 0.8.3
* new: harmonic pressure source

### 0.8.2
* fix: compatibility with jupyter notebook

### 0.8.1
* chg: mpltools to mplutils

### 0.8.0
* Changes in ofdlib2 & fdgrid

### 0.7.2
* minor changes

### 0.7.1
* fdgrid 0.6.4 changes : plot_obstacles -> plot_subdomains

### 0.7.0
* Curvilinear mesh official support

### 0.6.1
* ... with readme...
* Start to adapt code for adaptative/curvi. meshes

### 0.6.0
* Minor changes
* Fix no action specified issue.
* argparse + entry point

### 0.5.2
* Some minor changes about structure

### 0.5.1
* Probes added

### 0.5.0
* PML support

### 0.4.7
* Update for fdgrid 0.5.0
* Bugfixes

### 0.4.6
* Obstacles checking added

### 0.4.5
* Fix periodic boundaries

### 0.4.4
* Adapted to ofdlib2 v0.7

### 0.4.3
* README updated

### 0.4.2
* Check source location

### 0.4.1
* Fix capture + h5py deprecation warining

### 0.4.0
* Dispatch viscous flux
* minor changes
* Changes in config

### 0.3.0
* Dispatch capture + Minor changes

### 0.2.2-dev0
* dilatation + minor changes

### 0.2.1-dev0
* Minor changes. Patches removed.

### 0.2.0-dev0
* cin & sfilter simplified

### 0.1.3-dev0
* Filtering simplified

### 0.1.2-dev0
* minor changes

### 0.1.1-dev0
* Start patches
* Convert some functions to 'in-place function'
* Adaptated to ofdlib2 changes
* Upgraded to ofdlib2
* Viscosity basics
* First commit
