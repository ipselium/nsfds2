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
* ofdlib2 >= 0.9.3
* fdgrid >= 0.7.2
* mplutils >= 0.3.0

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

See `-h` option for further help :

```
nsfds2 solve -h
nsfds2 movie -h
nsfds2 show -h
```


## Custom Use

**nsfds2** can also be used as a classical Python package. The following
example gives the general philosophy :


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


## Config file

By default, ***nsfds2*** create the config file `~/.nfds2/nsfds2.conf`. This
file contains simulation paremeters that are used by the solver.


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
probes = True|False 			# Show probes
pml = True|False 			# Show PML

[save]
save = True|False 			# Activate save
path = results 				# Path to data file
filename = tmp 				# Data filename
compression = None|lzf 			# Activate compression
probes = [[128, 128], [128, 192]] # Probe locations. Must be list of lists
```
