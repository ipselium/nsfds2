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


## Sample Config file

```
[configuration]
timings = True
quiet = False

[simulation]
nt = 500
ns = 10
cfl = 0.5

[geometry]
mesh = regular
file = None
geoname = helmholtz_double
curvname = curvz
bc = ARRA
nx = 256
nz = 256
ix0 = 128
iz0 = 0
dx = 1e-4
dz = 1e-4

[PML]
beta = 0.0
alpha = 4.0
sigmax = 20
sigmaz = 20
npml = 15

[source]
type = pulse
ixs = 64
izs = 128
s0 = 1e6
f0 = 60000

[eulerian fluxes]
stencil = 11

[filtering]
filter = True
stencil = 11
stength = 0.75

[viscous fluxes]
viscosity = True
stencil = 7

[shock capture]
shock capture = True
stencil = 7
method = pressure

[figures]
figures = True

[save]
save = True
path = results
filename = tmp
compression = lzf
only p = False
probes = False
probes_locations = [[128, 128], [128, 192]]
```

## Changelog

### 0.8.2
* new: harmonic pressure source

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
