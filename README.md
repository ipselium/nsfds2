# nsfds2 : 2D Navier-Stokes Finite Differences Solver

nsfds2 is a 2D Navier-Stokes Solver that uses finite difference method.

# Dependencies

* numpy
* cython
* matplotlib
* h5py
* progressbar33
* ofdlib
* fdgrid

# Config file

```
[configuration]
timings = True
quiet = False

[simulation]
nt = 500
ns = 10
cfl = 0.5

[geometry]
file = None
name = square
bc = RRRR
nx = 256
nz = 256
ix0 = 0
iz0 = 0
dx = 1
dz = 1

[PML]
beta = 0.0
alpha = 4.0
sigmax = 20
sigmaz = 20
npml = 15

[source]
type = pulse
ixs = 32
izs = 128
s0 = 1e6

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

[probes]

[figures]
figures = True

[save]
save = True
path = results
filename = tmp
compression = lzf
only p = False
probes = True
probes_locations = [[128, 128], [128, 192]]
```

