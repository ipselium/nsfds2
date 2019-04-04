# nsfds2 : 2D Navier-Stokes Finite Differences Solver


# Dependencies

* numpy
* matplotlib
* cython
* h5py
* progressbar33
* ofdlib
* fdgrid

# Config file

```
[simulation]
nt = 150
ns = 10
cfl = 0.5
npml = 15

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

[source]
type = pulse
ixs = 32
izs = 32
s0 = 1e3

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
probes = False

[figures]
figures = True

[save]
save = True
filename = tmp
compression = lzf
only p = False
view = p
```

