# nsdf2 : 2D Navier-Stokes Finite Differences Solver


# Dependencies

pip install h5py==2.8.0rc1
pip install progressbar33
pip install ofdlib
pip install fdgrid

# Config file

```
[simulation]
viscosity = True
probes = True
nt = 150
ns = 10
nx = 256
nz = 256
ix0 = 0
iz0 = 0
dx = 1
dz = 1
cfl = 0.5
bc = RRRR
npml = 15
stencil = 11

[filtering]
filter = True

[shock capture]
shock capture = True
method = pressure

[source]
type = pulse
ixs = 64
izs = 64
s0 = 1e6

[probes]
probes = False

[save]
save = True
filename = tmp
compression = lzf
only p = False
view = p
```

