==========
Custom Use
==========

Introduction
============

**nsfds2** can also be used as a classical Python package. It provides
several main objects to perform numerical simulations :

- `init` package povides:

        - :py:class:`nsfds2.init.config.CfgSetup`: Parses all simulation
          parameters from *~/.nsfds2/nsfds2.conf*
        - :py:class:`nsfds2.init.fields.Fields`: Initialize all
          simulation parameters and fields

- `fdtd` module provides :py:class:`nsfds2.fdtd.FDTD` that runs the simulation
- `utils` package provides in particular :py:mod:`nsfds2.utils.graphics` module
  useful for post-treatment


 The following example gives the general philosophy to use **nsfds2**::

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


hdf5 files
==========

If `save` option is selected, hdf5 are created. They contain the following variables.

+-------------------+---------------------------------------------------+
| var               | variable                                          |
+===================+===================================================+
| rho               | density                                           |
+-------------------+---------------------------------------------------+
| rhou              | product of density and x component of velocity    |
+-------------------+---------------------------------------------------+
| rhov              | product of density and z component of velocity    |
+-------------------+---------------------------------------------------+
| rhoe              | product of density and energy                     |
+-------------------+---------------------------------------------------+
| x                 | x-grid                                            |
+-------------------+---------------------------------------------------+
| z                 | z-grid                                            |
+-------------------+---------------------------------------------------+
| probe_location    | coordinates of probes                             |
+-------------------+---------------------------------------------------+
| probe_value       | pressure at probe locations                       |
+-------------------+---------------------------------------------------+

For `rho`, `rhou`, `rhov`, and `rhoe` quantities, there are two distinct
variables. One with the `_init` suffix which represents the initial value of the
quantity (1d array) and one with the `_it` suffix which gathers the field at
each iteration (2d array).

To access the acoustic pressure, one can use:: 

    import numpy as np
    from ofdlib2 import fdtd
    ...
    p  = np.empty_like(rho) 
    fdtd.p(p, rho, rhou, rhov, rhoe, gamma)

Here `p` needs to be initialize first. Then, it is used as an input argument of
`fdtd.p` 
