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

