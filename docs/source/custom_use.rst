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
   cfg = CfgSetup()    # or cfg = CfgSetup('path_to_configfile.conf')

   # Define the mesh
   msh = Mesh((cfg.nx, cfg.nz), (cfg.dx, cfg.dz), origin=(cfg.ix0, cfg.iz0), obstacles=cfg.obstacles)

   # Init acoustic fields
   fld = Fields(msh, cfg)

   # Create simulation
   fdtd = FDTD(msh, fld, cfg)
   fdtd.run()

   # Figures
   plt.figure()
   plt.imshow(fld.p)
   plt.show()


One can also use the `pickle` module to load cfg instance from the .cfg file::


   import pickle
   with open('file.cfg', 'rb') as cfg_file:
       cfg = pickle.load(cfg_file)


hdf5 files
==========

If `save` option is selected, hdf5 are created. They contain the following variables.

+-------------------+--------------------------------------------------------+
| var               | variable                                               |
+===================+========================================================+
| rho_init          | initial density                                        |
+-------------------+--------------------------------------------------------+
| rhou_init         | initial product of density and x component of velocity |
+-------------------+--------------------------------------------------------+
| rhov_init         | initial product of density and z component of velocity |
+-------------------+--------------------------------------------------------+
| rhoe_init         | initial product of density and energy                  |
+-------------------+--------------------------------------------------------+
| rho_itX           | density                                                |
+-------------------+--------------------------------------------------------+
| rhou_itX          | product of density and x component of velocity         |
+-------------------+--------------------------------------------------------+
| rhov_itX          | product of density and z component of velocity         |
+-------------------+--------------------------------------------------------+
| rhoe_itX          | product of density and energy                          |
+-------------------+--------------------------------------------------------+
| x                 | x-grid                                                 |
+-------------------+--------------------------------------------------------+
| z                 | z-grid                                                 |
+-------------------+--------------------------------------------------------+
| probe_locations   | coordinates of probes                                  |
+-------------------+--------------------------------------------------------+
| probe_values      | pressure at probe locations                            |
+-------------------+--------------------------------------------------------+
| obstacles         | coordinates of obstacles                               |
+-------------------+--------------------------------------------------------+

The `rho_itX`, `rhou_itX`, `rhov_itX`, and `rhoe_itX` quantities gather the
field at each time iteration (2d array).

To access the acoustic pressure, one can use:: 

    import numpy as np
    from ofdlib2 import fdtd
    ...
    p  = np.empty_like(rho) 
    fdtd.p(p, rho, rhou, rhov, rhoe, gamma)

The variable `p` needs to be initialize first. Then, it is used as an input
argument of `fdtd.p`

Note that in curvilinear coordinates, there are following sets of coordinates :

+-------------------+--------------------------------------------------------+
| var               | variable                                               |
+===================+========================================================+
| xn                | numerical x-grid                                       |
+-------------------+--------------------------------------------------------+
| zn                | numerical z-grid                                       |
+-------------------+--------------------------------------------------------+
| xp                | physical x-grid                                        |
+-------------------+--------------------------------------------------------+
| zp                | physical z-grid                                        |
+-------------------+--------------------------------------------------------+
| J                 | Jacobian matrix of transformation                      |
+-------------------+--------------------------------------------------------+

cfg pickle files
================

A cfg file is also created for each simulation gathering its configuration. It
contains the configuration object relative to the simulation:: 

    import pickle
    
    with open(filename, 'rb') as _file:
        cfg = pickle.load(_file)

This object contains in particular :

- cfg.dx, cfg.dz, cfg.dt : spatial and time steps
- cfg.nx, cfg.nz, cfg.nt, cfg.ns : Number of points (spatial and temporal)
- cfg.p0, cfg.rho0, cfg.T0, cfg.c0, cfg.gamma, cfg.prandtl, cfg.mu : Thermophysical parameters
- ... and many other parameters.
