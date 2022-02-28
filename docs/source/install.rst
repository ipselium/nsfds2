Introducing **nsfds2**
======================

**nsfds2** is 2D Navier-Stokes Solver that uses finite difference method. In particular, **nsfds2** is specialized in acoustic simulations.

**nsfds2** is still in developpement. It is still full of bugs and comes with **ABSOLUTELY NO WARRANTY**.


Dependencies
------------

:python: >= 3.6
:numpy: >= 1.1
:matplotlib: >= 3.0
:h5py: >= 2.8
:progressbar33: >= 2.4
:ofdlib2: >= 0.10.0
:fdgrid: >= 0.10.0
:mplutils: >= 0.3.0

**Important:** To create animations using `nsfds2 make movie`, you also need to have **ffmpeg** installed on your system.


Installation
------------

Clone the repo at https://github.com/ipselium/nsfds2 and::

   python setup.py install

or

::

   pip install nsfds2


**Note:** To compile *ofdlib2*, OS X users may recquire :

::

   xcode-select --install
