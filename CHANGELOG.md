# Changelog

## 0.9.8-dev

new: probe spectrograms

## 0.9.7

* chg: update docs

## 0.9.6

* fix: bad ofdlib dependency

## 0.9.5

* new: wave source support
* new: modules utils.sounds
* new: argument parser has now make subcommand for 'movie' and 'sound'

* fix: fdgrid dependency
* fix: '~' now interpreted as home in config file
* fix: display probes in grid/pgrid
* fix: fdgrid v0.8.4 needed

## 0.9.4

* fix: moving bc with curvilinear mesh

## 0.9.3

* fix: Matplotlib frontend for OSX - 2nd try !

## 0.9.2

* new: bc_profiles available
* chg: update examples

* fix: Probes with curvilinear mesh
* fix: ofdlib 0.9.4 compatibility
* fix: Matplotlib frontend for OSX

## 0.9.1

* fix: Check config file version at startup

## 0.9.0

* new: Wall source support with fdgrid >=0.8
* doc: improve doc and provide sample configuration files

## 0.8.8

* doc: add example using custom geofile

## 0.8.7

* doc: update README

## 0.8.6

* add: '~' is now understood in config file
* add: nsfds2 show pgrid available

* fix: update config examples for nsfds=0.8.5
* fix: Bug loading curvilinear function


## 0.8.5

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

## 0.8.4
* new: white noise pressure source

## 0.8.3
* new: harmonic pressure source

## 0.8.2
* fix: compatibility with jupyter notebook

## 0.8.1
* chg: mpltools to mplutils

## 0.8.0
* Changes in ofdlib2 & fdgrid

## 0.7.2
* minor changes

## 0.7.1
* fdgrid 0.6.4 changes : plot_obstacles -> plot_subdomains

## 0.7.0
* Curvilinear mesh official support

## 0.6.1
* ... with readme...
* Start to adapt code for adaptative/curvi. meshes

## 0.6.0
* Minor changes
* Fix no action specified issue.
* argparse + entry point

## 0.5.2
* Some minor changes about structure

## 0.5.1
* Probes added

## 0.5.0
* PML support

## 0.4.7
* Update for fdgrid 0.5.0
* Bugfixes

## 0.4.6
* Obstacles checking added

## 0.4.5
* Fix periodic boundaries

## 0.4.4
* Adapted to ofdlib2 v0.7

## 0.4.3
* README updated

## 0.4.2
* Check source location

## 0.4.1
* Fix capture + h5py deprecation warining

## 0.4.0
* Dispatch viscous flux
* minor changes
* Changes in config

## 0.3.0
* Dispatch capture + Minor changes

## 0.2.2-dev0
* dilatation + minor changes

## 0.2.1-dev0
* Minor changes. Patches removed.

## 0.2.0-dev0
* cin & sfilter simplified

## 0.1.3-dev0
* Filtering simplified

## 0.1.2-dev0
* minor changes

## 0.1.1-dev0
* Start patches
* Convert some functions to 'in-place function'
* Adaptated to ofdlib2 changes
* Upgraded to ofdlib2
* Viscosity basics
* First commit
