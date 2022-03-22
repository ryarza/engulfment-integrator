Integrator for the equation of motion of an object embedded in a stellar envelope, subject to hydrodynamical and gravitational forces.

Automatically generated documentation and call graphs from Doxygen are available at `Link text <../../../doxygen/build/html/index.html>`_


Getting started
===============

Prerequisites
-------------

The requirements are:

  * GNU scientific library, which provides the ODE integrator.
  * HDF5, for I/O.

Compiling and testing
---------------------

You can run basic automated tests (Keplerian orbits) with::

  mkdir build
  cd build
  cmake ..
  make
  ctest

To compile the code, just run ``make`` on the top directory.

Integrating an inspiral
-----------------------

The ``extras`` folder has a sample ``input.txt`` that will integrate the trajectory of a Jupiter-like planet inside the envelope of a model of the Sun evolved to ten solar radii. You might need to edit ``input.txt`` so that it points to the actual path of the model on your system. After that, run the code with::

  engulfment-integrator input.txt

The code will produce a file called ``output.txt``, containing important variables as a function of time, and a file called ``scalars.txt``, containing CGS values for the simulation units as well as values for some input parameters. The columns in ``output.txt`` are

========== =======
Column     Meaning
========== =======
time       time
r          radial coordinate
theta      angular coordinate
rdot       rate of change of the radial coordinate
thetadot   rate of change of the angular coordinate
menc       stellar mass enclosed within ``r``
rho        stellar density at ``r``
q          mass ratio between the engulfed object and the enclosed mass at ``r``
eps_rho    dimensionless density gradient (number of density scale heights per accretion radius)
rp_over_ra ratio between geometrical and gravitational radii of the engulfed object
cp         ram pressure drag coefficient
cg         gravitational drag coefficient
f_ram      ram pressure drag force
f_grav     gravitational drag force
========== =======

All outputs are in simulation units.

Contents
========
.. toctree::
   :maxdepth: 2

   files

Indices and tables
==================

* :ref:`genindex`
