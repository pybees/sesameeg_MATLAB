.. SESAMEEG documentation master file, created by
   sphinx-quickstart on Mon Jan 18 14:44:12 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SESAMEEG
========

This is a MATLAB implementation of the Bayesian multi--dipole localization method SESAME [1]_
(SEquential Semi-Analytic Montecarlo Estimator) for the automatic estimation of
brain source currents from MEEG data, either in the time domain and in the frequency domain [2]_.

A mathematical description of the algorithm is available in the :doc:`documentation <explanation_algorithm>`.

The algorithm takes in input a sourcespace, a leadfield and a data time series, and outputs a posterior probability map for source locations, the estimated number of dipoles, their locations and their amplitudes.

:scpt:`script_run_SESAME` contains an example script that launches the algorithm
:func:`inverse_SESAME` contains the Monte Carlo algorithm 
:func:`inverse_SESAME_viewer` contains code for visualizing the output.


Bug reports
===========

Use the `github issue tracker <https://github.com/pybees/sesameeg/issues>`_ to report bugs.


Authors of the code
-------------------
| Alberto Sorrentino <sorrentino@dima.unige.it>
| Alessandro Viani <viani@dima.unige.it>
| Gianvittorio Luria <luria@dima.unige.it>,
| Sara Sommariva <sommariva@dima.unige.it>,



Cite our work
=============

If you use this code in your project, please consider citing our work:

.. [1] S. Sommariva and A. Sorrentino, `Sequential Monte Carlo samplers for semi-linear inverse problems and application to Magnetoencephalography <https://doi.org/10.1088/0266-5611/30/11/114020>`_. Inverse Problems, 30 114020 (2014).
.. [2] G. Luria et al., `Bayesian multi-dipole modelling in the frequency domain <https://doi.org/10.1016/j.jneumeth.2018.11.007>`_. J. Neurosci. Meth., 312 (2019) 27–36.

.. toctree::
   :maxdepth: 3
   :hidden:
   
   examples/index
   api
   explanation_algorithm
   
