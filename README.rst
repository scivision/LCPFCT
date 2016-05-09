.. image:: https://travis-ci.org/scienceopen/pyLCPFCT.svg?branch=master
    :target: https://travis-ci.org/scienceopen/pyLCPFCT
.. image:: https://codeclimate.com/github/scienceopen/pyLCPFCT/badges/gpa.svg
   :target: https://codeclimate.com/github/scienceopen/pyLCPFCT

==============
python-lcpfct
==============
NRL Flux-Corrected Transport algorithm for Solving Generalized Continuity Equations--now in Python with Examples!


.. image:: http://blogs.bu.edu/mhirsch/files/2015/03/fast2d.gif
   :alt: fancy output plot animated

demonstrates use of Fortran code called from Python. In this case, using Fortran code as a Python module
is about 50 times faster than the very slow writing Fortran output to a text file, and parsing the text
in Python or Matlab.

Install
=======
For the Python wrapping Fortran::

    python setup.py develop

(optional) to use just Fortran alone::
    
    cd bin
    cmake ..
    make


Examples
========

2-D explosion simulation
---------------------------
::

    ./runfast2d.py


1-D shock simulation
---------------------
::

    ./runshock.py

References
==========

http://www.nrl.navy.mil/lcp/LCPFCT
