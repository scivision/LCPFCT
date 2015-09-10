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
.. code:: bash

 $ git clone --depth 1 https://github.com/scienceopen/pyLCPFCT
 $ conda install --file requirements.txt 

Examples
========

2-D explosion simulation
---------------------------

.. code:: bash

 $ f2py -m fast2d -c fast2d.f gasdyn.f lcpfct.f 
 $ python plot_fast2d.py


1-D shock simulation
---------------------

.. code:: bash

  $ f2py -m shock -c shock.f gasdyn.f lcpfct.f 
  $ python plotshock.py

References
==========

http://www.nrl.navy.mil/lcp/LCPFCT
