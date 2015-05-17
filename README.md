[![Code Climate](https://codeclimate.com/github/scienceopen/pyLCPFCT/badges/gpa.svg)](https://codeclimate.com/github/scienceopen/pyLCPFCT)

# python-lcpfct
NRL Flux-Corrected Transport algorithm for Solving Generalized Continuity Equations--now in Python with Examples!


![fancy plot](http://blogs.bu.edu/mhirsch/files/2015/03/fast2d.gif)

demonstrates use of Fortran code called from Python. In this case, using Fortran code as a Python module
is about 50 times faster than the very slow writing Fortran output to a text file, and parsing the text
in Python or Matlab.

Installation:
------------
```
git clone --depth 1 https://github.com/scienceopen/pyLCPFCT
pip install -r requirements.txt 
```

2-D explosion simulation
---------------------------
```
f2py -m fast2d -c fast2d.f gasdyn.f lcpfct.f 
```
plot by:
```
python plot_fast2d.py
```

 1-D shock simulation
---------------------
```
f2py -m shock -c shock.f gasdyn.f lcpfct.f 
```
which creates a shock.so file to be used in Python, to see an example type
``` 
python plotshock.py 
``` 
and you'll see an animation of the 1-D shock.

These examples are from http://www.nrl.navy.mil/lcp/LCPFCT
