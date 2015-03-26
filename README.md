# python-lcpfct
NRL Flux-corrected Algorithm for Solving Generalized Continuity Equations--now in Python with Examples!

demonstrates use of Fortran code called from Python. In this case, using Fortran code as a Python module
is about 50 times faster than the loathsome writing Fortran output to a text file, and parsing the text
in Python or Matlab.

Prereqs:
--------
``` pip install -r requirements.txt ```

Fortran Compile:
----------------
```
f2py3 -m shock -h shock.pyf shock.f lcpfct.f gasdyn.f
f2py3 -c shock.pyf shock.f gasdyn.f lcpfct.f
```

That's all! By typing ``` python3 plotshock.py ``` you'll see an animation of the 1-D shock, which is an 
example from http://www.nrl.navy.mil/lcp/LCPFCT
