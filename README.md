# python-lcpfct
NRL Flux-Corrected Transport algorithm for Solving Generalized Continuity Equations--now in Python with Examples!


![alt fancy plot](http://blogs.bu.edu/mhirsch/files/2015/03/fast2d.gif)

demonstrates use of Fortran code called from Python. In this case, using Fortran code as a Python module
is about 50 times faster than the very slow writing Fortran output to a text file, and parsing the text
in Python or Matlab.

Prereqs:
--------
``` pip install -r requirements.txt ```

## 2-D explosion simulation
```
f2py3 -m fast2d -h fast2d.pyf fast2d.f lcpfct.f gasdyn.f
f2py3 -c fast2d.pyf fast2d.f gasdyn.f lcpfct.f
```
plot by:
```
python plot_fast2d.py
```

## 1-D shock simulation
```
f2py3 -m shock -h shock.pyf shock.f lcpfct.f gasdyn.f
f2py3 -c shock.pyf shock.f gasdyn.f lcpfct.f
```
type
``` 
python3 plotshock.py 
``` 
and you'll see an animation of the 1-D shock.

These examples are from http://www.nrl.navy.mil/lcp/LCPFCT
