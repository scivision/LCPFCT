# LCPFCT

[![image](https://zenodo.org/badge/32956122.svg)](https://zenodo.org/badge/latestdoi/32956122)
[![Actions Status](https://github.com/scivision/LCPFCT/workflows/ci/badge.svg)](https://github.com/scivision/LCPFCT/actions)


NRL Flux-Corrected Transport algorithm for Solving Generalized
Continuity Equations--now in Python with Examples!

![fancy output plot animated](tests/fast2d.gif)

demonstrates use of Fortran code called from Python. In this case, using
Fortran code as a Python module is about 50 times faster than the very
slow writing Fortran output to a text file, and parsing the text in
Python or Matlab.

## Install

For the Python wrapping Fortran:

    pip install -e .

(optional) to use just Fortran alone:

    cmake -B build
    cmake --build build --parallel

## Examples

### 2-D explosion

    ./runfast2d.py

### 1-D shock

    ./runshock.py

## References

<https://www.nrl.navy.mil/lcp/LCPFCT>
