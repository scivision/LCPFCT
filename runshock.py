#!/usr/bin/env python3
"""
Demonstration of NRL LCPFCT code in Python
Michael Hirsch

Old-fashioned way:
gfortran runshock.f shock.f lcpfct.f gasdyn.f -o shock
./shock  #writes fort.11 file
python runshock.py

New way-- 50+ times faster!:
f2py3 -m shock -h shock.pyf shock.f lcpfct.f gasdyn.f
f2py3 -c shock.pyf shock.f gasdyn.f lcpfct.f
python runshock.py

"""
import xarray
import numpy as np
from matplotlib.pyplot import draw, pause, figure
import shock  # fortran code needs f2py3 first as noted in comments

nx = 50


def runshock() -> xarray.Dataset:

    arr = shock.shock(nx)
    x = arr[:nx, 5]
    arr = arr[:, :5].reshape((-1, nx, 5), order='C')
    dt = 0.05
    time_sec = np.arange(0, 10.1, dt)

    ds = xarray.Dataset(
        {
            'Density': (('time', 'x'), arr[..., 0]),
            'Temperature': (('time', 'x'), arr[..., 1]),
            'Pressure': (('time', 'x'), arr[..., 2]),
            'Velocity': (('time', 'x'), arr[..., 3]),
            'Energy': (('time', 'x'), arr[..., 4]),
        },
        coords={'time': time_sec, 'x': x},
    )

    return ds


def plotshock(ds: xarray.Dataset):
    fg = figure()
    ax = fg.subplots(5, 1, sharex=True)
    fg.subplots_adjust(hspace=0.05)
    fg.suptitle('1-D shock')
    ax[-1].set_xlabel('x displacement')

    for i, k in enumerate(ds):
        ax[i].set_ylabel(k)

    for j in range(len(ds.time)):
        for i, k in enumerate(ds):
            ax[i].plot(ds[k][j, :])
        draw()
        pause(0.5)


def main():
    fdata = runshock()

    plotshock(fdata)


if __name__ == '__main__':
    main()
