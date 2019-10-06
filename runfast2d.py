#!/usr/bin/env python3
"""
Demonstration of NRL LCPFCT code in Python
"""
from matplotlib.pyplot import draw, pause, figure
import matplotlib.animation as anim
import fast2d
import xarray
from pathlib import Path
import numpy as np
import argparse

maxtstep = 801
WRITER = 'ffmpeg'
FPS = 10
CODEC = 'ffv1'
DPI = 72


def runfast2d():

    rho, vr, vz, erg = fast2d.fast2d()  # fortran to c order axes

    i = range(
        rho.shape[2]
    )  # time steps are dynamic, need to pass them out from Forran in the future
    ds = xarray.Dataset(
        {
            "density": (('time', 'r', 'z'), rho.transpose(2, 1, 0)),
            'velocity_r': (('time', 'r', 'z'), vr.transpose(2, 1, 0)),
            'velocity_z': (('time', 'r', 'z'), vz.transpose(2, 1, 0)),
            'pressure': (('time', 'r', 'z'), erg.transpose(2, 1, 0)),
        },
        coords={
            'time_index': i,
            'r': np.arange(0, 64 * 0.1, 0.1),
            'z': np.arange(0, 64 * 0.1, 0.1),
        },
    )
    return ds


def plotfast2d(ds: xarray.Dataset, outfile: Path):
    fg = figure(figsize=(8, 8))
    ax = fg.subplots(2, 2)
    ax = ax.ravel()
    ht = fg.suptitle('ti=0')

    cmap = (None, None, 'bwr', 'bwr')
    outfile = Path(outfile).expanduser()
    print('writing', outfile)
    Writer = anim.writers[WRITER]
    w = Writer(fps=FPS, codec=CODEC)
    with w.saving(fg, str(outfile), DPI):
        for i in ds.time_index:
            for j, k in enumerate(('density', 'pressure', 'velocity_r', 'velocity_z')):

                ht.set_text(f'ti={i.item()}')
                ax[j].cla()
                ax[j].imshow(ds[k][i, ...], origin='lower', cmap=cmap[j])
                ax[j].set_title(k)

            draw()
            pause(0.01)

            w.grab_frame()


def main():

    p = argparse.ArgumentParser()
    p.add_argument(
        'outfn', help='output movie file to write', nargs='?', default='fast2d.avi'
    )
    p = p.parse_args()

    ds = runfast2d()

    plotfast2d(ds, p.outfn)


if __name__ == '__main__':
    main()
