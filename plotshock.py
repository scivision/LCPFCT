#!/usr/bin/env python3
"""
Demonstration of NRL LCPFCT code in Python
Michael Hirsch

Old-fashioned way:
gfortran -O2 runshock.f shock.f lcpfct.f gasdyn.f -o shock
./shock  #writes fort.11 file
python plotshock.py

New way:
f2py3 -m shock -h shock.pyf shock.f lcpfct.f gasdyn.f
f2py3 -c shock.pyf shock.f gasdyn.f lcpfct.f

from shock import shock
>>> dir(shock)
['__doc__', '__file__', '__loader__', '__name__', '__package__', '__spec__', '__version__', 'arrays', 'conserve', 'fct_grid', 'fct_misc', 'fct_ndex', 'fct_scrh', 'fct_velo', 'gasdyn', 'lcpfct', 'makegrid', 'old_grid', 'residiff', 'set_grid', 'shock', 'sources', 'velocity', 'zerodiff', 'zeroflux']

"""
from __future__ import division
from pandas import Panel
from io import StringIO
from os.path import expanduser
from numpy import loadtxt, arange, asarray
from matplotlib.pyplot import draw, pause,subplots
from time import time

nx = 50

def runshock():
    from shock import shock #fortran code needs f2py3 first as noted in comments
    return shock(nx)

def readshock(fn):
    fn = expanduser(fn)

#%% get first header
    with open(fn,'r') as f:
        hd = f.readline().split() #column names for Panel
        hd2 = f.readline().split()
    if hd2[0] != '1':
        print('** warning, I appear to not be reading the header correctly')

    nx = int(hd2[-4])
    dt = float(hd2[-1])

#%% read data frames
    d = []; lines=''
    with open(fn,'r') as f:
        while True:
            lines=''
            #get to next record
            line = f.readline()
            if not line:
                break
            if line[0] != '1':
                continue
            #read the record
            for i in range(nx):
                lines += f.readline()
            d.append(loadtxt(StringIO(lines), usecols=(1,2,3,4,5,6)))
#%% setup output
    nt = len(d)
    t = arange(0,nt*dt,dt)

    dat = asarray(d)
    return Panel(dat[...,:5], t, dat[0,:,5], hd[1:-1])

def plotshock(dat):
    fg,ax = subplots(5,1,num=1,sharex=True)
    fg.subplots_adjust(hspace=0.05)
    ht = fg.suptitle('1-D shock')
    for i,j in enumerate(dat.minor_axis.values):
        ax[i].set_ylabel(j)

    for t,da in dat.iteritems():
        ht.set_text('1-D shock: t={:.3f} sec.'.format(t))
        for i,(j,d) in enumerate(da.iteritems()):
            ax[i].plot(d,label=j)
        draw()
        pause(0.5)



if __name__ == '__main__':
    from numpy.testing import assert_allclose
    try:
        tic = time()
        rhon = runshock()
        forttime = time()-tic
        print('fortran took {:0.3e} seconds'.format(forttime))
    except Exception as e:
        print('*** could not run shock via f2py, trying to read text output from previous run')
        print(str(e))


    tic = time()
    data = readshock('fort.11')
    print('reading and parsing text file took {:0.1f} times longer'.format(25*(time()-tic)/forttime))

    #assert_allclose(rhon[:50,:5],data[0,:50,:],rtol=1e-4,atol=1)
    #plotshock(data)