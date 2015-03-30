#!/usr/bin/env python3
"""
Demonstration of NRL LCPFCT code in Python
Michael Hirsch

f2py3 -m fast2d -h fast2d.pyf fast2d.f lcpfct.f gasdyn.f
f2py3 -c fast2d.pyf fast2d.f gasdyn.f lcpfct.f
python plot_fast2d.py
"""
from __future__ import division
from matplotlib.pyplot import draw, pause,subplots, show
maxtstep = 801

def runfast2d():
    from fast2d import fast2d #fortran code needs f2py3 first as noted in comments
    rho,vr,vz,erg = fast2d() #fortran to c order axes

    return rho.transpose(2,1,0), vr.transpose(2,1,0), vz.transpose(2,1,0), erg.transpose(2,1,0)

def plotfast2d(rho,vr,vz,erg):
    fg,ax = subplots(2,2)
    ax = ax.ravel()
    ht = fg.suptitle('ti=0')

    for i,(r,p,v,z) in enumerate(zip(rho,erg,vr,vz)):
        ht.set_text('ti='+str(i))
        ax[0].cla(); ax[0].imshow(r,origin='lower'); ax[0].set_title('density')

        ax[1].cla(); ax[1].imshow(p,origin='lower'); ax[1].set_title('pressure')

        ax[2].cla(); ax[2].imshow(v,cmap='bwr',origin='lower'); ax[2].set_title('vr')

        ax[3].cla(); ax[3].imshow(z,cmap='bwr',origin='lower'); ax[3].set_title('vz')

        draw()
        pause(0.01)
        #fg.savefig('/tmp/{:03d}.png'.format(i),dpi=150,bbox_inches='tight')


if __name__ == '__main__':
    rho,vr,vz,erg = runfast2d()



    plotfast2d(rho,vr,vz,erg)
    show()
