import os
import pylab as pl
from scipy.integrate import odeint
from scipy.integrate import ode
import numpy
from mpl_toolkits.mplot3d import Axes3D
import random



class surfCentreLineApx(object):
    def __init__(self,coef,k,w,drgCoef,dt):
        self.dt = dt 
        self.coef = coef
        self.k = k
        self.w = w
        self.drg = drgCoef
        self.sol = pl.array([]) 
    
    def set_sol(self,sol):
        self.sol=sol

    # just make normal functions to try to pass into odeint function. Should be much faster
    def f(self,xarr,t):
        temp = 0.0
    
        for i in range(2):
            temp+=pl.sin(self.k*xarr[2]-i*pl.pi)*pl.cos(self.w*t-i*pl.pi)/(pl.cosh(self.k*xarr[3])-pl.cos(self.k*xarr[2]-i*pl.pi)) 
        temp = temp*self.coef
        temp -= self.drg*xarr[0]
        x1dot = temp
        x2dot = 0.0
        x3dot = xarr[0]
        x4dot = 0.0
        return [x1dot,x2dot,x3dot,x4dot]
#**********************************************************************************************
#**********************************************************************************************
def get_poin(sol,dt):
    poin = pl.array([])
    for i in range(len(sol)):
        if(((i*dt)%(2.0*pl.pi))<=dt):
            poin = pl.append(poin,sol[i,:])
    poin = poin.reshape(-1,4)
    # flip order of array
    #poin = poin[::-1,:]
    return poin
#**********************************************************************************************
#**********************************************************************************************
def plot_sol(sol,dt):
    
    poin = get_poin(sol,dt)
    poin[:,2]=poin[:,2]%(2*pl.pi)

    first_fig = pl.figure()
    first_ax = first_fig.add_subplot(111)
    first_ax.scatter(poin[:,2],poin[:,0],color="Black",s=.1)
    first_ax.set_xlabel("$x_1$",fontsize=30)
    first_ax.set_ylabel("$x_2$",fontsize=30)
    first_ax.set_xlim([0,2*pl.pi])
    #first_ax.set_ylim([-1.3,1.3])
    first_fig.tight_layout()
    first_fig.savefig("plot.png")
    os.system("open plot.png")
    
#**********************************************************************************************
#**********************************************************************************************
def main():
    
    dt = .05 
    # total number of sterations to perform inorder to get cycles rioht
    totIter = 1000000
    totTime = totIter*dt
    first_time = pl.arange(0.0,totTime/4.0,dt)
    time = pl.arange(0.0,totTime,dt)
    
    
    surf = 1.0
    #coef = .3
    coef = 1.7
    k = 1.0
    w = 1.0
    damp = .1
    g = .1

    # how many cells is till periodicity use x = n*pi/k (n must be even #) modNum = 2*pl.pi/k
    modNum = 2.0*pl.pi
    
    # some random initial conditions
    initx = 4.549
    inity = 1.0
    initvx = -.2037
    initvy = 0.0

    # now find initial condition in stange atractor by running far in time. print this to use it for
    # later
    x0 = pl.array([initvx,initvy,initx,inity])
    
    apx = surfCentreLineApx(coef,k,w,damp,dt)
    first_sol = odeint(apx.f,x0,first_time)
    first_sol[:,2] = first_sol[:,2]%(2*pl.pi)
    sol = odeint(apx.f,first_sol[-1,:],time)
    plot_sol(sol,dt)

   
if __name__ == '__main__':
    main()
