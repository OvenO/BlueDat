import os
import pylab as pl
from scipy.integrate import odeint


class surfCentreLineApx(object):
    def __init__(self,coef,k,w,drgCoef):
        self.coef = coef
        self.k = k
        self.w = w
        self.drg = drgCoef
    
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


def main():
    
    dt = .2 
    # total number of iterations to perform
    totIter = 500000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    
    surf = 1.0
    coef = .2571
    k = 1.0
    w = 1.0
    damp = .1
    g = .1
# how many cells is till periodicity use x = n*pi/k (n must be even #) modNum = 2*pl.pi/k
    modNum = 2.0*pl.pi
    
    # initial conditions
    initx = pl.pi
    inity = 1.0
    initvx = 0.18
    initvy = 0.0
# initial conditions vector
    # set up: [xdot,ydot,x,y]
    x0 = pl.array([initvx,initvy,initx,inity])
    apx = surfCentreLineApx(coef,k,w,damp)
    sol = odeint(apx.f,x0,time)
    
    if "PositionTimeImgs" not in os.listdir("."):
        os.mkdir("PositionTimeImgs")
    os.chdir("PositionTimeImgs")
    
    sol[:,2] = sol[:,2]%modNum
    fig1 = pl.figure()
    ax1 = fig1.add_subplot(111)
    last_quarter = int(len(time)/8.0)
    ax1.plot(time[-last_quarter:],sol[-last_quarter:,2],color="Black")
    ax1.set_xlabel("$t$",fontsize=36)
    ax1.set_ylabel("$x$",fontsize=36)
    pl.savefig("t_vs_x.png")

    os.system("open t_vs_x.png")


if __name__ == '__main__':
    main()
