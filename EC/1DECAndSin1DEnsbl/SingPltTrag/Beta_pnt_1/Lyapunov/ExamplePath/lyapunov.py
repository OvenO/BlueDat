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
# Here we will find a second point that is some small distace (epsilon) from the x0 point.
# We need this function because we want the distace to be epsilon but the position about x0 to be
# random
def get_other_init(x0,epsilon):
    # we will use a change of coordinates to get the location of the particle relative to x0. First
    # we just find some random point a distace epsilon from the origin.
    theta = random.random()*2.0*pl.pi
    x = epsilon*pl.cos(theta)
    vx = epsilon*pl.sin(theta)
    
    # change of coordinates
    return x0 + pl.array([vx,0,x,0])

#**********************************************************************************************
#**********************************************************************************************
def get_poin(sol,dt):
    poin = pl.array([])
    for i in range(len(sol)):
        if(((i*dt)%(2.0*pl.pi))<dt):
            poin = pl.append(poin,sol[i,:])
    return poin.reshape(-1,4)
#**********************************************************************************************
#**********************************************************************************************
def main():
    
    # initial purturbation size
    epsilon = .01
    # The num of driving cycles before starting two particles back at epsilon distance from
    # eachother
    cycles = 1.0
    # number of cycle-sets to go through
    num = 2
    
    dt = .001 
    # total number of iterations to perform inorder to get cycles right
    totIter = cycles*2.0*pl.pi/dt
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    
    surf = 1.0
    coef = .3
    k = 1.0
    w = 1.0
    damp = .1
    g = .1

    # how many cells is till periodicity use x = n*pi/k (n must be even #) modNum = 2*pl.pi/k
    modNum = 2.0*pl.pi
    
    # initial conditions
    initx = 1.0
    inity = 1.0
    initvx = -.00
    initvy = 0.0

    data = open("data.txt","w")
    
    apx = surfCentreLineApx(coef,k,w,damp,dt)
    
    x0 = pl.array([initvx,initvy,initx,inity])
    x_other = get_other_init(x0,epsilon)
    
    print(x0)
    print(x_other)

    # two arrays to store solution data
    arr_orig = pl.array([])
    arr_othr = pl.array([])
    for i in range(num):
        sol = odeint(apx.f,x0,time)
        sol_other = odeint(apx.f,x_other,time)

        #sol[:,2]=sol[:,2]%(2*pl.pi)
        #sol_other[:,2]=sol_other[:,2]%(2*pl.pi)


        arr_orig = pl.append(arr_orig,sol)
        arr_othr = pl.append(arr_othr,sol_other)

# need to make an aray to apend and store these to for the final stuff. Allso need to put in the
        # stuff to calculate the LE.... or should I make a data file and make a seperate plot
        # proframe?? prablably the later of the two
        poin = get_poin(sol,dt)
        poin_other = get_poin(sol_other,dt)

        #arr_orig = pl.append(arr_orig,poin)
        #arr_othr = pl.append(arr_othr,poin_other)

        x0 = sol[-1,:]
        x_other = get_other_init(x0,epsilon)
    
    arr_orig = arr_orig.reshape(-1,4)
    arr_othr = arr_othr.reshape(-1,4)
    
    # plot of PC sections (red Blue)
    fig = pl.figure()
    ax = fig.add_subplot(111)
    #ax.scatter([0.0,pl.pi,2.0*pl.pi],[0.0,0.0,0.0],color="Red")
    ax.scatter(arr_orig[:,2],arr_orig[:,0],color="Red",s=.1)
    ax.scatter(arr_othr[:,2],arr_othr[:,0],color="Blue",s=.1)
    ax.set_xlabel("$x_1$",fontsize=25)
    ax.set_ylabel("$x_2$",fontsize=25)
    #ax.set_xlim([0,2*pl.pi])
    #ax.set_ylim([-1.3,1.3])
    fig.tight_layout()
    fig.savefig("sep.png")
    os.system("open sep.png")
    
if __name__ == '__main__':
    main()
