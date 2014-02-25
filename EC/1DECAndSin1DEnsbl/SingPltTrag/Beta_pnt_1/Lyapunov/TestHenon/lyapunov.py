import os
import pylab as pl
from scipy.integrate import odeint
from scipy.integrate import ode
import numpy
from mpl_toolkits.mplot3d import Axes3D
import random


#**********************************************************************************************
#**********************************************************************************************
# Here we will find a second point that is some small distace (epsilon) from the x0 point.
# We need this function because we want the distace to be epsilon but the position about x0 to be
# random
def get_first_init(x0,epsilon):
    # we will use a change of coordinates to get the location of the particle relative to x0. First
    # we just find some random point a distace epsilon from the origin.
    theta = random.random()*2.0*pl.pi
    x = epsilon*pl.cos(theta)
    vx = epsilon*pl.sin(theta)
    
    # change of coordinates
    return x0 + pl.array([vx,0,x,0])

#**********************************************************************************************
#**********************************************************************************************
# this function needs both the final positions of the two trajectories and the variable eplilon wich
# is the starting/reset DISTANCE of the two poinnts. jsee lab book #2 pg 59 for a
# more detailed explination. "renomalize" will retun a point that is on the line contecting the
# final points of the trajectorys and will this new point will be a distance epsilon from the
# unpurtubed trjectory (x_unpurt). x_putub is final position of purtubed trajectorie
def renormalize(x,y,a,b,epsilon):
    final_dist = pl.sqrt((a-x)**2 + (b-y)**2)
    # the new renormalized vx (see lab book #2 pg 61)
    xnew = x+(epsilon/final_dist)*(a-x)
    ynew = y+(epsilon/final_dist)*(b-y)

    return xnew,ynew,final_dist
#**********************************************************************************************
#**********************************************************************************************
def henon(xn,yn):
    a = 1.4
    b = 0.3
    xnp1 = 1 - a*xn**2 + yn
    ynp1 = b*xn

    return xnp1,ynp1

#**********************************************************************************************
#**********************************************************************************************
def main():
    
    # initial purturbation size
    epsilon = 1.0e-12
    
    # first throw away (to get initial conditions in atractor)
    throw_away_1 = 200
    # second throw away (to make sure puturbed trajectory is aligned)
    throw_away_2 = 200

    # number of cycle-sets to go through after throw_aways
    num = 6000
    
    # initial conditions
    xn = 1.0
    yn = 1.0
    
    
    for i in range(throw_away_1):
        xn,yn = henon(xn,yn)
    
    an = xn + epsilon
    bn = yn
    
    for i in range(throw_away_2):
        xn,yn = henon(xn,yn)
        an,bn = henon(an,bn)

        # renormalize the purtubed trajectory
        an,bn,dist = renormalize(xn,yn,an,bn,epsilon)

    # two arrays to store solution data
    xn_ar = pl.array([])
    yn_ar = pl.array([])
    an_ar = pl.array([])
    bn_ar = pl.array([])

    # final distance data
    darr = pl.array([])
    for i in range(num):
        xn,yn = henon(xn,yn)
        an,bn = henon(an,bn)

        xn_ar = pl.append(xn_ar,xn)
        yn_ar = pl.append(yn_ar,yn)
        an_ar = pl.append(an_ar,an)
        bn_ar = pl.append(bn_ar,bn)

        # renormalize the purtubed trajectory
        an,bn,dist = renormalize(xn,yn,an,bn,epsilon)
        
        darr = pl.append(darr,dist)
    
    
    eps_arr = pl.zeros(len(darr))+ epsilon
    le = pl.zeros(len(darr))+ pl.log(abs(darr/epsilon))

    le_avg = 0.0
    for i,j in enumerate(le):
        le_avg+=j
    le_avg = le_avg/len(le)

    print "le is:"
    print(le_avg)

    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.scatter(pl.arange(len(le)),le,s = .1)
    #ax.set_xlabel("$x$",fontsize=25)
    #ax.set_ylabel("$y$",fontsize=25)
    #ax.set_xlim([0,2*pl.pi])
    #ax.set_ylim([-1.3,1.3])
    fig.tight_layout()
    fig.savefig("le.png")
    os.system("open le.png")
    
    fig1 = pl.figure()
    ax1 = fig1.add_subplot(111)
    ax1.scatter(xn_ar,yn_ar,s = .1)
    #ax1.set_xlabel("$x$",fontsize=25)
    #ax1.set_ylabel("$y$",fontsize=25)
    #ax1.set_xlim([0,2*pl.pi])
    #ax1.set_ylim([-1.3,1.3])
    fig1.tight_layout()
    fig1.savefig("henon.png")
    os.system("open henon.png")
    
if __name__ == '__main__':
    main()
