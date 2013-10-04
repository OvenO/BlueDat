import os
import pylab as pl
from scipy.integrate import odeint
from scipy.integrate import ode
import numpy
from mpl_toolkits.mplot3d import Axes3D
import random

# SEE sprott.physics.wisc.edu/chaos/lornzle.htm for comparison

class lorenz(object):
    def __init__(self,p,r,b,dt):
        self.dt = dt 
        self.p = p
        self.r = r
        self.b = b
        self.sol = pl.array([]) 
    

    # just make normal functions to try to pass into odeint function. Should be much faster
    def f(self,xarr,t):
        x1dot = self.p*(xarr[1]-xarr[0])
        x2dot = -xarr[0]*xarr[2]+self.r*xarr[0]-xarr[1]
        x3dot = xarr[0]*xarr[1]-self.b*xarr[2]
    
        return [x1dot,x2dot,x3dot]

#**********************************************************************************************
#**********************************************************************************************
# function just finds the magnitude of the distance of the two input points
def distance(x1,x2):
    x_sqrd = (x1[0]-x2[0])**2
    y_sqrd = (x1[1]-x2[1])**2
    z_sqrd = (x1[2]-x2[2])**2
    return pl.sqrt(x_sqrd + y_sqrd + z_sqrd)
#**********************************************************************************************
#**********************************************************************************************
# this function needs both the final positions of the two trajectories and the variable eplilon wich
# is the starting/reset DISTANCE of the two poinnts. jsee lab book #2 pg 59 for a
# more detailed explination. "renomalize" will retun a point that is on the line contecting the
# final points of the trajectorys and will this new point will be a distance epsilon from the
# unpurtubed trjectory (x_unpurt). x_putub is final position of purtubed trajectorie
def renormalize(x_unpurt,x_puturb,epsilon):
    final_dist = distance(x_unpurt,x_puturb)
    print "final_dist is:"
    print final_dist
    xnew = pl.array([0.0,0.0,0.0])
    # the new renormalized vx (see lab book #2 pg 61)
    # THIS MIGHT BE WRONG!!!!! CHECK LATER FOR ERRORS IF NUMBER COMES OUT WRONG
    xnew[0] = x_unpurt[0]+(epsilon/final_dist)*(x_puturb[0]-x_unpurt[0])
    xnew[1] = x_unpurt[1]+(epsilon/final_dist)*(x_puturb[1]-x_unpurt[1])
    xnew[2] = x_unpurt[2]+(epsilon/final_dist)*(x_puturb[2]-x_unpurt[2])

    #figr = pl.figure()
    #axr = figr.add_subplot(111,projection="3d")
    #axr.plot([x_unpurt[0],x_puturb[0]],[x_unpurt[1],x_puturb[1]],[x_unpurt[2],x_puturb[2]],color="Black")
    #axr.scatter(x_unpurt[0],x_unpurt[1],x_unpurt[2],color="Red")
    #axr.scatter(xnew[0],xnew[1],xnew[2],s=10)
    #axr.set_xlabel("$x_1$",fontsize=25)
    #axr.set_ylabel("$x_2$",fontsize=25)
    #axr.set_zlabel("$x_3$",fontsize=25)
    ##axr.set_ylim([-1.3,1.3])
    #figr.tight_layout()
    #figr.savefig("checkrenormalize.png",dpi=300)
    #os.system("open checkrenormalize.png")
 
    return xnew
#**********************************************************************************************
#**********************************************************************************************
def main():
    
    # initial purturbation size
    epsilon = 1.0e-10
    #epsilon = 1.0e-12
    #epsilon = .0005
    # number of cycle-sets to throw away (even after transients) before keeping any information
    throw_away = 500
    # number of cycle-sets to go through after transiense
    #num = 3000
    num = 1000
    # time period of 1 cycle (because there is no driving frequency reference we are really just
    # guessing at this.
    step = .001
    

    dt = .001 
    totIter = step/dt
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    
    #time array for the first run to get initial condition in strange atractor
    first_time = pl.arange(0.0,.5/dt,dt)
    
    p = 10.0
    r = 28.0
    b = 8.0/3.0

    x0 = pl.array([.5,.5,.5])
    
    apx = lorenz(p,r,b,dt)
    first_sol = odeint(apx.f,x0,first_time)

    # plot of atractor
    #fig = pl.figure()
    #ax = fig.add_subplot(111,projection="3d")
    #ax.plot(first_sol[:,0],first_sol[:,1],first_sol[:,2],color="Black")
    #ax.set_xlabel("$x_1$",fontsize=25)
    #ax.set_ylabel("$x_2$",fontsize=25)
    #ax.set_zlabel("$x_3$",fontsize=25)
    ##ax.set_ylim([-1.3,1.3])
    #fig.tight_layout()
    #fig.savefig("sep.png",dpi=300)
    #os.system("open sep.png")
    
    # this should now be a poin in the atractor 
    x0 = first_sol[-1,:]
    print("x0 from first_sol is:")
    print(x0)

    x_other = x0+pl.array([pl.sqrt(epsilon),pl.sqrt(epsilon),0.0])

    # two arrays to store solution data
    arr_orig = pl.array([])
    arr_othr = pl.array([])
    final_dist_arr = pl.array([])

    for i in range(num + throw_away):

        sol = odeint(apx.f,x0,time)
        sol_other = odeint(apx.f,x_other,time)
        print sol
            
        if i>throw_away:
            dist = distance(sol[-1,:],sol_other[-1,:])
            final_dist_arr = pl.append(final_dist_arr,dist) 

        x0 = sol[-1,:]
        x_other = renormalize(x0,sol_other[-1,:],epsilon)
        print i
    
    arr_orig = arr_orig.reshape(-1,3)
    arr_othr = arr_othr.reshape(-1,3)
   
    #NOW CALCULATE LE
    eps_arr = pl.zeros(len(final_dist_arr))+epsilon
    le = pl.zeros(len(final_dist_arr)) + pl.log(abs(final_dist_arr/epsilon))
    print "le is:"
    print le[-1]
    print max(le)

    pl.scatter(pl.arange(len(le)),le,s=0.1)
    pl.savefig("le.png")
    os.system("open le.png")
    
    # plot of PC sections (red Blue)
    #fig = pl.figure()
    #ax = fig.add_subplot(111)
    ##ax.scatter([0.0,pl.pi,2.0*pl.pi],[0.0,0.0,0.0],color="Red")
    ##ax.scatter(arr_orig[:,2],arr_orig[:,0],color="Red",s=.1)
    ##ax.scatter(arr_othr[:,2],arr_othr[:,0],color="Blue",s=.1)
    #ax.scatter(full_orig[:,2],full_orig[:,0],color="Red",s=.1)
    #ax.scatter(full_othr[:,2],full_othr[:,0],color="Blue",s=.1)
    #ax.set_xlabel("$x_1$",fontsize=25)
    #ax.set_ylabel("$x_2$",fontsize=25)
    ##ax.set_xlim([0,2*pl.pi])
    ##ax.set_ylim([-1.3,1.3])
    #fig.tight_layout()
    #fig.savefig("sep.png")
    #os.system("open sep.png")
    
if __name__ == '__main__':
    main()
