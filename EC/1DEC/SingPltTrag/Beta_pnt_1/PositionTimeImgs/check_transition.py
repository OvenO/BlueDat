import os
import pylab as pl
from scipy.integrate import odeint
from scipy.integrate import ode
import numpy



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

    # define a funciton that grabs the matrix elements of the jacobian, set_sol must have already
    # been done for hhis to work
    def J(self,which_M,t):
        # to get the solution at a particular time we need the index that is assosiated witht that
        # time. we get this by taking the time value wanted and deviding by dt. In order for this to
        # work with single (non array values) of time we need a self,dt to be defined.
        x1 = self.sol[int(t/self.dt+.5),2]
        y  = self.sol[int(t/self.dt+.5),3]
        #print(t/self.dt)
        # define the matrix elements of the time dependent jacobian
        M11 = 0.0
        M12 = 1.0
        #M21 = self.coef*pl.cos(x1)*pl.cosh(y)*pl.cos(t)*(pl.cos(2.0*x1)+pl.cosh(2.0*y)-2.0)/(pl.cos(x1)**2-pl.cosh(y)**2)**2
        M21 = -2.0*self.coef*pl.cos(x1)*pl.cosh(y)*pl.cos(t)*(pl.cos(x1)**2-pl.cosh(y)**2+2.0*pl.sin(x1)**2)/(pl.cos(x1)**2-pl.cosh(y)**2)**2
        M22 = -self.drg

        if (which_M == "M11"):      
            return M11
        if (which_M == "M12"):      
            return M12
        if (which_M == "M21"):      
            return M21
        if (which_M == "M22"):      
            return M22

    def mw(self,warr,t):
        # to get the solution at a particular time we need the index that is assosiated witht that
        # time. we get this by taking the time value wanted and deviding by dt. In order for this to
        # work with single (non array values) of time we need a self,dt to be defined.

        dotW11 = warr[2]
        dotW12 = warr[3]
        dotW21 = warr[0]*self.J("M21",t)+warr[2]*(-self.drg)
        dotW22 = warr[1]*self.J("M21",t)+warr[3]*(-self.drg)
        return [dotW11,dotW12,dotW21,dotW22]
# functions looks to see weather or not the curent point is in the threshold radius of the first
# point
# returns True if NOT in threshhold radius
# returns False if we found our guy
def not_close(first_pnt,curnt_pnt,thresh):
    rf = pl.array([first_pnt[0] , first_pnt[2]])
    rs = pl.array([curnt_pnt[0] , curnt_pnt[2]])
    diff = rf-rs
    r = pl.sqrt(diff[0]**2+diff[1]**2)
    print("r is: "+str(r))

    if (r>thresh):
        return True
    else:
        return False


# find a single loop of the limit cycle. Might be periodoc over more than one cycle
# returns the solution of just that loop AND the periodicity of the loop
# takes a threshhold number. If it cant find a solution where the begining and end of the
# trajectroy lye within this threshold value than it quits and prints an error
#thresh is a distance in the phase plane
def find_one_full_closed(sol,thresh,dt):
    not_found = False
    # work our way backwards from last time value to find last period

    # first find last %2*pi position
    loc = len(sol[:,2])
    while ((loc*dt)%(2*pl.pi)>dt):
        loc-=1
    first_loc = loc 
    first_pnt = sol[first_loc,:]
    loc-=1
    # now find the next point where the orbit closes (going backward) 
    # orbits should have trajectories in multiples of 2*pi so only check those
    while ((loc*dt)%(2*pl.pi)>dt):
        loc-=1

    curnt_pnt = sol[loc,:]

    # for "slow" trajectories the point after the first may be within the threshold value. This is
    # not bad as it means the time step is definetely small enough but is messes up the next loop.
    # To fix this problem we will subtract more than one from loc. Not to much though otherwise we
    # risk crossing some 2*pi barier...though probably not. (left origonal loc-=1 for comparison).
    #loc -= 1
    # increas by pi/4
    loc -= int(pl.pi/4.0/dt)
    while (not_close(first_pnt,curnt_pnt,thresh)):
        if (loc == 0):
            print("Point in threshold not found!!")
            not_found = True
            #raise Exception("Point in threshold not found!!")
            break
        while ((loc*dt)%(2*pl.pi)>dt):
            loc-=1
        curnt_pnt = sol[loc,:]
        secnd_loc = loc
        loc-=1
    

    secnd_pnt = curnt_pnt

    if not_found:
        final = find_one_full_closed(sol,thresh+.003,dt)
    else:
        final = sol[secnd_loc:first_loc+1,:]
    
    return final
def main():

    dt = .002 
    # total number of iterations to perform
    totIter = 50000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    
    surf = 1.0
    coef = .202
    k = 1.0
    w = 1.0
    damp = .1
    g = .1

    # how many cells is till periodicity use x = n*pi/k (n must be even #) modNum = 2*pl.pi/k
    modNum = 2.0*pl.pi
    
    # initial conditions
    initx = pl.pi+.5
    inity = 1.0
    initvx = 0
    initvy = 0.0
    
    A = coef
    A_max = .275
    A_step = .00001
    
    count = 0

    # make arrays to keep eigen values. There willl be two eigen values so lets hve two seperate
    # arrays for them
    eigs1 = pl.array([])
    eigs2 = pl.array([])
    
    first_loop = pl.array([])
#7    while A < A_max:
    x0 = pl.array([initvx,initvy,initx,inity])
    for i in range(2):
        # initial conditions vector
        # set up: [xdot,ydot,x,y]
        apx = surfCentreLineApx(A,k,w,damp,dt)
        sol = odeint(apx.f,x0,time)
        
        sol[:,2]=sol[:,2]%(2*pl.pi)
        # find a single loop of the limit cycle. Might be periodoc over more than one cycle
        # returns the solution of just that loop AND the periodicity of the loop
        # takes a threshhold number. If it cant find a solution where the begining and end of the
        # trajectroy lye within this threshold value than it quits and prints an error
        #thresh is distance in the phase place
        thresh = .005
        #thresh = .005
        loop = find_one_full_closed(sol,thresh,dt)

        if "stop" in loop:
            break

        loop_t = pl.arange(0.0,(len(loop))*dt,dt)
        print(loop_t[-1])
        if i == 0:
            first_loop = loop[:,2]

        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.scatter([0.0,pl.pi,2.0*pl.pi],[0.0,0.0,0.0],color="Red")
        ax.plot(loop[:,2],loop[:,0])
        ax.set_xlabel("$x$",fontsize=25)
        ax.set_ylabel("$\dot{x}$",fontsize=25)
        fig.tight_layout()
        fig.savefig("LoopImgs/"+str(A)+".png")
        #os.system("open LoopImgs/" +str(A)+".png")
       

        apx.set_sol(loop)
    
        w0 = pl.array([1.0,0.0,0.0,1.0])
        w_of_t = odeint(apx.mw,w0,loop_t,hmax=dt,hmin=dt)
        #w_of_t = odeint(apx.mw,w0,loop_t)
        
        # print the period of the orbit we are working on
        print("q: " + str(loop_t[-1]/(2.0*pl.pi)))
    
        # make the matrix form of w_of_t
        matrix = w_of_t[-1,:].reshape(2,2)
        
        # use linalg to get the eigen values of the W(t=q) where q is the period time of the orbit
        vals,vect = numpy.linalg.eig(matrix) 
        
        eigs1 = pl.append(eigs1,vals[0])
        eigs2 = pl.append(eigs2,vals[1])



        count+=1
        x0 = loop[-1,:]
        A += A_step
        print(A)


    theta = pl.arange(0,10,.05)
    A_arr = pl.arange(coef,A,A_step)

    fig1 = pl.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(loop_t,first_loop)
    ax1.set_xlabel("time",fontsize=25)
    ax1.set_ylabel("firstloop",fontsize=25)
    fig1.tight_layout()
    fig1.savefig("first_loop.png")
    os.system("open first_loop.png")

    fig2 = pl.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(pl.arange(0,len(sol),1),sol[:,2])
    ax2.set_xlabel("time",fontsize=25)
    ax2.set_ylabel("sol x",fontsize=25)
    fig2.tight_layout()
    fig2.savefig("next_sol.png")
    os.system("open next_sol.png")


if __name__ == '__main__':
    main()
