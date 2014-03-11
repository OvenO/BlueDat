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

        dotW11 = warr[0]*self.J("M11",t)+warr[2]*self.J("M12",t)
        dotW12 = warr[1]*self.J("M11",t)+warr[3]*self.J("M12",t)
        dotW21 = warr[0]*self.J("M21",t)+warr[2]*self.J("M22",t)
        dotW22 = warr[1]*self.J("M21",t)+warr[3]*self.J("M22",t)
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
    loc -= int(pl.pi/4.0/dt)
    while (not_close(first_pnt,curnt_pnt,thresh)):
        if (loc == 0):
            print("Point in threshold not found!!")
            not_found =True
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

    dt = .001
    # total number of iterations to perform
    totIter = 500000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    
    surf = 1.0
    coef = .21
    k = 1.0
    w = 1.0
    damp = .1
    g = .1
# how many cells is till periodicity use x = n*pi/k (n must be even #) modNum = 2*pl.pi/k
    modNum = 2.0*pl.pi
    
    # initial conditions
    initx = 3.4653
    inity = 1.0
    initvx = -0.0242
    initvy = 0.0
# initial conditions vector
    # set up: [xdot,ydot,x,y]
    x0 = pl.array([initvx,initvy,initx,inity])
    apx = surfCentreLineApx(coef,k,w,damp,dt)
    sol = odeint(apx.f,x0,time)
    
    sol[:,2]=sol[:,2]%(2*pl.pi)
    # find a single loop of the limit cycle. Might be periodoc over more than one cycle
    # returns the solution of just that loop AND the periodicity of the loop
    # takes a threshhold number. If it cant find a solution where the begining and end of the
    # trajectroy lye within this threshold value than it quits and prints an error
    #thresh is distance in the phase place
    thresh = .005
    loop = find_one_full_closed(sol,thresh,dt)
    loop_t = pl.arange(0.0,(len(loop))*dt,dt)
    apx.set_sol(loop)

    w0 = pl.array([1.0,0.0,0.0,1.0])
    w_of_t = odeint(apx.mw,w0,loop_t,hmax=dt,hmin=dt)
    #w_of_t = odeint(apx.mw,w0,loop_t)
    print("len w_of_t: " + str(len(w_of_t)))

    # make the matrix form of w_of_t
    matrix = w_of_t[-1,:].reshape(2,2)
    
    # use linalg to get the eigen values of the W(t=q) where q is the period time of the orbit
    print(numpy.linalg.eig(matrix))
    # test the w solution matrix by refinding the sololution with it
    initial_x  = loop[0,2] 
    initial_vx = loop[0,0]
    ut = w_of_t[:,0]*initial_x + w_of_t[:,1]*initial_vx
    vt = w_of_t[:,2]*initial_x + w_of_t[:,3]*initial_vx

    fig4 = pl.figure()
    ax4 = fig4.add_subplot(111)
    #ax1.plot(sol[-int(.5*len(sol[:,0])):,2],sol[-int(.5*len(sol[:,0])):,0])
    ax4.plot(ut,vt)
    ax4.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
    #ax1.set_title("Time Slice")
    #ax1.legend(loc = "best")
    ax4.set_xlabel("$x$",fontsize="30")
    ax4.set_ylabel("$\dot{x}$",fontsize="30")
    fig4.savefig("try_w_sol.png")
    
    os.system("open try_w_sol.png")

    fig2,ax2 = pl.subplots(2,sharex=True)
    #ax2[0].plot(loop_t,ut+loop[:,2])
    #ax2[1].plot(loop_t,vt+loop[:,0])
    ax2[0].plot(loop_t,ut)
    ax2[1].plot(loop_t,vt)
    ax2[1].set_xlabel("$t$")
    fig2.savefig("test_sol_m.png")
    os.system("open test_sol_m.png")


    fig3,axarr = pl.subplots(4,sharex=True)
    axarr[0].plot(pl.arange(0,len(w_of_t)),w_of_t[:,0])
    axarr[1].plot(pl.arange(0,len(w_of_t)),w_of_t[:,1])
    axarr[2].plot(pl.arange(0,len(w_of_t)),w_of_t[:,2])
    axarr[3].plot(pl.arange(0,len(w_of_t)),w_of_t[:,3])

    axarr[0].set_ylabel("0")
    axarr[1].set_ylabel("1")
    axarr[2].set_ylabel("2")
    axarr[3].set_ylabel("3")
    fig3.savefig("w_of_t.png")
    os.system("open w_of_t.png")


    print("loopj_t[-1]: " +str(loop_t[-1]))
    print("loop_t")
    print(len(loop_t)) 
    print(len(loop))
    
    fig1 = pl.figure()
    ax1 = fig1.add_subplot(111)
    #ax1.plot(sol[-int(.5*len(sol[:,0])):,2],sol[-int(.5*len(sol[:,0])):,0])
    ax1.plot(loop[:,2],loop[:,0])
    ax1.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
    #ax1.set_title("Time Slice")
    #ax1.legend(loc = "best")
    ax1.set_xlabel("$x$",fontsize="30")
    ax1.set_ylabel("$\dot{x}$",fontsize="30")
    fig1.savefig("loop_plot.png")
    
    os.system("open loop_plot.png")    


    print("w matrix")
    print([w_of_t[-1,0],w_of_t[-1,1]])
    print([w_of_t[-1,2],w_of_t[-1,3]])

    #for a in range(len(sol[:,0])):
    #    sol[a,2] = sol[a,2]%modNum
    #    sol[a,3] = sol[a,3]+.063
    #    fig1 = pl.figure()
    #    ax1 = fig1.add_subplot(111)
    #    ax1.scatter(sol[a,2],abs(sol[a,3]-1.0)+1.0,marker = "o", s=100,color="Grey")
    #    surf_arr = pl.arange(0.0,10,.2)
    #    pl.axis([0.0,2.0*pl.pi,-.5,3.0])
    #    ax1.fill_between(surf_arr,pl.zeros([len(surf_arr)])+1,pl.zeros([len(surf_arr)])-.5,where=None,color="Black")
    #    ax1.set_xlabel("$x$",fontsize=36)
    #    ax1.set_ylabel("$y$",fontsize=36)
    #    ax1.tick_params(axis="both",labelsize=15)
    #    ax1.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],s=70,color = "White", marker="o",label="Electrodes")
    #    pl.savefig(str(a)+".png")

   
    dat_file = open("data.txt","w") 
    for i in range(len(sol[:,2])):
        dat_file.write(str(sol[i,0])+" "+str(sol[i,1])+" "+str(sol[i,2])+" "+str(sol[i,3])+"\n")
    dat_file.close()

        
        
    # make text file with all extra information
    outFile = open("info.dat","w")
    outFile.write("Info \n coefficient: " + str(coef) \
            + "\nwave number: " +str(k)\
            + "\nomega: " + str(w)\
            + "\ndamping: " + str(damp)\
            + "\ng: " + str(g)\
            + "\ntime step: " + str(dt)\
            + "\ntotal time: " + str(dt*totIter)\
            + "\ntotal iterations: " + str(totIter)\
            + "\nInitial Conditions: \n" +
            "initial x: " +str(initx) \
            +"\ninitial y: " +str(inity) \
            +"\ninitial vx: " +str(initvx)\
            +"\ninitial vy: " +str(initvy) )
    outFile.close()


if __name__ == '__main__':
    main()
