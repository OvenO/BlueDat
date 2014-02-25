import os
import pylab as pl
from scipy.integrate import odeint
from scipy.integrate import ode
import numpy

# this is just to look at the stability of the pi fixed point

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
        #x1 = self.sol[int(t/self.dt+.5),2]
        #y  = self.sol[int(t/self.dt+.5),3]
        x1 = pl.pi
        y = 1.0
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

        #dotW11 = warr[0]*self.J("M11",t)+warr[2]*self.J("M12",t)
        #dotW12 = warr[1]*self.J("M11",t)+warr[3]*self.J("M12",t)
        #dotW21 = warr[0]*self.J("M21",t)+warr[2]*self.J("M22",t)
        #dotW22 = warr[1]*self.J("M21",t)+warr[3]*self.J("M22",t)
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
            raise Exception("Point in threshold not found!!")
        while ((loc*dt)%(2*pl.pi)>dt):
            loc-=1
        curnt_pnt = sol[loc,:]
        secnd_loc = loc
        loc-=1

    secnd_pnt = curnt_pnt
     
    # just error checking for now...but secnt_pnt and first_pnt should be about the same so print
    # and see:
    return sol[secnd_loc:first_loc+1,:]
def main():
    is_transparent = False 
    # this variable just exsits so we dont print the A value of the bifurca
    # once.
    found_bif = False
    dt = .001 
    # total number of iterations to perform
    totIter = 50000
    totTime = totIter*dt

    time = pl.arange(0.0,totTime,dt)

    # initial conditions
    initx = pl.pi
    inity = 1.0
    initvx = 0.0
    initvy = 0.0
    
    surf = 1.0
    coef = .1
    k = 1.0
    w = 1.0
    damp = .5
    g = .1

    # how many cells is till periodicity use x = n*pi/k (n must be even #) modNum = 2*pl.pi/k
    modNum = 2.0*pl.pi
   
    A = coef
    A_max = .30
    A_step = .001
    
    count = 0

    # make arrays to keep eigen values. There willl be two eigen values so lets hve two seperate
    # arrays for them
    eigs1 = pl.array([])
    eigs2 = pl.array([])

    while A < A_max:
        # initial conditions vector
        # set up: [xdot,ydot,x,y]
        x0 = pl.array([initvx,initvy,initx,inity])
        apx = surfCentreLineApx(A,k,w,damp,dt)
        #sol = odeint(apx.f,x0,time)
        
        #sol[:,2]=sol[:,2]%(2*pl.pi)
        # find a single loop of the limit cycle. Might be periodoc over more than one cycle
        # returns the solution of just that loop AND the periodicity of the loop
        # takes a threshhold number. If it cant find a solution where the begining and end of the
        # trajectroy lye within this threshold value than it quits and prints an error
        #thresh is distance in the phase place
        thresh = .005
        #loop = find_one_full_closed(sol,thresh,dt)
        
        #loop for the 0 fixed point is just 0,0,0,0,0,0,0

        loop_t = pl.arange(0.0,2.0*pl.pi,dt)
         
        w0 = pl.array([1.0,0.0,0.0,1.0])
        w_of_t = odeint(apx.mw,w0,loop_t,hmax=dt,hmin=dt)
        #w_of_t = odeint(apx.mw,w0,loop_t)
        print("len w_of_t: " + str(len(w_of_t)))
    
        # make the matrix form of w_of_t
        matrix = w_of_t[-1,:].reshape(2,2)

        # print the determinante of the matrix. should be constant value (.2846)
        print("determinant W")
        print(matrix[0,0]*matrix[1,1]-matrix[1,0]*matrix[0,1])
        
        # use linalg to get the eigen values of the W(t=q) where q is the period time of the orbit
        vals,vect = numpy.linalg.eig(matrix) 
       
        if((abs(vals[0])>=1.0) & (not found_bif)):
            print("this is the bifurcation point (l1)")
            print(A)
            found_bif = True
        if((abs(vals[1])>=1.0) & (not found_bif)):
            print("this is the bifurcation point (l2)")
            print(A)
            found_bif = True
        eigs1 = pl.append(eigs1,vals[0])
        eigs2 = pl.append(eigs2,vals[1])

        count+=1
#        x0 = loop[-1,:]
        A += A_step
        print(A)


    theta = pl.arange(0,10,.01)


    fig1 = pl.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(pl.cos(theta),pl.sin(theta),color = "Black")
    ax1.plot([k.real for k in eigs1],[l.imag for l in eigs1],color = "Black")
    ax1.set_xlabel("Re[$\lambda_1$]",fontsize=25)
    ax1.set_ylabel("Im[$\lambda_1$]",fontsize=25)
    fig1.tight_layout()
    fig1.savefig("pi_eig1.png",dpi=300,transparent=is_transparent)
    os.system("open pi_eig1.png")

    fig2 = pl.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(pl.cos(theta),pl.sin(theta),color = "Black")
    ax2.plot([k.real for k in eigs2],[l.imag for l in eigs2],color = "Black")
    ax2.set_xlabel("Re[$\lambda_2$]",fontsize=25)
    ax2.set_ylabel("Im[$\lambda_2$]",fontsize=25)
    fig2.tight_layout()
    fig2.savefig("pi_eig2.png",dpi=300,transparent=is_transparent)
    os.system("open pi_eig2.png")
    
    A_arr = pl.arange(coef,A_max,A_step)
    print("length of As and eigs")
    print(len(A_arr))
    print(len([k.real for k in eigs1]))
    while len(A_arr)>len([k.real for k in eigs1]):
        A_arr = A_arr[:-1]
    while len(A_arr)<len([k.real for k in eigs1]):
        A_arr = pl.append(A_arr,A_arr[-1]+A_step)
    print(len(A_arr))
    print(len([k.real for k in eigs1]))

    fig3, ax3 = pl.subplots(2,sharex=True)
    ax3[0].plot(A_arr,[k.real for k in eigs1],color = "Black")
    ax3[1].plot(A_arr,[k.imag for k in eigs1],color = "Black")
    ax3[0].set_ylabel("Re[$\lambda_1$]",fontsize=25)
    ax3[1].set_ylabel("Im[$\lambda_1$]",fontsize=25)
    ax3[1].set_xlabel("$A$",fontsize=25)
    
    fig3.tight_layout()
    fig3.savefig("pi_A_vs_eig1.png",dpi=300,transparent=is_transparent)
    os.system("open pi_A_vs_eig1.png")

    fig4, ax4 = pl.subplots(2,sharex=True)
    ax4[0].plot(A_arr,[k.real for k in eigs2],color = "Black")
    ax4[1].plot(A_arr,[k.imag for k in eigs2],color = "Black")
    ax4[0].set_ylabel("Re[$\lambda_2$]",fontsize=25)
    ax4[1].set_ylabel("Im[$\lambda_2$]",fontsize=25)
    ax4[1].set_xlabel("$A$",fontsize=25)
    fig4.tight_layout()
    fig4.savefig("pi_A_vs_eig2.png",dpi=300,transparent=is_transparent)
    os.system("open pi_A_vs_eig2.png")

    fig5, ax5 = pl.subplots(2,sharex=True)
    ax5[0].plot(A_arr,abs(eigs1),color = "Black")
    ax5[1].plot(A_arr,abs(eigs2),color = "Black")
    ax5[0].set_ylabel("$\lambda_1$",fontsize=25)
    ax5[1].set_ylabel("$\lambda_2$",fontsize=25)
    ax5[1].set_xlabel("$A$",fontsize=25)
    fig5.tight_layout()
    fig5.savefig("pi_A_vs_mag_eigs.png",dpi=300,transparent=is_transparent)
    os.system("open pi_A_vs_mag_eigs.png")


    eig_file = open("pi_data.txt","w") 
    eig_file.write("eig1   eig2   A\n")
    for i in range(len(eigs1)):
        eig_file.write(str(eigs1[i])+" "+str(eigs2[i])+" "+str(A_arr[i])+"\n")
    eig_file.close()

        
        
    ## make text file with all extra information
    #outFile = open("info.dat","w")
    #outFile.write("Info \n coefficient: " + str(coef) \
    #        + "\nwave number: " +str(k)\
    #        + "\nomega: " + str(w)\
    #        + "\ndamping: " + str(damp)\
    #        + "\ng: " + str(g)\
    #        + "\ntime step: " + str(dt)\
    #        + "\ntotal time: " + str(dt*totIter)\
    #        + "\ntotal iterations: " + str(totIter)\
    #        + "\nInitial Conditions: \n" +
    #        "initial x: " +str(initx) \
    #        +"\ninitial y: " +str(inity) \
    #        +"\ninitial vx: " +str(initvx)\
    #        +"\ninitial vy: " +str(initvy) )
    #outFile.close()


if __name__ == '__main__':
    main()
