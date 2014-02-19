#!/usr/bin/python
import pylab as pl
from scipy.integrate import odeint
from scipy.integrate import ode
import time as thetime
import argparse
from datetime import datetime
#from mpl_toolkits.mplot3d import Axes3D
import random
import argparse
import o_funcs as of
import os 
import sys
sys.path.append(os.path.expanduser("~")+"/datasphere/ECproject")
import ECclass as ec

#import imp
#foo = imp.load_source('module.name', '/path/to/file.py')
#foo.MyClass()




#**********************************************************************************************
#**********************************************************************************************
# Here we will find a second point that is some small distace (epsilon) from the x0 point.
# We need this function because we want the distace to be epsilon but the position about x0 to be
# random
# Particle variable is the particle we want to purturb
# N is the number of particles 
def get_first_init(x0,epsilon,particle,N):
    x_new = pl.copy(x0)
    print('getting the first initial condition')
    print('fiducial initial: '+str(x0))
    # multi particle array layout [nth particle v, (n-1)th particle v , ..., 0th v, nth particle x, x, ... , 0th particle x]
    
    # we will use a change of coordinates to get the location of the particle relative to x0. First
    # we just find some random point a distace epsilon from the origin.
    theta = random.random()*2.0*pl.pi
    x = epsilon*pl.cos(theta)
    vx = epsilon*pl.sin(theta)

    x_new[particle] = x0[particle]+vx
    x_new[N+particle] = x0[N+particle]+x

    print('randomly purturbed initial: ' + str(x0))
    
    return x_new

#**********************************************************************************************
#**********************************************************************************************
# function just finds the magnitude of the distance of the two input points. This needs to be a
# little intelegent to make sure the periodic boundary conditions are included corectly. Lets have
# it return a boolean value specifying weather it was closer across the boundary or not.
def distance(x1,x2,particle,N):
    # see weather the particle is closer in the normal way or across the seam of the bounday.
    x_normal =  abs(x1[N+particle]-x2[N+particle])
    print('x_normal: ' + str(x_normal))
    x_period =  abs(2.0*pl.pi - x_normal)
    print('x_periodi: ' + str(x_period))
    
    if x_normal < x_period:
        #print('x_normal < x_period (distance NOT across boundry is smaller)')
        across_seam = False
        x_sqrd = x_normal**2
    if x_normal > x_period:
        #print('x_normal > x_period (distance across boundry is smaller)')
        across_seam = True
        x_sqrd = x_period**2
    
    vx_sqrd = (x1[particle]-x2[particle])**2

    dist = pl.sqrt(vx_sqrd + x_sqrd)
    
    #if dist>2.0*pl.pi:
    #    dist -= 2.0*pl.pi
    
    print('returning across_seam = ' +str(across_seam))
    return dist, across_seam
#**********************************************************************************************
#**********************************************************************************************
     # gram schmidt orthogonalize ONLY the trajectories that have grown beyond
     # first_epsilon to appropriat distances from fudicial. Set
     # all others (ones that dont grow) back to their futicial conterpart positions. All
     # together these will be iniitial conditinos for next run.

# see lab book #2 pg 59 for a
# more detailed explination. "renomalize" will retun a point that is on the line contecting the
# final points of the trajectorys and will this new point will be a distance epsilon from the
# unpurtubed trjectory (x_unpurt). x_putub is final position of purtubed trajectorie

# For a detailed description of what we are doing in the "renormalize" see notebook 2 pg 117
# (115,116 too if all particles might be chaotic but this is not neccassary)

def renormalize(x_unpurt,x_purt,total_epsilon,N):

    xnew = pl.copy(x_unpurt)
    
    # first determine which trajectories get to small
    # this variable is the amount by which a trajectory must grow otherwise it is going to get set
    # to fiducial position and prbably will not grow after that 
    check = total_epsilon/pl.sqrt(N)/2.0 
    
    # keep track of the index of particles that are "chaotic"
    chaotic = pl.array([])
    # keep track of their final distances too
    chaotic_dists = pl.array([])
    # now check
    for l in range(N):
        # what is the distace between the lth fiducial particle adn the perturbed?
        lth_distance,across_seam = distance(x_unpurt,x_purt,l,N)
        # if the distace is shorter across the seam of the boundary (across_seam = True) then set
        # the purtubed particle position to +- 2pi so the points are next to eachother. if the
        # purturbed particle position is less than the unpurturbed add 2pi. if the otherway around
        # minus 2pi
        if across_seam:
            if x_unpurt[l+N] < x_purt[l+N]:
                x_purt[l+N] = x_purt[l+N]-2.0*pl.pi
            if x_unpurt[l+N] > x_purt[l+N]:
                x_purt[l+N] = x_purt[l+N]+2.0*pl.pi
            # check new distance
            temp_dist,temp_across = distance(x_unpurt,x_purt,l,N)
            print('fixed a distance in renormalize. New dist is: ' +str(temp_dist))
        # after all this we need to re modulus the system --> at the end of this function

        if lth_distance < check:
            # set velocity to fiducial
            xnew[l] = x_unpurt[l]
            # set position to fiducial
            xnew[N+l] = x_unpurt[N+l]
        else:
            chaotic = pl.append(chaotic,l)
            chaotic_dists = pl.append(chaotic_dists,lth_distance)

    print('number of chaotic is: ' + str(len(chaotic)))

    # If there are no chaotic paticles than it is still useful to measure the
    # negative LE. To do this --> if there are no chaotic particle --> make a randome set of
    # purturbed initial conditions just like whats done in for the very first set of initial
    # conditions
    if len(chaotic)==0:
        for k in range(N):
            xnew = get_first_init(xnew,total_epsilon/pl.sqrt(N),k,N)

    else:
        # The best we can do is split the total epsilon up evenly between the chaotic particles.
        to_perturb = total_epsilon/pl.sqrt(float(len(chaotic)))
        # where len(chaotic) is the number of chaotic particles

        # enumerate chaotic for the indicies of particles to renormalize to distance to_perturb from
        # fiducial
        for n,m in enumerate(chaotic):
            # particle m will be chaotic_dists[n] from fudical
            # velocity
            xnew[m] = x_unpurt[m] + (to_perturb/chaotic_dists[n])*(x_purt[m]-x_unpurt[m])
            # position
            xnew[N+m] = x_unpurt[N+m] + (to_perturb/chaotic_dists[n])*(x_purt[N+m]-x_unpurt[N+m])

    # because of seam checking we need to re modulus the system
    xnew[N:]=xnew[N:]%(2.0*pl.pi)
    return xnew

#**********************************************************************************************
#**********************************************************************************************
# NOT GOOD FOR NEW VERSION
def plot_first_sol(sol,dt):
    
    poin = get_poin(sol,dt)
    poin[:,2]=poin[:,2]%(2*pl.pi)

    first_fig = pl.figure()
    first_ax = first_fig.add_subplot(111)
    first_ax.scatter(poin[:,2],poin[:,0],color="Black",s=.1)
    first_ax.set_xlabel("$x_1$",fontsize=25)
    first_ax.set_ylabel("$x_2$",fontsize=25)
    #first_ax.set_xlim([0,2*pl.pi])
    #first_ax.set_ylim([-1.3,1.3])
    first_fig.savefig("first_plot.png")
    os.system("open first_plot.png")
    
#**********************************************************************************************
#**********************************************************************************************
def main():

    # Structure:
    # 1) take multi particle run of interest -> this is fiducial trajectory.
    # 1a) slice the data
    # 2) throw awway transients
    # 3) take PC of fiducal as initial conditions for a seperate run 
    # 4) purturb ALL the partiicls epsilon/N
    # 5) Run puturbed version of system for a period
    # 6) measure All the distance between the fiducial trajectories and the new trajectories after
    # the one period. Store the sum of these distances and each of them individualy. FOR ALL PARTICLES. 
    # 7) move all particles back along line to fiducial particls to get expansion to get iniitial conditinos for next run.
    # 8) do this many times and average the total separations to find LE

    # 1)
    # getting the data
    parser = argparse.ArgumentParser()
    # d is for directory
    parser.add_argument('-d',action='store',dest = 'd',type = str, required = False, default = './')
    # f is for file
    parser.add_argument('-f',action='store',dest = 'f',type = str, required = False)
    # n is for the particle of interest
    #parser.add_argument('-n',action='store',dest = 'n',type = int, required = True)
    # eplsilon should be passable. THis is total epsilon
    parser.add_argument('-e',action='store',dest = 'e',type = float, required = False,default = 1e-8)


    inargs = parser.parse_args()
    d = inargs.d
    f = inargs.f
    #particle = inargs.n

    # each run's distances are sumed and stored
    sum_final_dist_sqrd_arr = pl.array([])
    # individual particles distances from fidicual are stored in a matrix
    final_dist_sqrd_arr = pl.array([])

    # get system info
    qq,dt,beta,A_sweep_str,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim = of.get_system_info()
    print("A from get_system_info() is: " +str(A_sweep_str))
    
    file_object = open(f,'r')
    first_line = file_object.readline()
    print('first line of data file: ' +str(first_line))
    A = pl.zeros(N)+float(first_line.split()[-1])
    print("A from working file is: " + str(A))

    data = pl.genfromtxt(file_object)
    print('shape of fiducial data before doing anything: ' +str(pl.shape(data)))

    # 1a)
    # first lets slice the data so that we have poincare sections of minimum field potential. This
    # should make it so slight errors in initial conditions of unpurturbed particles are less
    # prevelent. Why? Because small error in position when \Phi not zero means an error in evergy as
    # well. This will happen with the velocities (KE) but we can minimizi it by using zerro
    # potential poincare sections.

    # data for values of t=pi(2*n + 1/2) (zero potential poincare seciton)
    new_data = pl.array([])
    checked = True
    for i in range(len(data)):
        # This is getting values of time that are at makimum potentials!!! WRONG
        # check_time = i*dt%(pl.pi*2.0)
        # Right
        check_time = (i*dt+pl.pi/2.0)%(pl.pi*2.0)
        if check_time < dt and check_time > 0.0:
            new_data = pl.append(new_data,data[i,:])
            if checked:
                first_i = i
                checked=False


    sliced_data = new_data.reshape(-1,2*N)
    print('shape of sliced fiducial data before doing anything: ' +str(pl.shape(data)))

    sliced_data[:,N:] = sliced_data[:,N:]%(2.0*pl.pi)

    # for full data set just throw away the first few points so it lines up with the sliced data
    # startwise
    data = data[first_i:,:]

    # SETTING SOME OTHER VARIABLS NOW
    # the perturb distance for INDIVIDUAL particles is different from THE SYSTEM PURTERBATION
    total_epsilon = inargs.e
    first_epsilon = total_epsilon/pl.sqrt(N)
    print('total_epsilon: ' + str(total_epsilon))
    print('first_epsilon: ' + str(first_epsilon))
    # good
    #epsilon = 1.0e-10
    # works 
    # epsilon = 1.0e-10
    # 1.0e-13

    # period variable should probably just be kept at the actual period (2*pl.pi) but the purpose of
    # this variable is to alow us to changhe the lenght of time we wate before we colect the
    # distance information and renormalize. This is also nesssasary in the final calculation of the
    # LE becae LE = 1/period * ln(rm/r0).
    period = 2.0*pl.pi
    print('period is: ' + str(period))


    # time array of one period. BE CAREFULL -> must start at pi/2 and go to 
    t = pl.arange(3.0*pl.pi/2.0,3.0*pl.pi/2.0+2.0*pl.pi,dt)
    
    # 2)
    # throw awway transients
    # fiducial run is probably not very long becasue of the need for a high time resolution so lets
    # throw awway the whole first half
    sliced_data = sliced_data[(len(sliced_data[:,0])/2):,:]
    print('shape of sliced fiducial data after getting rid of transients: ' +str(pl.shape(sliced_data)))

    # 3) get our initial conditions for the purturbed run
    fiducial_start = sliced_data[0,:]
    print('First fiducial: '+str(fiducial_start))
    # 4) purturb all particles by first_epsilon
    init = pl.copy(fiducial_start)
    for i in range(N):
        init = get_first_init(init,first_epsilon,i,N)
    print('first purturbed initial conditions: '+str(init))
    # 5) Run puturbed version of system for a period
    elec = ec.Sin1D(qq,A,beta,x_num_cell)

    watch = pl.array([])
    watch_le = 0.0

    # In order to watch whats happening we are going to add some plotting stuff that we can put in
    # "if" or just comment out later
    os.mkdir('WatchingLE')
    #for i in range(1,len(sliced_data)-1):
    for i in range(1,50):
        print(i)

        purt_sol = odeint(elec.f,init,t)
        
        #print('init: ' + str(init))
        #print('sliced_data[i-1,:] = ' + str(sliced_data[i-1,:]))
        
        ## Lets see what happens when we run the fudicial trajectory again with
        #test_sol = odeint(elec.f,sliced_data[i-1,:],t)
        #test_sol[:,N:(2*N)] = test_sol[:,N:(2*N)]%(2.0*pl.pi)

        # make sure particle are in the right modulus space
        #purt_sol[:,N:(2*N)] = purt_sol[:,N:(2*N)]%(2.0*pl.pi)

        # 6) measure the distance between the fiducial trajectores and the new trajectories after the one period
        fiducial_end = sliced_data[i,:]
        purt_end = purt_sol[-1,:]
        #purt_end = test_sol[-1,:]

        first_fig = pl.figure()
        first_ax = first_fig.add_subplot(111)
        for gamma in range(N):
            first_ax.scatter(purt_sol[:,gamma+N],purt_sol[:,gamma],color="Red",s=5)
            first_ax.scatter(data[(int(2.0*pl.pi/.001)*(i-1)):(int(2.0*pl.pi/.001)*i),gamma+N],data[(int(2.0*pl.pi/.001)*(i-1)):(int(2.0*pl.pi/.001)*i),gamma],color="Blue",s=5)
            if gamma == 1:
                first_ax.annotate('start',xy=(purt_sol[0,gamma+N],purt_sol[0,gamma]),xytext=(pl.pi/2,-1.5),arrowprops=dict(facecolor='black',shrink=0.05))
                first_ax.annotate('stop' ,xy=(purt_sol[-1,gamma+N],purt_sol[-1,gamma]),xytext=(pl.pi/2,1.5)    ,arrowprops=dict(facecolor='black',shrink=0.05))
        first_ax.set_xlim([0.0,2.0*pl.pi])
        first_ax.set_ylim([-2.0,2.0])
        #first_ax.set_xlabel("$x_1$",fontsize=25)
        #first_ax.set_ylabel("$x_2$",fontsize=25)
        #first_ax.set_xlim([0,2*pl.pi])
        #first_ax.set_ylim([-1.3,1.3])
        first_fig.savefig('WatchingLE/'+str(i)+".png")
        pl.close(first_fig)


        # (for more than notes here see notebood pg 115)
        # this is for summing the distances. We need to be carful summing the distances. If we think
        # about the system as a "single particle" with 2*N*Dim degrees of freedom than the puturbed
        # distance d will be d = sqrt(x1^2 + x2^2 + x3^3 + ... xN^2). This is different than
        # sqrt(position_1^2 + vel_1^2) + sqrt(position_2^2 + vel_2^2) + ... sqrt(position_N^2 + vel_N^2)
        # So...
        # Distance Method 1) 
        # We will store the distances in phase space between individual particles of the purturbed and
        # fiducial run to know wich particles are stable. Then we can find the LE of those individual
        # particles.
        # Distance Method 2) 
        # We will also store the purturbed "distance" of the whole system from the fiducial run also
        # and calculate the LE of this.
        to_sum_sqrd_dist = 0.0
        print_str_dists = ''
        for j in range(N):
            final_dist,across_seam = distance(fiducial_end,purt_end,j,N)
            print_str_dists += str(final_dist)+' '
            # For Distance Method 2 must square final dist to get xi^2 + vi^2
            to_sum_sqrd_dist += final_dist**2
            # For Distance Method 1 just append the distance between fiducial and purtubed particl i.
            # (will reshape this arrasy later)... SQUARED? --> See notebood pg 115.
            final_dist_sqrd_arr = pl.append(final_dist_sqrd_arr,final_dist**2)


        print(print_str_dists)
        # append the summed distance
        print('to_sum_sqrd_dist: '+str(to_sum_sqrd_dist))
        sum_final_dist_sqrd_arr = pl.append(sum_final_dist_sqrd_arr,to_sum_sqrd_dist)
        
        # 7) gram schmidt orthogonalize ONLY the trajectories that have grown beyond
        # first_epsilon to appropriat distances from fudicial. Set
        # all others (ones that dont grow) back to their futicial conterpart positions. All
        # together these will be iniitial conditinos for next run.
        init = renormalize(fiducial_end,purt_end,total_epsilon,N)
        print('renormalized... new initial conditions are: ' + str(init))
        print('compair above to fiducial final position ->: ' + str(fiducial_end))
        print('totat_epsilon: ' + str(total_epsilon))
        print('one minus the other sqrd sqrted (should be total_epsilon): '+
                str(pl.sqrt(((init-fiducial_end)**2).sum())))

        watch_le += pl.log(abs(sum_final_dist_sqrd_arr[-1]/total_epsilon))
        cur_avg = watch_le/i/period
        watch = pl.append(watch,cur_avg)

    # reshape final_dist_arr
    final_dist_sqrd_arr = final_dist_sqrd_arr.reshape(-1,N)

    eps_arr = pl.zeros(len(sum_final_dist_sqrd_arr))+total_epsilon
    le = pl.log(pl.sqrt(abs(sum_final_dist_sqrd_arr))/total_epsilon)/period

    print('mean LE (LE is): ' +str(le.mean()))
    print('standard deviation LE: ' +str(le.std()))

    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.scatter(pl.arange(len(watch)),watch,s=.1)
    #ax.set_xlabel("$x_1$",fontsize=25)
    #ax.set_ylabel("$x_2$",fontsize=25)
    #ax.set_xlim([0,2*pl.pi])
    #ax.set_ylim([-1.3,1.3])
    #fig.tight_layout()
    fig.savefig("convergence_LE.png")
    os.system("open convergence_LE.png")
   



    
#OLD PROGRAM FOR REFFERENCE
#    # initial purturbation size
#    # try
#    epsilon = 1.0e-7
#    # good
#    #epsilon = 1.0e-10
#    # works 
#    # epsilon = 1.0e-10
#    # 1.0e-13
#    
#    # period variable should probably just be kept at the actual period (2*pl.pi) but the purpose of
#    # this variable is to alow us to changhe the lenght of time we wate before we colect the
#    # distance information and renormalize. This is also nesssasary in the final calculation of the
#    # LE becae LE = 1/period * ln(rm/r0).
#    period = 2.0*pl.pi
#    print('period is: ' + str(period))
#
#    print('epsilon is: ' + str(epsilon)+'\n')
#    # works
#
#    # Throw away transients of 
#    throw_away_1 = 100
#    print('throw_away_1 is: ' +str(throw_away_1)+'\n') 
#
#    # number of cycle-sets to go through
#    num = 10000
#    print('number of cycle-sets to go through is: ' + str(num))
#    # works
#    #num = 120000
#    #num = 60000
#    #num = 16000
#    #num = 30000
#    #num = 8000
#    
#    dt = .0001 
#    print('dt is: ' + str(dt))
#    # total number of sterations to perform inorder to get cycles rioht
#    totTime = period
#    time = pl.arange(0.0,totTime,dt)
#    
#    #time array for the first run to get initial condition in strange atractor
#    first_time = pl.arange(0.0,throw_away_1*2.0*pl.pi,dt)
#    print('len(first_time): ' + str(len(first_time)))
#    
#    surf = 1.0
#    coef = 1.7
#    k = 1.0
#    w = 1.0
#    damp = .1
#    g = .1
#
#    # how many cells is till periodicity use x = n*pi/k (n must be even #) modNum = 2*pl.pi/k
#    modNum = 2.0*pl.pi
#    
#    # some random initial conditions
#    initx = 3.0
#    inity = 1.0
#    initvx = .1600
#    initvy = 0.0
#
#    # now find initial condition in stange atractor by running far in time. print this to use it for
#    # later
#    x0 = pl.array([initvx,initvy,initx,inity])
#    
#    apx = surfCentreLineApx(coef,k,w,damp,dt)
#    first_sol = odeint(apx.f,x0,first_time)
#    #plot_first_sol(first_sol,dt)
#    # now reset x0
#    first_sol[:,2]=first_sol[:,2]%(2*pl.pi)
#    x0 = first_sol[-1,:]
#    
#    x_other = get_first_init(x0,epsilon)
#    
#    print("x0 from first_sol is:")
#    print(x0)
#    print("x_other (purturbed is:")
#    print x_other
#
#    # two arrays to store solution data
#    arr_orig = pl.array([])
#    arr_othr = pl.array([])
#
#    #full_orig = pl.array([])
#    #full_othr = pl.array([])
#    
#    # array to keep distance after driving cycle 
#    darr = pl.array([]) 
#    watch = pl.array([])
#    watch_le = 0.0
#
#    for i in range(num + throw_away_2):
#        sol = odeint(apx.f,x0,time)
#        sol_other = odeint(apx.f,x_other,time)
#
#        sol[:,2]=sol[:,2]%(2*pl.pi)
#        sol_other[:,2]=sol_other[:,2]%(2*pl.pi)
#
#        #arr_orig = pl.append(arr_orig,sol[-1,:])
#        #arr_othr = pl.append(arr_othr,sol_other[-1,:])
#    
#        if i> throw_away_2: 
#            #get the new distance
#            darr = pl.append(darr,distance(sol[-1,:],sol_other[-1,:]))
#            
#            watch_le += pl.log(abs(darr[-1]/epsilon))
#            cur_avg = watch_le/(i-throw_away_2+1)/period
#            watch = pl.append(watch,cur_avg)
#
#
#        #full_orig = pl.append(full_orig,sol)
#        #full_othr = pl.append(full_othr,sol_other)
#        
#        # This is to see that our orthogonalization is working and that we are in the chaotic
#        # atractor
#        #if i> throw_away_2:
#        #    poin = get_poin(sol,dt)
#        #    poin_other = get_poin(sol_other,dt)
#
#        #    fig = pl.figure()
#        #    ax = fig.add_subplot(111)
#        #    ax.scatter(sol[-1,2],sol[-1,0],s=1.0,color="Red")
#        #    ax.scatter(sol_other[-1,2],sol[-1,0],s=.5,color="Blue")
#        #    #ax.set_xlabel("$x_1$",fontsize=25)
#        #    #ax.set_ylabel("$x_2$",fontsize=25)
#        #    #ax.set_xlim([0,2*pl.pi])
#        #    #ax.set_ylim([-1.3,1.3])
#        #    fig.savefig("LyapImgs/"+str(i)+".png")
#
#        #    fig2 = pl.figure()
#        #    ax2 = fig2.add_subplot(111)
#        #    ax2.scatter([0,pl.pi,2*pl.pi],[0,0,0],s=10.0)
#        #    ax2.plot(sol[:,2],sol[:,0])
#        #    fig2.savefig("LyapImgs/"+str(i+1000)+".png")
#
#
#        x0 = sol[-1,:]
#        x_other = renormalize(x0,sol_other[-1,:],epsilon)
#    
#    #arr_orig = arr_orig.reshape(-1,4)
#    #arr_othr = arr_othr.reshape(-1,4)
#    
#    #full_orig = full_orig.reshape(-1,4)
#    #full_othr = full_othr.reshape(-1,4)
#
#    eps_arr = pl.zeros(len(darr))+epsilon
#    le = pl.zeros(len(darr))+pl.log(abs(darr/epsilon))/period
#
#    le_avg = 0.0
#    for i,j in enumerate(le):
#        le_avg += j
#    le_avg = le_avg/len(le)
#
#    print("le is")
#    print(le_avg)
#
#    fig = pl.figure()
#    ax = fig.add_subplot(111)
#    ax.scatter(pl.arange(len(watch)),watch,s=.1)
#    #ax.set_xlabel("$x_1$",fontsize=25)
#    #ax.set_ylabel("$x_2$",fontsize=25)
#    #ax.set_xlim([0,2*pl.pi])
#    #ax.set_ylim([-1.3,1.3])
#    #fig.tight_layout()
#    fig.savefig("convergence.png")
#    os.system("open convergence.png")
# 
#
#    ## plot of PC sections (red Blue)
#    #fig = pl.figure()
#    #ax = fig.add_subplot(111)
#    ##ax.scatter([0.0,pl.pi,2.0*pl.pi],[0.0,0.0,0.0],color="Red")
#    ##ax.scatter(arr_orig[:,2],arr_orig[:,0],color="Red",s=.1)
#    ##ax.scatter(arr_othr[:,2],arr_othr[:,0],color="Blue",s=.1)
#    #ax.scatter(full_orig[:,2],full_orig[:,0],color="Red",s=.1)
#    #ax.scatter(full_othr[:,2],full_othr[:,0],color="Blue",s=.1)
#    #ax.set_xlabel("$x_1$",fontsize=25)
#    #ax.set_ylabel("$x_2$",fontsize=25)
#    ##ax.set_xlim([0,2*pl.pi])
#    ##ax.set_ylim([-1.3,1.3])
#    #fig.tight_layout()
#    #fig.savefig("sep.png")
#    #os.system("open sep.png")
    
if __name__ == '__main__':
    main()
