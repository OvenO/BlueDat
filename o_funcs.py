#!/usr/bin/python
import os
import scipy as sp
# An assortment of usful functions for data analysis

#***********************************************************************************************
#***********************************************************************************************


#***********************************************************************************************
#***********************************************************************************************
#Big ol' class for system parameters so a given analysis doesn't have to wory about variable passing
class anal_run():

    def __init__(self):

        self.list_dir = []
        self.together = False

    def get_info(self):
        # need to go into one of the directories to get the system info if we are doing 'together'.
        if self.together: 
            we_not_in = True
            where_in_list = 0
            list_dir = os.listdir('.')
            while we_not_in:
                if 'Block_' in list_dir[where_in_list]:
                    os.chdir(os.listdir('.')[where_in_list])
                    we_not_in = False
                where_in_list += 1

        self.qq,\
            self.dt,\
            self.beta,\
            self.A,\
            self.cycles,\
            self.N,\
            self.x_num_cell,\
            self.y_num_cell,\
            self.order,\
            self.sweep_str,\
            self.Dim,\
            self.O_mega,\
            self.system = get_system_info()

    # returns list of  poindat files if we are doing normal run otherwise...
    # The purpos of this function is to grab all the path names of all the files in different block
    # directories that all have the same N. We need to be slightly carful using this because several
    # things are based off of this assuption that N is the sam and the number of point of the sweep
    # variable are all the same.
    # Also this assums that all the directories in the given directory (d) is stuff you want in the
    # plot.
    # Function returns a list of all the paths.
    def set_list_dir(self):
            
        # need to get back out for later on
        if self.together: 
            os.chdir('..')

            print('going deep')
    
            # get the list of all the Block_... directories that contain the fiels with the data
            all_files = os.listdir(".")
            all_block_dir = []
            # filter out non Block* files
            for i,j in enumerate(all_files):
                if 'Block_' in j:
                    all_block_dir.append(j)
                    
            # This is needs to be a little clever to avoid going into each one.
            # by assuming that they are all have the name number of divisions of the sweep parameter we can
            # just count the number of data fiels in one Block directory and then we know then names of
            # rest.
            print('counting sweep parameter divisions from directory: ' +all_block_dir[0])
            # count the number of poindat fiels
            num_runs = 0
            for i,j in enumerate(os.listdir(all_block_dir[0])):
                if "poindat" in j:
                    num_runs+=1
            print('num_runs is: '+str(num_runs))
    
            # initialize the list to store the final path names is
            list_dir = [] 
    
            for a,b in enumerate(all_block_dir):
                for l in range(num_runs):
                    list_dir.append(b+'/'+str(l)+'poindat.txt')

        else:
            all_dir = os.listdir(".")
            list_dir = [] 

            num_runs = 0
            for i,j in enumerate(all_dir):
                if ("poindat" in j) and ('.gz' not in j):
                    list_dir.append(j)
                    num_runs+=1
            print('num_runs is: '+str(num_runs))
            self.num_runs = num_runs
       
        self.list_dir = list_dir



#***********************************************************************************************
#***********************************************************************************************
# This returns sliced data that is sliced at the zero potential poincare section (zpps)
def get_zpps(sol,Dim,N,dt):
    # slice the data so we only have data for values of t=pi(2*n + 1/2)
    new_data = sp.array([])
    for i in range(len(sol)):
        # This is getting values of time that are at maximum potentials!!! WRONG
        # check_time = i*dt%(pl.pi*2.0)
        check_time = (i*dt+sp.pi/2.0)%(sp.pi*2.0)
        if check_time < dt and check_time > 0.0:
            new_data = sp.append(new_data,sol[i,:])
    
    return new_data.reshape(-1,Dim*2*N)


#***********************************************************************************************
#***********************************************************************************************
# This returns sliced data that is sliced at the t=0 (Max potential) poincare section (mpps)
def get_mpps(sol,Dim,N,dt):
    # slice the data so we only have data for values of t=pi*2*n)
    new_data = sp.array([])
    for i in range(len(sol)):
        # This is getting values of time that are at maximum potentials!!! WRONG
        # check_time = i*dt%(pl.pi*2.0)
        check_time = (i*dt)%(sp.pi*2.0)
        if check_time < dt and check_time > 0.0:
            new_data = sp.append(new_data,sol[i,:])
    
    return new_data.reshape(-1,Dim*2*N)


#***********************************************************************************************
#***********************************************************************************************

# This function taks a multiparticle solution as an input. It returns a solution of the same form
# but where the trajectory with the smallest average squared velocity is at the centre of the plot.
# We are only going to be looking at the last quarter of the trajectory. There is no reason to look
# at more than that.
# For now this function is limited in that it only works for runs in a single unit cell.
# Also this function is limited in that it only works for 1D.
def center_orbit(sol,N,Dim):
    how_much = 20
    #sol[:,N:(2*N)]=sol[:,N:(2*N)]%(2.0*sp.pi)

    print('looking for orbit for with smallest velocity amplitude. This only works for 1D systems')

    # make an array of the mean velocity squared values. The array will be orderd as the solution is
    # ordered
    mean_vel_arr = sp.array([])
    # we also need an array of mean_positions to figure out where these things are
    mean_pos_arr = sp.array([])
    count = 0
    while count < N:
        mean_vel_arr = sp.append(mean_vel_arr,(sol[(-len(sol)/how_much):,count]**2).mean())
        print('mean vel appended: ' +str((sol[(-len(sol)/how_much):,count]**2).mean()))

        # so this is a little trickey... 
        # We also need the standard deveation becuase if the paticle is oscilating at the boundary
        # so it is right around 0 AND 2pi then its mean position is pi. We can use the standard
        # deviation to tell weather or not it is actualy at pi or at the boundary. The standard
        # deveation should certinaly be less than pi/2 unless there is only 1 particle.
        cur_mean_pos = (sol[(-len(sol)/how_much):,N+count]).mean()
        cur_mean_std = (sol[(-len(sol)/how_much):,N+count]).std()
        if (cur_mean_std > 3.0):
            print('foud particle oscilating at boundary')
            print('standard deviation: ' +str(cur_mean_std))
            cur_mean_pos = 0.0
        mean_pos_arr = sp.append(mean_pos_arr,cur_mean_pos)
        print('mean pos appended: ' +str(cur_mean_pos))

        count+=1

    print('mean pos arr: ' +str(mean_pos_arr))
    print('mean vel arr: ' +str(mean_vel_arr))

    # which particle is the one with the smallest v^2 trajectory? the_one will be the index of this
    # particle
    the_one = sp.argmin(mean_vel_arr)
    print('orbit with smallest velocity amplitued: '+str(the_one))
    print('mean vel of the_one: ' +str(mean_vel_arr[the_one]))
    print('mean pos of the_one: ' +str(mean_pos_arr[the_one]))

    # Now we need to shift everything to get it into the center. 
    # there are a few ways to do this. We are going to try this one but it might not be the best one

    shift = sp.pi-mean_pos_arr[the_one]
    # now shift everything by the right amount
    sol[:,N:] = (sol[:,N:]+shift)%(2.0*sp.pi)

    return sol

#    # first devide the system (lenght 2*pi) by N. Then shift system by some integer amount of
#    # (2*pi)/N. that gets the particle of interest at the center
#    divs = 2.0*sp.pi/N
#    # We can find the integer to shif by with by making an array of divs going up. subract from this
#    # array the mean position of the particle -> index of the minimum value is going to be the
#    # integer to shift by
#
#    divs_arr = sp.arange(0,2.0*sp.pi,divs)
#    dist_arr = divs_arr - mean_pos_arr[the_one]
#    print('dist_arr: '+str(dist_arr))
#    int_dist = sp.argmin(dist_arr)
#    
#    print('int_dist: ' +str(int_dist))
#    print('divs: ' +str(divs))
#    print('int_dist*divs: ' + str(int_dist*divs))
#    shift = sp.pi-int_dist*divs
#    print('shif distance: ' +str(shift))
#
#    # now shift everything by the right amount
#    sol[:,N:] = (sol[:,N:]+shift)%(2.0*sp.pi)
#
#    return sol

#***********************************************************************************************
#***********************************************************************************************
# This is going to run a floquet stability analysis.
# 2 different aproaches:
#   1) find full loop of trajectory numericaly (kinda a crap shoot for dense large period orbits)
#   then run floquet analysis using the loop we found. continue for increasing A using small.
#       a) must use small chages of A to keep the right orbit
#       b) must use small dt so the time value of output (to become input) is as close to n*2*pi as
#       posible
#   2) find fixed points in PC section. fit line to pradict where unstable fixed point is. do
#   stability analysis of the period of the stable orbit on the unstable fixed point.

#***********************************************************************************************
#***********************************************************************************************
# this needs the dimesion of the system in question so it can interpret the file corectly
def get_system_info():

    # define a variable that can contain the type of system
    system = 'not defined'

    print('getting system info')
    # depending on the variable we want we need a number that controles how we slice the data
    # lines. With the dimension of the system we can slice everything right. 
    # in order to get info files there must be poindat.txt files in the same directory. Grab one at
    # random to get the shape of the data SO WE CAN FIND THE DIMESTION
    for i,j in enumerate(os.listdir('.')):
        print('grabing file to find dimesion of system')
        if 'poindat.txt' in j:
            print('found file: '+j)
            the_file = open(j,'r')
            # get rid of first line
            the_file.readline()
            data = sp.genfromtxt(the_file) 
            the_file.close()
            break 

    # get the  system info
    info_f = open('info.txt','r')
    l = info_f.readlines()

    looking_for=['qq','dt:','beta','A ','cycles','particle number','x_num_cell','y_num_cell','order','O_mega']
    values =    [ 0.0, 0.0, 0.0  ,0.0, 0.0    , 0.0             , 0.0        , 0.0        , 0.0   ,0.0]
    # make sure to type cast all of these right later (above)

    print('going into enumeration of info file')
    for i,j in enumerate(l):
        for a,b in enumerate(looking_for):
            if b in j:
                if 'O_mega' in j:
                    system = 'OurHMF'
                if 'sweep' in j:
                    values[a] = 'sweep'
                else:
                    # in the case where values[a] is a string -> keep string because its a note
                    if type(values[a]) == str: continue
                    else: 
                        values[a] = float(j.split()[-1])

    qq         = values[0]
    dt         = values[1]
    beta       = values[2]
    A          = values[3]
    cycles     = values[4]
    N          = values[5]
    x_num_cell = values[6]
    y_num_cell = values[7]
    order      = values[8]
    O_mega     = values[9]

    # type cast 
    if type(qq        ) != str: qq         = float(qq        )
    else: sweep_str = r'$q_i q_j$'
    if type(dt        ) != str: dt         = float(dt        )
    else: sweep_str = r'$\delta t$'
    if type(beta      ) != str: beta       = float(beta      )
    else: sweep_str = r'$\beta$'
    if type(A         ) != str: A          = float(A         )
    else: sweep_str = r'$A$'
    if type(cycles    ) != str: cycles     = int(  cycles    )
    else: sweep_str = '$Number of Cycles$'
    if type(N         ) != str: N          = int(  N         )
    else: sweep_str = r'$N$'
    if type(O_mega    ) != str: O_mega     = float(O_mega    )
    else: sweep_str = r'$\Omega$'
    if type(x_num_cell) != str: x_num_cell = int(  x_num_cell)
    if type(y_num_cell) != str: y_num_cell = int(  y_num_cell)
    if type(order     ) != str: order      = int(  order     )
    else: sweep_str = r'$O$'

    print('qq is: ' +str(qq))
    print('dt is: '+str(dt))
    print('beta is: '+str(beta))
    print('A is: '+str(A))
    print('cycles is: '+str(cycles))
    print('N is: ' + str(N))
    print('O_mega is: ' + str(O_mega))
    print('x_num_cell is: '+str(x_num_cell))
    print('y_num_cell is: '+str(y_num_cell))
    print('order is: '+str(order))
    print('sweep_str is: '+str(sweep_str))

    # need N to calculate Dim so this needs to be down here
    Dim = (sp.shape(data)[1])/(2*N)
    print('dimension: '+str(Dim))

    # take care of the fact that 1D sytems num_cell IS NOT x_num_cell and therfore not found
    if Dim ==1:
        y_num_cell = ' No y '
        order = 'polygamma'
        for i,j in enumerate(l):
            if 'num_cell' in j:
                print('got num_cell 1D: ' +str(int(float(j.split()[-1]))))
                x_num_cell = int(float(j.split()[-1]))

    print('returning informatin')
    
    return qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim,O_mega,system


