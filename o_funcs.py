#!/usr/bin/python
import os
import scipy as sp
# An assortment of usful functions for data analysis

# This is going to run a floquet stability analysis.
# 2 different aproaches:
#   1) find full loop of trajectory numericaly (kinda a crap shoot for dense large period orbits)
#   then run floquet analysis using the loop we found. continue for increasing A using small.
#       a) must use small chages of A to keep the right orbit
#       b) must use small dt so the time value of output (to become input) is as close to n*2*pi as
#       posible
#   2) find fixed points in PC section. fit line to pradict where unstable fixed point is. do
#   stability analysis of the period of the stable orbit on the unstable fixed point.

# this needs the dimesion of the system in question so it can interpret the file corectly
def get_system_info():
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

    looking_for=['qq','dt','beta','A','cycles','particle number','x_num_cell','y_num_cell','order']
    values =    [ 0.0, 0.0, 0.0  ,0.0, 0.0    , 0.0             , 0.0        , 0.0        , 0.0   ]
    # make sure to type cast all of these right later (above)

    print('going into enumeration of info file')
    for i,j in enumerate(l):
        for a,b in enumerate(looking_for):
            if b in j:
                if 'sweep' in j:
                    values[a] = 'sweep'
                else:
                    # in the case where values[a] is a string -> keep string because its a note
                    if type(values[a]) == str: continue
                    else: values[a] = float(j.split()[-1])

    qq         = values[0]
    dt         = values[1]
    beta       = values[2]
    A          = values[3]
    cycles     = values[4]
    N          = values[5]
    x_num_cell = values[6]
    y_num_cell = values[7]
    order      = values[8]
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
    if type(x_num_cell) != str: x_num_cell = int(  x_num_cell)
    if type(y_num_cell) != str: y_num_cell = int(  y_num_cell)
    if type(order     ) != str: order      = int(  order     )
    else: sweep_str = r'$O$'

    # need N to calculate Dim so this needs to be down here
    Dim = (sp.shape(data)[1])/(2*N)
    print('dimension: '+str(Dim))

    if Dim ==1:
        y_num_cell = ' No y '
        order = 'polygamma'
        for i,j in enumerate(l):
            if 'num_cell' in j:
                print('got num_cell 1D: ' +str(int(float(j.split()[-1]))))
                x_num_cell = int(float(j.split()[-1]))

    print('returning informatin')
    return qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim

