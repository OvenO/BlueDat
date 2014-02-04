import os
import pylab as pl
import argparse
from datetime import datetime
import time as thetime
import paramiko
import random
import shutil

class Sin2D(object):
    # N -> number of particles
    # number_of -> number of blocks (usualy 500). i.e. how many steps
    # start -> parameter starting poin. 
    # stop -> parameter stoping point (step size will be calculated). 
    # cycles-> runtime cycles. 
    # num_cell -> number of unit cells to make simulation length over
    # qq -> intercharge force
    # d -> length of system as determined by the num_cell
    # A -> this is the interaction amplitude with the EC fild. It is an array so
    # that each particle may have its own interaction with the field
    def __init__(self):
        self.block_dir = ''
        self.var = 'A'
        self.script_dir= '/users/o/m/omyers/datasphere/ECproject/Sin2D/'
        self.number_of = 300
        self.start     = 0.0
        self.stop      = 2.7
        self.dt        = 0.05
        self.cycles    = 150
        self.N = 4
        self.qq = 1.0
        self.beta = .6
        self.x_num_cell = 1.0
        self.y_num_cell = 1.0
        self.x_periodic = True
        self.order = 2
        self.dx = self.x_num_cell * 2.0 * pl.pi
        self.dy = self.y_num_cell * 2.0 * pl.pi
        #self.A = pl.zeros(2*self.N) +.5
        self.A = 1.34

        # Full trajectory or just Poincare sections
        self.sliced = False
        #self.sliced = True

    def init_block(self):

        # find the amount we want to increase the parameter by
        inc_param = (self.stop-self.start)/self.number_of
    
        cur_var = self.start 
    
        get_make_dir(self)
        
        # make info file
        info_file = open(self.block_dir + '/info.txt','w')
        info_file.write('dir: '+str(self.block_dir)+'\nsurf: 1.0') 
        info_file.write('\ndt: '+str(self.dt))
        if self.var == 'beta':
            info_file.write('\nbeta (damping): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\nbeta (damping): '+str(self.beta))
        if self.var == 'cycles': 
            info_file.write('\ncycles (run time): sweep variable '+str(self.start)+'-'+str(self.stop))
        else: 
            info_file.write('\ncycles (run time): '+str(self.cycles))
        if self.var == 'N': 
            info_file.write('\nN (particle number): sweep variable '+str(self.start)+'-'+str(self.stop))
        else: 
            info_file.write('\nN (particle number): '+str(self.N))
        if self.var == 'qq': 
            info_file.write('\nqq (particle interaction strength): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\nqq (particle interaction strength): '+str(self.qq))
        if self.var == 'x_num_cell':
            info_file.write('\nx_num_cell (system length in x): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\nx_num_cell (system length in x): '+str(self.x_num_cell))
        if self.var == 'y_num_cell':
            info_file.write('\ny_num_cell (system length): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\ny_num_cell (system length): '+str(self.y_num_cell))
        if self.var == 'x_periodic':
            info_file.write('\nx_periodic (peridic in x?): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\nx_periodic (peridic in x?): '+str(self.x_periodic))
        if self.var == 'order':
            info_file.write('\norder (order of periodic wrap force): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\norder (order of periodic wrap force): '+str(self.order))
        if self.var == 'A': 
            info_file.write('\nA (interaction amplitude): sweep variable '+str(self.start)+'-'+str(self.stop))
        else: 
            if type(self.A)==float: 
                info_file.write('\nA (interaction amplitude): '+str(self.A))
            else: 
                info_file.write('\nA (interaction amplitude): '+str(min(self.A))+'-'+str(max(self.A)))
        info_file.close()

        for l in range(self.number_of):
            self.x0 = pl.zeros([4*self.N])
            
            for i,j in enumerate(self.x0):
                if i in range(2*self.N,3*self.N):
                    self.x0[i] = random.random()*self.dx
                if i in range(3*self.N,4*self.N):
                    self.x0[i] = random.random()*self.dy

            #make the file
            cur_poin_file = open(self.block_dir+'/'+str(l)+'poindat.txt','w')
    
            #write the first few lines of the file
            cur_poin_file.write(str(self.var)+' --> ' +str(cur_var)+'\n')
            
            # this just prints the numbers with one space imbetween. the .replace gets rid of the \n
            # but i'm still not really sure why they are there in the first place.
            cur_poin_file.write(str(self.x0)[1:-1].replace('\n',''))
            cur_poin_file.write('\n')
            
            cur_poin_file.close()
            cur_var += inc_param

        # now that all the file are made make an "Orig" directory in the block directory to store
        # the origonal blocks of initial conditions. This way if something goes wrong we can pull
        # them out and run again.py again by hand
        os.mkdir(self.block_dir + '/Orig')
        os.system('cp ' + self.block_dir +'/*.txt ' + self.block_dir + '/Orig/')
        os.system('scp -r ./' + self.block_dir + ' omyers@bluemoon-user2.uvm.edu:/users/o/m/omyers/Data/EC/2DBlock/Old/')
        print 'self.block_dir'
        print self.block_dir
        os.system('say files copied to cluster')  

class One_Particle_Ensbl_Sin_1D(object):
    # number_of -> number of blocks (usualy 500). start -> parameter starting poin. stop ->
    # parameter stoping point (step size will be calculated). cycles-> runtime cycles. x_,vx_dim->
    # dimensions of block of initial conditions in phase space. x_,vx_num number of points in the x and
    # vx. The <>_corner variables define the lower left corner of the block in phase space
    def __init__(self):
        self.block_dir = ''
        self.var = 'A'
        self.script_dir= '/users/o/m/omyers/datasphere/ECproject/One_Particle_Sin_1D/'
        self.number_of = 300
        self.start     = 0.7
        self.stop      = 1.375
        self.dt        = 0.05
        self.cycles    = 100
        self.beta = .6
        self.num_cell = 1.0
        self.d = self.num_cell * 2.0 * pl.pi
        self.A = .5
        self.x_dim     = 6.28
        self.vx_dim    = 6.0
        self.vx_corner = -3.0
        self.x_corner  = 0.0
        self.vy_corner = 0.0
        self.y_corner  = 0.0

        self.x_num     = 60
        self.vx_num    = 30
        # Full trajectory or just Poincare sections
        #self.sliced = False
        self.sliced = True


    def init_block(self):
        if (self.x_num == 0) or (self.x_dim==0.0):
            increment_x = 0.0
            if(self.x_num == self.x_dim):
                print('there is a descrepency between number and dimesions of block --> check to_run file')
        else:
            increment_x = self.x_dim/self.x_num 
            print('increment x: ' + str(increment_x))

        if (self.vx_num == 0) or (self.vx_dim==0.0):
            increment_vx = 0.0
            if(self.vx_num == self.vx_dim):
                print('there is a descrepency between number and dimesions of block --> check to_run file')
        else:
            increment_vx = self.vx_dim/self.vx_num 



        # find the amount we want to increase the parameter by
        inc_param = (self.stop-self.start)/self.number_of
    
        cur_coef = self.start + inc_param
    
        get_make_dir(self)
        
        # make info file
        info_file = open(self.block_dir + '/info.txt','w')
        info_file.write('dir: '+str(self.block_dir)+'\nsurf: 1.0') 
        info_file.write('\ndt: '+str(self.dt))
        if self.var == 'beta':
            info_file.write('\nbeta (damping): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\nbeta (damping): '+str(self.beta))
        if self.var == 'cycles': 
            info_file.write('\ncycles (run time): sweep variable '+str(self.start)+'-'+str(self.stop))
        else: 
            info_file.write('\ncycles (run time): '+str(self.cycles))
        if self.var == 'num_cell':
            info_file.write('\nnum_cell (system length): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\nnum_cell (system length): '+str(self.num_cell))
        if self.var == 'A': 
            info_file.write('\nA (interaction amplitude): sweep variable '+str(self.start)+'-'+str(self.stop))
        else: 
            if type(A)==float: 
                info_file.write('\nA (interaction amplitude): '+str(self.A))
            else: 
                print('A NOT A FLOAT')
                #info_file.write('\nA (interaction amplitude): '+str(min(self.A))+'-'+str(max(self.A)))
        info_file.close()

        print('--dir is: '+str(self.block_dir))

        for l in range(self.number_of):
            #make the file
            cur_poin_file = open(self.block_dir+'/'+str(l)+'poindat.txt','w')
    
            #write the first few lines of the file
            cur_poin_file.write('coef --> ' +str(cur_coef)+'\n')
            cur_poin_file.write('ZONE I=something DATAPACKING=POINT\n')
            
            for kapa in range(self.vx_num):
                for alpha in range(self.x_num+1):
                    # write vx
                    cur_poin_file.write(str(self.vx_corner + kapa*increment_vx)+ "  ")
                    # vy
                    cur_poin_file.write(str(self.vy_corner)+"  ")
                    # x
                    cur_poin_file.write(str(self.x_corner + alpha*increment_x)+"  ")
                    # y
                    cur_poin_file.write(str(self.y_corner))
                    # write an enter
                    cur_poin_file.write('\n')
            
            cur_poin_file.close()
            cur_coef += inc_param

        # now that all the file are made make an "Orig" directory in the block directory to store
        # the origonal blocks of initial conditions. This way if something goes wrong we can pull
        # them out and run again.py again by hand
        os.system('mkdir ' + self.block_dir + '/Orig')
        os.system('cp ' + self.block_dir +'/*.txt ' + self.block_dir + '/Orig/')
        os.system('scp -r ./' + self.block_dir + ' omyers@bluemoon-user2.uvm.edu:/users/o/m/omyers/Data/EC/2DBlock/Old/')
        print self.block_dir

class One_Particle_Ensbl_Sin_2D(object):
    # number_of -> number of blocks (usualy 500). start -> parameter starting poin. stop ->
    # parameter stoping point (step size will be calculated). cycles-> runtime cycles. x_,vx_dim->
    # dimensions of block of initial conditions in phase space. x_,vx_num number of points in the x and
    # vx. The <>_corner variables define the lower left corner of the block in phase space
    def __init__(self):
        self.block_dir = ''
        self.var = 'A'
        self.script_dir= '/users/o/m/omyers/datasphere/ECproject/One_Particle_Sin_2D/'
        self.number_of = 30
        self.start     = 0.0
        self.stop      = 1.5
        self.dt        = 0.05
        self.cycles    = 10
        self.beta = .6
        self.x_num_cell = 1.0
        self.y_num_cell = 1.0
        self.dx = self.x_num_cell * 2.0 * pl.pi
        self.dy = self.y_num_cell * 2.0 * pl.pi
        self.A = .5
        self.x_dim     = 6.28
        self.y_dim    = 6.28
        self.vx_corner = 0.0
        self.x_corner  = 0.0
        self.vy_corner = 0.0
        self.y_corner  = 0.0

        self.x_num     = 40
        self.y_num    = 40
        # Full trajectory or just Poincare sections
        self.sliced = False


    def init_block(self):
        if (self.x_num == 0) or (self.x_dim==0.0):
            increment_x = 0.0
            if(self.x_num == self.x_dim):
                print('there is a descrepency between number and dimesions of block --> check to_run file')
        else:
            increment_x = self.x_dim/self.x_num 
            print('increment x: ' + str(increment_x))

        if (self.y_num == 0) or (self.y_dim==0.0):
            increment_y = 0.0
            if(self.y_num == self.y_dim):
                print('there is a descrepency between number and dimesions of block --> check to_run file')
        else:
            increment_y = self.y_dim/self.y_num 



        # find the amount we want to increase the parameter by
        inc_param = (self.stop-self.start)/self.number_of
    
        cur_coef = self.start + inc_param
    
        get_make_dir(self)
        
        # make info file
        info_file = open(self.block_dir + '/info.txt','w')
        info_file.write('dir: '+str(self.block_dir)+'\nsurf: 1.0') 
        info_file.write('\ndt: '+str(self.dt))
        if self.var == 'beta':
            info_file.write('\nbeta (damping): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\nbeta (damping): '+str(self.beta))
        if self.var == 'cycles': 
            info_file.write('\ncycles (run time): sweep variable '+str(self.start)+'-'+str(self.stop))
        else: 
            info_file.write('\ncycles (run time): '+str(self.cycles))
        if self.var == 'x_num_cell':
            info_file.write('\nx_num_cell (system length): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\nx_num_cell (system length): '+str(self.x_num_cell))
        if self.var == 'y_num_cell':
            info_file.write('\ny_num_cell (system length): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\ny_num_cell (system length): '+str(self.y_num_cell))


        if self.var == 'A': 
            info_file.write('\nA (interaction amplitude): sweep variable '+str(self.start)+'-'+str(self.stop))
        else: 
            if type(A)==float: 
                info_file.write('\nA (interaction amplitude): '+str(self.A))
            else: 
                print('A NOT A FLOAT')
                #info_file.write('\nA (interaction amplitude): '+str(min(self.A))+'-'+str(max(self.A)))
        info_file.close()

        print('--dir is: '+str(self.block_dir))

        for l in range(self.number_of):
            #make the file
            cur_poin_file = open(self.block_dir+'/'+str(l)+'poindat.txt','w')
    
            #write the first few lines of the file
            cur_poin_file.write('coef --> ' +str(cur_coef)+'\n')
            cur_poin_file.write('ZONE I=something DATAPACKING=POINT\n')
            
            for kapa in range(self.y_num+1):
                for alpha in range(self.x_num+1):
                    # write vx
                    cur_poin_file.write(str(self.y_corner)+ "  ")
                    # vy
                    cur_poin_file.write(str(self.vy_corner)+"  ")
                    # x
                    cur_poin_file.write(str(self.x_corner + alpha*increment_x)+"  ")
                    # y
                    cur_poin_file.write(str(self.y_corner + kapa*increment_y)+ "  ")
                    # write an enter
                    cur_poin_file.write('\n')
            
            cur_poin_file.close()
            cur_coef += inc_param

        # now that all the file are made make an "Orig" directory in the block directory to store
        # the origonal blocks of initial conditions. This way if something goes wrong we can pull
        # them out and run again.py again by hand
        os.system('mkdir ' + self.block_dir + '/Orig')
        os.system('cp ' + self.block_dir +'/*.txt ' + self.block_dir + '/Orig/')
        os.system('scp -r ./' + self.block_dir + ' omyers@bluemoon-user2.uvm.edu:/users/o/m/omyers/Data/EC/2DBlock/Old/')
        print self.block_dir


class Sin1D(object):
    # N -> number of particles
    # number_of -> number of blocks (usualy 500). i.e. how many steps
    # start -> parameter starting poin. 
    # stop -> parameter stoping point (step size will be calculated). 
    # cycles-> runtime cycles. 
    # num_cell -> number of unit cells to make simulation length over
    # qq -> intercharge force
    # d -> length of system as determined by the num_cell
    # A -> this is the interaction amplitude with the EC fild. It is an array so
    # that each particle may have its own interaction with the field
    def __init__(self):
        self.block_dir = ''
        self.var = 'A'
        self.script_dir= '/users/o/m/omyers/datasphere/ECproject/Sin1D/'
        self.number_of = 200
        self.start     = 1.0
        self.stop      = 1.4
        self.dt        = 0.001
        self.cycles    = 400
        self.N = 4
        self.qq = 1.0
        self.beta = .6
        self.num_cell = 1.0
        self.d = self.num_cell * 2.0 * pl.pi
        self.A = pl.zeros(2*self.N) +.5

        # Full trajectory or just Poincare sections
        self.sliced = False
        #self.sliced = True

    def init_block(self):

        # find the amount we want to increase the parameter by
        inc_param = (self.stop-self.start)/self.number_of
    
        cur_var = self.start
    
        get_make_dir(self)
        
        # make info file
        info_file = open(self.block_dir + '/info.txt','w')
        info_file.write('dir: '+str(self.block_dir)+'\nsurf: 1.0') 
        info_file.write('\ndt: '+str(self.dt))
        if self.var == 'beta':
            info_file.write('\nbeta (damping): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\nbeta (damping): '+str(self.beta))
        if self.var == 'cycles': 
            info_file.write('\ncycles (run time): sweep variable '+str(self.start)+'-'+str(self.stop))
        else: 
            info_file.write('\ncycles (run time): '+str(self.cycles))
        if self.var == 'N': 
            info_file.write('\nN (particle number): sweep variable '+str(self.start)+'-'+str(self.stop))
        else: 
            info_file.write('\nN (particle number): '+str(self.N))
        if self.var == 'qq': 
            info_file.write('\nqq (particle interaction strength): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\nqq (particle interaction strength): '+str(self.qq))
        if self.var == 'num_cell':
            info_file.write('\nnum_cell (system length): sweep variable '+str(self.start)+'-'+str(self.stop))
        else:
            info_file.write('\nnum_cell (system length): '+str(self.num_cell))
        if self.var == 'A': 
            info_file.write('\nA (interaction amplitude): sweep variable '+str(self.start)+'-'+str(self.stop))
        else: 
            if type(A)==float: 
                info_file.write('\nA (interaction amplitude): '+str(self.A))
            else: 
                info_file.write('\nA (interaction amplitude): '+str(min(self.A))+'-'+str(max(self.A)))
        info_file.close()

        for l in range(self.number_of):
            self.x0 = pl.zeros([2*self.N])
            
            for i,j in enumerate(self.x0):
                if i in range(self.N,2*self.N):
                    print(i)
                    self.x0[i] = random.random()*self.d
                    continue

            #make the file
            cur_poin_file = open(self.block_dir+'/'+str(l)+'poindat.txt','w')
    
            #write the first few lines of the file
            cur_poin_file.write(str(self.var)+' --> ' +str(cur_var)+'\n')
            
            # this just prints the numbers with one space imbetween. the .replace gets rid of the \n
            # but i'm still not really sure why they are there in the first place.
            cur_poin_file.write(str(self.x0)[1:-1].replace('\n',''))
            cur_poin_file.write('\n')
            
            cur_poin_file.close()
            cur_var += inc_param

        # now that all the file are made make an "Orig" directory in the block directory to store
        # the origonal blocks of initial conditions. This way if something goes wrong we can pull
        # them out and run again.py again by hand
        os.mkdir(self.block_dir + '/Orig')
        os.system('cp ' + self.block_dir +'/*.txt ' + self.block_dir + '/Orig/')
        os.system('scp -r ./' + self.block_dir + ' omyers@bluemoon-user2.uvm.edu:/users/o/m/omyers/Data/EC/2DBlock/Old/')
        print 'self.block_dir'
        print self.block_dir

class MB1DEC(object):
    # N -> number of particles
    # number_of -> number of blocks (usualy 500). i.e. how many steps
    # start -> parameter starting poin. 
    # stop -> parameter stoping point (step size will be calculated). 
    # cycles-> runtime cycles. 
    # num_cell -> number of unit cells to make simulation length over
    # qq -> intercharge force
    # d -> length of system as determined by the num_cell
    # A -> this is the interaction amplitude with the EC fild. It is an array so
    # that each particle may have its own interaction with the field
    def __init__(self):
        self.block_dir = ''
        self.var = 'A'
        self.script_dir= '/users/o/m/omyers/datasphere/ECproject/MB1DEC/'
        self.number_of = 300
        self.start     = 0.5
        self.stop      = 2.0
        self.dt        = 0.05
        self.cycles    = 120
        self.N = 20
        self.qq = .01
        self.beta = .6
        self.num_cell = 3.0
        self.d = self.num_cell * 2.0 * pl.pi
        self.A = pl.zeros(2*self.N) + .5

        # Full trajectory or just Poincare sections
        #self.sliced = False
        self.sliced = True

    
    def init_block(self):

        # find the amount we want to increase the parameter by
        inc_param = (self.stop-self.start)/self.number_of
    
        cur_var = self.start + inc_param
    
        get_make_dir(self)
        
        # make info file
        info_file = open(self.block_dir + '/info.txt','w')
        info_file.write('dir: '+str(self.block_dir)+'\nsurf: 1.0') 
        info_file.write('\ndt: '+str(self.dt))
        if self.var == 'beta':
            info_file.write('\nbeta (damping): sweep variable ')
        else:
            info_file.write('\nbeta (damping): '+str(self.beta))
        if self.var == 'cycles': 
            info_file.write('\ncycles (run time): sweep variable ')
        else: 
            info_file.write('\ncycles (run time): '+str(self.cycles))
        if self.var == 'N': 
            info_file.write('\nN (particle number): sweep variable ')
        else: 
            info_file.write('\nN (particle number): '+str(self.N))
        if self.var == 'qq': 
            info_file.write('\nqq (particle interaction strength): sweep variable ')
        else:
            info_file.write('\nqq (particle interaction strength): '+str(self.qq))
        if self.var == 'num_cell':
            info_file.write('\nnum_cell (system length): sweep variable ')
        else:
            info_file.write('\nnum_cell (system length): '+str(self.num_cell))
        if self.var == 'A': 
            info_file.write('\nA (interaction amplitude): sweep variable ')
        else: 
            if type(A)==float: 
                info_file.write('\nA (interaction amplitude): '+str(self.A))
            else: 
                info_file.write('\nA (interaction amplitude): '+str(min(self.A))+'-'+str(max(self.A)))
        info_file.close()

        for l in range(self.number_of):
            self.x0 = pl.zeros([2*self.N])
            
            for i,j in enumerate(self.x0):
                if i in range(self.N,2*self.N):
                    print(i)
                    self.x0[i] = random.random()*self.d
                    continue

            #make the file
            cur_poin_file = open(self.block_dir+'/'+str(l)+'poindat.txt','w')
    
            #write the first few lines of the file
            cur_poin_file.write(str(self.var)+' --> ' +str(cur_var)+'\n')
            
            # this just prints the numbers with one space imbetween. the .replace gets rid of the \n
            # but i'm still not really sure why they are there in the first place.
            cur_poin_file.write(str(self.x0)[1:-1].replace('\n',''))
            cur_poin_file.write('\n')
            
            cur_poin_file.close()
            cur_var += inc_param

        # now that all the file are made make an "Orig" directory in the block directory to store
        # the origonal blocks of initial conditions. This way if something goes wrong we can pull
        # them out and run again.py again by hand
        os.mkdir(self.block_dir + '/Orig')
        os.system('cp ' + self.block_dir +'/*.txt ' + self.block_dir + '/Orig/')
        os.system('scp -r ./' + self.block_dir + ' omyers@bluemoon-user2.uvm.edu:/users/o/m/omyers/Data/EC/2DBlock/Old/')
        print 'self.block_dir'
        print self.block_dir


   
#*****************************************************************************                         
#*****************************************************************************                         
class One_D_EC(object):
    # number_of -> number of blocks (usualy 500). start -> parameter starting poin. stop ->
    # parameter stoping point (step size will be calculated). cycles-> runtime cycles. x_,vx_dim->
    # dimensions of block of initial conditions in phase space. x_,vx_num number of points in the x and
    # vx. The <>_corner variables define the lower left corner of the block in phase space
    def __init__(self,the_file):
        lines = the_file.readlines()
        new_lines = []
        for i,j in enumerate(lines):
            if j=='\n':
                continue
            if '#' in j.split()[0][0]:
                continue
            new_lines.append(j)

        self.block_dir = ''
        self.script_dir= '/users/o/m/omyers/datasphere/ECproject/1DEC/'
        self.number_of = int(new_lines[0].split()[0] )
        self.start     = float(new_lines[1].split()[0] )
        self.stop      = float(new_lines[2].split()[0] )
        self.k         = float(new_lines[3].split()[0] )
        self.w         = float(new_lines[4].split()[0] )
        self.beta      = float(new_lines[5].split()[0] )
        self.dt        = float(new_lines[6].split()[0] )
        self.cycles    = float(new_lines[7].split()[0] )
        self.x_dim     = float(new_lines[8].split()[0] )
        self.vx_dim    = float(new_lines[9].split()[0] )
        self.vx_corner = float(new_lines[10].split()[0])
        self.vy_corner = float(new_lines[11].split()[0])
        self.x_corner  = float(new_lines[12].split()[0])
        self.y_corner  = float(new_lines[13].split()[0])
        self.x_num     = int(new_lines[14].split()[0]) 
        self.vx_num    = int(new_lines[15].split()[0]) 
        if 'Full' in new_lines[16].split()[0]:
            self.sliced = False
        if 'Poin' in new_lines[16].split()[0]:
            self.sliced = True

    def init_block(self):
        if (self.x_num == 0) or (self.x_dim==0.0):
            increment_x = 0.0
            if(self.x_num == self.x_dim):
                print('there is a descrepency between number and dimesions of block --> check to_run file')
        else:
            increment_x = self.x_dim/self.x_num 
            print('increment x: ' + str(increment_x))

        if (self.vx_num == 0) or (self.vx_dim==0.0):
            increment_vx = 0.0
            if(self.vx_num == self.vx_dim):
                print('there is a descrepency between number and dimesions of block --> check to_run file')
        else:
            increment_vx = self.vx_dim/self.vx_num 



        # find the amount we want to increase the parameter by
        inc_param = (self.stop-self.start)/self.number_of
    
        cur_coef = self.start + inc_param
    
        get_make_dir(self)
        
        # make info file
        info_file = open(self.block_dir + '/info.txt','w')
        info_file.write('--dir\n'+str(self.block_dir)+'\n--surf \n1.0 \n--k \n'+str(self.k)\
                +'\n--w\n'+str(self.w)+'\n--beta\n'+str(self.beta)+'\n--g\n0.1\n--dt\n'\
                +str(self.dt)+'\n\n')
        info_file.close()

        print('--dir is ( I think there is a problem with this): '+str(self.block_dir))

        for l in range(self.number_of):
            #make the file
            cur_poin_file = open(self.block_dir+'/'+str(l)+'poindat.txt','w')
    
            #write the first few lines of the file
            cur_poin_file.write('coef --> ' +str(cur_coef)+'\n')
            cur_poin_file.write('ZONE I=something DATAPACKING=POINT\n')
            
            for kapa in range(self.vx_num):
                for alpha in range(self.x_num+1):
                    # write vx
                    cur_poin_file.write(str(self.vx_corner + kapa*increment_vx)+ "  ")
                    # vy
                    cur_poin_file.write(str(self.vy_corner)+"  ")
                    # x
                    cur_poin_file.write(str(self.x_corner + alpha*increment_x)+"  ")
                    # y
                    cur_poin_file.write(str(self.y_corner))
                    # write an enter
                    cur_poin_file.write('\n')
            
            cur_poin_file.close()
            cur_coef += inc_param

        # now that all the file are made make an "Orig" directory in the block directory to store
        # the origonal blocks of initial conditions. This way if something goes wrong we can pull
        # them out and run again.py again by hand
        os.system('mkdir ' + self.block_dir + '/Orig')
        os.system('cp ' + self.block_dir +'/*.txt ' + self.block_dir + '/Orig/')
        os.system('scp -r ./' + self.block_dir + ' omyers@bluemoon-user2.uvm.edu:/users/o/m/omyers/Data/EC/2DBlock/Old/')
        print self.block_dir

   
#*****************************************************************************                         
#*****************************************************************************                         

class Twin_EC(object):
    # number_of -> number of blocks (usualy 500). start -> parameter starting poin. stop ->
    # parameter stoping point (step size will be calculated). cycles-> runtime cycles. x_,vx_dim->
    # dimensions of block of initial conditions in phase space. x_,vx_num number of points in the x and
    # vx. The <>_corner variables define the lower left corner of the block in phase space
    def __init__(self,the_file):
        lines = the_file.readlines()
        new_lines = []
        for i,j in enumerate(lines):
            if j=='\n':
                continue
            if '#' in j.split()[0][0]:
                continue
            new_lines.append(j)

        self.block_dir = ''
        self.script_dir= '/users/o/m/omyers/datasphere/ECproject/TwinEC/'
        self.number_of = int(new_lines[0].split()[0] )
        self.start     = float(new_lines[1].split()[0] )
        self.stop      = float(new_lines[2].split()[0] )
        # dimesionless parameters (d,x,y,,beta)
        self.beta      = float(new_lines[3].split()[0] )
        self.d         = float(new_lines[4].split()[0] )
        self.dt        = float(new_lines[5].split()[0] )
        self.cycles    = float(new_lines[6].split()[0] )
        self.x_dim     = float(new_lines[7].split()[0] )
        self.y_dim     = float(new_lines[8].split()[0] )
        self.vx_dim    = float(new_lines[9].split()[0] )
        self.vy_dim    = float(new_lines[10].split()[0] )
        self.vx_corner = float(new_lines[11].split()[0])
        self.vy_corner = float(new_lines[12].split()[0])
        self.x_corner  = float(new_lines[13].split()[0])
        self.y_corner  = float(new_lines[14].split()[0])
        self.x_num     = int(new_lines[15].split()[0]) 
        self.y_num     = int(new_lines[16].split()[0]) 
        self.vx_num    = int(new_lines[17].split()[0]) 
        self.vy_num    = int(new_lines[18].split()[0]) 
        # this is an important variable that we need to be carefull with. THis is how we tell the
        # program that we want full trajectoies (LOOOTS of memory) or Poincare sections. It is a
        # buliean value 'sliced' as in 'time sliced' i.e. poincare sections. 
        if 'Full' in new_lines[19].split()[0]:
            self.sliced = False
        if 'Poin' in new_lines[19].split()[0]:
            self.sliced = True

   
    def init_block(self):

        # find the amount we want to increase the parameter by
        if (self.x_num == 0) or (self.x_dim==0.0):
            increment_x = 0.0
            if(self.x_num != self.x_dim):
                print('there is a descrepency between number and dimesions of block --> check to_run file')
        else:
            increment_x = self.x_dim/self.x_num 

        if (self.y_num == 0) or (self.y_dim==0.0):
            increment_y = 0.0
            if(self.y_num != self.y_dim):
                print('there is a descrepency between number and dimesions of block --> check to_run file')
        else:
            increment_y = self.x_dim/self.y_num 

        if (self.vx_num == 0) or (self.vx_dim==0.0):
            increment_vx = 0.0
            if(self.vx_num != self.vx_dim):
                print('there is a descrepency between number and dimesions of block --> check to_run file')
        else:
            increment_vx = self.x_dim/self.vx_num 

        if (self.vy_num == 0) or (self.vy_dim==0.0):
            increment_vy = 0.0
            if(self.vy_num != self.vy_dim):
                print('there is a descrepency between number and dimesions of block --> check to_run file')
        else:
            increment_vy = self.x_dim/self.vy_num 


        inc_param = (self.stop-self.start)/self.number_of
    
        cur_coef = self.start + inc_param
    
        get_make_dir(self)
        
        # make info file
        info_file = open(self.block_dir + '/info.txt','w')
        info_file.write('--dir\n'+str(self.block_dir)+'\n--d \n'+str(self.d)\
                +'\n--beta\n'+str(self.beta)+'\n--dt\n'\
                +str(self.dt)+'\n\n')
        info_file.close()
    
        for l in range(self.number_of):
            #make the file
            cur_poin_file = open(self.block_dir+'/'+str(l)+'poindat.txt','w')
    
            #write the first few lines of the file
            cur_poin_file.write('coef --> ' +str(cur_coef)+'\n')
            cur_poin_file.write('ZONE I=something DATAPACKING=POINT\n')
            
            for kapa in range(self.y_num+1):
                for alpha in range(self.x_num+1):
                    for gama in range(self.vx_num+1):
                        for nu in range(self.vy_num+1):
                            # write vx
                            cur_poin_file.write(str(self.vx_corner+ gama*increment_vx)+"  ")
                            # vy
                            cur_poin_file.write(str(self.vy_corner + nu*increment_vy)+"  ")
                            # x
                            cur_poin_file.write(str(self.x_corner + alpha*increment_x)+"  ")
                            # y
                            cur_poin_file.write(str(self.y_corner + kapa*increment_y)+"  ")
                            # write an enter
                            cur_poin_file.write('\n')
            
            cur_poin_file.close()
            cur_coef += inc_param

        # now that all the file are made make an "Orig" directory in the block directory to store
        # the origonal blocks of initial conditions. This way if something goes wrong we can pull
        # them out and run again.py again by hand
        os.system('mkdir ' + self.block_dir + '/Orig')
        os.system('cp ' + self.block_dir +'/*.txt ' + self.block_dir + '/Orig/')
        os.system('scp -r ./' + self.block_dir + ' omyers@bluemoon-user2.uvm.edu:/users/o/m/omyers/Data/EC/4DBlock/Old/')
        print self.block_dir

#**********Functions**********************************************************                         
def we_done_yet(ssh,job_id,to_run_object):
    stdin,stdout,stderr = ssh.exec_command('ls | grep '+str(job_id))
    out_lines = stdout.readlines()
    err_lines = stderr.readlines()

    number_compleated = len(out_lines)
    print('This is for checkin stuff in the we_done_yet() function:')
    print('number_compleated is: ' +str(number_compleated))
    print('to_run_object.number_of is: ' +str(to_run_object.number_of))

    if number_compleated == to_run_object.number_of:
        return True
    else:
        return False

def get_make_dir(to__object):
    # use time to name directories
    timestr = thetime.asctime()
    newtimestr = "Block_"
    for a,b in enumerate(timestr.split()):
        newtimestr+=b+"_"
    newtimestr = newtimestr.replace(":","")
    print(newtimestr)
    os.system("mkdir "+newtimestr)
    to__object.block_dir = newtimestr


#*****************************************************************************                         

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--keyword',action='store',dest = 'keyword',type = str, required = True)
    inargs = parser.parse_args() 
    key_word = inargs.keyword

    to_run_file = open('to_run_'+key_word+'.txt','r')
     
    if key_word == '1DEC':
        # pass the file so it can get all of the variables
        to_run_object = One_D_EC(to_run_file) 

    if key_word == 'TwinEC':
        to_run_object = Twin_EC(to_run_file)

    if key_word == 'MB1DEC':
        to_run_object = MB1DEC()

    if key_word == 'Sin1D':
        to_run_object = Sin1D()
    if key_word == 'Sin2D':
        to_run_object = Sin2D()
    if key_word == 'Ensbl_Sin1D':
        to_run_object = One_Particle_Ensbl_Sin_1D()
    if key_word == 'Ensbl_Sin2D':
        to_run_object = One_Particle_Ensbl_Sin_2D()

    
    # this initializes the block on this computer and then scp it to cluster
    to_run_object.init_block()

    # now that block is on the cluster we need to run again.py on it
    # first get totIter from the number of cycles to give to again.script
    totIter = int(to_run_object.cycles*2.0*pl.pi/to_run_object.dt)
    # We also need to give again.script the block name from to_run_object.block_dir.
    # sleep for a couple of min so transfer can compleet
    thetime.sleep(200)
    print('totIter: '+str(totIter))

    # The way to do this is definetely with paramiko
    # make paramiko object. call it ssh (why not)
    ssh = paramiko.SSHClient()
    # now we are going to make it so we can log into zoo. This is only ok to do becuase we are
    # behind their firewall while we are on campus anyway.
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    # now connect!
    ssh.connect('bluemoon-user2.uvm.edu',username='omyers',password='376.re.1873.oven')
    a,b,c = ssh.exec_command('pwd')
    print('working directory after login: ' + str(b.readlines()))
    print('to_run_object.block_dir is: '+str(to_run_object.block_dir))

    print('cp command: '+ 'cp '+to_run_object.script_dir+'again.script again.script')
    stdin,stdout,stderr = ssh.exec_command('cp '+to_run_object.script_dir+'again.script again.script')
    out_lines = stdout.readlines()
    err_lines = stderr.readlines()
    print('stdout lines after copy: ' + str(out_lines))
    print('stderr lines after copy: ' + str(err_lines))

    #sleep_secs = 100
    #print('sleeping for ... '+str(sleep_secs)+' seconds')
    #thetime.sleep(sleep_secs)

    # variable into qsub works through the following format:
    # qsub -v param1=val,param2=val,... script.sh
    # FIN is the totatl number of blocks we are running. if we are running 500 blocks FIN = 499
    # -q shortq  
    print('qsub command: '+'qsub -v DIR='+to_run_object.block_dir+ ',FIN='+str(to_run_object.number_of-1)+ ',TOTITER='+str(totIter)+ ',SLICED='+str(to_run_object.sliced)+ ' -t 0-' +str(to_run_object.number_of-1)+ ' again.script')

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('bluemoon-user2.uvm.edu',username='omyers',password='376.re.1873.oven')

    stdin,stdout,stderr = ssh.exec_command('qsub -v DIR='+to_run_object.block_dir+
            ',FIN='+str(to_run_object.number_of-1)+
            ',TOTITER='+str(totIter)+
            ',SLICED='+str(to_run_object.sliced)+
            ' -t 0-' +str(to_run_object.number_of-1)+ 
            ' again.script')
    out_lines = stdout.readlines()
    err_lines = stderr.readlines()

    print('length of job id out_lines: ' + str(len(out_lines)))
    if len(out_lines) == 0:
        shutil.rmtree(to_run_object.block_dir)

    job_id = out_lines[0][:7]
    print('job id: ' + str(job_id))
    if len(job_id)==0:
        os.system('say job did not not not run on cluster')
    else: os.system('say job is successfully running on cluster')

    print('stdout lines: ' + str(out_lines))
    print('stderr lines: ' + str(err_lines))
    #stdin,stdout,stderr = ssh.exec_command('rm tempagain.script')
    
    ## listen and wait for all out put to be done.
    ## for testing lets just wait 5min
    #while not we_done_yet(ssh,job_id,to_run_object):
    #    #thetime.sleep(300)
    #    thetime.sleep(10)
    #
    #print('WE MADE IT YAY')

    # a couple of audible bells to let me know its done
    print('\a')

if __name__ == '__main__':
    main()
