import os
import pylab as pl
import argparse
from datetime import datetime
import time as thetime
import paramiko

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

    def get_make_dir(self):
        # use time to name directories
        timestr = thetime.asctime()
        newtimestr = "Block_"
        for a,b in enumerate(timestr.split()):
            newtimestr+=b+"_"
        newtimestr = newtimestr.replace(":","")
        print(newtimestr)
        os.system("mkdir "+newtimestr)
        self.block_dir = newtimestr
    
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
    
        self.get_make_dir()
        
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

    def we_done_yet(self,ssh,job_id):
        stdin,stdout,stderr = ssh.exec_command('ls | grep '+str(job_id))
        out_lines = stdout.readlines()
        err_lines = stderr.readlines()

        number_compleated = len(out_lines)
        print('This is for checkin stuff in the we_done_yet() function:')
        print('number_compleated is: ' +str(number_compleated))
        print('self.number_of is: ' +str(self.number_of))

        if number_compleated == self.number_of:
            return True
        else:
            return False

   
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


    def get_make_dir(self):
        # use time to name directories
        timestr = thetime.asctime()
        newtimestr = "Block_"
        for a,b in enumerate(timestr.split()):
            newtimestr+=b+"_"
        newtimestr = newtimestr.replace(":","")
        print(newtimestr)
        os.system("mkdir "+newtimestr)
        self.block_dir = newtimestr
    
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
    
        self.get_make_dir()
        
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

    def we_done_yet(self,ssh,job_id):
        stdin,stdout,stderr = ssh.exec_command('ls | grep '+str(job_id))
        out_lines = stdout.readlines()
        err_lines = stderr.readlines()

        number_compleated = len(out_lines)
        print('This is for checkin stuff in the we_done_yet() function:')
        print('number_compleated is: ' +str(number_compleated))
        print('self.number_of is: ' +str(self.number_of))

        if number_compleated == self.number_of:
            return True
        else:
            return False

#*****************************************************************************                         
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
    # variable into qsub works through the following format:
    # qsub -v param1=val,param2=val,... script.sh
    # FIN is the totatl number of blocks we are running. if we are running 500 blocks FIN = 499
    print('qsub command: '+'qsub -v DIR='+to_run_object.block_dir+ ',FIN='+str(to_run_object.number_of-1)+ ',TOTITER='+str(totIter)+ ',SLICED='+str(to_run_object.sliced)+ ' -t 0-' +str(to_run_object.number_of-1)+ 'again.script')
    stdin,stdout,stderr = ssh.exec_command('qsub -v DIR='+to_run_object.block_dir+
            ',FIN='+str(to_run_object.number_of-1)+
            ',TOTITER='+str(totIter)+
            ',SLICED='+str(to_run_object.sliced)+
            ' -t 0-' +str(to_run_object.number_of-1)+ 
            ' again.script')
    out_lines = stdout.readlines()
    err_lines = stderr.readlines()

    job_id = out_lines[0][:7]

    print('stdout lines: ' + str(out_lines))
    print('stderr lines: ' + str(err_lines))
    #stdin,stdout,stderr = ssh.exec_command('rm tempagain.script')
    
    # listen and wait for all out put to be done.
    # for testing lets just wait 5min
    while not to_run_object.we_done_yet(ssh,job_id):
        thetime.sleep(300)
    
    # remove the origonal files
    stdin,stdout,stderr = ssh.exec_command('rm \
            /users/o/m/omyers/Data/EC/2DBlock/Old/'+to_run_object.block_dir+'/*poindat.txt')
    # unzip the new ones. 
    # stdin,stdout,stderr = ssh.exec_command('gunzip /users/o/m/omyers/Data/EC/2DBlock/'+to_run_object.block_dir/+'*poindat.txt.gz')
    
    print('WE MADE IT YAY')

if __name__ == '__main__':
    main()
