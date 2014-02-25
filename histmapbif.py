#!/usr/bin/python
import pylab as pl
import os
import argparse
import scipy as np
import o_funcs as of
import subprocess

# The purpos of this function is to grab all the path names of all the files in different block
# directories that all have the same N. We need to be slightly carful using this because several
# things are based off of this assuption that N is the sam and the number of point of the sweep
# variable are all the same.
# Also this assums that all the directories in the given directory (d) is stuff you want in the
# plot.
# Function returns a list of all the paths.
def go_deep():
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
    
    return list_dir

def main():

    min_count = 1
    parser = argparse.ArgumentParser()
    # d is for directory
    parser.add_argument('-d',action='store',dest = 'd',type = str, required = False, default = './')
    # number of PC sections we want to use (going back from the last)
    parser.add_argument('-n',action='store',dest = 'n',type = int, required = False)
    # project onto which axis?
    parser.add_argument('-p',action='store',dest = 'p',type = str, required = False, default = 'x')
    # grid size
    parser.add_argument('-g',action='store',dest = 'g',type = int, required = False, default = 300)
    # This is going to be a unique option to collect AAAALLLLLL the data of several individual
    # bifurcation runs (everything the same except random initial conditions) into one so that the
    # diagram does not look so messy. The argument being passed is a directory --> no default
    # option.
    parser.add_argument('--together',action='store',dest='together',type = bool, required = False, default=False)
    inargs = parser.parse_args()
    d = inargs.d
    projection = inargs.p
    grid_size = inargs.g 
    together = inargs.together


    num_pc = inargs.n

    print('type num_pc: ' + str(type(num_pc)))

    os.chdir(d)

    # need to go into one of the directories to get the system info if we are doing 'together'.
    if together: os.chdir(os.listdir('.')[0])

    # get the  system info
    # note: if Dim = 1 y_num_cell -> ' No y ' and order -> 'polygamma'
    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim = of.get_system_info()
    print('Dim is: ' +str(Dim))
    print('x_num_cell: ' + str(x_num_cell))

    # need to get back out for later on
    if together: os.chdir('..')

    # depending on the projection type we will determine "ax_num" variable. This number is an easy
    # way to controle whicn axis we are projecting onto by the way we slce the data. form is this:
    # data[time,ax_num*N:(ax_num + 1)*N] 
    # for ax_num = 0 ---> vx data
    # for ax_num = 1 ---> vy data
    # for ax_num = 2 ---> x data
    # for ax_num = 3 ---> y data
    # lets also define the string for the axis lable depending on the projection
    if projection == 'x':
        ax_str = r'$x$'
        slice_num = Dim
    if projection == 'y':
        if Dim == 1: 
            print('No y projectino in 1D')
            quit()
        slice_num = Dim+1
        ax_str = r'$y$'
    if projection == 'vx':
        ax_str = r'$\dot{x}$'
        slice_num = 0
    if projection == 'vy':
        if Dim == 1: 
            print('No vy projectino in 1D')
            quit()
        slice_num = Dim-1
        ax_str = r'$\dot{y}$'
    
    if together:
        list_dir = go_deep()
    else:
        all_dir = os.listdir(".")
        list_dir = [] 

        num_runs = 0
        for i,j in enumerate(all_dir):
            if "poindat" in j:
                list_dir.append(j)
                num_runs+=1
        print('num_runs is: '+str(num_runs))
            
   
    build = pl.zeros(N)+1.0
    # this array keeps track of the variable values for each particle
    var_arr = pl.array([])
    all_data = pl.array([])
    # the data to be ploted. projection type will determine which data gets stored in this array
    to_plot = pl.array([])

    for i,j in enumerate(list_dir):
        print('cur_file should be: ' + str(j))
        cur_file = open(j,"r")
        # get the varible for the plot
        var = float(cur_file.readline().split()[-1])
        data = np.genfromtxt(cur_file)
        print('len(data): '+str(len(data)))
        print('shape(data): ' + str(pl.shape(data)))
        # modulus data by system length
        data[:,Dim*N:(Dim+1)*N] = data[:,Dim*N:(Dim+1)*N]%(x_num_cell*2.0*pl.pi)
        if Dim ==2:
            print('in Dim==2 section')
            data[:,(Dim+1)*N:(Dim+2)*N] = data[:,(Dim+1)*N:(Dim+2)*N]%(y_num_cell*2.0*pl.pi)
        #print('shape(data) after % x/y_num_cell: ' + str(pl.shape(data)))
        # if we want more poincare section then count back from the last checking to see if we are
        # at the t%2pl.pi=0 point and add data accordingly
        if num_pc != None:
            # indexing variable
            counting = 0
            # keep track of how many poincare sections we have saved
            passed_pc = 0
            while passed_pc <= num_pc:
                if (counting*dt)%(2.0*pl.pi)<=dt:
                    print('passed_pc: '+str(passed_pc))
                    to_plot = pl.append(to_plot,data[-counting,slice_num*N:(slice_num+1)*N])
                    print('making to_plot')
                    #print(to_plot)
                    var_arr = pl.append(var_arr,build*var)
                    passed_pc += 1
                counting +=1
        else:
            print('in else')
            #print('to be added to to_plot: '+ str(data[-1,slice_num*N:(slice_num+1)*N]))
            to_plot = pl.append(to_plot,data[-1,slice_num*N:(slice_num+1)*N])
            print('making to_plot')
            #print(to_plot)
            #all_data = pl.append(all_data,data)
            var_arr = pl.append(var_arr,build*var)
    
        cur_file.close()

    print('len(to_plot): '+str(len(to_plot)))
    print(to_plot)
    print('len(var_arr): ' + str(len(var_arr)))
    print(var_arr)

    pl.hexbin(var_arr,to_plot,gridsize=grid_size,bins='log',mincnt=min_count,edgecolor='none') 
    #cmap=pl.cm.YlOrRd_r
    #pl.hexbin(A_arr,x,gridsize=400,cmap=pl.cm.Greys,bins='log',mincnt=2,edgecolor='none') 
    #pl.hexbin(A_arr,x,bins='log') 

    pl.colorbar()
    pl.xlabel(sweep_str,fontsize=30)
    pl.ylabel(ax_str,fontsize=30)
    #pl.tight_layout()
    # No more over writing images
    number = 0
    save_str = 'paper_'+projection+'_bif_hist'+str(number)+'.png'
    while save_str in os.listdir('.'):
        number +=1
        save_str = 'paper_'+projection+'_bif_hist'+str(number)+'.png'

    pl.savefig(save_str,dpi=200)
    os.system('open paper_'+projection+'_bif_hist'+str(number)+'.png')

    subprocess.call('say Finished the bifurcation histogram',shell = True)

    #ax.axis([0.0,2*pl.pi,-1.8,1.8])


if __name__ == "__main__":
    main()