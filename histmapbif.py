#!/usr/bin/python
import pylab as pl
import os
import argparse
import scipy as np
import o_funcs as of
import subprocess


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

    # ancl --> anal class
    ancl = of.anal_run()
    ancl.together = together
    ancl.get_info()
    ancl.set_list_dir()

    # get the  system info
    # note: if Dim = 1 y_num_cell -> ' No y ' and order -> 'polygamma'
    print('Dim is: ' +str(ancl.Dim))
    print('x_num_cell: ' + str(ancl.x_num_cell))

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
        slice_num = ancl.Dim
    if projection == 'y':
        if ancl.Dim == 1: 
            print('No y projectino in 1D')
            quit()
        slice_num = ancl.Dim+1
        ax_str = r'$y$'
    if projection == 'vx':
        ax_str = r'$\dot{x}$'
        slice_num = 0
    if projection == 'vy':
        if Dim == 1: 
            print('No vy projectino in 1D')
            quit()
        slice_num = ancl.Dim-1
        ax_str = r'$\dot{y}$'
    
   
    build = pl.zeros(ancl.N)+1.0
    # this array keeps track of the variable values for each particle
    var_arr = pl.array([])
    all_data = pl.array([])
    # the data to be ploted. projection type will determine which data gets stored in this array
    to_plot = pl.array([])

    for i,j in enumerate(ancl.list_dir):
        print('cur_file should be: ' + str(j))
        cur_file = open(j,"r")
        # get the varible for the plot
        var = float(cur_file.readline().split()[-1])
        data = np.genfromtxt(cur_file)
        print('len(data): '+str(len(data)))
        print('shape(data): ' + str(pl.shape(data)))
        # modulus data by system length
        data[:,ancl.Dim*ancl.N:(ancl.Dim+1)*ancl.N] = data[:,ancl.Dim*ancl.N:(ancl.Dim+1)*ancl.N]%(ancl.x_num_cell*2.0*pl.pi)
        if ancl.Dim ==2:
            print('in Dim==2 section')
            data[:,(ancl.Dim+1)*ancl.N:(ancl.Dim+2)*ancl.N] = data[:,(ancl.Dim+1)*ancl.N:(ancl.Dim+2)*ancl.N]%(ancl.y_num_cell*2.0*pl.pi)
        #print('shape(data) after % x/y_num_cell: ' + str(pl.shape(data)))
        # if we want more poincare section then count back from the last checking to see if we are
        # at the t%2pl.pi=0 point and add data accordingly
        if num_pc != None:
            # indexing variable
            counting = 0
            # keep track of how many poincare sections we have saved
            passed_pc = 0
            while passed_pc <= num_pc:
                if (counting*ancl.dt)%(2.0*pl.pi)<=ancl.dt:
                    print('passed_pc: '+str(passed_pc))
                    to_plot = pl.append(to_plot,data[-counting,slice_num*ancl.N:(slice_num+1)*ancl.N])
                    print('making to_plot')
                    #print(to_plot)
                    var_arr = pl.append(var_arr,build*var)
                    passed_pc += 1
                counting +=1
        else:
            print('in else')
            #print('to be added to to_plot: '+ str(data[-1,slice_num*N:(slice_num+1)*N]))
            to_plot = pl.append(to_plot,data[-1,slice_num*ancl.N:(slice_num+1)*ancl.N])
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
    pl.xlabel(ancl.sweep_str,fontsize=30)
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
