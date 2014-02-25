import pylab as pl
import os
import argparse
import scipy as np


def main():
    grid_size = 400
    min_count = 1
    parser = argparse.ArgumentParser()
    # d is for directory
    parser.add_argument('-d',action='store',dest = 'd',type = str, required = True)
    # number of PC sections we want to use (going back from the last)
    parser.add_argument('-n',action='store',dest = 'n',type = int, required = False)
    inargs = parser.parse_args()
    d = inargs.d
    num_pc = inargs.n

    print('type num_pc: ' + str(type(num_pc)))


    os.chdir(d)

    all_dir = os.listdir(".")
    list_dir = [] 

    num_runs = 0
    for i,j in enumerate(all_dir):
        if "poindat" in j:
            list_dir.append(j)
            num_runs+=1
    print('num_runs is: '+str(num_runs))
            
    # get the  system info
    info_f = open('info.txt','r')
    l = info_f.readlines()

    for i,j in enumerate(l):
        if 'cell' in j:
            num_cell = float(j.split()[-1])
            print('num_cell: ' + str(num_cell))
        if 'qq' in j:
            qq = float(j.split()[-1])
            print('qq: ' + str(qq))
        if 'dt' in j:
            dt = float(j.split()[-1])
            print('dt: ' + str(dt))
        if 'beta' in j:
            beta = float(j.split()[-1])
            print('beta: ' + str(beta))
        if 'A:' in j:
            A = float(j.split()[-1])
            print('A: ' + str(A))
        if 'cycles' in j:
            cycles = int(j.split()[-1])
            print('cycles: ' + str(cycles))
        if 'particle number' in j:
            N = int(j.split()[-1])
            print('N: ' + str(N))

   
    build = pl.zeros(N)+1.0
    # this array keeps track of the variable values for each particle
    var_arr = pl.array([])
    all_data = pl.array([])
    x = pl.array([])

    for i,j in enumerate(list_dir):
        print('cur_file should be: ' + str(j))
        cur_file = open(j,"r")
        # get the varible for the plot
        var = float(cur_file.readline().split()[-1])
        data = np.genfromtxt(cur_file)
        #print(pl.shape(data))
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
                    x = pl.append(x,data[-counting,N:(2*N)]%(num_cell*2.0*pl.pi))
                    #all_data = pl.append(all_data,data)
                    var_arr = pl.append(var_arr,build*var)
                    passed_pc += 1
                counting +=1
        else:
            x = pl.append(x,data[-1,N:(2*N)]%(num_cell*2.0*pl.pi))
            #all_data = pl.append(all_data,data)
            var_arr = pl.append(var_arr,build*var)
    
        cur_file.close()

    # forms all_data into the following: [var run, time, degrees of freedom]
    #all_data = all_data.reshape(num_runs,-1,2*N)
    # modulus by the length of the system
    #all_data[:,:,N:2*N] = all_data[:,:,N:2*N]%(num_cell*2.0*pl.pi)
        
    pl.hexbin(var_arr,x,gridsize=grid_size,bins="log",mincnt=min_count,edgecolor="none") 
    #pl.hexbin(var_arr,x,gridsize=grid_size,bins="log",mincnt=2,edgecolor="none") 

    #cmap=pl.cm.YlOrRd_r
    #pl.hexbin(A_arr,x,gridsize=400,cmap=pl.cm.Greys,bins="log",mincnt=2,edgecolor="none") 
    #pl.hexbin(A_arr,x,bins="log") 

    pl.colorbar()
    pl.xlabel("$A$",fontsize=30)
    pl.ylabel("$x$",fontsize=30)
    pl.tight_layout()
    pl.savefig("paper_bif_hist.png",dpi=300)
    os.system("open paper_bif_hist.png")

    #ax.set_xlabel("$A$",fontsize =  25 )
    #ax.set_ylabel("$x$", fontsize = 25 )
    #ax.axis([0.0,2*pl.pi,-1.8,1.8])
    # j[:-11] is the number infrot of poindat.txt
    #fig.savefig("bif_hist.png")


if __name__ == "__main__":
    main()
