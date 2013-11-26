import pylab as pl
import os
import argparse
import plot_one as po

# autocorrelation fuction. tau is the comonly used variable for the length of time one is looking
# for correlations over. tau is going to be just an integer. If a specific value of time is desired
# it needs to be turned into an integer via deviding by dt. here is only an integer
def acf(tau,data):
    # data has already been passed such that is is a 1d array of time serries form. Be carefule with
    # the input then!
    return sum(data[:-tau]*data[tau:])/(len(data)-tau)

def make_acf_tau_plot(f,v):

    cur_poin_num = int(f[:f.find('p')])
    if (str(cur_poin_num)+'RunImages') not in os.listdir('.'):
        os.mkdir(str(cur_poin_num)+'RunImages')

    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str  = po.get_system_info()
    to_save_dir = str(cur_poin_num)+'RunImages/TauACF_Var_'+v
    os.mkdir(to_save_dir)
    x_lbl = r'tau'
    y_lbl = 'Autocorrelation'

    work_file = open(f,'r')
    sweep_var = float(work_file.readline().split()[-1])
    data = pl.genfromtxt(work_file)
    work_file.close()


#    # depending on the variable we want we need a number that controles how we slice the data
#    # lines. With the dimension of the system we can slice everything right. Can get the dimension with the
#    # following line.
#    Dim = (pl.shape(data)[1])/(2*N)
#    # can use the dimension to slice everything right
#    if v == 'x':
#        input_arr = data[:,Dim*N:(Dim+1)*N]
#    if v == 'vx':
#        input_arr = data[:,0:Dim*N]
#    # Will only get into these if asking for 2D stuff anyway
#    if v == 'y':
#        input_arr = data[:,(Dim+1)*N:(Dim+2)*N]
#    if v == 'vy':
#        input_arr = data[:,(Dim-1)*N:(Dim)*N]

   

    # make em for every particle in the simulation
    for a in range(N):

        # depending on the variable we want we need a number that controles how we slice the data
        # lines. With the dimension of the system we can slice everything right. Can get the dimension with the
        # following line.
        Dim = (pl.shape(data)[1])/(2*N)
        # can use the dimension to slice everything right
        if v == 'x':
            input_arr = data[:,Dim*N+a]
        if v == 'vx':
            input_arr = data[:,a]
        # Will only get into these if asking for 2D stuff anyway
        if v == 'y':
            input_arr = data[:,(Dim+1)*N+a]
        if v == 'vy':
            input_arr = data[:,(Dim-1)*N+a]

        # lets try this for tau from dt to 5 cycles
        tau_arr = range(int(5*2.0*pl.pi/dt))
        output_arr = pl.array([])
        for i,j in enumerate(tau_arr):
            output_arr = pl.append(output_arr,acf(j,input_arr))
        
        
        fig = pl.figure()
        ax = fig.add_subplot(111)
        #ax.set_xlim(x_rng)
        #ax.set_ylim(y_rng)
        ax.set_xlabel(x_lbl,fontsize=30)
        ax.set_ylabel(y_lbl,fontsize=30)
        ax.scatter(tau_arr,output_arr,c='k')
        fig.tight_layout()
        fig.savefig(to_save_dir+'/%(number)04d.png'%{'number':a})
        pl.close(fig)

def make_acf_sweep_plot(v):

    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str  = get_system_info()
    to_save_dir = 'SweepACF_Var_'+v
    os.mkdir(to_save_dir)
    x_lbl = sweep_str
    y_lbl = 'Autocorrelation'

def main():

    # autocorrelation function as function of sweep parameter
    # we need two differnt casses:
    #   1) acf as a function of choice of tau for a single file
    #   2) acf with a set tau as a function of the sweep

    parser = argparse.ArgumentParser()
    # d is for directory
    parser.add_argument('-d',action='store',dest = 'd',type = str, required = True)
    # f is for file this is needed for make_acf_tau_plot()
    parser.add_argument('-f',action='store',dest = 'f',type = str, required = False)
    # plot type
    parser.add_argument('-t',action='store',dest = 't',type = str,required = True)
    # which variable are we looking for corrilations is
    parser.add_argument('-v',action='store',dest = 'v',type = str,required = False,default='x')

    inargs = parser.parse_args()
    d = inargs.d
    f = inargs.f
    v = inargs.v
    plot_type = inargs.t


    os.chdir(d)
    # This is also a single particle operation. Lets make a directory and have an individual image
    # for each run.

    if plot_type == 'acf_tau':
        make_acf_tau_plot(f,v)
    if plot_type == 'acf_sweep':
        make_acf_sweep_plot(v)



if __name__ == '__main__':
    main()

