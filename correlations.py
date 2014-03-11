#!/usr/bin/python
import pylab as pl
import os
import argparse
import o_funcs as of
import scipy.integrate as scint

# autocorrelation fuction. tau is the comonly used variable for the length of time one is looking
# for correlations over. tau is going to be just an integer. If a specific value of time is desired
# it needs to be turned into an integer via deviding by dt. here is only an integer
def acf(tau,data):
    # data has already been passed such that is is a 1d array of time serries form. Be carefule with
    # the input then!
    return sum(data[:-tau]*data[tau:])/(len(data)-tau)
    print('len of auto correlation interval: ' +str(len(data[:-tau]*data[tau:])))
    #return sum(data[:-tau]*data[tau:])
    
    #return sum(data[t]*data[t+tau:])/len(data[tau:])

# Big Note: Becasue we cannot keep track of which particle is which (if there can even be a
# conection) we need to integrate over average acf_tau where the average is over the N particles.
def make_int_c_sqrd_plot(v):

    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()
    
    y_lbl = r'$\int C^2(\tau)$'
    x_lbl = sweep_str

    sweep_var_arr = pl.array([])
    int_c_sqrd_arr = pl.array([])
    # loop over all files
    for i,j in enumerate(os.listdir('.')):
        if 'poindat.txt' not in j:
            continue
        
        work_file = open(j,'r')
        sweep_var_arr = pl.append(sweep_var_arr,float(work_file.readline().split()[-1]))
        data = pl.genfromtxt(work_file)
        work_file.close()

        average_out = pl.array([])

        # average over every particle in the simulation
        for a in range(N):
            
            if v == 'x':
                input_arr = data[:,Dim*N+a]
            if v == 'vx':
                input_arr = data[:,a]
            # Will only get into these if asking for 2D stuff anyway
            if v == 'y':
                if Dim==1:
                    print('No y in 1D')
                    quit()
                input_arr = data[:,(Dim+1)*N+a]
            if v == 'vy':
                if Dim==1:
                    print('No y in 1D')
                    quit()
                input_arr = data[:,(Dim-1)*N+a]

            print('shape of input_arr: ' + str(pl.shape(input_arr)))

            # lets try this for t from 0 to 100 cycles
            #t_arr = pl.arange(0,int(100.0*2.0*pl.pi),int(2.0*pl.pi))
            tau_arr = pl.arange(1,int(cycles*2.0*pl.pi),1)
            output_arr = pl.array([])
            for i,j in enumerate(tau_arr):
                cur = acf(j,input_arr)
                #print('current run of acf is (should be one number): ' +str(cur))
                
                output_arr = pl.append(output_arr,cur)

            # for average acf plot
            if a == 0:
                average_out = pl.append(average_out,output_arr)
            else:
                average_out += output_arr

        average_out = average_out/a
        print('shape of average_out (should be 1D array): ' + str(pl.shape(average_out)))

        # Romberg Integration integrate averages**2 to get int c^2 number
        #average_out = average_out[:257]
        #print('shape of :2**13+1 averag out: ' +str(len(average_out)))
        #int_c_sqrd = scint.romb(average_out**2,show=True)

        # simpson integration
        int_c_sqrd = scint.simps(average_out, x=None, dx=1, axis=-1, even='avg')
        int_c_sqrd_arr = pl.append(int_c_sqrd_arr, int_c_sqrd)
        print('int_c_sqrd_arr (should be number): ' + str(int_c_sqrd_arr))

         
            
    fig = pl.figure()
    ax = fig.add_subplot(111)
    #ax.set_xlim(x_rng)
    #ax.set_ylim(y_rng)
    ax.set_xlabel(x_lbl,fontsize=30)
    ax.set_ylabel(y_lbl,fontsize=30)
    ax.scatter(sweep_var_arr,int_c_sqrd_arr,c='k')
    fig.tight_layout()
    fig.savefig('int_c_sqrd_vs_sweep.png')
    pl.close(fig)

    print('\a')
    print('\a')
    os.system('open int_c_sqrd_vs_sweep.png')

def working_on_making_int_c_sqrd_treating_in_correct_way():
    # When v = 'all' it means we are looking at the magnitude of the velocities as seen from
    # looking at the system as a vector with 2*N*D degrees of freedom. eg. if there are 2 in
    # a 1D potential particles we will be looking for auto correlations in sqrt(vx_1^2 +vx_2^2)
    #if v == 'all':


def make_acf_tau_plot(f,v):

    cur_poin_num = int(f[:f.find('p')])
    if (str(cur_poin_num)+'RunImages') not in os.listdir('.'):
        os.mkdir(str(cur_poin_num)+'RunImages')

    to_save_dir = str(cur_poin_num)+'RunImages/TauACF_Var_'+v
    os.mkdir(to_save_dir)
    x_lbl = r'$\tau$'
    y_lbl = 'Autocorrelation'

    work_file = open(f,'r')
    sweep_var = float(work_file.readline().split()[-1])
    data = pl.genfromtxt(work_file)
    work_file.close()


    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()

    # lets also make a plot of the averages of all the autocorolations for each particle ... so far
    # they are all pretty simular looking anyway. I'm currious if I can use the averages for the
    # c_sqrd plot (see c_sqrd).
    average_out = pl.array([])

    # make em for every particle in the simulation
    for a in range(N):

        if v == 'x':
            input_arr = data[:,Dim*N+a]
        if v == 'vx':
            input_arr = data[:,a]
        # Will only get into these if asking for 2D stuff anyway
        if v == 'y':
            if Dim==1:
                print('No y in 1D')
                quit()
            input_arr = data[:,(Dim+1)*N+a]
        if v == 'vy':
            if Dim==1:
                print('No y in 1D')
                quit()
            input_arr = data[:,(Dim-1)*N+a]

        print('shape of input_arr: ' + str(pl.shape(input_arr)))

        # lets try this for t from 0 to 100 cycles
        #t_arr = pl.arange(0,int(100.0*2.0*pl.pi),int(2.0*pl.pi))
        tau_arr = pl.arange(1,int(100.0*2.0*pl.pi),1)
        output_arr = pl.array([])
        for i,j in enumerate(tau_arr):
            cur = acf(j,input_arr)
            print('current run of acf is (should be one number): ' +str(cur))
            output_arr = pl.append(output_arr,cur)

        # for average acf plot
        if a ==0:
            average_out = pl.append(average_out,output_arr)
        else:
            average_out += output_arr
     
        
        fig = pl.figure()
        ax = fig.add_subplot(111)
        #ax.set_xlim(x_rng)
        #ax.set_ylim(y_rng)
        ax.set_xlabel(x_lbl,fontsize=30)
        ax.set_ylabel(y_lbl,fontsize=30)
        ax.scatter(tau_arr*dt,output_arr,c='k')
        fig.tight_layout()
        fig.savefig(to_save_dir+'/%(number)04d.png'%{'number':a})
        pl.close(fig)

    average_out = average_out/a
    fig=pl.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(x_lbl,fontsize=30)
    ax.set_ylabel(r'$\langle$'+y_lbl+r'$\rangle$',fontsize=30)
    ax.scatter(tau_arr*dt,average_out,c='k')
    fig.tight_layout()
    fig.savefig(to_save_dir+'/average_all.png')
    pl.close(fig)

    print('\a')
    print('\a')
    os.system('open '+to_save_dir+'/*.png')

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
    parser.add_argument('-d',action='store',dest = 'd',type = str, required = False, default = './')
    # f is for file this is needed for make_acf_tau_plot()
    parser.add_argument('-f',action='store',dest = 'f',type = str, required = False)
    # plot type
    parser.add_argument('-t',action='store',dest = 't',type = str,required = True)
    # which variable are we looking for corrilations is
    parser.add_argument('-v',action='store',dest = 'v',type = str,required = False,default='all')

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
    if plot_type == 'int_c_sqrd':
        make_int_c_sqrd_plot(v)


if __name__ == '__main__':
    main()

