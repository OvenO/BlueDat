#!/usr/bin/python
import pylab as pl
import o_funcs as of
import argparse 
import os

def main():
    parser = argparse.ArgumentParser()
    # f is for file
    parser.add_argument('-f',action='store',dest = 'f',type = str, required = False)
    # plot type
    parser.add_argument('-t',action='store',dest = 't',type = str,required = True)
    # plot sub type
    parser.add_argument('--st',action='store',dest = 'st',type = str,required = False)
    # plot every image or everyother or ...
    # to plot evey image skip = 1. Skip cannot be less than 1 or else you get devide by zero
    # error in the moddulus.
    parser.add_argument('-s',action='store',dest = 's',type = int,required = False,default = 2)

    inargs = parser.parse_args()
    f = inargs.f
    plot_type = inargs.t
    print('plot_type is: ' + str(plot_type))
    plot_sub_type = inargs.st
    skip = inargs.s

    # ancl --> anal class
    ancl = of.anal_run()
    ancl.get_info()
    ancl.set_list_dir()

    # how many cycles do we want to get trow away as transients
    # i.e. what do we think transients are
    how_many_get_rid = 50 

    if plot_type == 'one':
        working_file = open(f,"r")

        # get the varible for the plot
        var = float(working_file.readline().split()[-1])
        sol = pl.genfromtxt(working_file)
        working_file.close()

        # puts the file number ahead of the RunImages folder to prevent overwriting
        f_num_str = f[:f.find('p')]

        print('making directory: ' + f_num_str+'RunImages/PhaseSpace')
        os.mkdir(f_num_str+'RunImages')
        os.mkdir(f_num_str+'RunImages/PhaseSpace')

        for i in range(len(sol)):
            print('in make_one loop')
            if i%skip!=0:
                continue
            print(i)
            print(sol[i,:ancl.N])

            r_fig = pl.figure()
            r_ax = r_fig.add_subplot(111)
            r_ax.set_xlim([-pl.pi/ancl.N,2.0*pl.pi-pl.pi/ancl.N])
            r_ax.set_ylim([-1,1])
            r_ax.set_xlabel(r'$\phi$',fontsize=30)
            r_ax.set_ylabel(r'$\dot{\phi}$',fontsize=30)
            # Need to plot phi. not theta
            for j in range(ancl.N):
                # these two lines are just to make the phi more readable
                theta_j = sol[i,ancl.N+j]
                theta_dot_j = sol[i,j]

                
                phi_j = (pl.pi/ancl.N)*(j*2+pl.sin(theta_j)) 
                phi_dot_j = (pl.pi/ancl.N)*pl.cos(theta_j)*theta_dot_j

                r_ax.scatter(phi_j,phi_dot_j)

            r_fig.tight_layout()
            r_fig.savefig(f_num_str+'RunImages/PhaseSpace/%(number)04d.png'%{'number':i})
            pl.close(r_fig)
        
    # plot the velocity distribution for each file
    if plot_type == 'veldist':
        to_save_dir = 'VelDistMovie_neg100-neg50'
        os.mkdir(to_save_dir)

        cycle_int = int(2.0*pl.pi/ancl.dt)
        for i,j in enumerate(ancl.list_dir):
            working_file = open(j,'r')
            cur_sweep_var = float(working_file.readline().split()[-1])
            cur_data = pl.genfromtxt(working_file)
            working_file.close()
            
            # get rid of transient cycles
            cur_data = cur_data[-100*cycle_int:-50*cycle_int,:]
            # get one array of all the particles velocity maginudes
            theta_dot_arr = pl.array([])
            theta_dot_arr = pl.append(theta_dot_arr,cur_data[:,:ancl.N])
    
            fig = pl.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel(r'$\dot{\theta}$',fontsize=30)
            ax.set_ylabel('Count',fontsize = 25)
            ax.hist(theta_dot_arr,bins = 20)
            fig.tight_layout()
            fig.savefig(to_save_dir+'/%(number)04d.png'%{'number':i})
            pl.close(fig)
    # x magnitism
    # This may not be the right order parameter but I'm curious to see what this looks like. The
    # details of why this is an order parameter in the HMF method are in "Nonequilibrium statistical
    # mechanics of systems with long-range interactions" Physics reports, 2014.
    if plot_type == 'mx':
        #os.mkdir('MxMovie')
        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Cycle',fontsize=30)
        ax.set_ylabel(r'$M_x$',fontsize = 25)
        for i,j in enumerate(ancl.list_dir):
            working_file = open(j,'r')
            cur_sweep_var = float(working_file.readline().split()[-1])
            cur_data = pl.genfromtxt(working_file)
            working_file.close()
            
            # get rid of transient cycles
            cur_data = cur_data[how_many_get_rid:,:]
            
            # want to plot the magnitism as a function of time for like the last 10 cycles
            # make an array of the average values at a given time for differnt values of time
            theta_arr = pl.array([])
            for a in range(int((10*2*pl.pi)/ancl.dt)):
                theta_arr = pl.append(theta_arr,pl.cos(cur_data[-(a+1),ancl.N:]).sum()/ancl.N)
    
            ax.plot(pl.linspace(-10,0,int((10*2*pl.pi)/ancl.dt)),theta_arr)
        fig.tight_layout()
        #fig.savefig('MxMovie/%(number)04d.png'%{'number':i})
        fig.savefig('mx.png')
        pl.close(fig)


if __name__ == '__main__':
    main()
