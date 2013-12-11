#!/usr/bin/python
import pylab as pl
import os
import argparse
import o_funcs as of
import scipy.constants as constants

def diffusion_coef(f):

#    cur_poin_num = int(f[:f.find('p')])
#    if (str(cur_poin_num)+'RunImages') not in os.listdir('.'):
#        os.mkdir(str(cur_poin_num)+'RunImages')

    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()

    print('dimension: '+str(Dim))

#    to_save_dir = str(cur_poin_num)+'RunImages/TauACF_Var_'+v
#    os.mkdir(to_save_dir)

    work_file = open(f,'r')
    print('working file is: ' +str(work_file))
    sweep_var = float(work_file.readline().split()[-1])
    print('sweep variable is: ' + str(sweep_var))
    data = pl.genfromtxt(work_file)
    print('shape of data is: ' + str(pl.shape(data)))
    work_file.close()

    # depending on the variable we want we need a number that controles how we slice the data
    # lines. With the dimension of the system we can slice everything right. Can get the dimension with the
    # following line.
  
    # For including the average of all. NOT DONE YET HERE
    diffusion_arr = pl.array([])
    # make em for every particle in the simulation
    for a in range(N):

        # can use the dimension to slice everything right
        cur_x = data[:,Dim*N+a]
        if Dim ==2:
            cur_y = data[:,(Dim+1)*N+a]
            distance_arr = pl.sqrt(cur_x**2 + cur_y**2)
        else:
            distance_arr = cur_x

        #print('shape of distance array: ' + str(pl.shape(input_arr)))
        
        # This equation generaly holds and for now this is how we are going to calculate the
        # diffusion coefficien:
        # <x^2>=q*D*t where q is numerical constant that depends on dimesionality. q=2*Dim. D is the
        # diffusion coefficient, and t is the time for which <x^2 is calculated>
        # this and some usefull equations from
        # http://www.life.illinois.edu/crofts/bioph354/diffusion1.html
        # So we can find D
        dist_arr_sqrd = distance_arr**2
        mean_sqrd_dist = dist_arr_sqrd.mean()
        total_time = len(distance_arr)*dt
        Diff_coef = mean_sqrd_dist/(2.0*Dim*total_time)
        
        print('Diffusion coefficient for particle ' + str(a) + ' = ' + str(Diff_coef))

        # Lets also try to print a tempature from the equation D = kT/f. T -> Absolute tempature.
        # k -> boltzman constat. f -> frictional constant (beta)
        Temp = Diff_coef * beta 
        
        print('Tempature from Diffusion coef = ' +str(Temp))
        

    print('\a')
    print('\a')

def kenetic_energy():
    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()

    averages= pl.array([])
    var_arr = pl.array([])
    for i,j in enumerate(os.listdir('.')):
        cur_file = open(j,'r')
        cur_sweep_var = cur_file.readline()
        cur_data=pl.genfromtxt(cur_file)
        cur_file.close()

        var_arr = pl.append(var_arr,cur_sweep_var)

        if Dim==1:
            mag_vel_arr_sqrd = cur_data[:,:N]**2
        if Dim==2:
            mag_vel_arr_sqrd = cur_data[:,:N]**2+cur_data[:,N:2*N]**2
        
        avg_mag_vel_sqrd = mag_vel_arr_sqrd.mean()
        averages = pl.append(averages,avg_mag_vel_sqrd)
        
    fig = pl.figure()
    ax = fig.add_subplot(111)
    pl.scatter(var_arr,averages)
    ax.set_xlabel(sweep_str,fontsize=30)
    ax.set_ylabel(r'$\langle v^2 \rangle$',fontsize=30)
    fig.tight_layout()
    fig.savefig('v_sqrd_avg.png',dpi=300)
    pl.close(fig)


def main():

    parser = argparse.ArgumentParser()
    # d is for directory
    parser.add_argument('-d',action='store',dest = 'd',type = str, required = True)
    # f is for file this is needed for make_acf_tau_plot()
    parser.add_argument('-f',action='store',dest = 'f',type = str, required = False)
    # plot type
    parser.add_argument('-t',action='store',dest = 't',type = str,required = True)

    inargs = parser.parse_args()
    d = inargs.d
    f = inargs.f
    plot_type = inargs.t


    os.chdir(d)
    print('changed directory to'+ d)
    # This is also a single particle operation. Lets make a directory and have an individual image
    # for each run.

    if plot_type == 'diffusion':
        print('calling diffusion_coef(f)')
        diffusion_coef(f)
    if plot_type == 'ke':
        kenetic_energy()
        print('calling kenetic_energy()')

if __name__ == '__main__':
    main()

#        fig = pl.figure()
#        ax = fig.add_subplot(111)
#        #ax.set_xlim(x_rng)
#        #ax.set_ylim(y_rng)
#        ax.set_xlabel(x_lbl,fontsize=30)
#        ax.set_ylabel(y_lbl,fontsize=30)
#        ax.scatter(tau_arr,output_arr,c='k')
#        fig.tight_layout()
#        fig.savefig(to_save_dir+'/%(number)04d.png'%{'number':a})
#        pl.close(fig)

