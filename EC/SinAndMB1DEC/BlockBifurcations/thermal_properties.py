import pylab as pl
import os
import argparse
import plot_one as po
import scipy.constants as constants

def diffusion_coef(f):
    # depending on the variable we want we need a number that controles how we slice the data
    # lines. With the dimension of the system we can slice everything right. Can get the dimension with the
    # following line.
    Dim = (pl.shape(data)[1])/(2*N)
    print('dimension: '+str(Dim))


#    cur_poin_num = int(f[:f.find('p')])
#    if (str(cur_poin_num)+'RunImages') not in os.listdir('.'):
#        os.mkdir(str(cur_poin_num)+'RunImages')

    if Dim == 2:
        qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str  = po.get_system_info()
    if Dim == 1:



#    to_save_dir = str(cur_poin_num)+'RunImages/TauACF_Var_'+v
#    os.mkdir(to_save_dir)

    work_file = open(f,'r')
    sweep_var = float(work_file.readline().split()[-1])
    data = pl.genfromtxt(work_file)
    work_file.close()

  

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
    # This is also a single particle operation. Lets make a directory and have an individual image
    # for each run.

    if plot_type == 'diffusion':
        diffusion_coef(f)

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

