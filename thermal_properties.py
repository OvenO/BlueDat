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
    # how much of soluton do we want to use? 1 -> all 0 -> none
    how_much = .2

    averages= pl.array([])
    var_arr = pl.array([])
    # See paper j. chem phys., vol 120, No 1, 1 Jan 2004
    energy_stuff = pl.array([])
    for i,j in enumerate(os.listdir('.')):
        if 'poindat.txt' not in j:
            continue
        cur_file = open(j,'r')
        cur_sweep_var = float(cur_file.readline().split()[-1])
        cur_data=pl.genfromtxt(cur_file)
        cur_file.close()

        var_arr = pl.append(var_arr,cur_sweep_var)

        if Dim==1:
            mag_vel_arr_sqrd = cur_data[int(-how_much*len(cur_data)):,:N]**2
        if Dim==2:
            mag_vel_arr_sqrd = cur_data[int(-how_much*len(cur_data)):,:N]**2+cur_data[int(-how_much*len(cur_data)):,N:2*N]**2
        if i==0: print('shape of mag_vel_arr_sqrd'+str(pl.shape(mag_vel_arr_sqrd)))
        
        to_sum = 0.0
        to_sum_avg_en = 0.0
        for a in range(len(mag_vel_arr_sqrd[0,:])):
            for b in range(len(mag_vel_arr_sqrd[0,:])):
                first = (mag_vel_arr_sqrd[:,a]*mag_vel_arr_sqrd[:,b]).mean()
                second = (mag_vel_arr_sqrd[:,a].mean())*(mag_vel_arr_sqrd[:,b].mean())
                
                to_sum += first - second
                to_sum_avg_en += second
        
        #print('to_sum: ' +str(to_sum))
        #print('to_sum_avg_en: ' +str(to_sum_avg_en))
        energy_stuff = pl.append(energy_stuff,to_sum/to_sum_avg_en)

        # this is more or less wrong. you could find the varience and the eaverage but for the
        # enerygy you need to sum the KE of each particle. look as second varible.
        #avg_mag_vel_sqrd = mag_vel_arr_sqrd.mean()
        
        averages = pl.append(averages,pl.sqrt(to_sum_avg_en))
        
    fig = pl.figure()
    ax = fig.add_subplot(111)
    pl.scatter(var_arr,averages/2)
    #ax.set_xlim([.5,1.5])
    ax.set_xlabel(sweep_str,fontsize=30)
    ax.set_ylabel(r'$\langle v^2 \rangle / 2$',fontsize=30)
    fig.tight_layout()
    fig.savefig('v_sqrd_avg.png',dpi=300)
    pl.close(fig)

    fig = pl.figure()
    ax = fig.add_subplot(111)
    # deided by 4 comes from KE -> (1/2)**2
    pl.scatter(var_arr,energy_stuff/4)
    #ax.set_xlim([.5,1.5])
    ax.set_xlabel(sweep_str,fontsize=30)
    ax.set_ylabel(r'$\frac{\langle KE^2 \rangle - \langle KE \rangle ^ 2}{\langle KE \rangle ^ 2}$',fontsize=30)
    fig.tight_layout()
    fig.savefig('energy_stuff.png',dpi=300)
    pl.close(fig)

def frequency(f):
    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()


# We know that the energy of the system when t=n pi/2  is only KE because the potential U(x,t) is
# flat at those times. We need to slice at t=pi(2*n + 1/2) in order to actualy et PC sections though. 
def zero_potential_pc():
    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()

    data_file_name = 'sliced_energy_data.txt'
    if data_file_name in os.listdir('.'):
        data_file = open(data_file_name,'r')
        # first line is labels
        labels = data_file.readline()
        plotting_data = pl.genfromtxt(data_file)
        #first column sweep variables
        var_arr = plotting_data[:,0]
        # evergy_stuff is next coulumn
        energy_stuff_1 = plotting_data[:,1]
        energy_stuff_2 = plotting_data[:,2]
        averages_1 = plotting_data[:,3]
        averages_2 = plotting_data[:,4]
        # std_arr is only for averages_2 
        std_arr = plotting_data[:,5]

    else:
        data_file = open(data_file_name,'w')
        data_file.write('sweep_var   energy_stuff_1   energy_stuff_2   averages_1   averages_2 standard_dev\n')

        # how much of soluton do we want to use? 1 -> all, 0 -> none
        # this can be bigger than in the unsliced ones becasue transients are more or less gone after
        # several PCs anyway.
        how_much = .8

        # See paper j. chem phys., vol 120, No 1, 1 Jan 2004
        # energy_stuff_1 and averages_1 do not asume that the velocites are independently
        # distributed. These use above paper eqn 14
        # energy_sruff_2 and averages_2 assume this is a thermal system and use above paper eqn 15
        energy_stuff_1 = pl.array([])
        energy_stuff_2 = pl.array([])
        var_arr = pl.array([])
        averages_1 = pl.array([])
        averages_2 = pl.array([])
        std_arr = pl.array([])
        for i,j in enumerate(os.listdir('.')):
            next_line = ''
            if 'poindat.txt' not in j:
                continue
            cur_file = open(j,'r')
            cur_sweep_var = float(cur_file.readline().split()[-1])
            cur_data=pl.genfromtxt(cur_file)
            cur_file.close()

            next_line += str(cur_sweep_var)+ '   '
            var_arr = pl.append(var_arr,cur_sweep_var)

            # slice the data so we only have data for values of t=pi(2*n + 1/2)
            new_data = pl.array([])
            for i in range(len(cur_data)):
                check_time = i*dt%(pl.pi*2.0)
                if check_time < dt and check_time > 0.0:
                    new_data = pl.append(new_data,cur_data[i,:])

            cur_data = new_data.reshape(-1,Dim*2*N)

            if Dim==1:
                mag_vel_arr_sqrd = cur_data[int(-how_much*len(cur_data)):,:N]**2
            if Dim==2:
                mag_vel_arr_sqrd = cur_data[int(-how_much*len(cur_data)):,:N]**2+cur_data[int(-how_much*len(cur_data)):,N:2*N]**2
            if i==0: print('shape of mag_vel_arr_sqrd'+str(pl.shape(mag_vel_arr_sqrd)))
            
            cur_en_stuff_2 = N*((mag_vel_arr_sqrd**2).mean() - (mag_vel_arr_sqrd.mean())**2)/((mag_vel_arr_sqrd.mean())**2)/4
            energy_stuff_2 = pl.append(energy_stuff_2, cur_en_stuff_2)

            cur_av_2 = pl.sqrt(N*mag_vel_arr_sqrd.mean()**2/4)
            averages_2 = pl.append(averages_2,cur_av_2)
            cur_std = pl.sqrt((N*mag_vel_arr_sqrd**2/4).std())
            std_arr = pl.append(std_arr,cur_std)
            
            to_sum = 0.0
            to_sum_avg_en = 0.0
            for a in range(len(mag_vel_arr_sqrd[0,:])):
                for b in range(len(mag_vel_arr_sqrd[0,:])):
                    first = (mag_vel_arr_sqrd[:,a]*mag_vel_arr_sqrd[:,b]).mean()
                    second = (mag_vel_arr_sqrd[:,a].mean())*(mag_vel_arr_sqrd[:,b].mean())
                    
                    to_sum += first - second
                    to_sum_avg_en += second
            
            #print('to_sum: ' +str(to_sum))
            #print('to_sum_avg_en: ' +str(to_sum_avg_en))
            cur_en_stuff_1 = to_sum/to_sum_avg_en/4
            energy_stuff_1 = pl.append(energy_stuff_1,cur_en_stuff_1)

            
            cur_av_1 = pl.sqrt(to_sum_avg_en/4)
            averages_1 = pl.append(averages_1,cur_av_1)

            next_line += str(cur_en_stuff_1)+'   '+str(cur_en_stuff_2)+'   '+str(cur_av_1)+ \
                    '   '+str(cur_av_2)+'   '+str(cur_std)+'\n'
            data_file.write(next_line)
        
    fig = pl.figure()
    ax = fig.add_subplot(111)
    # form of errorbar(x,y,xerr=xerr_arr,yerr=yerr_arr)
    pl.scatter(var_arr,averages_1,c='r')
    pl.errorbar(var_arr,averages_2,yerr=std_arr,c='b',ls='none',fmt='o')
    ax.set_xlabel(sweep_str,fontsize=30)
    ax.set_ylabel(r'$ \langle E \rangle $',fontsize=30)
    fig.tight_layout()
    fig.savefig('sliced_E_avg.png',dpi=300)
    pl.close(fig)

    fig = pl.figure()
    ax = fig.add_subplot(111)
    # deided by 4 comes from KE -> (1/2)**2
    pl.scatter(var_arr,energy_stuff_1,c='r')
    #pl.scatter(var_arr,energy_stuff_2,c='b')
    #ax.set_xlim([0.0,1.6])
    #ax.set_ylim([0.0,.04])
    #ax.set_xlim([.5,1.5])
    ax.set_xlabel(sweep_str,fontsize=30)
    ax.set_ylabel(r'$\frac{\langle E^2 \rangle - \langle E \rangle ^ 2}{\langle E \rangle ^ 2}$',fontsize=30)
    fig.tight_layout()
    fig.savefig('sliced_energy_stuff_1.png',dpi=300)
    pl.close(fig)
    
    fig = pl.figure()
    ax = fig.add_subplot(111)
    # deided by 4 comes from KE -> (1/2)**2
    pl.scatter(var_arr,energy_stuff_2,c='b')
    ax.set_xlabel(sweep_str,fontsize=30)
    ax.set_ylabel(r'$\frac{\langle E^2 \rangle - \langle E \rangle ^ 2}{\langle E \rangle ^ 2}$',fontsize=30)
    fig.tight_layout()
    fig.savefig('sliced_energy_stuff_2.png',dpi=300)
    pl.close(fig)

    print('\a')

def main():

    parser = argparse.ArgumentParser()
    # d is for directory
    parser.add_argument('-d',action='store',dest = 'd',type = str, required = False, default = './')
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
        print('calling kenetic_energy()')
        kenetic_energy()
    if plot_type == 'frequency':
        print('calling frequency')
        frequency(f)
    # We know that the energy of the system when t=n pi/2  is only KE because the potential U(x,t) is
    # flat at those times. We are going to aslo plot these results as a way fo observing total
    # energy
    if plot_type == 'sliced_E':
        zero_potential_pc()

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

