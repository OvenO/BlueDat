#!/usr/bin/python
import pylab as pl
import os
import argparse
import o_funcs as of
import scipy.constants as constants

# sliced specific heat. using known zero potential of slices to find full energy variance and
# using A as tempature in calculation of specific heat. (pg 254 coputational physics book).
# C = (delta E)^2/kb*T^2
def ssheat():

    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()

    print('sweep must be over A!! sweep_str is: '+sweep_str)

    data_file_name = 'ssheat_data.txt'
    if data_file_name in os.listdir('.'):
        data_file = open(data_file_name,'r')
        # first line is labels
        labels = data_file.readline()
        plotting_data = pl.genfromtxt(data_file)
        #first column sweep variables
        var_arr = plotting_data[:,0]
        # evergy_stuff is next coulumn
        delta_E_sqrd = plotting_data[:,1]
        s_heat = plotting_data[:,2]

    else:
        data_file = open(data_file_name,'w')
        data_file.write('sweep_var   delta_E_sqrd   specific_heat \n')

        # how much of soluton do we want to use? 1 -> all, 0 -> none
        # this can be bigger than in the unsliced ones becasue transients are more or less gone after
        # several PCs anyway.
        how_much = .8

        delta_E_sqrd = pl.array([])
        s_heat = pl.array([])
        var_arr = pl.array([])
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
                # This is getting values of time that are at makimum potentials!!! WRONG
                # check_time = i*dt%(pl.pi*2.0)
                check_time = (i*dt+pl.pi/2.0)%(pl.pi*2.0)
                if check_time < dt and check_time > 0.0:
                    new_data = pl.append(new_data,cur_data[i,:])

            cur_data = new_data.reshape(-1,Dim*2*N)

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
            
            cur_en_stuff = to_sum/4
            delta_E_sqrd = pl.append(delta_E_sqrd,cur_en_stuff)
            
            cur_s_heat = cur_en_stuff/(N*cur_sweep_var**2)
            s_heat = pl.append(s_heat,cur_s_heat)
            
            next_line += str(cur_en_stuff)+'   '+str(cur_s_heat)+ '\n'

            data_file.write(next_line)
        
    fig = pl.figure()
    ax = fig.add_subplot(111)
    # form of errorbar(x,y,xerr=xerr_arr,yerr=yerr_arr)
    pl.scatter(var_arr,s_heat,c='k')
    #pl.errorbar(var_arr,averages_2,yerr=std_arr,c='b',ls='none',fmt='o')
    ax.set_xlabel(sweep_str,fontsize=30)
    ax.set_ylabel('Specific heat per particle',fontsize=20)
    #ax.set_ylim([0.0,0.2])
    fig.tight_layout()
    fig.savefig('specific_heat_per_particle.png',dpi=300)
    pl.close(fig)


    print('\a')
    os.system('say finnished calculating specific heat')

# sliced is boolean

def sheat_vs_tempature():

    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()

    if 'ssheat_data.txt' not in os.listdir('.'):
        print('Need specific heat data')
        os.system('say Need specific heat data')
    if 'temp_granular_sliced.txt' not in os.listdir('.'):
        print('Need granular tempature data')
        os.system('say Need granular tempature data')

    tempature_file = open('temp_granular_sliced.txt','r')
    sheat_file = open('ssheat_data.txt','r')

    # first line is labels
    tempature_labels = tempature_file.readline()
    tempature_plotting_data = pl.genfromtxt(tempature_file)
    tempature_arr = tempature_plotting_data[:,1]

    # first line is labels
    sheat_labels = sheat_file.readline()
    sheat_plotting_data = pl.genfromtxt(sheat_file)
    #first column sweep variables
    var_arr = sheat_plotting_data[:,0]
    # evergy_stuff is next coulumn
    delta_E_sqrd = sheat_plotting_data[:,1]
    s_heat_arr = sheat_plotting_data[:,2]

    fig = pl.figure()
    ax = fig.add_subplot(111)
    # form of errorbar(x,y,xerr=xerr_arr,yerr=yerr_arr)
    pl.scatter(tempature_arr,s_heat_arr,c='k')
    #pl.errorbar(var_arr,averages_2,yerr=std_arr,c='b',ls='none',fmt='o')
    ax.set_xlabel(r'T_g',fontsize=30)
    ax.set_ylabel('Specific heat per particle',fontsize=20)
    fig.tight_layout()
    fig.savefig('T_vs_s_heat.png',dpi=300)
    pl.close(fig)


    print('\a')
    os.system('say finnished plotting tempature against specific heat')



def temp_granular(sliced):

    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()

    if sliced: 
        data_file_name = 'temp_granular_sliced.txt'
        save_str = 'granular_tempature_sliced.png'
    else: 
        data_file_name = 'temp_granular_not_sliced.txt'
        save_str = 'granular_tempature_not_sliced.png'

    if data_file_name in os.listdir('.'):
        data_file = open(data_file_name,'r')
        # first line is labels
        labels = data_file.readline()
        plotting_data = pl.genfromtxt(data_file)
        #first column sweep variables
        var_arr = plotting_data[:,0]
        # tempature is next coulumn
        temp_arr = plotting_data[:,1]

    else:
        data_file = open(data_file_name,'w')
        data_file.write('sweep_var   granular_tempature\n')

        # how much of soluton do we want to use? 1 -> all, 0 -> none
        # this can be bigger than in the unsliced ones becasue transients are more or less gone after
        # several PCs anyway.
        how_much = .6

        temp_arr = pl.array([])
        var_arr = pl.array([])
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

            if sliced:
                # slice the data so we only have data for values of t=pi(2*n + 1/2)
                new_data = pl.array([])
                for i in range(len(cur_data)):
                    # This is getting values of time that are at makimum potentials!!! WRONG
                    # check_time = i*dt%(pl.pi*2.0)
                    check_time = (i*dt+pl.pi/2.0)%(pl.pi*2.0)
                    if check_time < dt and check_time > 0.0:
                        new_data = pl.append(new_data,cur_data[i,:])

                cur_data = new_data.reshape(-1,Dim*2*N)

            if Dim==1:
                mag_vel_arr_sqrd = cur_data[int(-how_much*len(cur_data)):,:N]**2
            if Dim==2:
                mag_vel_arr_sqrd = cur_data[int(-how_much*len(cur_data)):,:N]**2+cur_data[int(-how_much*len(cur_data)):,N:2*N]**2
            if i==0: print('shape of mag_vel_arr_sqrd'+str(pl.shape(mag_vel_arr_sqrd)))
            
            cur_temp = mag_vel_arr_sqrd.sum()/N     
            temp_arr = pl.append(temp_arr,cur_temp)
                   
            next_line += str(cur_temp)+ '\n'

            data_file.write(next_line)

        
    fig = pl.figure()
    ax = fig.add_subplot(111)
    # form of errorbar(x,y,xerr=xerr_arr,yerr=yerr_arr)
    pl.scatter(var_arr,temp_arr,c='k')
    #pl.errorbar(var_arr,averages_2,yerr=std_arr,c='b',ls='none',fmt='o')
    ax.set_xlabel(sweep_str,fontsize=30)
    ax.set_ylabel(r'$T_g$',fontsize=30)
    fig.tight_layout()
    fig.savefig(save_str,dpi=300)
    pl.close(fig)

    data_file.close()


    print('\a')
    os.system('say finnished calculating granular tempature')


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

#def kenetic_energy():
#    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()
#    # how much of soluton do we want to use? 1 -> all 0 -> none
#    how_much = .2
#
#    averages= pl.array([])
#    var_arr = pl.array([])
#    # See paper j. chem phys., vol 120, No 1, 1 Jan 2004
#    energy_stuff = pl.array([])
#    for i,j in enumerate(os.listdir('.')):
#        if 'poindat.txt' not in j:
#            continue
#        cur_file = open(j,'r')
#        cur_sweep_var = float(cur_file.readline().split()[-1])
#        cur_data=pl.genfromtxt(cur_file)
#        cur_file.close()
#
#        var_arr = pl.append(var_arr,cur_sweep_var)
#
#        if Dim==1:
#            mag_vel_arr_sqrd = cur_data[int(-how_much*len(cur_data)):,:N]**2
#        if Dim==2:
#            mag_vel_arr_sqrd = cur_data[int(-how_much*len(cur_data)):,:N]**2+cur_data[int(-how_much*len(cur_data)):,N:2*N]**2
#        if i==0: print('shape of mag_vel_arr_sqrd'+str(pl.shape(mag_vel_arr_sqrd)))
#        
#        to_sum = 0.0
#        to_sum_avg_en = 0.0
#        for a in range(len(mag_vel_arr_sqrd[0,:])):
#            for b in range(len(mag_vel_arr_sqrd[0,:])):
#                first = (mag_vel_arr_sqrd[:,a]*mag_vel_arr_sqrd[:,b]).mean()
#                second = (mag_vel_arr_sqrd[:,a].mean())*(mag_vel_arr_sqrd[:,b].mean())
#                
#                to_sum += first - second
#                to_sum_avg_en += second
#        
#        #print('to_sum: ' +str(to_sum))
#        #print('to_sum_avg_en: ' +str(to_sum_avg_en))
#        energy_stuff = pl.append(energy_stuff,to_sum/to_sum_avg_en)
#
#        # this is more or less wrong. you could find the varience and the eaverage but for the
#        # enerygy you need to sum the KE of each particle. look as second varible.
#        #avg_mag_vel_sqrd = mag_vel_arr_sqrd.mean()
#        
#        averages = pl.append(averages,pl.sqrt(to_sum_avg_en))
#        
#    fig = pl.figure()
#    ax = fig.add_subplot(111)
#    pl.scatter(var_arr,averages/2)
#    #ax.set_xlim([.5,1.5])
#    ax.set_xlabel(sweep_str,fontsize=30)
#    ax.set_ylabel(r'$\langle v^2 \rangle / 2$',fontsize=30)
#    fig.tight_layout()
#    fig.savefig('v_sqrd_avg.png',dpi=300)
#    pl.close(fig)
#
#    fig = pl.figure()
#    ax = fig.add_subplot(111)
#    # deided by 4 comes from KE -> (1/2)**2
#    pl.scatter(var_arr,energy_stuff/4)
#    #ax.set_xlim([.5,1.5])
#    ax.set_xlabel(sweep_str,fontsize=30)
#    ax.set_ylabel(r'$\frac{\langle KE^2 \rangle - \langle KE \rangle ^ 2}{\langle KE \rangle ^ 2}$',fontsize=30)
#    fig.tight_layout()
#    fig.savefig('energy_stuff.png',dpi=300)
#    pl.close(fig)

def frequency(f):
    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()


# We know that the energy of the system when t=n pi/2  is only KE because the potential U(x,t) is
# flat at those times. We need to slice at t=pi(2*n + 1/2) in order to actualy et PC sections though. 
def energy_fluctuation(keyword):

    avg_save_str = 'avg_energy_stuff'
    en_stuff_save_str = 'energy_stuff_'

    if keyword == 'slice':
        average_y_lbl = r"$ \langle E' \rangle $"
        en_stuff_y_lbl = r"$\frac{\langle E' ^2 \rangle - \langle E' \rangle ^ 2}{\langle E' \rangle ^ 2}$"
        avg_save_str = 'sliced_E_'+avg_save_str
        en_stuff_save_str = 'sliced_E_' +en_stuff_save_str
    if keyword == 'ke':
        average_y_lbl = r'$ \langle KE \rangle $'
        en_stuff_y_lbl = r'$\frac{\langle KE^2 \rangle - \langle KE \rangle ^ 2}{\langle KE \rangle ^ 2}$'

    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()

    if keyword == 'slice': data_file_name = 'sliced_energy_data.txt'
    if keyword == 'ke': data_file_name = 'ke_energy_data.txt'
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

            if keyword == 'slice':
                # slice the data so we only have data for values of t=pi(2*n + 1/2)
                new_data = pl.array([])
                for i in range(len(cur_data)):
                    # This is getting values of time that are at makimum potentials!!! WRONG
                    # check_time = i*dt%(pl.pi*2.0)
                    # This is right
                    check_time = (i*dt+pl.pi/2.0)%(pl.pi*2.0)
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

            
            cur_av_1 = pl.sqrt(to_sum_avg_en/4/N)
            averages_1 = pl.append(averages_1,cur_av_1)

            next_line += str(cur_en_stuff_1)+'   '+str(cur_en_stuff_2)+'   '+str(cur_av_1)+ \
                    '   '+str(cur_av_2)+'   '+str(cur_std)+'\n'
            data_file.write(next_line)
        
    fig = pl.figure()
    ax = fig.add_subplot(111)
    # form of errorbar(x,y,xerr=xerr_arr,yerr=yerr_arr)
    pl.scatter(var_arr,averages_1,c='k')
    #pl.errorbar(var_arr,averages_2,yerr=std_arr,c='b',ls='none',fmt='o')
    ax.set_xlabel(sweep_str,fontsize=30)
    ax.set_ylabel(average_y_lbl,fontsize=30)
    fig.tight_layout()
    fig.savefig(avg_save_str+'.png',dpi=300)
    pl.close(fig)

    fig = pl.figure()
    ax = fig.add_subplot(111)
    # deided by 4 comes from KE -> (1/2)**2
    pl.scatter(var_arr,energy_stuff_1,c='r')
    #pl.scatter(var_arr,energy_stuff_2,c='b')
    #ax.set_xlim([0.0,1.6])
    #ax.set_ylim([0.0,.04])
    #ax.set_xlim([.5,1.5])
    #ax.set_ylim([0,.03])
    ax.set_xlabel(sweep_str,fontsize=30)
    ax.set_ylabel(en_stuff_y_lbl,fontsize=30)
    fig.tight_layout()
    fig.savefig(en_stuff_save_str+'1.png',dpi=300)
    pl.close(fig)
    
    fig = pl.figure()
    ax = fig.add_subplot(111)
    # deided by 4 comes from KE -> (1/2)**2
    pl.scatter(var_arr,energy_stuff_2,c='b')
    ax.set_xlabel(sweep_str,fontsize=30)
    ax.set_ylabel(en_stuff_y_lbl,fontsize=30)
    #ax.set_xlim([.5,1.5])
    #ax.set_ylim([0,.03])
    fig.tight_layout()
    fig.savefig(en_stuff_save_str+'2.png',dpi=300)
    pl.close(fig)

    print('\a')
    
def spatio_temporal():

    os.mkdir('SpatioTemporalVels')

    print('RIGHT NOW THIS IS ONLY FOR VX!!!!!!!')
    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str,Dim  = of.get_system_info()

    p_arr = pl.arange(0,N)
    

    # How many cycles do we want to look at?
    how_many = 10

    var_arr = pl.array([])
    for i,j in enumerate(os.listdir('.')):
        if 'poindat.txt' not in j:
            continue
        print('working on file ' + j)
        poin_num = int(j[:j.find('p')])
        cur_file = open(j,'r')
        cur_sweep_var = float(cur_file.readline().split()[-1])
        cur_data=pl.genfromtxt(cur_file)
        cur_file.close()

        var_arr = pl.append(var_arr,cur_sweep_var)
        
        count = 0
        grid = cur_data[-int(how_many*2.0*pl.pi/dt):,:N]

        # in 1D because particles never cross eachother we can order them in the images to mathch
        # their physical order.
        grid_ordered = pl.zeros(pl.shape(grid))
        # can just use the initial conditions to figure out where each is
        init_x = cur_data[0,N:2*N]
        sorted_x = sorted(init_x)
        for a,alpha in enumerate(sorted_x):
            for b,beta in enumerate(init_x):
                if alpha == beta:
                    grid_ordered[:,a]=grid[:,b]
        
    
        print('shape of grid_ordered: ' + str(pl.shape(grid_ordered)))
        
        fig = pl.figure()
        ax = fig.add_subplot(111)
        # form of errorbar(x,y,xerr=xerr_arr,yerr=yerr_arr)
        ax.imshow(grid_ordered,interpolation="nearest", aspect='auto')
        ax.set_xlabel('Particle',fontsize=30)
        #ax.set_aspect('equal')
        ax.set_ylabel(r'$ t $',fontsize=30)
        fig.tight_layout()
        fig.savefig('SpatioTemporalVels/%(number)04d.png'%{'number':poin_num})
        pl.close(fig)



        # slice the data so we only have data for values of t=pi(2*n + 1/2)

#        new_data = pl.array([])
#        for i in range(len(cur_data)):
#            check_time = i*dt%(pl.pi*2.0)
#            if check_time < dt and check_time > 0.0:
#                new_data = pl.append(new_data,cur_data[i,:])
#
#        cur_data = new_data.reshape(-1,Dim*2*N)
#
#        if Dim==1:
#            mag_vel_arr_sqrd = cur_data[int(-how_much*len(cur_data)):,:N]**2
#        if Dim==2:
#            mag_vel_arr_sqrd = cur_data[int(-how_much*len(cur_data)):,:N]**2+cur_data[int(-how_much*len(cur_data)):,N:2*N]**2
#        if i==0: print('shape of mag_vel_arr_sqrd'+str(pl.shape(mag_vel_arr_sqrd)))
#        
#        cur_en_stuff_2 = N*((mag_vel_arr_sqrd**2).mean() - (mag_vel_arr_sqrd.mean())**2)/((mag_vel_arr_sqrd.mean())**2)/4
#        energy_stuff_2 = pl.append(energy_stuff_2, cur_en_stuff_2)
#
#        cur_av_2 = pl.sqrt(N*mag_vel_arr_sqrd.mean()**2/4)
#        averages_2 = pl.append(averages_2,cur_av_2)
#        cur_std = pl.sqrt((N*mag_vel_arr_sqrd**2/4).std())
#        std_arr = pl.append(std_arr,cur_std)
#        
#        to_sum = 0.0
#        to_sum_avg_en = 0.0
#        for a in range(len(mag_vel_arr_sqrd[0,:])):
#            for b in range(len(mag_vel_arr_sqrd[0,:])):
#                first = (mag_vel_arr_sqrd[:,a]*mag_vel_arr_sqrd[:,b]).mean()
#                second = (mag_vel_arr_sqrd[:,a].mean())*(mag_vel_arr_sqrd[:,b].mean())
#                
#                to_sum += first - second
#                to_sum_avg_en += second
#        
#        #print('to_sum: ' +str(to_sum))
#        #print('to_sum_avg_en: ' +str(to_sum_avg_en))
#        cur_en_stuff_1 = to_sum/to_sum_avg_en/4
#        energy_stuff_1 = pl.append(energy_stuff_1,cur_en_stuff_1)
#
#        
#        cur_av_1 = pl.sqrt(to_sum_avg_en/4)
#        averages_1 = pl.append(averages_1,cur_av_1)
#
#        next_line += str(cur_en_stuff_1)+'   '+str(cur_en_stuff_2)+'   '+str(cur_av_1)+ \
#                '   '+str(cur_av_2)+'   '+str(cur_std)+'\n'
#        data_file.write(next_line)
#        
#    fig = pl.figure()
#    ax = fig.add_subplot(111)
#    # form of errorbar(x,y,xerr=xerr_arr,yerr=yerr_arr)
#    pl.scatter(var_arr,averages_1,c='r')
#    pl.errorbar(var_arr,averages_2,yerr=std_arr,c='b',ls='none',fmt='o')
#    ax.set_xlabel(sweep_str,fontsize=30)
#    ax.set_ylabel(r'$ \langle E \rangle $',fontsize=30)
#    fig.tight_layout()
#    fig.savefig('sliced_E_avg.png',dpi=300)
#    pl.close(fig)
#
#    fig = pl.figure()
#    ax = fig.add_subplot(111)
#    # deided by 4 comes from KE -> (1/2)**2
#    pl.scatter(var_arr,energy_stuff_1,c='r')
#    #pl.scatter(var_arr,energy_stuff_2,c='b')
#    #ax.set_xlim([0.0,1.6])
#    #ax.set_ylim([0.0,.04])
#    #ax.set_xlim([.5,1.5])
#    ax.set_xlabel(sweep_str,fontsize=30)
#    ax.set_ylabel(r'$\frac{\langle E^2 \rangle - \langle E \rangle ^ 2}{\langle E \rangle ^ 2}$',fontsize=30)
#    fig.tight_layout()
#    fig.savefig('sliced_energy_stuff_1.png',dpi=300)
#    pl.close(fig)
#    
#    fig = pl.figure()
#    ax = fig.add_subplot(111)
#    # deided by 4 comes from KE -> (1/2)**2
#    pl.scatter(var_arr,energy_stuff_2,c='b')
#    ax.set_xlabel(sweep_str,fontsize=30)
#    ax.set_ylabel(r'$\frac{\langle E^2 \rangle - \langle E \rangle ^ 2}{\langle E \rangle ^ 2}$',fontsize=30)
#    fig.tight_layout()
#    fig.savefig('sliced_energy_stuff_2.png',dpi=300)
#    pl.close(fig)
#
#    print('\a')

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
        print('calling energy_fluctuation(ke)')
        energy_fluctuation('ke')
    if plot_type == 'frequency':
        print('calling frequency')
        frequency(f)
    # We know that the energy of the system when t=n pi/2  is only KE because the potential U(x,t) is
    # flat at those times. We are going to aslo plot these results as a way fo observing total
    # energy
    if plot_type == 'sliced_E':
        print('calling energy_fluctuation(slice)')
        energy_fluctuation('slice')
    if plot_type =='st':
        print('calling spatio_temporal')
        spatio_temporal()
    # sliced specific heat. using known zero potential of slices to find full energy variance and
    # using A as tempature in calculation of specific heat. (pg 254 coputational physics book).
    # C = (delta E)^2/kb*T^2
    if plot_type == 'ssheat':
        print('calling ssheat()')
        ssheat()
    if plot_type == 'temp_sliced':
        temp_granular(True)
    if plot_type == 'temp_not_sliced':
        temp_granular(False)
    if plot_type == 'sheat_v_temp':
        sheat_vs_tempature()

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
#        fig.savefig(to_save_dir+'/%(number)04d.png'%{'number':cur_poin_num})
#        pl.close(fig)

