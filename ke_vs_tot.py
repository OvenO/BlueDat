import pylab as pl
import os


def period_energies(ancl,keyword):

    # how many idecies do we want to use? One Period 
    num = int(2.0*pl.pi/ancl.dt)

    # make directory for images
    os.mkdir('EnergyLastPerCompair')
    
    for i,j in enumerate(ancl.list_dir):
        print('working with file ' + str(j))
        if 'poindat.txt' not in j:
            continue
        cur_file = open(j,'r')
        cur_sweep_var = float(cur_file.readline().split()[-1])
        cur_data=pl.genfromtxt(cur_file)
        cur_file.close()

        var_arr = pl.append(var_arr,cur_sweep_var)

        mag_vel_arr_sqrd = cur_data[-num:,:ancl.N]**2
        positions = cur_data[-num:,ancl.N:]

        #if ancl.Dim==2:
        #    mag_vel_arr_sqrd = cur_data[-int(2.0*pl.pi/ancl.dt):,:ancl.N]**2+cur_data[-int(2.0*pl.pi/ancl.dt):,ancl.N:2*ancl.N]**2
        if i==0: print('shape of mag_vel_arr_sqrd'+str(pl.shape(mag_vel_arr_sqrd)))
        
        print('len(mag_vel_arr_sqrd): ' +str(len(mag_vel_arr_sqrd)))
        
        ke = pl.array([]) #kinetic energy
        pot = pl.array([]) #potential from time dependent potential
        inter = pl.array([]) #interparticle potential

        # sum the potential energies of each particle and sum the kenetic energies of each particle
        for a in range(num): 
            ke = pl.append(ke,.5*mag_vel_arr_sqrd[a,:].sum())
            pot = pl.append(pot,pl.cos(positions[a,:])*pl.cos(ancl.dt*a).sum())

            inter = HEERE

        
        fig = pl.figure()
        ax = fig.add_subplot(111)
        # form of errorbar(x,y,xerr=xerr_arr,yerr=yerr_arr)
        pl.scatter(var_arr,averages_1,c='k')
        #pl.errorbar(var_arr,averages_2,yerr=std_arr,c='b',ls='none',fmt='o')
        ax.set_xlabel(ancl.sweep_str,fontsize=30)
        ax.set_ylabel(average_y_lbl,fontsize=30)
        fig.tight_layout()
        fig.savefig(avg_save_str+'.png',dpi=300)
        pl.close(fig)

        fig = pl.figure()
        ax = fig.add_subplot(111)
        pl.scatter(var_arr,energy_stuff_1,c='r')
        #pl.scatter(var_arr,energy_stuff_2,c='b')
        #ax.set_xlim([0.0,1.6])
        #ax.set_ylim([0.0,.04])
        #ax.set_xlim([.5,1.0])
        #ax.set_ylim([0,.24])
        ax.set_xlabel(ancl.sweep_str,fontsize=30)
        ax.set_ylabel(en_stuff_y_lbl,fontsize=30)
        fig.tight_layout()
        fig.savefig(en_stuff_save_str+'1.png',dpi=300)
        pl.close(fig)
    
#    fig = pl.figure()
#    ax = fig.add_subplot(111)
#    pl.scatter(var_arr,energy_stuff_2,c='b')
#    ax.set_xlabel(ancl.sweep_str,fontsize=30)
#    ax.set_ylabel(en_stuff_y_lbl,fontsize=30)
#    #ax.set_xlim([.5,1.5])
#    #ax.set_ylim([0,.03])
#    fig.tight_layout()
#    fig.savefig(en_stuff_save_str+'2.png',dpi=300)
#    pl.close(fig)

    print('\a')



def main():

if __name__ == '__main__':
    main()
