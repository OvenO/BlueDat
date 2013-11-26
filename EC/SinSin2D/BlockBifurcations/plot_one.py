import pylab as pl
import os
import argparse
import numpy as np

# Still crossing over program from 1D to 2D. 
# working features:
# 1) movie of real trajectories for single file (make_real)
# 2) movie of last section of trajectories of all files (lm)
# 3) t=0 poincare section (tzpc)


def get_system_info():
    # get the  system info
    info_f = open('info.txt','r')
    l = info_f.readlines()

    looking_for=['qq','dt','beta','A','cycles','particle number','x_num_cell','y_num_cell','order']
    values =    [ 0.0, 0.0, 0.0  ,0.0, 0.0    , 0.0             , 0.0        , 0.0        , 0.0   ]
    # make sure to type cast all of these right later (above)

    for i,j in enumerate(l):
        for a,b in enumerate(looking_for):
            if b in j:
                if 'sweep' in j:
                    values[a] = 'sweep'
                else: 
                    values[a] = float(j.split()[-1])
                
    qq         = values[0]               
    dt         = values[1]
    beta       = values[2]
    A          = values[3]
    cycles     = values[4]
    N          = values[5]
    x_num_cell = values[6]
    y_num_cell = values[7]
    order      = values[8]
    # type cast 
    if type(qq        ) != str: qq         = float(qq        )
    else: sweep_str = r'$q_i q_j$'
    if type(dt        ) != str: dt         = float(dt        )
    else: sweep_str = r'$\delta t$'
    if type(beta      ) != str: beta       = float(beta      )
    else: sweep_str = r'$\beta$'
    if type(A         ) != str: A          = float(A         )
    else: sweep_str = r'$A$'
    if type(cycles    ) != str: cycles     = int(  cycles    )
    else: sweep_str = '$Number of Cycles$'
    if type(N         ) != str: N          = int(  N         )
    else: sweep_str = r'$N$'
    if type(x_num_cell) != str: x_num_cell = int(  x_num_cell)
    if type(y_num_cell) != str: y_num_cell = int(  y_num_cell)
    if type(order     ) != str: order      = int(  order     )
    else: sweep_str = r'$O$'

    return qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str

def set_projection(projection,x_num_cell,y_num_cell):
    # when defining the projection string: <horizontal axis><vertical axis>. eg xvx -> vx phase
    # space with x as horizontal axis.
    # depending on the projection type we need a number that controles how we slice the data
    # lines. An array catered to the slicing (see plot to make more sence of it). For now
    # projection will not be anything pathological (i.e only velocity as verticasl
    # axis in x,y phase projections
    if projection == 'xy':
        ax_num = pl.array([2,3,3,4])
        # make strings for plots acordingly
        x_lbl = r'$x$'
        y_lbl = r'$y$'
        # set axis limits (ranges) acordingly
        x_rng = [0.0,2.0*pl.pi*x_num_cell]
        y_rng = [0.0,2.0*pl.pi*y_num_cell]
        
    if projection == 'xvx':
        ax_num = pl.array([2,3,0,1])
        ax_num = 1
        x_lbl = r'$x$'
        y_lbl = r'$\dot{x}$'
        x_rng = [0.0,2.0*pl.pi*x_num_cell]
        y_rng = [-2.5,2.5]
    if projection == 'yvy':
        ax_num = pl.array([3,4,1,2])
        ax_num = 2
        x_lbl = r'$y$'
        y_lbl = r'$\dot{y}$'
        x_rng = [0.0,2.0*pl.pi*y_num_cell]
        y_rng = [-2.5,2.5]
    if projection == 'vxvy':
        ax_num = pl.array([0,1,1,2])
        ax_num = 3
        x_lbl = r'$\dot{x}$'
        y_lbl = r'$\dot{y}$'
        x_rng = [-2.5,2.5]
        y_rng = [-2.5,2.5]

    return x_lbl,y_lbl,x_rng,y_rng,ax_num 

def main():
    # 'one option with a plane projection makes a time movies of mparticles motions on that plane
    make_one_movie = False
    # last movie (lm) Define projection type
    make_lm = False
    # t=0 poincare section
    make_tzpc = False
    # v=0 poincare section
    make_vzpc = False
    # plot final veloceties
    make_fv = False
    # plot final velocities for ALL OF THEM. (final velocity movie)
    make_fvm = False
    # plot velocity averages and standard deviatinos for now lets do averages over last half of run
    # times
    make_averages = False
    # plot autocorrelation values as variable is changed
    make_acf = False

    parser = argparse.ArgumentParser()
    # d is for directory
    parser.add_argument('-d',action='store',dest = 'd',type = str, required = True)
    # f is for file
    parser.add_argument('-f',action='store',dest = 'f',type = str, required = False)
    # plot type
    parser.add_argument('-t',action='store',dest = 't',type = str,required = True)
    # projection on to which axis
    parser.add_argument('-p',action='store',dest = 'p',type = str,required = False,default='x')
    inargs = parser.parse_args()
    d = inargs.d
    f = inargs.f
    plot_type = inargs.t
    projection = inargs.p




    # plot type
    if plot_type ==  'one':
        make_one_movie =  True
    if plot_type == 'lm':
        make_lm = True

    if plot_type =='fv':
        make_fv = True
        want_bins = 20
    if plot_type =='fvm':
        make_fvm = True
        want_bins = 20
    # plot velocity averages and standard deviatinos 
    if plot_type == 'averages':
        make_averages = True
        # for now lets do averages over last half of run times
        avg_over = .5

    # ultimitaly we want to work backwards from the final time so that if we want to we can cut off
    # transients but for now we are just going to look at the whole solution because
    # otherwise we are not going to have very many data points
    # throw away first 5th
    throw_away = 5
    if plot_type == 'vzpc':
        make_vzpc = True
    if plot_type == 'tzpc':
        make_tzpc = True
        


    os.chdir(d)

    qq,dt,beta,A,cycles,N,x_num_cell,y_num_cell,order,sweep_str = get_system_info()

    x_lbl,ylbl,x_rng,y_rng,ax_num = set_projection(projection,x_num_cell,y_num_cell)
    
    if make_lm or make_fvm or make_averages:
        count = 0
        var_arr = pl.array([])
        if make_lm:
            os.mkdir('LastMovie_'+projection)
        if make_fvm:
            os.mkdir('LastVelMovie')
            # lets see what the velocity disrobution is when the potential is zero
            # go_back =int(pl.pi/2.0/dt)
            go_back = -1
        if make_averages:
            avges_fig = pl.figure()
            avges_ax = avges_fig.add_subplot(111)
            avges_ax.set_xlabel(r'$A$',fontsize = 25)
            #avges_ax.set_ylabel(r'$\langle v \rangle \pm \sigma$')
            avges_ax.set_ylabel(r'$ \sigma ^2$', fontsize = 25)
            avg_vels = pl.array([])
            std_vels = pl.array([])
        for i,j in enumerate(os.listdir('.')):
            if 'poindat.txt' in j:
                cur_poin_num = int(j[:j.find('p')])
                cur_file = open(j,"r")
                # get the varible for the plot
                var = float(cur_file.readline().split()[-1])
                print(var)
                var_arr = pl.append(var_arr,var)

                sol = np.genfromtxt(cur_file)
                cur_file.close()

                # modulus by the lengths of the system
                sol[:,2*N:3*N] = sol[:,2*N:3*N]%(x_num_cell*2.0*pl.pi)
                sol[:,3*N:4*N] = sol[:,3*N:4*N]%(y_num_cell*2.0*pl.pi)
                
                if make_averages:
                    temp_v = pl.array([])
                    temp_v = pl.append(temp_v,sol[int(len(sol)*avg_over):,:N])
                    avg_vels = pl.append(avg_vels,temp_v.mean())
                    std_vels = pl.append(std_vels,temp_v.std())
                
                if make_lm:
                    lm_fig = pl.figure()
                    lm_ax = lm_fig.add_subplot(111)
                    lm_ax.set_xlabel(x_lbl,fontsize=30)
                    lm_ax.set_ylabel(y_lbl,fontsize=30)
                    # how much of the last part of the solution do you want to plot
                    amount = len(sol)-int(len(sol)/6.0)
                    lm_ax.scatter(sol[amount:,(ax_num[0])*N:(ax_num[1])*N],sol[amount:,(ax_num[2])*N:(ax_num[3])*N],c='k')
                    lm_fig.tight_layout()
                    lm_fig.savefig('LastMovie_'+projection+'/%(number)04d.png'%{'number':cur_poin_num})
                    pl.close(lm_fig)
                    
                if make_fvm:
                    fv_fig = pl.figure()
                    fv_ax = fv_fig.add_subplot(111)
                    fv_ax.hist(sol[-go_back,:N],bins=want_bins)
                    fv_ax.set_xlabel(r'$v$'      ,fontsize=30)
                    #fv_ax.set_ylabel(r'$$',fontsize=30)
                    #fv_ax.scatter(sol[i,N:2*N]%d,sol[i,:N],c='b')
                    fv_fig.tight_layout()
                    fv_fig.savefig('LastVelMovie/%(number)04d.png'%{'number':cur_poin_num})
                    pl.close(fv_fig)

                count+=1
        if make_averages:
            # plot the mean values
            #avges_ax.scatter(var_arr,avg_vels,c='b'         ,s=1)
            #avges_ax.scatter(var_arr,avg_vels+std_vels,c='r',s=1)
            #avges_ax.scatter(var_arr,avg_vels-std_vels,c='r',s=1)
            avges_ax.scatter(var_arr,std_vels**2,c='r',s=1)
            avges_fig.tight_layout()
            avges_fig.savefig('average_vels.png',dpi=300)

   
    else:
        cur_file = open(f,"r")
        # get the varible for the plot
        var = float(cur_file.readline().split()[-1])
        sol = np.genfromtxt(cur_file)
        cur_file.close()

        # modulus by the length of the system
        sol[:,2*N:3*N] = sol[:,2*N:3*N]%(x_num_cell*2.0*pl.pi)
        sol[:,3*N:4*N] = sol[:,3*N:4*N]%(y_num_cell*2.0*pl.pi)

        # puts the file number ahead of the RunImages folder to prevent overwriting
        f_num_str = f[:f.find('p')]
        if (f_num_str+'RunImages') not in os.listdir('.'):
            os.mkdir(f_num_str+'RunImages')
        if make_one_movie:
            os.mkdir(f_num_str+'RunImages/Space_'+projection)
        if make_tzpc:
            # see if the data is already in poincare section form
            already_pc = False
            if cycles >= len(sol)-10:
                already_pc = True
            # ultimitaly we want to work backwards from the final time so that if we want to we can cut off
            # transients but for now we are just going to look at the whole solution because
            # otherwise we are not going to have very many data points
            tz_fig = pl.figure()
            tz_ax = tz_fig.add_subplot(111)
            tz_ax.set_xlabel(x_lbl,fontsize=30)
            tz_ax.set_ylabel(y_lbl,fontsize=30)
            for i in range(len(sol)/throw_away,len(sol)):
                if(((i*dt)%(2.0*pl.pi))<dt) or already_pc:
                    #tz_ax.set_ylim([-2.0,2.0])
                    tz_ax.scatter(sol[i,ax_num[0]*N:ax_num[1]*N],sol[i,ax_num[2]*N:ax_num[3]*N],c='b',s=1)
            tz_fig.tight_layout()
            tz_fig.savefig(f_num_str+'RunImages/'+f_num_str+'tzpc.png')
            pl.close(tz_fig)
        if make_vzpc:
            # see if the data is already in t=0 poincare section form
            # if it is we cannot do this so exit
            if cycles >= len(sol)-10:
                print('Already t=0 sliced data -> can not make v=0 poincare section')
                raise SystemExit
            vz_fig = pl.figure()
            vz_ax = vz_fig.add_subplot(111)
            #vz_ax.set_xlim([0,d])
            #vz_ax.set_ylim([-2.0,2.0])
            vz_ax.set_xlabel(r'$x$'      ,fontsize=30)
            vz_ax.set_ylabel(r'$t$',fontsize=30)
            for i in range(len(sol)/throw_away,len(sol)-1):
                for j in range(N):
                    # if the velocoty changes sign then scatter plot it
                    if (sol[i,j]>=0.0 and sol[i+1,j]<=0.0) or (sol[i,j]<=0.0 and sol[i+1,j]>=0.0):
                        print('i is: '+str(i))
                        print('numbers are:')
                        print(sol[i,j])
                        print(sol[i+1,j])
                        vz_ax.scatter(sol[i,N+j]%d,(i*dt)%(2.0*pl.pi),c='b',s=1)
            vz_fig.tight_layout()
            vz_fig.savefig(f_num_str+'RunImages/'+f_num_str+'vzpc.png')
            pl.close(vz_fig)


        
        if make_one_movie:
            for i in range(len(sol)):
                print(i)
                print(sol[i,:N])
                # scatter for different As solutions
                #ax.scatter(sol[i,2*N:(2*N+N/2)],sol[i,3*N:(3*N+N/2)]%d,c='r')
                #ax.scatter(sol[i,(2*N+N/2):3*N],sol[i,(3*N+N/2):4*N]%d,c='b')
                if make_one_movie:
                    r_fig = pl.figure()
                    r_ax = r_fig.add_subplot(111)
                    r_ax.set_xlim(x_rng)
                    r_ax.set_ylim(y_rng)
                    r_ax.set_xlabel(x_lbl,fontsize=30)
                    r_ax.set_ylabel(x_lbl,fontsize=30)
                    r_ax.scatter(sol[i,ax_num[0]*N:ax_num[1]*N],sol[i,ax_num[2]*N:ax_num[3]*N],c='b',s=1)
                    r_fig.tight_layout()
                    r_fig.savefig(f_num_str+'RunImages/Space_'+projection+'/%(number)04d.png'%{'number':i})
                    pl.close(r_fig)

        if make_fv:
            # lets see what the velocity disrobution is when the potential is zero
            # go_back =int(pl.pi/2.0/dt)
            go_back = -1
            print('number of data points going into histogram: ' +str(len(sol[-1,:N])))
            fv_fig = pl.figure()
            fv_ax = fv_fig.add_subplot(111)
            pl.hist(sol[-go_back,:N],bins=want_bins)
            fv_ax.set_xlabel(r'$v$'      ,fontsize=30)
            #fv_ax.set_ylabel(r'$$',fontsize=30)
            #fv_ax.scatter(sol[i,N:2*N]%d,sol[i,:N],c='b')
            fv_fig.tight_layout()
            fv_fig.savefig(f_num_str+'RunImages/final_vel_his.png',dpi=300)
            pl.close(fv_fig)

        
    #pl.colorbar()
    #pl.xlabel("$A$",fontsize=30)
    #pl.ylabel("$x_1$",fontsize=30)
    #pl.savefig("paper_bif_hist.png",dpi=300)
    #os.system("open paper_bif_hist.png")

    #ax.set_xlabel("$A$",fontsize =  25 )
    #ax.set_ylabel("$x$", fontsize = 25 )
    #ax.axis([0.0,2*pl.pi,-1.8,1.8])
    # j[:-11] is the number infrot of poindat.txt
    #fig.savefig("bif_hist.png")


if __name__ == "__main__":
    main()
