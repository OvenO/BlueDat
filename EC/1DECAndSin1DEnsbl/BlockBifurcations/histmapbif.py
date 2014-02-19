import pylab as pl
import os
import argparse
import numpy as np


def main():
    # paper or presentation format
    paper = True

    # how many time slice (ts) sections do we want to include in the images?
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', action = 'store', dest = 'n', type = int, default = 1)
    parser.add_argument('-p', action = 'store', dest = 'p',type = int, default = 1)
    inargs = parser.parse_args()
    
    num_ts = inargs.n
    b_pts = inargs.p

    # opacity of points 0-1
        
    all_dir = os.listdir(".")
    list_dir = [] 

    for i,j in enumerate(all_dir):
        if "poindat" in j:
            list_dir.append(j)
    
    #fig = pl.figure()
    #ax = fig.add_subplot(111)
    x = pl.array([])
    y = pl.array([])
    vx = pl.array([])
    vy = pl.array([])
    
    # this array keeps track of the coefficient values for each particle
    A_arr = pl.array([])
    build = pl.zeros([num_ts*b_pts])
    build+=1.0

    for i,j in enumerate(list_dir):
        location = -1

        cur_file = open(j,"r")

        # get the coeficient for the plot
        coef = float(cur_file.readline().split()[-1])
        print(i)
        
        data = np.genfromtxt(cur_file,comments="Z")
        x = pl.append(x,data[-b_pts*num_ts:,2]%(2.0*pl.pi))
        A_arr = pl.append(A_arr,build*coef)
        
    if (paper):
        pl.hexbin(A_arr,x,gridsize=300,bins="log",mincnt=1,edgecolor="none") 
        #cmap=pl.cm.YlOrRd_r
    else:
        pl.hexbin(A_arr,x,gridsize=300,cmap=pl.cm.Greys,bins="log",mincnt=1,edgecolor="none") 
    #pl.hexbin(A_arr,x,bins="log") 


    pl.colorbar()
    pl.xlabel("$A$",fontsize=30)
    pl.ylabel("$x_1$",fontsize=30)
    if paper:
        pl.savefig("paper_bif_hist.png",dpi=600)
    else:
        pl.savefig("pres_bif_hist.png",dpi=600,transparent = True)
    
    if paper:
        os.system("open paper_bif_hist.png")
    else:
        os.system("open pres_bif_hist.png")

    #ax.set_xlabel("$A$",fontsize =  25 )
    #ax.set_ylabel("$x$", fontsize = 25 )
    #ax.axis([0.0,2*pl.pi,-1.8,1.8])
    # j[:-11] is the number infrot of poindat.txt
    #fig.savefig("bif_hist.png")


if __name__ == "__main__":
    main()
