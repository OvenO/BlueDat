import pylab as pl
import os
import argparse
import numpy as np


def main():
    # paper or presentation format
    paper = True
    # define variable to get rid of first time slices
    # ( is in actual number of coef or beta that we want to start from)
    cutoff = .2

    # how many time slice (ts) sections do we want to include in the images?
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', action = 'store', dest = 'n', type = int, default = 1)
    parser.add_argument('-p', action = 'store', dest = 'p',type = int, default = 1830)
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
        
        if coef>cutoff:
            data = np.genfromtxt(cur_file,comments="Z")
            x = pl.append(x,data[-b_pts*num_ts:,2])
            y = pl.append(y,data[-b_pts*num_ts:,3])
            A_arr = pl.append(A_arr,build*coef)
        
    y = 1.0+abs(1.0-y)
    x = x%(2.0*pl.pi)
    A_arr = A_arr[100:]
    y= y[100:]
    if (paper):
        pl.hexbin(A_arr,y,gridsize=200,bins="log",mincnt=2,edgecolor='none') 
        #H, xe, ye = pl.histogram2d(A_arr,y,bins=190) 
        #cmap=pl.cm.YlOrRd_r
    #else:
    #    pl.hexbin(A_arr,y,gridsize=400,cmap=pl.cm.Greys,bins="log",mincnt=2) 
    #pl.hexbin(A_arr,x,bins="log") 
    
    #H = pl.log(H)
    #pl.imshow(H,interpolation = "nearest")
    pl.colorbar()
    pl.xlabel("$beta$",fontsize=30)
    pl.ylabel("$x$",fontsize=30)
    #pl.ylim([1.0,30])
    if paper:
        pl.savefig("y_paper_bif_hist.png",dpi=600)
    else:
        pl.savefig("y_pres_bif_hist.png",dpi=600,transparent = True)
    
    if paper:
        os.system("open y_paper_bif_hist.png")
    else:
        os.system("open y_pres_bif_hist.png")

    #ax.set_xlabel("$A$",fontsize =  25 )
    #ax.set_ylabel("$x$", fontsize = 25 )
    #ax.axis([0.0,200,1.0,80])
    # j[:-11] is the number infrot of poindat.txt
    #fig.savefig("bif_hist.png")


if __name__ == "__main__":
    main()
