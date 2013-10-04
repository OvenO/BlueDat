import pylab as pl
import os
import argparse
import numpy as np


def main():
    # how many time slice (ts) sections do we want to include in the images?
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', action = 'store', dest = 'n', type = int, default = 1)
    parser.add_argument('-p', action = 'store', dest = 'p',type = int, default = 1800)
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
        print(coef)
        
        data = np.genfromtxt(cur_file,comments="Z")
        x = pl.append(x,data[-b_pts*num_ts:,2])
        A_arr = pl.append(A_arr,build*coef)
        

        #ax.scatter(pl.zeros([len(x)])+coef,x,marker='s',s=5,alpha = opacity,edgecolors = "None")
        #ax.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Grey", marker = "o", label = "Electrodes")

    
    H,xedges,yedges = pl.histogram2d(x,A_arr,bins=(20,20))
    print(xedges)
    extent = [yedges[0],yedges[-1],0.0,2*pl.pi]
    print(yedges[0]) 
    print(yedges[-1]) 
    

    pl.imshow(H,extent = extent,aspect="auto", interpolation = "nearest",cmap=pl.cm.Greys)
    pl.colorbar()
    pl.xlabel("$A$",fontsize=30)
    pl.ylabel("$x$",fontsize=30)
    pl.savefig("bif_hist.png")
    
    os.system("open bif_hist.png")
    #ax.set_xlabel("$A$",fontsize =  25 )
    #ax.set_ylabel("$x$", fontsize = 25 )
    #ax.axis([0.0,2*pl.pi,-1.8,1.8])
    # j[:-11] is the number infrot of poindat.txt
    #fig.savefig("bif_hist.png")


if __name__ == "__main__":
    main()
