import pylab as pl
import os
import numpy as np
def main():
    f = open("final_position.txt","r")
    data = pl.genfromtxt(f,comments = "L")
    
    # need to get every other
    x = pl.array([])
    y = pl.array([])
    for i,j in enumerate(data[:-7,2]):
        if i%2 == 0:
            x = pl.append(x,data[i,4])
            y = pl.append(y,j)
    
    print(x)
    print(y)
    fit = np.polyfit(x,y,2)

    print(fit)
    
    #fited = fit[0]+fit[1]*x + fit[2]*x**2
    fited = np.poly1d(fit)
    print(fited)

    pl.plot(pl.append(x,[.262,.264,.266]),fited(pl.append(x,[.262,.264,.266])),color="black")
    pl.scatter(x,y,color = "black")
    pl.xlabel("$A$",fontsize="30")
    pl.ylabel("$x$",fontsize="30")
    pl.savefig("fin_pts.png",transparent=True,dpi=300)
    
    os.system("open fin_pts.png")

if __name__ == "__main__":
    main()
