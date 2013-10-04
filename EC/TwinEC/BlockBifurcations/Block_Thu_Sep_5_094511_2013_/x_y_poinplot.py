import pylab as pl
import os


def main():
    which = 136
    thedata = open(str(which)+"poindat.txt","r")
    os.mkdir(str(which)+"PoinImages")
    os.chdir(str(which)+"PoinImages")

    # The first 3 lines are just stuff for the tecplot program so lets just get rid of those.
    for i in range(2):
        thedata.readline()
 
    # need to parse and colect into appropriate arrays. use split fuction built into string object.
    curline = thedata.readline()
    parsedline = curline.split()

    i = 0
   
    while(curline!= ""):
        # need figure object so we dont end up with overlaping images
        fig = pl.figure()
        ax = fig.add_subplot(111)
        # Now we get to the data part.
        # first lets define our arrays so we can collect the data appropriately.
        # We could make one matrix but lets use individual arrays for the sake of readability.
        x = pl.array([])
        y = pl.array([])
        vx = pl.array([])
        vy = pl.array([])
        while(parsedline[0]!="ZONE"):
            x = pl.append(x,float(parsedline[2]))%(2.0*pl.pi)
            y = pl.append(y,float(parsedline[3]))
            vx = pl.append(vx,float(parsedline[0]))
            vy = pl.append(vy,float(parsedline[1]))

            curline = thedata.readline()
            parsedline = curline.split()


        ax.scatter(x,y,marker='o',s=.1,label="Particle")
        ax.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", marker = "o", label = "Electrodes")
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        #ax.axis([0.0,2*pl.pi,-6.0,6.0])
        fig.savefig(str(i)+".png")
    
        i += 1
        curline = thedata.readline()
        parsedline = curline.split()

if __name__ == "__main__":
    main()
