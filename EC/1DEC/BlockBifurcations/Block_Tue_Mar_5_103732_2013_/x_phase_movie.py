import pylab as pl
import os
import argparse


def main():
    # how many time slice (ts) sections do we want to include in the images?
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', action = 'store', dest = 'n', type = int, default = 1)
    inargs = parser.parse_args()
    
    num_ts = inargs.n

        
    os.mkdir("Images")
    
    all_dir = os.listdir(".")
    list_dir = [] 

    for i,j in enumerate(all_dir):
        if "poindat" in j:
            list_dir.append(j)
    
    print(list_dir)

    for i,j in enumerate(list_dir):
        location = -1
        fig = pl.figure()
        ax = fig.add_subplot(111)

        cur_file = open(j,"r")
        lines = cur_file.readlines()
        
        cur_line = lines[location]
        parsed_line = cur_line.split()

        x = pl.array([])
        y = pl.array([])
        vx= pl.array([])
        vy= pl.array([])
        for a in range(num_ts):
            while (parsed_line[0] != "ZONE"):

                x = pl.append(x,float(parsed_line[2]))
                y = pl.append(y,float(parsed_line[3]))
                vx = pl.append(vx,float(parsed_line[0]))
                vy = pl.append(vy,float(parsed_line[1]))
                
                location -=1
                cur_line = lines[location]
                parsed_line = cur_line.split()

            location -=1                   
            cur_line = lines[location]
            parsed_line = cur_line.split()

        ax.scatter(x,vx,marker='o',s=2,label="Particle")
        #ax.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Grey", marker = "o", label = "Electrodes")
        ax.set_xlabel("$x$",fontsize = 30)
        ax.set_ylabel("$\dot{x}$", fontsize = 30)
        #ax.axis([0.0,2*pl.pi,-1.8,1.8])
        # j[:-11] is the number infrot of poindat.txt
        fig.savefig(j[:-11]+".png")
    
    os.system("mv *.png Images")

if __name__ == "__main__":
    main()
