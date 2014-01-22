from visual import *
import argparse
import time
import os

def main():
    
    #parser = argparse.ArgumentParser(description="Tell the program to collect video images?")
    #parser.add_argument('v', dest = "video", action="store",type=string,

    video = False

    file = open("realvisual.txt","r")

    curline = file.readline()
    curline = file.readline()
    
    curnum = curline.split()

    ball = sphere(pos=(float(curnum[0]),float(curnum[1]),0),color = color.red, radius =.1)

    floor = box(pos=(0,1.0,0),length = 6.5, height = 0.5, width = 3.0, color=color.blue)

    electrode1 = cylinder(pos=(-3.14159,0,-1.5),axis=(0,0,3),radius = .1)
    electrode2 = cylinder(pos=(-0,0,-1.5),axis=(0,0,3),radius = .1)
    electrode3 = cylinder(pos=(3.14159,0,-1.5),axis=(0,0,3),radius = .1)
    
    if(video):
        time.sleep(10) 

    count = 1
    while(curline != ""):
        
        rate(100)
            
        if (video):
            os.system("screencapture "+str(count)+".png")
            count += 1


        curnum = curline.split()
        ball.pos = [float(curnum[0])-3.14159,float(curnum[1]),0.0]

        curline = file.readline()



if __name__ == "__main__":
    main()
