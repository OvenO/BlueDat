import os
import scipy as sp
import argparse
def main():

    parser = argparse.ArgumentParser()
    # d is for directory
    parser.add_argument('-d',action='store',dest = 'd',type = str, required = True)
    inargs = parser.parse_args()
    d = inargs.d
 
    os.chdir(d)

    if 'DidntFinish' not in os.listdir('.'):
        os.mkdir('DidntFinish')

    for i,j in enumerate(os.listdir('.')):
        if 'poindat.txt' in j:
            f = open(j,'r')
            first = f.readline()
            second = f.readline()
            third = f.readline()
            fourth = f.readline()
            f.close()
            if fourth == '':
                os.system('mv '+j+' DidntFinish/')
                print('moving '+j)

if __name__ == '__main__':
    main()
