#!/usr/bin/python
import os
import scipy as sp
import argparse
import time as thetime
import shutil

def move_files():
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

def remove_zeros():
    for i,j in enumerate(os.listdir('.')):
        if 'poindat.txt' in j:
            # difine this so we only work on the file if we need to
            work_on_file = True
            print('working file: '+str(j))
            f = open(j,'r')
            # get rid of first line (string)
            var_line = f.readline()
            # get all the data
            data = sp.genfromtxt(f)
            f.close()

            print('shape of data after genfromtxt: ' + str(sp.shape(data)))

            if len(sp.shape(data))==1:
                print('This file did not run: ' + j)
                print('moving to failed_run_'+j[:j.find('p')]+'.txt')
                shutil.move(j,'failed_run_'+j[:j.find('p')]+'.txt')
                continue 
                

            # if the sum of the squares is zero the line is all zerros
            count = 1
            while count<len(data):
                line_sqrd = data[-count,:]**2
                if line_sqrd.sum() != 0.0:
                    # if the line_sqrn isnt zero and count =1 there are no zero lines so move on.
                    # i.e dont work_on_file.
                    if count == 1:
                        work_on_file = False

                        
                    break
                count+=1

            # gets rid of actual data line unless you do this
            count -=1

            if work_on_file:
                # write file 
                print('writing to file: ' + str(j))
                wf = open(j,'w')
                wf.write(var_line)

                # write lines ommiting zerro lines
                print('shape of data array: '+str(sp.shape(data)))
                for a in range(len(data[:-count,0])):
                    # 1:-1 string slice gets rid of the [ ] brackets in the string
                    wf.write(str(data[a,:])[1:-1].replace('\n',''))
                    wf.write('\n')
                #print('last line after corections: ' +str(data[:-(count+1),:]))
                
                wf.close()

def main():

    parser = argparse.ArgumentParser()
    # d is for directory
    parser.add_argument('-d',action='store',dest = 'd',type = str, required = True)
    # what are we doing? t = type
    parser.add_argument('-t',action='store',dest = 't',type = str, required = True)

    inargs = parser.parse_args()
    d = inargs.d
    do_what = inargs.t
 
    os.chdir(d)

    if do_what == 'mv':
        move_files()
    if do_what == 'zeros':
        remove_zeros()

    # a couple of audible bells to let me know its done
    for i in range(5):
        thetime.sleep(.5)
        print('\a')

if __name__ == '__main__':
    main()
