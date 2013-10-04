This directory is going to contain the highest fuction of the analysis programs in "datashpere". The
purpose of the programs in directory "Head" is to make the implementation of the analysis sofware
easy. In short the program head.py will take the input variables and feed them to the VACC (or
whichever computer cluster you would like).

Keep in mind that at first head.py will not be very sofisticated and therefore error checking will
need to be done by hand. As time goes on I will make the program capable of trobleshooting errors by
reading the .o files and hadeling the probplem apropriately. 

Also you may want to run head.py on a different machine (like the linux one) so that if your main
computer needs to be restarted or put to sleep nothing will be inturupted.

head.py will need to read a file that contains which system to work on and the governing parameters
-> make text file "to_run_<key_word>.txt" that contains the parameters and runtime for the
<key_word> system. Fill in <key_word> with the system of interest as the to_run files will be very
different from eachother. The posible key_words right now(kida) are:
    a) 1DEC -> one-dimesional electric curtain
    b) 2DEC -> two-dimesional electric curtain
    c) 1DPP -> one-dimesional periodic potential
    c) 2DPP -> two-dimesional periodic potential
This file will be structured like this:
<number> discription of the number

Here is the phsudo code for head.py:

1) input argument (argparse) tells us wich system and therefor which to_run file to read

2) open appropriate "to_run_<key_word>.txt"

3) parse parameters and asign parameters to varables

4) initialize a block(s) (localy) -> store in this directory
    a) this will reqire putting the initabloc.py code in head 
    b) include a Orig file to keep initial conditions incase something goes wrong

5) move the block directory to the right "Old" directory on cluster

6) remotly tell cluster to run the corect again.py. Need to put the block name in again.script
    a) for this restructure again.py so that it can take the runtime as an argparse argument. NEED
    TO PASS THIS ARGUMENT FROM again.script !!!

7) wait for all again.py to finish
    a) this needs to be done by checking for compleation files in regular intervals of time.
    b) NOTE: I think we should make this check the .o fiels to see if there was an "exceeded
    walltime limit" for any of them
8) after everything is finised run takelast.py. program will need to modify or pass the block name to
takelast.script

9) wait for takelast to finish.

10) make directory in the right FINAL place for block info

12) scp Last stuff to the above directory

13) Notify Owen when done
