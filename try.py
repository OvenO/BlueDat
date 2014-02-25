import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a',action='store',dest='a',type=bool,required=False)
    inargs = parser.parse_args()
    to_print = inargs.a 

    print(to_print)
    print('type(to_print): '+ str( type(to_print)))
    print('to_print == None: ' + str(to_print==None))


if __name__ == '__main__':
    main()
