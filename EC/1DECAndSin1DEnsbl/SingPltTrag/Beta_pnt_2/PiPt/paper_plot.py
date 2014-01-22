import os
import pylab as pl
from scipy.integrate import odeint
from scipy.integrate import ode
import numpy

def main():

    is_transparent = False 
    
    f = open("pi_data.txt","r")
    
    # this is a little different than normal becase of the complex data for the floquet stability
    # multipliers. When we use the "dtype" option we get a single array of tuples so slicing is a
    # little more awkward has to look like data[#][#] to get a single value NOT data[#,#].
    data = pl.genfromtxt(f,comments="e",dtype="complex,complex,float")
   
    eigs1 = pl.array([])
    eigs2 = pl.array([])
    A = pl.array([])
    
    for i,j in enumerate(data):
        eigs1 = pl.append(eigs1,j[0])
        eigs2 = pl.append(eigs2,j[1])
        A = pl.append(A,j[2])

    fig1, ax1 = pl.subplots(2,2,sharex=True)
    ax1[0,0].plot(A,[k.real for k in eigs1],color = "Black")
    ax1[1,0].plot(A,[k.imag for k in eigs1],color = "Black")
    ax1[0,1].plot(A,[k.real for k in eigs2],color = "Black")
    ax1[1,1].plot(A,[k.imag for k in eigs2],color = "Black")

    ax1[0,0].set_ylabel("Re[$\lambda_1$]",fontsize=25)
    ax1[1,0].set_ylabel("Im[$\lambda_1$]",fontsize=25)
    ax1[0,1].set_ylabel("Re[$\lambda_2$]",fontsize=25)
    ax1[1,1].set_ylabel("Im[$\lambda_2$]",fontsize=25)
    ax1[1,0].set_xlabel("$A$",fontsize=25)
    ax1[1,1].set_xlabel("$A$",fontsize=25)
    fig1.tight_layout()
    fig1.savefig("paper_A_vs_eigs.png",dpi=300,transparent=is_transparent)
    os.system("open paper_A_vs_eigs.png")

   
   
    #fig3, ax3 = pl.subplots(2,sharex=True)
    #ax3[0].plot(A_arr,[k.real for k in eigs1],color = "Black")
    #ax3[1].plot(A_arr,[k.imag for k in eigs1],color = "Black")
    #ax3[0].set_ylabel("Re[$\lambda_1$]",fontsize=25)
    #ax3[1].set_ylabel("Im[$\lambda_1$]",fontsize=25)
    #ax3[1].set_xlabel("$A$",fontsize=25)
    #fig3.tight_layout()
    #fig3.savefig("pi_A_vs_eig1.png",dpi=300,transparent=is_transparent)
    #os.system("open pi_A_vs_eig1.png")

    #fig4, ax4 = pl.subplots(2,sharex=True)
    #ax4[0].plot(A_arr,[k.real for k in eigs2],color = "Black")
    #ax4[1].plot(A_arr,[k.imag for k in eigs2],color = "Black")
    #ax4[0].set_ylabel("Re[$\lambda_2$]",fontsize=25)
    #ax4[1].set_ylabel("Im[$\lambda_2$]",fontsize=25)
    #ax4[1].set_xlabel("$A$",fontsize=25)
    #fig4.tight_layout()
    #fig4.savefig("pi_A_vs_eig2.png",dpi=300,transparent=is_transparent)
    #os.system("open pi_A_vs_eig2.png")

    #fig5, ax5 = pl.subplots(2,sharex=True)
    #ax5[0].plot(A_arr,abs(eigs1),color = "Black")
    #ax5[1].plot(A_arr,abs(eigs2),color = "Black")
    #ax5[0].set_ylabel("$\lambda_1$",fontsize=25)
    #ax5[1].set_ylabel("$\lambda_2$",fontsize=25)
    #ax5[1].set_xlabel("$A$",fontsize=25)
    #fig5.tight_layout()
    #fig5.savefig("pi_A_vs_mag_eigs.png",dpi=300,transparent=is_transparent)
    #os.system("open pi_A_vs_mag_eigs.png")


    #fig1 = pl.figure()
    #ax1 = fig1.add_subplot(111)
    #ax1.plot(pl.cos(theta),pl.sin(theta),color = "Black")
    #ax1.plot([k.real for k in eigs1],[l.imag for l in eigs1],color = "Black")
    #ax1.set_xlabel("Re[$\lambda_1$]",fontsize=25)
    #ax1.set_ylabel("Im[$\lambda_1$]",fontsize=25)
    #fig1.tight_layout()
    #fig1.savefig("pi_eig1.png",dpi=300,transparent=is_transparent)
    #os.system("open pi_eig1.png")

    #fig2 = pl.figure()
    #ax2 = fig2.add_subplot(111)
    #ax2.plot(pl.cos(theta),pl.sin(theta),color = "Black")
    #ax2.plot([k.real for k in eigs2],[l.imag for l in eigs2],color = "Black")
    #ax2.set_xlabel("Re[$\lambda_2$]",fontsize=25)
    #ax2.set_ylabel("Im[$\lambda_2$]",fontsize=25)
    #fig2.tight_layout()
    #fig2.savefig("pi_eig2.png",dpi=300,transparent=is_transparent)
    #os.system("open pi_eig2.png")
 

if __name__ == '__main__':
    main()
