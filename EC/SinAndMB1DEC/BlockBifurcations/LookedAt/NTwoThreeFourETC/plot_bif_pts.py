import pylab as pl

def main():
    x_N3 = pl.array([1.23,1.235])
    x_N4 = pl.array([0.98,0.985])
    x_N5 = pl.array([1.836,1.84275])
    x_N6 = pl.array([1.520,1.535])
    x_N7 = pl.array([2.555625,2.564375])
    x_N8 = pl.array([4.186666667,4.2])
    x_N9 = pl.array([3.68125,3.6875])

    #split into even and odd, [N=1,N=3,...]
    odd_bif_pts = pl.array([.75365,x_N3.mean(),x_N5.mean(),x_N7.mean(),x_N9.mean()])
    even_bif_pts= pl.array([.615,x_N4.mean(),x_N6.mean(),x_N8.mean()])

    # error arrays
    odd_err = pl.array([5e-5,x_N3.std(),x_N5.std(),x_N7.std(),x_N9.std()])
    even_err= pl.array([.005,x_N4.std(),x_N6.std(),x_N8.std()])

    # N arrays
    N_o = [1,3,5,7,9]
    N_e = [2,4,6,8]

    fig = pl.figure()
    ax = fig.add_subplot(111)

    #ax.plot([1,2,3,4],a_arr)
    ax.errorbar(N_o,odd_bif_pts,yerr=odd_err,linestyle = "-",marker="^",color="Blue")
    ax.errorbar(N_e,even_bif_pts,yerr=even_err,linestyle = "-",marker="s",color="Red")
    ax.set_xlabel(r'$N$',fontsize=25)
    ax.set_ylabel(r'$A_{c1}$',fontsize=25)
    fig.savefig('multi_bif_vals.eps')


if __name__ == "__main__":
    main()
