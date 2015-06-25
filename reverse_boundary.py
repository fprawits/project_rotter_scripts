#!/usr/bin/env python2.7

#infile = "rdmboundary.dat"
#outfile = "rev_" + infile

def reverse_boundary(infile="rdmboundary.dat"):
    """
    Reverse the rdm boundary for injection from right and update parameters.xml
    
    infile:     file to read the random boundary from; reversed boundary will
                be written to outfile = "rev_" + infile
    """

    import numpy as np
    outfile = "rev_" + infile

    with open(infile) as f:
        print "found " + infile 

    print "reversing boundary..."
    x,y = np.loadtxt(infile, unpack=True)
    y_rev = y[::-1]
    
    np.savetxt(outfile, np.transpose([x,y_rev]), fmt="%.16f")
    print "reversed boundary saved to " + outfile

    with open("parameters.xml", "r") as params_file:
        params = "".join(params_file.readlines())
        print "updating parameters.xml..."
        params = params.replace(infile, outfile)

    with open("parameters.xml", "w") as params_file:
        params_file.write(params)
        print "parameters.xml updated"

    return

if __name__ == "__main__":
    reverse_boundary()
