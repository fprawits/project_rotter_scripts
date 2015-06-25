#!/usr/bin/env python2.7

import numpy as np
from useful_functions import get_parameter_value

import argparse
from argparse import ArgumentDefaultsHelpFormatter as Default_Formatter
parser=argparse.ArgumentParser(
        description='Change the height of a notch of the boundary and save '
                    'changes to a new file', formatter_class=Default_Formatter)

parser.add_argument('-f', '--infile',
                    default="rdmboundary.dat",
                    type=str,
                    help="file containing the random boundary to be modified")
parser.add_argument('-i', '--idx',
                    required=True,
                    type=int,
                    help='index of the notch to be modified; first notch bears'
                    ' index 0. If a negative index is given, counting is done '
                    'backwards, as for usual python iterables')
parser.add_argument('-m', '--modlist', 
                    metavar='INT',
                    nargs='+',
                    type=int, 
                    default=[-1,0,1],
                    help="for each integer a seperate modified boundary will "
                    "be generated, with the specified notch's height modified "
                    "by dy times that integer") 
parser.add_argument('-n',
                    type=int,
                    help='number of notches in cavity') 
parser.add_argument('-d', '--doc',
                    action='store_true',
                    default=False,
                    help='save essential input parametes to '
                    'modified_boundary.doc') 


def change_notch(
                infile="rdmboundary.dat",
                idx=None, 
                modlist=None, 
                n=None, 
                doc=False
                ):
    """Change the height of a notch of the boundary and save changes to new file

    infile: file containing the random boundary to be modified
    idx: index of the notch to be modified, the first notch beares index 0
    modlist: list of integers, specifying how often dy is be applied to notch
    n: number of notches
    """

    with open("parameters.xml") as params_file:
        print "Reading parameters..."
        params = params_file.read()
        pphw = get_parameter_value(params, "points_per_halfwave")
        modes = get_parameter_value(params, "modes")
        dy = 1.0/(modes*pphw+1)
        print """
        modes = {MODES}
        pphw = {PPHW}
        dy = {DY}
        """.format(MODES=modes, PPHW=pphw, DY=dy)

    if doc:
        doc_file = "modified_boundary.doc"
        print "Creating " + doc_file
        with open(doc_file, "w") as f:
            doc = """
            index of modified notch:    {INDEX}
            number of notches:          {N}
            modlist:                    {MODLIST}
            dy:                         {DY}
            """.format(INDEX=idx, N=n, MODLIST=modlist, DY=dy)
            f.write(doc)


    print "Reading boundary from {}...".format(infile)
    x, boundary = np.loadtxt(infile, unpack=True)
    # Find the start and end of the cavity, because boundary contains leads
    cavity_idx = boundary.nonzero()[0]
    start = cavity_idx[0]
    end = cavity_idx[-1] + 1

    cavity = boundary[start:end].reshape( (n,-1) )

    # Make sure modlist consists of integers, and change it such that it can be
    # used conveniently for iteration
    modlist = np.array(modlist, dtype=int)
    modlist[1:] = np.diff(modlist)
    
    for i, factor in enumerate(modlist):
        print "modifying height of notch {}...".format(idx)
        cavity[idx] += factor*dy
        ofile = "rdmboundary_{}.dat".format(i)
        print "saving modified boundary to {}".format(ofile)
        np.savetxt(ofile, np.transpose( (x,boundary) ), fmt="%g \t %g")
    

if __name__ == "__main__":
    args = vars(parser.parse_args())
    change_notch(**args)
