#!/usr/bin/env python2.7
# TO DO:    better implementation of plot option
#           maybe implement lead joining with numpy vectors
##############################################################################
########################        rdmstep.py      ##############################
##############################################################################

"""Generate the one sided surface disorder"""

import sys, random, argparse
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentDefaultsHelpFormatter as Default_Formatter

parser=argparse.ArgumentParser(description='Generate the one sided surface ' 
                               'disorder', formatter_class=Default_Formatter)
parser.add_argument('delta', type=float, \
                    help='range of random notch height in boundary')
parser.add_argument('length', type=float, help='length of individual notch') 
parser.add_argument('notches', type=int, help='number of notches in boundary ' 
                    '(Note: a distance of jnotches*length without surface '
                    'disorder will be added left and right in order to join '
                    'the lead properly!)')
parser.add_argument('ppn', type=int, help='number of points per notch; '\
                    'the actual points per notch are chosen s. t. ppn/length '\
                    'is divisible without remainder')
parser.add_argument('-p', '--plot', action='store_true', default=False,
                    help='plot the potential (implemented unefficiently)')
parser.add_argument('-o', '--outfile', type=str, metavar='filename.ending',
                    default='rdmboundary.dat', help='name of output file')
parser.add_argument('-s', '--seed', nargs='?', type=int, default=1, const=None,
                    help='seed for the random number generator, '
                    'invoking this option without a value will result in '
                    'pseudo-random behaviour with no fixed seed')
parser.add_argument('-d', '--doc', action="store_true", default=False,
                    help="Create a boundary.doc file containing all arguments") 
 
                       
def generate_boundary(delta=None, length=None, notches=None, ppn=None,
                      plot=False, outfile="rdmboundary.dat", seed=None,
                      doc = True, **kwargs):
    """
    Generate the one sided surface disorder
    
    delta: range of random notch height in boundary 
    length: length of individual notch
    notches: number of notches in boundary
    ppn: number of points per notch the actual points per notch are chosen 
         s. t. ppn/length is divisible without remainder

    options:
    plot: plot data (implemented unefficiently)
    outfile: name of output file
    seed: seed for the random number generator (default: None)
    doc: Create a boundary.doc file containing all arguments (default: True)
    kwargs: dump variable for convenient usage with iterate.py
    """

# length of the cavity
    L = notches*length
# resetting ppn such that ppl is definitely an integer
    ppn = int(int(ppn/length) * length)
# resolution in x-direction
    dx = length/ppn
# length of the pre- and post-leads, used to join cavity properly
# if value is changed from 1.0, length_lead/dx might not be integer any longer
    length_lead = 1.0
# points per lead
    ppl = int(round(length_lead/dx))
# total number of points in x direction
    N_file = notches*ppn + 2*ppl + 1
# total length of the system including pre- and post-leads:
    L_prime = L + 2*length_lead

# for the purpose of debugging, function should be called with fixed seed
    random.seed(seed)
    profile = open(outfile, 'w')

# add no surface disorder of length length_lead (ppl*dx) to the left,
# in order to join lead properly
    x = 0.0
    y = 0.0
    for i in range(ppl):
        profile.write('%(x)g \t %(y)g\n' %vars()) 
        x += dx


    for i in range(notches):
        y = random.uniform(-delta/2, delta/2)       # get new step height
        for j in range(ppn):
            profile.write('%g \t %g\n' %(x, y))
            x += dx


# add no surface disorder of length length_lead (ppl*dx) to the right,
# in order to join lead properly
    y = 0.0
    for i in range(ppl+1):
        profile.write('%(x)g \t %(y)g\n' %vars()) 
        x += dx

    profile.close()

# options:
    if plot == True:
        x,y = np.loadtxt(outfile, unpack = True)
        plt.plot(x,y)
        plt.show()

    if doc:
        doc_file = open("boundary.doc", "w")
        doc_file.write(
                """
                input parameters:
                ----------------
                delta:      {}
                length:     {}
                notches:    {}
                ppn:        {}

                options:
                -------
                seed:       {}
                outfile:    {}

                calculation parameters:
                ----------------------
                L:          {}
                L':         {}
                dx:         {}
                N_file:     {}
                """.format(delta, length, notches, ppn, seed, outfile,\
                           L, L_prime, dx, N_file))
        doc_file.write("\n")
        doc_file.close()
    return

###############################################################################

if __name__ == '__main__':

    args=vars(parser.parse_args())
    generate_boundary(**args)

