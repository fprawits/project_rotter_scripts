#!/usr/bin/env python2.7

import numpy as np

def step(x,a,d):
    """Return 1 if x lies in the intervall [a,a+d) and 0 else"""

    return np.logical_and(x >= a, x < a+d).astype(int)


def generate_boundary2(n, lenth, modes, pphw, outfile = "rdmboundary.dat", seed = None,
                       doc = True):
    """Generate a boundary with random notch heights, consisting of n notches

    n:          integer number of notches
    modes:      number of open modes, therefore defining scattering energy
                usually this is a float value
    pphw:       points per halfwave (accuracy parameter)
    outfile:    filename of stored boundary (default: 'rdmboundary.dat')
    seed:       seed for the random number generator (default: None)
    """

    length = 0.2
    offset = 1.0
    L = n*length
    L_prime = L + 2*offset
    W = 1.0
    dx = W/(modes*pphw + 1)
    r_nx = int(L_prime/dx)
    x = np.linspace(0, L_prime, r_nx)

    np.random.seed(seed)
    delta = 0.3
    height = (np.random.random(n)*2 - 1) * delta/2
    y = np.array([step(x, i*length + offset, length) * height[i] for i in range(n)])
    y = y.sum(axis=0)
    
    np.savetxt(outfile, zip(x,y))

    if doc:
        doc_file = open("boundary.doc", "w")
        doc_file.write(
                """#boundary generated with generate_boundary2.py
                input parameters:
                ----------------
                notches:    {}
                length:     {}
                modes:      {}
                pphw:       {}
                seed:       {}


                calculation parameters:
                ----------------------
                delta:      {}
                length:     {}
                L:          {}
                offset:     {}
                dx:         {}
                N_file:     {}
                L':         {}
                """.format(n, modes, pphw, seed, delta, L, L_prime, dx, N_file))
        doc_file.close()


if __name__ == '__main__':

    generate_boundary2(n=5, modes=10.5, pphw=20, seed=1)
