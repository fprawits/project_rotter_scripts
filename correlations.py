#!/usr/bin/env python2.7

import numpy as np
import os
import glob
import re
import Q_Matrix

from rename_images import rename_images

from useful_functions import update
from useful_functions import run_code
from useful_functions import get_parameter_value
from useful_functions import convert_to_complex
from useful_functions import backup_files
from useful_functions import remove_files

import argparse
from argparse import ArgumentDefaultsHelpFormatter as help_formatter

parser=argparse.ArgumentParser(
    description="""Calculate the correlations <psi_B|psi_B+delta_B> of Q_B and
    transmission eigenstates to a reference B""",
    formatter_class=help_formatter) 

parser.add_argument('-B', '--ref-B',
                    type=float,
                    required=True,
                    help="""reference value for magentic field used to
                    calculate correlations""")
parser.add_argument('-d', '--delta-B',
                    type=float,
                    default=0.5,
                    help="""stepsize for magnetic field scan""")
parser.add_argument('-n',
                    type=int,
                    default=5,
                    help="""number of steps in one direction, the total number
                    of points is given by 2*n + 1""")
parser.add_argument('-p', '--multi-core',
                    metavar='#cores',
                    type=int,
                    nargs='?',
                    default=1,
                    const=4,
                    help="""when running the code locally, this option enables
                    parallel computing. If no argument is provided, all 4 cores
                    will be used!""")
parser.add_argument('-b', '--backup',
                    action="store_true",
                    default=False,
                    help="""Prevent overriding of data, by backing up all files
                    to directory bak/ with the current B-value prepended to 
                    their filename""")

args=vars(parser.parse_args())


# Get parameters
with open("parameters.xml") as params_file:
    print "Reading parameters..."
    params = params_file.read()
    pphw = get_parameter_value(params, "points_per_halfwave")
    modes = get_parameter_value(params, "modes")
    dx = 1.0/(modes*pphw+1)
    print """
    modes = {MODES}
    pphw = {PPHW}
    dx = {DX}
    """.format(MODES=modes, PPHW=pphw, DX=dx)


# Calculate wavefunction for reference B
# --------------------------------------
print "Starting calculation for reference B = {}".format(args["ref_B"])

with open("parameters.xml", "w") as params_file:
    print "Updating parameters..."
    params = update(params, Bfield=args["ref_B"])
    params_file.write(params)
    dB = get_parameter_value(params, "dB") * args["ref_B"]
    print """
    Bfield = {BFIELD}
    dB = {DERIVATIVE_STEPSIZE}
    """.format(BFIELD=args["ref_B"], DERIVATIVE_STEPSIZE=dB)

run_code("input.xml", opt=args["multi_core"])
Q_B = Q_Matrix.Time_Delay_Matrix(derivative_stepsize=dB)
Q_B.write_eigenstates()
Q_B.write_eigenvalues()
run_code("output.xml", opt=args["multi_core"])

cwd = os.getcwd()
ascii_files = glob.glob(cwd+"/Q_states.*.coeff.*.purewavefunction.ascii")
ascii_files.sort()
ascii_regex = r'Q_states\.(\d{,4})\.coeff\.(\d{,4})\.purewavefunction\.ascii'
ascii_pattern = re.compile(ascii_regex)

# Get file length to allocate psi_ref
with open(ascii_files[0]) as f:
    ascii_f_len = sum(1 for line in f)
    psi_ref = np.zeros( (int(modes),ascii_f_len), dtype=complex ) 

for f in ascii_files:
    print "Reading {FILE}...".format(FILE=f)
    mode = int(ascii_pattern.search(f).group(1))
    psi_ref[mode] = np.loadtxt(f, usecols=(1,), dtype=complex,
                    converters={1: convert_to_complex})

if args["backup"]:
    print "Renaming images..."
    rename_images()
    print "Backing up files..."
    backup_files(glob_arg="Q_states", prefix="B_ref_{}_".format(args["ref_B"]))


# Scan over all B values
# ----------------------
# the total number of B_values is 2*n + 1
B_values = np.array([ args["ref_B"] + i*args["delta_B"] 
                    for i in range(-args["n"], args["n"]+1) ])
correlation_table = np.zeros( (B_values.size, modes) )
with open("max_mode.log", "w") as log:
    log.write("# Bfield file_mode ref_mode\n")

#for i,B in enumerate(B_values):
for B, corr in zip(B_values, correlation_table):
    print "Starting calculation for B = {}".format(B)

    with open("parameters.xml", "w") as params_file:
        print "Updating parameters..."
        params = update(params, Bfield=B)
        params_file.write(params)
        dB = get_parameter_value(params, "dB") * B
        print """
        Bfield = {BFIELD}
        dB = {DERIVATIVE_STEPSIZE}
        """.format(BFIELD=B, DERIVATIVE_STEPSIZE=dB)
        
    if B == args["ref_B"]:
        print "REFERENZ"
#FIX THIS because vdot will flatten axis!
        # unbedingt noch dx**2 einbauen, vorerst weggelassen, um Vergleichbarkeit mit notebook zu erleichtern
        tmp = np.absolute(np.sum(psi_ref.conjugate() * psi_ref, axis=1)) # * dx**2
        #np.copyto(correlation_table[i], tmp)
        np.copyto(corr, tmp)
        continue

    
    # This is a little fix for the problem stated below, if some ascii_files are 
    # missing, the calculation is repeated until all ascii files have been
    # generated. The sets are used, since this is only compares hashes of the
    # elements contained in lists, and is thus the fastest way for unique elems
    current_output = []
    fix_flag = 0
    while set(current_output) != set(ascii_files):
        run_code("input.xml", opt=args["multi_core"])
        Q_B = Q_Matrix.Time_Delay_Matrix(derivative_stepsize=dB)
        Q_B.write_eigenstates()
        Q_B.write_eigenvalues()
        # Remove ascii_files from previous calculation, as solve_xml_mumps 
        # sometimes omits the output from one mode, in which case the old ascii
        # file would be used again.
        remove_files(ascii_files)
        run_code("output.xml", opt=args["multi_core"])
        current_output = glob.glob(cwd+"/Q_states.*.coeff.*.purewavefunction.ascii")
        if fix_flag: print "Repeating calc for B={} due to bug!".format(B)
        fix_flag += 1


    allowed_modes = range(int(modes))
    for f in ascii_files:
        print "Reading {FILE}...".format(FILE=f)
        mode = int(ascii_pattern.search(f).group(1))
        psi_B = np.zeros_like(psi_ref[0])
        psi_B = np.loadtxt(f, usecols=(1,), dtype=complex,
                        converters={1: convert_to_complex})
        #correlation_table[i,mode] = np.absolute(np.vdot(psi_ref[mode], psi_B)) #* dx**2
        
        # Since the ordering of the output can't be controlled well, the 
        # reference psi and the current psi get misaligned, that is the mode
        # number of the ascii file is different from the one in psi_ref.
        # Therefore the correlation with every reference psi is calculated and
        # the maximum value is used as the correlation. Once a mode from psi_ref
        # has been used for a correlation, it can't be used for the other ascii
        # files to prevent using the same reference mode twice.

        # old
        #corr[mode] = np.absolute(np.vdot(psi_ref[mode], psi_B)) #* dx**2

        tmp = np.absolute(np.sum(np.conjugate(psi_ref[allowed_modes]) * psi_B, axis=1)) # * dx**2
        imax = tmp.argmax()
        use_mode = allowed_modes[imax]
        corr[use_mode] = tmp[imax]
        allowed_modes.remove(use_mode)

        if use_mode != mode:
            with open("max_mode.log", "a") as log:
                log.write("{BFIELD} {FILEMODE} {REFMODE}\n".format(
                    BFIELD=B, FILEMODE=mode, REFMODE=use_mode) )


    if args["backup"]:
        print "Renaming images..."
        rename_images()
        print "Backing up files..."
        backup_files( glob_arg="Q_states", prefix="B_{}_".format(B) )


# Save results to file:
# --------------------
with open("correlation.dat", "a") as ofile:
    print "Saving results to {}".format(ofile)
    HEADER = "Bfield" + "  mode {}  "*int(modes)
    HEADER = HEADER.format(*range(int(modes)))
    np.savetxt( ofile, 
                np.c_[B_values, correlation_table],
                header = HEADER,
                fmt = "%.2f" + " %.8f"*int(modes) )
