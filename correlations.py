#!/usr/bin/env python2.7

import numpy as np
import os
import glob
import re
import Q_Matrix
import T_Matrix

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
parser.add_argument('-g', '--full-grid',
                    action="store_true",
                    default=False,
                    help="""Calculate the correlations for entire cavity (e.g.
                    use all grid points) instead of the output interface (e.g. 
                    only last y-slice)""")

args=vars(parser.parse_args())


# Get parameters
with open("parameters.xml") as params_file:
    print "Reading parameters..."
    params = params_file.read()
    pphw = get_parameter_value(params, "points_per_halfwave")
    modes = get_parameter_value(params, "modes")
    dx = 1.0/(modes*pphw+1)
    r_ny = int(pphw*modes + 1)
    print """
    modes = {MODES}
    pphw = {PPHW}
    dx = {DX}
    r_ny = {R_NY}
    """.format(MODES=modes, PPHW=pphw, DX=dx, R_NY=r_ny)


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

TT = T_Matrix.T_Matrix()
TT.write_eigenstates()
TT.write_eigenvalues()
run_code("T_output.xml", opt=args["multi_core"])

ascii_files = glob.glob("[Q,T]_states.*.coeff.*.purewavefunction.ascii")
ascii_files.sort()
ascii_regex = r'[Q,T]_states\.(\d{,4})\.coeff\.(\d{,4})\.purewavefunction\.ascii'
ascii_pattern = re.compile(ascii_regex)

# Get file length to allocate psi_ref
with open(ascii_files[0]) as f:
    ascii_f_len = int(sum(1 for line in f))

    if args["full_grid"]:
        grid = ascii_f_len
    else:
        grid = r_ny

    grid_delta = ascii_f_len - grid
    psi_ref = np.zeros( (int(modes), grid), dtype=complex ) 
    chi_ref = np.zeros_like(psi_ref)

    #template = np.zeros( (int(modes), grid), dtype=complex ) 
    #psi_ref = {"Q": np.zeros_like(template), "T": np.zeros_like(template)}


for f in ascii_files:
    print "Reading {FILE}...".format(FILE=f)
    mode = int(ascii_pattern.search(f).group(1))

    if f.startswith("Q"):
        psi_ref[mode] = np.loadtxt(f, usecols=(1,), skiprows=grid_delta, 
                        dtype=complex, converters={1: convert_to_complex})
    elif f.startswith("T"):
        chi_ref[mode] = np.loadtxt(f, usecols=(1,), skiprows=grid_delta,
                        dtype=complex, converters={1: convert_to_complex})

    #psi_ref[f[0]][mode] = np.loadtxt(f, usecols=(1,), skiprows=grid_delta, 
    #                        dtype=complex, converters={1: convert_to_complex})

if args["backup"]:
    print "Renaming images..."
    rename_images()
    print "Backing up files..."
    backup_files(glob_arg="[Q,T]_states", 
                 prefix="B_ref_{}_".format(args["ref_B"]))


# Scan over all B values
# ----------------------
# the total number of B_values is 2*n + 1
B_values = np.array([ args["ref_B"] + i*args["delta_B"] 
                    for i in range(-args["n"], args["n"]+1) ])
correlation_table = np.zeros( (B_values.size, modes) )
correlation_table_chi = np.zeros_like(correlation_table)
with open("max_mode.log", "w") as log:
    log.write("# state Bfield file_mode ref_mode\n")

#for i,B in enumerate(B_values):
for B, corr, corr_chi in zip(B_values, correlation_table, correlation_table_chi):
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

        # unbedingt noch dx einbauen, vorerst weggelassen, um Vergleichbarkeit mit notebook zu erleichtern
        # tmp = np.absolute(np.sum(psi_ref.conjugate() * psi_ref, axis=1)) # * dx
        # tmp_chi = np.absolute(np.sum(chi_ref.conjugate() * chi_ref, axis=1)) # * dx
        # np.copyto(corr, tmp)
        # np.copyto(corr_chi, tmp_chi)

        # use different correlation (cos alpha) which will give 1.0 for reference
        # additionally all dx will drop out
        tmp = np.ones_like(corr)
        np.copyto(corr, tmp)
        np.copyto(corr_chi, tmp)

        #np.copyto(correlation_table[i], tmp)
        continue

    
    # This is a little fix for the problem stated below, if some ascii_files are 
    # missing, the calculation is repeated until all ascii files have been
    # generated. The sets are used, since this is only compares hashes of the
    # elements contained in lists, and is thus the fastest way for unique elems
    current_output = []
    fix_flag = 0
    while set(current_output) != set(ascii_files):

        # Remove ascii_files from previous calculation, as solve_xml_mumps 
        # sometimes omits the output from one mode, in which case the old ascii
        # file would be used again.
        remove_files(ascii_files)

        run_code("input.xml", opt=args["multi_core"])

        Q_B = Q_Matrix.Time_Delay_Matrix(derivative_stepsize=dB)
        Q_B.write_eigenstates()
        Q_B.write_eigenvalues()
        run_code("output.xml", opt=args["multi_core"])

        TT = T_Matrix.T_Matrix()
        TT.write_eigenstates()
        TT.write_eigenvalues()
        run_code("T_output.xml", opt=args["multi_core"])

        current_output = glob.glob("[Q,T]_states.*.coeff.*.purewavefunction.ascii")
        if fix_flag: print "Repeating calc for B={} due to bug!".format(B)
        fix_flag += 1


    allowed_modes = range(int(modes))
    allowed_modes_chi = range(int(modes))
    for f in ascii_files:
        print "Reading {FILE}...".format(FILE=f)
        mode = int(ascii_pattern.search(f).group(1))
        psi_B = np.zeros_like(psi_ref[0])
        psi_B = np.loadtxt(f, usecols=(1,), dtype=complex, skiprows=grid_delta,
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

        if f.startswith("Q"):
            bra_ket = np.absolute(np.sum(
                        np.conjugate(psi_ref[allowed_modes]) * psi_B, axis=1
                      ))
            ref_norm = np.sqrt(np.absolute(np.sum(
                        np.conjugate(psi_ref[allowed_modes]) * psi_ref[allowed_modes], axis=1
                       )))
            current_norm = np.sqrt(np.absolute(
                            np.sum(np.conjugate(psi_B) * psi_B)
                           ))
            cos_alpha = bra_ket / (ref_norm * current_norm)
            imax = cos_alpha.argmax()
            use_mode = allowed_modes[imax]
            corr[use_mode] = cos_alpha[imax]
            allowed_modes.remove(use_mode)
        elif f.startswith("T"):
            bra_ket = np.absolute(np.sum(
                        np.conjugate(chi_ref[allowed_modes_chi]) * psi_B, axis=1
                      ))
            ref_norm = np.sqrt(np.absolute(np.sum(
                        np.conjugate(chi_ref[allowed_modes_chi]) * chi_ref[allowed_modes_chi], axis=1
                       )))
            current_norm = np.sqrt(np.absolute(np.sum(np.conjugate(psi_B) * psi_B)))
            cos_alpha = bra_ket / (ref_norm * current_norm)
            imax = cos_alpha.argmax()
            use_mode = allowed_modes_chi[imax]
            corr_chi[use_mode] = cos_alpha[imax]
            allowed_modes_chi.remove(use_mode)

        if use_mode != mode:
            with open("max_mode.log", "a") as log:
                log.write("{STATE} {BFIELD} {FILEMODE} {REFMODE}\n".format(
                    STATE=f[0], BFIELD=B, FILEMODE=mode, REFMODE=use_mode) )


    if args["backup"]:
        print "Renaming images..."
        rename_images()
        print "Backing up files..."
        backup_files( glob_arg="[Q,T]_states", prefix="B_{}_".format(B) )


# Save results to file:
# --------------------
with open("correlation.dat", "a") as Qfile, open("correlation_T.dat", "a") as Tfile:
    HEADER = "Bfield" + "  mode {}  "*int(modes)
    HEADER = HEADER.format(*range(int(modes)))
    print "Saving results to {}".format(Qfile)
    np.savetxt( Qfile, 
                np.c_[B_values, correlation_table],
                header = HEADER,
                fmt = "%.2f" + " %.8f"*int(modes) )
    print "Saving results to {}".format(Tfile)
    np.savetxt( Tfile, 
                np.c_[B_values, correlation_table_chi],
                header = HEADER,
                fmt = "%.2f" + " %.8f"*int(modes) )


# Document settings used in calculation:
# ------------------------------------
with open("correlation.doc", "w") as docfile:
    doc = """# Settings used by correlation.py
    reference B: {ref_B}
    delta B: {delta_B}
    n: {n}
    mutlti core: {multi_core}
    backup: {backup}
    full grid: {full_grid}
    """.format(**args)
    docfile.write(doc)


