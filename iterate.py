#!/usr/bin/env python2.7
"""Iterate over different lengths of the one-side disordered wire"""

# TO DO:
# print input for checking
# Add accepting a list of different L or notch values instead of min, max, step
# let the program iterate different variables, implement either by choices
# or childparsers
# implement a dryrun for testing purpose, that is create dir and input.xml,
# but dont start subSGE.py
# rewrite update function, s.t. it is updating both args and parameters, maybe
# even add a key-value pair "L" to args for that purpose
# rewrite clean_up s.t. shell=true isnt used anymore
# include a --log flag to save all output going to stdout by default in logfile


# Import modules
# --------------
import datetime
today = datetime.date.today().strftime('%Y%m%d')
import os, sys, subprocess
from generate_boundary import generate_boundary
import math


# Parsing of command line
# -----------------------
import argparse
from argparse import ArgumentDefaultsHelpFormatter as help_formatter

parser=argparse.ArgumentParser(description="Iterate a parameter, "
                               "make dir, generate boundary, start job",
                               formatter_class=help_formatter) 

# --> imperative arguments 
parser.add_argument('itermin', metavar='min', type=int,
                    help="smallest value of variable to be iterated")
parser.add_argument('itermax', metavar='max', type=int,
                    help="highest value of variable to be iterated")
parser.add_argument('iterstep', type=int, help="step size of iteration")

# --> options
parser.add_argument('-d', '--date', nargs='?', metavar='yyyymmdd', type=str,
                    default=today, const='', help='date prepended to directory,'
                    ' if no argument is provided no date will be prepended')
parser.add_argument('-p', '--multi-core', metavar='#cores', type=int, nargs='?',
                    default=1, const=4, help='when running the code locally, '
                    'this option enables parallel computing. If no argument '
                    'is provided, all 4 cores will be used!')
parser.add_argument('-c', '--clean-up', action='store_true', default=False,
                    help='remove irrelevant output-files')
parser.add_argument('-I', '--draw-image', action='store_true', default=False,
                    help='draw a contour plot for each mode, store as jpq')

# --> parameters for generate_boundary
parser.add_argument('--delta', type=float, default=0.3,
                    help='range of random notch height in boundary')
parser.add_argument('--length', type=float, default=0.2,
                    help='length of individual notch')
parser.add_argument('-s', '--seed', nargs='?', default=1, const=None,
                    help='set the seed for the random number generator, pseudo-'
                    'random behaviour if no argument is given')
parser.add_argument('-o', '--outfile', metavar='filename.ending', type=str,
                    help='random boundary filename', default='rdmboundary.dat')

# --> parameters for solve_xml_mumps
parser.add_argument('-i', '--infile', metavar='filename.xml', type=str, 
                    help='name of input file to create for greens code',\
                    default='input.xml') 
parser.add_argument("--modes", type=float, default=2.5,
                    help="floor(modes) is the number of transverse modes")
parser.add_argument("--pphw", type=int, default=50, dest="points_per_halfwave", 
                    help="points_per_halfwave, accuracy-parameter")
parser.add_argument("-W", type=float, default=1.0, help="width of the wire")
parser.add_argument('-B', '--magnetic', nargs='?', metavar='B-field', type=str,
                    default=None, const='dynamic', help='Turn on the B-field, '
                    'if no argument is provided, the B-field is calculated '
                    's.t. the cyclotron radius is kept at 3*W')


# --> gather variables in 'args' & 'parameters' dictionary
# ----> set the resolution of generate_boundary higher than input.xml
args = vars(parser.parse_args())
dy = args["W"]/(1+args["points_per_halfwave"]*args["modes"])
args["ppn"] = int(2*args["length"]/dy)
args["notches"] = args["itermin"]

parameters = {'"modes"><':
              '"modes"> {} <'.format(args["modes"]),
              '"points_per_halfwave"><': 
              '"points_per_halfwave"> {} <'.format(args["points_per_halfwave"]),
              '"do_single_slices"><':
              '"do_single_slices"> 2 <',
              '"W"><': '"W"> {} <'.format(args["W"]),
              '"L"><': '"L"> {} <'.format((args["notches"]+10)*args["length"]),
              '"N_file"><': 
              '"N_file"> {} <'.format(1+args["ppn"]*(args["notches"]+10)),
              '"rdm_boundary"><': 
              '"rdm_boundary"> {} <'.format(args["outfile"])}

# ----> include the B-field in parameters
if args["magnetic"]:
    try:
        args["magnetic"] = float(args["magnetic"])

    except:
        if args["magnetic"]!="dynamic":
            print "ERROR: None or float value needed for B-field"
            print "aborting calculation"
            sys.exit(1)

        print "B-field will be derived from r_c = 3*W"
        kF = args["modes"]*math.pi/args["W"]
        args["magnetic"] = -1.0 * kF/(3*args["W"])

    finally:
        parameters['"Bfield"><'] = '"Bfield"> {} <'.format(args["magnetic"])

# ----> include pictures in parameters
if args["draw_image"]:
    parameters['"pics_min"><'] = '"pics_min"> 0 <'
    parameters['"pics_step"><'] = '"pics_step"> 1 <'
    parameters['"pics_N"><'] = '"pics_N"> floor($modes) <'
 

# Get template for input.xml
# --------------------------
if args["magnetic"]:
    template_name= "template_mag_input.xml"
else:
    template_name= "template_input.xml"
    
script_loc = sys.path[0] + "/"
with open(script_loc + template_name, "r") as template_file:
    template_text = "".join(template_file.readlines())


# Code splitted in several functions to make Iteration LOOP more clear
# --------------------------------------------------------------------
def create_input(template, parameters, **kwargs):
    """Create input.xml for greens code and generate boundary"""

    input_file = open(kwargs["infile"], "w")
    head, body = template.split("<!--splitter-->")

    for old in parameters:                          # do replacements only on
        head = head.replace(old, parameters[old])   # head for efficiency

    input_file.write(head)
    input_file.write(body)
    input_file.close()
    print "{INPUT} created".format(INPUT=kwargs["infile"])
 
    generate_boundary(**kwargs)
    print "{BOUNDARY} created".format(BOUNDARY=kwargs["outfile"])


def run_code(in_file, opt=args["multi_core"]):
    """Start the calculation via solve_xml_mumps"""

    if os.environ.get('TMPDIR') and os.environ.get('NSLOTS'):
        print "running code on cluster..."
        print "$TMPDIR", os.environ.get('TMPDIR')
        print "$NSLOTS", os.environ.get('NSLOTS')
        cmd = ("mpirun -machinefile {TMPDIR}/machines -np {NSLOTS} "
               "solve_xml_mumps -i {INPUT}").format(INPUT=in_file, **os.environ)
    else:
        if opt == 1:
            print "running code locally..."
            cmd = "solve_xml_mumps -i {INPUT}".format(INPUT=in_file)
        elif opt > 1 and opt < 5:
            print "running code locally usind {} cores...".format(opt)
            cmd = ("mpirun -np {} solve_xml_mumps -i {}").format(opt, in_file)

    subprocess.call(cmd.split())


def run_postprocessor(L, infile="S_matrix.dat"):
    """Calculate scattering matrix for each L and gather results in one file"""

    cmd = "S_Matrix.py"
    subprocess.call(cmd.split())
    outfile = "Smat_amp_seed_{}.dat".format(args["seed"])
    with open(infile) as src, open(outfile, "a") as dest:
        text = src.readlines()
        if not os.stat(outfile).st_size > 0:
            dest.write(text[0].replace("#", "#\tL"))
        dest.write("\t{}".format(L) + text[1])

    cmd += " -p"
    subprocess.call(cmd.split())
    outfile = "Smat_prb_seed_{}.dat".format(args["seed"])
    with open(infile) as src, open(outfile, "a") as dest:
        text = src.readlines()
        if not os.stat(outfile).st_size > 0:
            dest.write(text[0].replace("#", "#\tL"))
        dest.write("\t{}\t".format(L) + text[1])


def clean_up(**kwargs):
    """Remove the unnecessary files created by solve_xml_mumps"""

    try:
        subprocess.call("exrm -f {DOC} {infile} {outfile} Smat_*".format(
                        DOC="boundary.doc", **kwargs), shell=True)
    except:
        print "could not do automatic clean up"
    else:
        print "cleaning up..."



# Set up directory
# ----------------
dirname = "seed" + "_{}".format(args["seed"]) 
topdir = os.getcwd()

# -->some options for the directory name
if args["magnetic"]:
    dirname = "B_" + dirname
if args["date"]:
    dirname = args["date"] + "_" + dirname

if not os.path.exists(dirname):
    print "creating directory:  {}".format(dirname)
    os.mkdir(dirname)
else:
    print "WARNING: directory already exists, results will be appended to "\
          "exisiting data"
os.chdir(dirname)


# Iteration LOOP
# --------------
while args["notches"] <= args['itermax']:

    # +10 due to lead joining
    L = (args["notches"] + 10) * args["length"] 
    N_file = (args["notches"] + 10) * args["ppn"] + 1

    # update parameters
    parameters['"L"><'] = '"L"> {} <'.format(L)
    parameters['"N_file"><'] = '"N_file"> {} <'.format(N_file)

    create_input(template_text, parameters, **args)
    run_code(args["infile"])
    run_postprocessor(L)

    args["notches"] += args["iterstep"]


# Finalizing
# ----------
if args["clean_up"]:
    clean_up(**args)
os.chdir(topdir)
print "FINISHED calculation for seed {}".format(args["seed"])

################################################################################
# unfinished:

# def update(parameters, varlist, **kwargs):
#     """
#     Update the key value pairs of dictionary,
#     non existing pairs will be added to dictionary
# 
#     arguments:
#     dictionary  -   the dictionary that is being updated
#     keys        -   the keys whose values will be updated,
#                     if a key-value pair doesn't exist in the
#                     dictionary it will be added
#     values      -   the values to the corresponding key
#     """
# 
#     L = (notches + 10) * length         # +10 due to lead joining
# 
#     for key in varlist:
#         parameters['"'+key+'"><'] ='"'+key+'"> {} <'.format(eval(key))

