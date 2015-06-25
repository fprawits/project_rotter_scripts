#!/usr/bin/env python2.7
# TO DO:
# exit program if no parameter is iterated
# Add accepting a list of different L values instead of min, max, step
# let the program iterate different variables, implement either by choices
# or childparsers
# implement a dryrun for testing purpose, that is create dir and input.xml,
# but dont start subSGE.py
# rewrite update function, s.t. it is updating both args and parameters, maybe
# even add a key-value pair "L" to args for that purpose
"""Iterate a parameter, make dir, generate boundary, start job"""

import datetime, argparse, os, sys
from argparse import ArgumentDefaultsHelpFormatter as help_formatter
from generate_boundary import generate_boundary

script_loc = sys.path[0] + "/"
today = datetime.date.today().strftime('%Y%m%d')

parser=argparse.ArgumentParser(description="Iterate a parameter, "
                               "make dir, generate boundary, start job",
                               formatter_class=help_formatter) 
# parser.add_argument('itervar', choices=['delta', 'l', 'notches', 'ppn'],
#                    help='choose which parameter shall be iterated')
parser.add_argument('itermin', metavar='min', type=int,
                    help="smallest value of variable to be iterated")
parser.add_argument('itermax', metavar='max', type=int,
                    help="highest value of variable to be iterated")
parser.add_argument('iterstep', type=int, help="step size of iteration")
parser.add_argument('--delta', type=float, default=0.3,
                    help='range of random notch height in boundary')
parser.add_argument('--length', type=float, default=0.2,
                    help='length of individual notch')
parser.add_argument('--notches', type=int, default=40,
                    help='number of notches in boundary')
parser.add_argument('--ppn', type=int, default=100,
                    help="total number of points in x direction")
parser.add_argument('-i', '--infile', metavar='filename.xml', type=str, 
                    help='name of input file used by greens code',\
                    default='input.xml') 
parser.add_argument('-o', '--outfile', metavar='filename.ending', type=str,
                    help='random boundary filename', default='rdmboundary.dat')
parser.add_argument('-d', '--date',nargs='?', metavar='yyyymmdd', type=str,
                    default=today, const='', help='date prepended to directory,'
                    ' if no argument is provided no date will be prepended')
parser.add_argument('-s', '--seed', nargs='?', default=1, const=None,
                    help='set the seed for the random number generator, pseudo-'
                    'random behaviour if no argument is given')
parser.add_argument('-S', '--showseed', action='store_true', default=False,
                    help="show the used seed in the directory name")
parser.add_argument("--modes", type=float, default=2.5,
                    help="floor(modes) is the number of transverse modes")
parser.add_argument("--pphw", type=int, default=50, dest="points_per_halfwave", 
                    help="points_per_halfwave, accuracy-parameter")
parser.add_argument("-W", type=float, default=1.0, help="width of the wire")

args = vars(parser.parse_args())
parameters = {'"modes"><': '"modes"> {} <'.format(args["modes"]),
              '"points_per_halfwave"><': '"points_per_halfwave"> {} <'.format(args["points_per_halfwave"]),
              '"do_single_slices"><': '"do_single_slices"> 2 <',
              '"W"><': '"W"> {} <'.format(args["W"]),
              '"L"><': '"L"> {} <'.format((args["notches"]+10)*args["length"]),
              '"N_file"><': '"N_file"> {} <'.format(1+args["ppn"]*(args["notches"]+10)),
              '"rdm_boundary"><': '"rdm_boundary"> {} <'.format(args["outfile"])}


itername = "L"
dirbase = "random_boundary"
topdir = os.getcwd()

# some options for the directory name
if args["date"]: 
    args["date"] += "_"
if args["showseed"]: 
    dirbase = dirbase + "_s_" + str(args["seed"]) 
dirbase = args["date"] + dirbase + "_{}_".format(itername)

# get the template for input.xml
template_file = open(script_loc + "template_input.xml", "r")
template_text = "".join(template_file.readlines())
template_file.close()



def set_up(dirname, template, parameters, **kwargs):
    """Create directory and set up all input files for calculation"""

    os.mkdir(dirname)
    os.chdir(dirname)
    input_file = open(kwargs["infile"], "w")
    head, body = template.split("<!--splitter-->")

    for old in parameters:                          # do replacements only on
        head = head.replace(old, parameters[old])   # head for efficiency

    input_file.write(head)
    input_file.write(body)
    input_file.close()
 
    generate_boundary(**kwargs)
    os.chdir(topdir)
    return



# Iteration LOOP
# --------------

# set the resolution of generate_boundary higher than input.xml
dy = args["W"]/(1+args["points_per_halfwave"]*args["modes"])
args["ppn"] = int(2*args["length"]/dy)
itervar = args['itermin']

while itervar <= args['itermax']:

    args["notches"] = itervar       # hier sollte ein itername her!!!
    L = (args["notches"] + 10) * args["length"] # +10 due to lead joining
    N_file = (args["notches"] + 10) * args["ppn"] + 1

   # update(parameters, ["L", "N_file"], **args)

    parameters['"L"><'] = '"L"> {} <'.format(L)
    parameters['"N_file"><'] = '"N_file"> {} <'.format(N_file)

    dirname = dirbase + str(L)     #str.format might be good, since L is float

    set_up(dirname, template_text, parameters, **args)

    itervar += args["iterstep"]

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

