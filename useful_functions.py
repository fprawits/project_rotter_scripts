"""Collection of often used functions
    
    update(text, **kwargs):
        Update parameters of an xml file with new values.

    run_code(infile, opt=args["multi_core"]):
        Start the calculation via solve_xml_mumps.
    
    conver_to_complex(s):
        Convert a string of the form (x,y) to a complex number z = x+1j*y.

    ascii_to_complex(string):
        Convert the ascii '(x,y)' strings as found in solve_xml_mumps output to
        a complex number z = x + i*y

    list_to_complex(t):
        Convert a 2 element tuple or list t to complex number.
"""

import ast
import re
import os
import glob
import sys
import shutil
import subprocess

def update(text, **kwargs):
    """Return text with new VALUES assigned to parameters (KEYS)

    text:       string containing all contents of xml file to be updated
    kwargs:     dictionary of KEY-VALUE pairs where KEY is the parameter
                in xml file that shall be updated to VALUE

    The update is performed via a REGEX-search on the string text, containing
    all the content of the xml file. Each KEY in kwargs corresponds to a
    parameter of the xml-file and will be updated to its assigned VALUE.
    """

    #import re
    #import sys

    #regular expression for arbitrary parameter (digit) value
    re_digit = r"""
            [-+]?                   #optional sign
            (
                \d+(\.\d*)? |       #number with or without decimal places or
                \.\d+               #number only with decimal places
            ) 
            ([eE][-+]?\d+)?         #number in floating point representation
            """
    
    for key, value in kwargs.iteritems():

        pattern = re.compile(r"""
                            <param \s name="{PARAM}"> \s* {DIGIT} \s* </param>
                            """.format(PARAM=key, DIGIT=re_digit), 
                            re.VERBOSE)
        replacement = r'<param name="{PARAM}"> {VALUE} </param>'.format(
                      PARAM=key, VALUE=value)

        text, cnt = pattern.subn(replacement, text)

        if not cnt:
            print """ERROR: Updating {PARAM} to {VALUE} in xml-file failed
            Aborting script execution""".format(PARAM=key, VALUE=value)
            sys.exit(1)

    return text


def get_parameter_value(text, param):
    
    #import re

    #regular expression for arbitrary parameter (digit) value
    re_digit = r"""
            [-+]?                   #optional sign
            (
                \d+(\.\d*)? |       #number with or without decimal places or
                \.\d+               #number only with decimal places
            ) 
            ([eE][-+]?\d+)?         #number in floating point representation
            """
    #regular expression for arbitrary binary operator
    re_op = r'[+\-\*/]'
    #regular expression for dereferenced variable
    re_var = r'\$\w+'

    pattern = re.compile(r"""
                        <param \s name="{PARAM}"> \s*
                        ({VAR} \s* {OP})* \s*
                        (?P<value>{DIGIT}) \s*
                        ({OP} \s* {VAR})* \s* </param>
                        """.format(
                                    PARAM=param,
                                    DIGIT=re_digit,
                                    OP=re_op,
                                    VAR=re_var), 
                        re.VERBOSE)

    value = float(pattern.search(text).group('value'))
    return value


def run_code(in_file, opt=1):
    """Start the calculation via solve_xml_mumps"""

    #import os
    #import subprocess

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
            print "running code locally using {} cores...".format(opt)
            cmd = ("mpirun -np {} solve_xml_mumps -i {}").format(opt, in_file)

    subprocess.call(cmd.split())


def convert_to_complex(s):
    """Convert a string of the form (x,y) to a complex number z = x+1j*y."""

    regex = re.compile(r'\(([^,\)]+),([^,\)]+)\)')
    x, y = map(float, regex.match(s).groups())

    return x + 1j*y


def asciii_to_complex(string):
    """Convert the ascii '(x,y)' strings as found in solve_xml_mumps output to
    a complex number z = x + i*y
    
    Note: since this function involves ast.literal_eval, it is usually faster
    to use the regex search from convert_to_complex.
    """
    
    x,y = ast.literal_eval(string)
    return complex(x,y)


def list_to_complex(t):
    """Convert a 2 element tuple or list t to complex number"""

    if len(t) == 2:
        try:
            z = complex(float(t[0]), float(t[1]))

        except:
            print "list_to_complex: bad input, list or tuple required"
            raise
    else:
        print "list_to_complex: Input list/tuple must consist of two elements"

    return z


def backup_files(file_list=None, prefix="bak_", bak_dir="bak", glob_arg=None):
    """Backup all files either found in file_list or containing the substring 
    glob_arg, by renaming them with prefix and moving them to bak_dir (can be
    cwd if bak_dir is '' or None)
    """
# Include handling of files with full path name

    if bak_dir:
        prefix = os.path.join(bak_dir,prefix)
        try:
            os.mkdir(bak_dir)
        except OSError:
            if not os.path.isdir(bak_dir):
                raise


    if not file_list:
        file_list = []
    
# The list comprehension used here to avoid adding another prefix to an already
# backed up file when using this function inside a loop is not ideal as the 
# initial glob_list created by glob.glob can grow very large for many loops.
# this can make this function very slow, therefore it might be better to move
# all files to a backup directory in order to keep the working directory clean
# and the initial glob list small in size
    if glob_arg:
        glob_list = glob.glob( "*{}*".format(glob_arg) )
        glob_list = [ f for f in glob_list if not f.startswith(prefix) ]
        file_list += glob_list

    # file_list should contain only unique entries
    file_list = set(file_list)
    
    for f in file_list:
        print "Moving {OLD} to {NEW}".format(OLD=f, NEW=prefix+f)
        shutil.move(f, prefix+f)


def remove_files(file_list):
    """Iterate through all files in file_list and try to remove them. If the 
    file names have no path to them, the current working directory is assumed
    as their path.
    """

    for f in file_list:
        try:
            print "Removing file {}".format(f)
            os.remove(f)

        except OSError:
            print(
                    "OSError while calling remove_files:\n"
                    "file {} probably not found ".format(f)
                 )

