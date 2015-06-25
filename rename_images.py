#!/usr/bin/env python2.7

# Currently only working on cwd, improve by adding a path as arg
def rename_images():
    """ Renames the *.jpg files created by solve_xml_mumps in the current
    working directory such that they are sorted in a more convenient way.
    """
    import glob
    import os
    import re
    import sys
    import shutil

    pic_prefix = "pic"
    pic_rename_marker = "layer0"
    pics = glob.glob("*{}*".format(pic_rename_marker)) 
    #pics = [l for l in os.listdir(path) if pic_rename_marker in l]
    
    if pics: 
        print "{} files for renaming found".format(len(pics))

    else:
        print """Warning: no file containing marker '{MARKER}' found 
        Renaming has probably already been performed
        Aborting execution
        """.format(MARKER=pic_rename_marker)
        sys.exit(10)
    
    # if the jobname contains a '.' then this will fail!
   # start = pics[0].find(".") + 1
   # end = pics[0].find(".", start)
   # jobname = pics[0][start:end] 
   # print "JOBNAME {}".format(jobname)

    regex = r"""
            {PREFIX}\.
            (?P<jobname>\w+)\.
            (?P<id>(\d+\.)?\d+)\.
            (?P<case>\w+)\.
            \d+\.
            {MARKER}\.
            (?P<suffix>\S+)
            """.format(PREFIX=pic_prefix, MARKER=pic_rename_marker)

    pattern = re.compile(regex, re.VERBOSE)

    # can be used to revert changes
    pics_new = []
    
    cnt_jpg = 0
    cnt_dat = 0

    print "Renaming following file(s) from OLD to NEW:"

    for pic in pics:
        
        match = pattern.search(pic)

        if pic.endswith("jpg"):
            cnt_jpg += 1
        elif pic.endswith("average_over-y.dat"):
            cnt_dat += 1
        else:
            print "Error while reading files, unknown file type included"
            sys.exit(1)


        new_name =  "{PREFIX}.{JOBNAME}.{CASE}.{ID}.{SUFFIX}".format(
                    PREFIX = pic_prefix, 
                    JOBNAME = match.group("jobname"), 
                    CASE = match.group("case"), 
                    ID = match.group("id"), 
                    SUFFIX = match.group("suffix"))

        print "{OLD}  ->  {NEW}".format(OLD=pic, NEW=new_name)
        shutil.move(pic, new_name) 
        pics_new.append(new_name)


    if len(pics) == len(pics_new):
        print "Renamed {COUNT} jpg file(s)".format(COUNT=cnt_jpg)
        print "Renamed {COUNT} dat file(s)".format(COUNT=cnt_dat)

    else:
        print "Error: Some files were not renamed"
        sys.exit(2)


if __name__ == "__main__":

    import os

    cwd = os.getcwd()
    print "Renaming images in {CWD}".format(CWD=cwd)
    rename_images()
