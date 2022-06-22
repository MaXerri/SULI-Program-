#! /bin/bash
#mario Xerri 6/15/22

import pandas as pd
import numpy as np 
import subprocess as sub
import os 
import subprocess 
import sys 

#just a random cell input and could be changed to prompt user input
cell = "C 10 10 8 90 90 78"

bash_script = subprocess.Popen(["/home/mxerri/SULI/ncdist6_14_22/ncdist/dc7unsrtdist_mat --info"],shell = True,executable = "/bin/bash",stdin = subprocess.PIPE,stdout = subprocess.PIPE)
bash_script.stdin.write(cell)
out = bash_script.communicate()[0]
bash_script.stdin.close()

#this index  of the list from the split function can be changed based on what do you want dc7unsrtdist_mat to output
primitive =(((out.split("\n")[1]).split(":")[1]).lstrip()).split()

print([float(x) for x in primitive])

#this is just what is written above as a function
def toPrimitive(cell):
        """ This converts the non-primitive cells and primitive cells to niggli reduced primitive G6 representation
            The code called was devleoped by Herbert Bernstein and Lawrence Andrews and can be found in https://github.com/yayahjb/ncdist

            Parameters: pandas DataFrame of a b c  alpha beta gamma spg_type spg z
        """

        cell_string = " ".join(str(item) for item in cell)
        bash_script = subprocess.Popen(["/home/mxerri/SULI/ncdist21J/ncdist/dc7unsrtdist_mat --info"],shell = True,executable = "/bin/bash",stdin = subprocess.PIPE,stdout = subprocess.PIPE)
        bash_script.stdin.write(cell_string)
        out = bash_script.communicate()[0]
        bash_script.stdin.close()
        primitive =(((out.split("\n")[1]).split(":")[1]).lstrip()).split()
        return([float(x) for x in primitive])





