#!/bin/bash
import pandas as pd
import numpy as np 
import subprocess as sub
import os 
import subprocess 
import sys 

cell = "C 10 10 8 90 90 78"

bash_script = subprocess.Popen(["/home/mxerri/SULI/ncdist6_14_22/ncdist/dc7unsrtdist_mat --info"],shell = True,executable = "/bin/bash",stdin = subprocess.PIPE,stdout = subprocess.PIPE)
bash_script.stdin.write(cell)
out = bash_script.communicate()[0]
bash_script.stdin.close()
primitive =(((out.split("\n")[0]).split(":")[1]).lstrip()).split()

#print(primitive)




