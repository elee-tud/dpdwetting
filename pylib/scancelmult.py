#!/home/elee/programs/Python-3.6.10/bin/python

import sys
import subprocess as sp


nstart=int(sys.argv[1])
nend=int(sys.argv[2])

for i in range(nstart, nend+1):
    sp.check_output("scancel "+str(i), shell=True)
