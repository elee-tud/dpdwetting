#!/usr/bin/python

import sys
import numpy as np
from fileio.fileread import file_to_array as fta
from iocontrol.options import get_options




inname, outname, skip=get_options(sys.argv, ['-i', '-o', '-s'], ['str', 'str', 'int'], ['file.in', 'file.out', 1])
indata, comment=fta(inname)
outdata=[]
for idx, data in enumerate(indata):
    if idx%skip==0:
        outdata.append(data)
header=''
for c in comment:
    header += c

np.savetxt(outname, np.array(outdata), fmt='%15.4f', delimiter='', newline='\n', header=header, footer='' ,comments='#')
print("Data is written in %s."%outname)



