#!/home/elee/programs/Python-3.6.10/bin/python

import numpy as np

def file_to_array(filename, rowlen=-1):
    comment_string=['#', '@', '\n']
    fp=open(filename, "r")
    data=[]
    comment=''
    for line in fp.readlines():
        if line[0] not in comment_string:
            linedata=[]
            if rowlen==-1:
                for d in line.split():
                    if d=='nan':
                        linedata.append(np.nan)
                    else:
                        linedata.append(float(d))
            else:
                for d in line.split()[:rowlen]:
                    if d=='nan':
                        linedata.append(np.nan)
                    else:
                        linedata.append(float(d))

            data.append(linedata)
        else:
            comment += line
    fp.close()

    return np.array(data), comment





