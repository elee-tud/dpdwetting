#!/home/elee/programs/Python-3.6.10/bin/python

from fileio.fileread import *
from scipy.ndimage import gaussian_filter1d
import sys
import operator
import numpy as np
from iocontrol.options import get_options as getopt

options=['-i', '-o', '-t']
types=['str', 'str', 'float']
defaults=['file.in', 'default', 1.0]
filename, outname, tau=getopt(sys.argv, options, types, defaults)

if outname=='default':
    outname=filename.replace('.out', '_gf_sig{:.1f}.out'.format(tau))

data,comment=file_to_array(filename)
dt=data[1,0]-data[0,0]
sigmat=tau/dt
new_data=np.zeros(data.shape)
new_data[:,0]=[t for t in data[:,0]]
for i in range(1,data.shape[1]):
    new_data[:,i]=gaussian_filter1d(data[:,i], sigma=sigmat)

outheader='#Data in {} after Gaussian filter with sigma={}\n'.format(filename, tau)
np.savetxt(outname, new_data, fmt='%1.5e', newline='\n', header=outheader)





