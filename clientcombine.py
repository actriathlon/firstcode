import numpy as np
import sys
import math


chosen = int(sys.argv[1])
z = sys.argv[2]

inname = sys.argv[3]
outname = sys.argv[4]

data = np.loadtxt(inname)

a = np.mean(data[:,1])
b = np.mean(data[:,2])
c = np.mean(data[:,3])
d = np.mean(data[:,4])
e = np.mean(data[:,5])
e1 = np.std(data[:,1])
e2 = np.std(data[:,2])
e3 = np.std(data[:,3])
e4 = np.std(data[:,4])
e5 = np.std(data[:,5])


o = open(outname,'a')

o.write(str(chosen))
o.write(' ')
o.write(str(a))
o.write(' ')
o.write(str(e1))
o.write(' ')
o.write(str(b))
o.write(' ')
o.write(str(e2))
o.write(' ')
o.write(str(c))
o.write(' ')
o.write(str(e3))
o.write(' ')
o.write(str(d))
o.write(' ')
o.write(str(e4))
o.write(' ')
o.write(str(e))
o.write(' ')
o.write(str(e5))
o.write('\n')
