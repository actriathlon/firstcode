import numpy as np
import sys


fname = sys.argv[1]
#data = np.loadtxt(fname,skiprows= 100)
eqtime = sys.argv[2]
column = sys.argv[3]
skips = sys.argv[4]
skip = int(skips)

data = np.loadtxt(fname,skiprows= skip)

#print chainlength
print eqtime
eq = int(eqtime)
c = int(column)

a = np.mean(data[eq:,c])
b = np.std(data[eq:,c])
print a, b


