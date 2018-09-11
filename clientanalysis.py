import numpy as np
import sys
import math

data = np.loadtxt('bonddate.dat')
cluster = int(sys.argv[1])
bulk = int(sys.argv[2])
choose = int(sys.argv[3])
ncl = int(sys.argv[4])
nbu = int(sys.argv[5])

if choose ==3:
	a = float(data[0]+data[1]+data[2])/cluster
	b = float(data[4]+data[5]+data[6])/bulk
	c =data[0]+data[1]+data[2]
	d = data[4]+data[5]+data[6]
	e = data[10]+data[11]+data[12] #check this
        g = float(data[0]+data[1]+data[2])/ncl
        h = float(data[4]+data[5]+data[6])/nbu

elif choose ==2:
        a = float(data[0]+data[1])/cluster
        b = float(data[3]+data[4])/bulk
	c = data[0]+data[1]
	d = data[3] + data[4]
	e = data[8]+data[9] 
        g = float(data[0]+data[1])/ncl
        h = float(data[3]+data[4])/nbu

elif choose ==1:
	a = float(data[0])/cluster
	b = float(data[2])/bulk
	c = data[0]
	d = data[2]
	e = data[9]   #add this value in
        g = float(data[0])/ncl
        h = float(data[2])/nbu


print a,b,a/b,c,d,e
print 'ratios',c/e,d/e

if choose == 3:
	f = open('../bindingaverages3.dat','a')
elif choose == 2:
	f = open('../bindingaverages2.dat','a')
elif choose == 1:
	f = open('../bindingaverages1.dat','a')


f.write(str(a))
f.write(' ')
f.write(str(b))
f.write(' ')
f.write(str(c))
f.write(' ')
f.write(str(d))
f.write(' ')
f.write(str(e))
f.write(' ')
f.write(str(g))
f.write(' ')
f.write(str(h))
f.write(' ')
f.write(str(cluster))
f.write(' ')
f.write(str(bulk))
f.write(' ')
f.write(str(ncl))
f.write(' ')
f.write(str(nbu))
f.write('\n')
