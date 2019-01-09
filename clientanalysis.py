import numpy as np
import sys
import math

data = np.loadtxt('bonddate.dat')
choose = int(sys.argv[1])
octname = sys.argv[2]

if choose ==4:
        a= data[19]
	b = data[20]
        c = data[18]

elif choose ==5:
	a= data[22]
        b = data[23]
        c = data[21]

elif choose ==3:

        a= data[16]
        b = data[17]
        c = data[15]



elif choose ==2:

        a= data[13]
        b = data[14]
        c = data[12]



elif choose ==1:
        a= data[10]
        b = data[11]
        c = data[9]


loc=np.loadtxt('location.dat')
d=loc[0]
e=loc[1]

f=open(octname,'a')

f.write(str(choose))
f.write(' ')
f.write(str(a))
f.write(' ')
f.write(str(abs(b)))
f.write(' ')
f.write(str(c))
f.write(' ')
f.write(str(d))
f.write(' ')
f.write(str(e))
f.write('\n')
