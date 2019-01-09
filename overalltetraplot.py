import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys


name = sys.argv[1]
name2 = sys.argv[2]


results = np.loadtxt('0.1/trend.dat')
z = np.arange(0,5,0.01)


a = (results[:,0])
b = (results[:,1])
c = (results[:,3])
d = (results[:,7])
e = (results[:,9])


results2 = np.loadtxt('0.2/trend.dat')


b2 = (results2[:,1])
c2 = (results2[:,3])
d2 = (results2[:,7])
e2 = (results2[:,9])


results3 = np.loadtxt('0.3/trend.dat')


b3 = (results3[:,1])
c3 = (results3[:,3])
d3 = (results3[:,7])
e3 = (results3[:,9])


results4 = np.loadtxt('0.4/trend.dat')


b4 = (results4[:,1])
c4 = (results4[:,3])
d4 = (results4[:,7])
e4 = (results4[:,9])


results5 = np.loadtxt('0.5/trend.dat')


b5 = (results5[:,1])
c5 = (results5[:,3])
d5 = (results5[:,7])
e5 = (results5[:,9])


results6 = np.loadtxt('0.6/trend.dat')


b6 = (results6[:,1])
c6 = (results6[:,3])
d6 = (results6[:,7])
e6 = (results6[:,9])


results7 = np.loadtxt('0.8/trend.dat')


b7 = (results7[:,1])
c7 = (results7[:,3])
d7 = (results7[:,7])
e7 = (results7[:,9])


results8 = np.loadtxt('0.01/trend.dat')


b8 = (results8[:,1])
c8 = (results8[:,3])
d8 = (results8[:,7])
e8 = (results8[:,9])


results9 = np.loadtxt('0.05/trend.dat')


b9 = (results9[:,1])
c9 = (results9[:,3])
d9 = (results9[:,7])
e9 = (results9[:,9])

results10 = np.loadtxt('0.15/trend.dat')


b10 = (results10[:,1])
c10 = (results10[:,3])
d10 = (results10[:,7])
e10 = (results10[:,9])


plt.figure(1)
#plt.subplot(221)
plt.title("Client Recruitment")
plt.xlabel("Valency of Client",fontsize=12)
plt.ylabel("Left vs Right",fontsize=12)


plt.plot(a,b8,'x',label="0.01")
plt.plot(a,b9,'x',label="0.05")
plt.plot(a,b,'x',label="0.1")
plt.plot(a,b10,'x',label="0.15")
plt.plot(a,b2,'x',label="0.2")
plt.plot(a,b3,'x',label="0.3")
plt.plot(a,b4,'x',label="0.4")
plt.plot(a,b5,'x',label="0.5")
plt.plot(a,b6,'x',label="0.6")
plt.plot(a,b7,'x',label="0.8")


leg = plt.legend(loc='best')
leg.get_frame().set_alpha(0.5)
plt.savefig(name)

#plt.subplot(222)

plt.figure(2)

plt.title("Client Location")
plt.xlabel("Valency of Client",fontsize=12)
plt.ylabel("Left vs Right",fontsize=12)


plt.plot(a,d8/e8,'x',label="0.01")
plt.plot(a,d9/e9,'x',label="0.05")
plt.plot(a,d/e,'x',label="0.1")
plt.plot(a,d10/e10,'x',label="0.15")
plt.plot(a,d2/e2,'x',label="0.2")
plt.plot(a,d3/e3,'x',label="0.3")
plt.plot(a,d4/e4,'x',label="0.4")
plt.plot(a,d5/e5,'x',label="0.5")
plt.plot(a,d6/e6,'x',label="0.6")
plt.plot(a,d7/e7,'x',label="0.8")


leg = plt.legend(loc='best')
leg.get_frame().set_alpha(0.5)

plt.savefig(name2)
#plt.show()       


