import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys


name = sys.argv[1]


results = np.loadtxt('0.1/trend.dat')
#z = np.arange(0,5,0.01)
#comb=[]

plt.figure(1)

for nm in range(5):

	comba=np.arange(13)
	print comba
	combb=np.arange(13)
	print combb
	combc=np.arange(13)
	print combc




	comba[0] = (results[nm,1])
	combb[0] = (results[nm,7])
	combc[0] = (results[nm,9])


	results2 = np.loadtxt('0.2/trend.dat')


	comba[1] = (results2[nm,1])
	combb[1] = (results2[nm,7])
	combc[1] = (results2[nm,9])


	results3 = np.loadtxt('0.3/trend.dat')


	comba[2] = (results3[nm,1])
	combb[2] = (results3[nm,7])
	combc[2] = (results3[nm,9])


	results4 = np.loadtxt('0.4/trend.dat')


	comba[3] = (results4[nm,1])
	combb[3] = (results4[nm,7])
	combc[3] = (results4[nm,9])


	results5 = np.loadtxt('0.5/trend.dat')


	comba[4] = (results5[nm,1])
	combb[4] = (results5[nm,7])
	combc[4] = (results5[nm,9])


	results6 = np.loadtxt('0.6/trend.dat')


	comba[5] = (results6[nm,1])
	combb[5] = (results6[nm,7])
	combc[5] = (results6[nm,9])


	results7 = np.loadtxt('0.8/trend.dat')


	comba[6] = (results7[nm,1])
	combb[6] = (results7[nm,7])
	combc[6] = (results7[nm,9])


	results8 = np.loadtxt('0.01/trend.dat')


	comba[7] = (results8[nm,1])
	combb[7] = (results8[nm,7])
	combc[7] = (results8[nm,9])

	results9 = np.loadtxt('0.05/trend.dat')


	comba[8] = (results9[nm,1])
	combb[8] = (results9[nm,7])
	combc[8] = (results9[nm,9])


	results10 = np.loadtxt('0.15/trend.dat')


	comba[9] = (results10[nm,1])
	combb[9] = (results10[nm,7])
	combc[9] = (results10[nm,9])


        results11 = np.loadtxt('0.25/trend.dat')


        comba[10] = (results11[nm,1])
        combb[10] = (results11[nm,7])
        combc[10] = (results11[nm,9])
 
	results12 = np.loadtxt('0.35/trend.dat')


        comba[11] = (results12[nm,1])
        combb[11] = (results12[nm,7])
        combc[11] = (results12[nm,9])
        
	results13 = np.loadtxt('0.7/trend.dat')


        comba[12] = (results13[nm,1])
        combb[12] = (results13[nm,7])
        combc[12] = (results13[nm,9])



	z=np.arange(0.1,1.4,0.1)
	z[6] = 0.8
	z[7] = 0.01
	z[8] = 0.05
	z[9] = 0.15
	z[10] = 0.25
	z[11] = 0.35
	z[12] = 0.7
	print z


	#print a
	tr='valency'+' '+str(nm+1)
	plt.plot(z,comba,'x',label=tr)



#np.append(comb(b))

#plt.figure(1)
#plt.subplot(221)
plt.title("Client Binding")
plt.xlabel("Valency of Client",fontsize=12)
plt.ylabel("Left vs Right",fontsize=12)
plt.ylim(0.0,None)
plt.xlim(0.0,0.8)


leg = plt.legend(loc='best')
leg.get_frame().set_alpha(0.5)
plt.savefig(name)

#plt.subplot(222)



