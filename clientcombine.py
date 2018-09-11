import numpy as np
import sys
import math

data = np.loadtxt('bindingaverages3.dat')
z = sys.argv[1]

cluster = np.mean(data[:,2])
bulk = np.mean(data[:,3])
unbound = np.mean(data[:,4])
sitescluster = np.mean(data[:,7])
sitesbulk = np.mean(data[:,8])
neighbourscluster = np.mean(data[:,9])
neighboursbulk = np.mean(data[:,10])

e1 = np.std(data[:,2])
e2 = np.std(data[:,3])
e3 = np.std(data[:,4])
e7 = np.std(data[:,7])
e8 = np.std(data[:,8])
e9 = np.std(data[:,9])
e10 = np.std(data[:,10])



a = cluster/sitescluster
b = bulk/sitesbulk
c = cluster/neighbourscluster
d = bulk/neighboursbulk
e = cluster/unbound
f = bulk/unbound


err1 = np.sqrt(((e1/unbound)**2)+(((cluster/(unbound*unbound))*e3)**2))
err2 = np.sqrt(((e2/unbound)**2)+(((bulk/(unbound*unbound))*e3)**2))
err3 = np.sqrt(((e1/sitescluster)**2)+(((cluster/(sitescluster*sitescluster))*e7)**2))
err4 = np.sqrt(((e2/sitesbulk)**2)+(((bulk/(sitesbulk*sitesbulk))*e8)**2))
err5 = np.sqrt(((e1/neighbourscluster)**2)+(((cluster/(neighbourscluster*neighbourscluster))*e9)**2))
err6 = np.sqrt(((e2/neighboursbulk)**2)+(((bulk/(neighboursbulk*neighboursbulk))*e10)**2))
err7 = np.sqrt(((err3/b)**2)+(((a/(b*b))*err4)**2))
err8 = np.sqrt(((err5/d)**2)+(((c/(d*d))*err6)**2))
#err8 = np.sqrt(((e1/()**2)+(((c/(d*d))*err6)**2))
err9 = np.sqrt(((err3/unbound)**2)+(((a/(unbound*unbound))*e3)**2))
err10 = np.sqrt(((err4/unbound)**2)+(((b/(unbound*unbound))*e3)**2))
err11 = np.sqrt(((err1/f)**2)+(((e/(f*f))*err2)**2))




o = open('../clientresults3.dat','a')

o.write(str(z))
o.write(' ')
o.write(str(cluster/unbound))
o.write(' ')
o.write(str(err1))
o.write(' ')
o.write(str(bulk/unbound))
o.write(' ')
o.write(str(err2))
o.write(' ')
o.write(str(a))
o.write(' ')
o.write(str(err3))
o.write(' ')
o.write(str(b))
o.write(' ')
o.write(str(err4))
o.write(' ')
o.write(str(c))
o.write(' ')
o.write(str(err5))
o.write(' ')
o.write(str(d))
o.write(' ')
o.write(str(err6))
o.write(' ')
o.write(str(a/b))
o.write(' ')
o.write(str(err7))
o.write(' ')
o.write(str(c/d))
o.write(' ')
o.write(str(err8))
o.write(' ')
o.write(str(a/unbound))
o.write(' ')
o.write(str(err9))
o.write(' ') 
o.write(str(b/unbound))
o.write(' ')
o.write(str(err10))
o.write(' ')
o.write(str(e/f))
o.write(' ')
o.write(str(err11))
o.write(' ')
o.write(str(cluster))
o.write(' ')
o.write(str(bulk))
o.write(' ')
o.write(str(unbound))
o.write(' ')
o.write('\n')
