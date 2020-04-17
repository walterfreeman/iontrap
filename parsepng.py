import numpy as np
import sys
import matplotlib.pyplot as plt

print ("Filename appears to be",sys.argv[1])
x=np.loadtxt(sys.argv[1], usecols = (2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34), delimiter=" ")

azmin=x[0,9]
qzmin=x[0,10]

azmax=x[-1,9]
qzmax=x[-1,10]

azindexmax=x[-1,14]
qzindexmax=x[-1,15]

print ("range of az:",azmin,azmax)
print ("range of qz:",qzmin,qzmax)
print ("index max of az:",azindexmax)
print ("index max of qz:",qzindexmax)

mesh=np.zeros((int(qzindexmax+1),int(azindexmax+1)))
xx=np.zeros((int(qzindexmax+1),int(azindexmax+1)))
yy=np.zeros((int(qzindexmax+1),int(azindexmax+1)))

for i in range(len(x)):
    azindex=int(x[i, 14])
    qzindex=int(x[i, 15])
    mesh[qzindex, azindex]=x[i, 16]
    xx[qzindex, azindex]=x[i, 9]
    yy[qzindex, azindex]=x[i, 10]

plt.pcolormesh(yy, xx, mesh)
plt.xlabel("q_z")
plt.ylabel("a_z")
ghz = f"%05.2f" % (x[-1,2]*1e-9)
filename = f"stability-%sGHz.png" % ghz
plt.title(f"Stability analysis for omega = %s GHz" % ghz)
plt.savefig(filename, bbox_inches='tight')
