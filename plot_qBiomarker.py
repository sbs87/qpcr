import bmath
import sys
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

file= bmath.readfile("newinfile.txt")
print file[4][0]
def read_agreement():
	dates=file[0][1:]#["W3D2","W7D4","W8D3","W9D4"]
	species=[]
	for i in range(len(file)):
		if(file[i][0]=="Species"):
			next
		else:
			species.append(file[i][0])
	#species=["L.crispatus","L.iners","G.vaginalis","P.acnes","L.jensenii"]
	X=np.array([(1,2,1,1),(1,2,1,1),(1,3,2,1),(1,1,1,3),(1,1,1,3)])
	#Species1(T1,T2,T3...Tn), Species2(T1,T2,T3...Tn)
	return dates, species, X


dates, species, X=read_agreement()


fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(X, cmap=cm.jet, interpolation='nearest')

numrows, numcols = X.shape
numdates=len(dates)
numspecies=len(species)
xticks(range(numdates),dates,rotation=0)
yticks(range(numspecies),species,rotation=0)

#def format_coord(x, y):
#    col = int(x+0.5)
#    row = int(y+0.5)
#    if col>=0 and col<numcols and row>=0 and row<numrows:
#        z = X[row,col]
#        return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, z)
#    else:
#        return 'x=%1.4f, y=%1.4f'%(x, y)

#ax.format_coord = format_coord
plt.show()