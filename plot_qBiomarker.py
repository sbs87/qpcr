import cutsom_utilities
import sys
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import color_maps as cmc

file= cutsom_utilities.readfile("qPCR_dCt_alltimes")

cmap3=cmc.discrete_cmap(4)

def read_agreement():
	dates=file[0][1:]#Cols
	ndates=len(dates)
	species=[]
	X=np.zeros([0,ndates])
	j=1 #Counter for only non zero counted species
	for i in range(len(file)):
		if(file[i][0]=="Species"):
			next
		else:
			catval=0
			for val in file[i][1:]:
			    catval=catval+double(val)
			if catval==2:
			    next
			else:
			    species.append(file[i][0])
			    X=np.insert(X,j-1,[file[i][1:]],axis=0)
			    j=j+1
	return dates, species, X

dates, species, X=read_agreement()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(X, cmap=cm.get_cmap("afmhot"), interpolation='nearest',aspect='auto')

numrows, numcols = X.shape
numdates=len(dates)

numspecies=len(species)
xticks(range(numdates),dates,rotation=0)
yticks(range(numspecies),species,rotation=0)

def format_coord(x, y):
    col = int(x+0.5)
    row = int(y+0.5)
    if col>=0 and col<numcols and row>=0 and row<numrows:
        z = X[row,col]
        return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, z)
    else:
        return 'x=%1.4f, y=%1.4f'%(x, y)

ax.format_coord = format_coord
plt.show()



#--------------------------------------------------------------
Y=([1,20,21],[40,3,8])
#Species1[Method1,Method2]
print Y
plot(Y)

show()