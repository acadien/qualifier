#!/usr/bin/python

import pylab as pl
import sys
#mine
from poscarsPlot import plotsimulation

data=open(sys.argv[1],"r").readlines()
if data[0][0]=="#":
    data.pop(0)
atoms=list()
msdNeighb=list()
msdIdeal=list()
energies=list()
volumes=list()
natoms=list()

basis=[[11.0,0,0],[0,11.0,0],[0,0,11.0]]
pos=0
while len(data)>1:
    natoms.append(int(data[pos+1]))
    energies.append(float(data[pos+3]))
    volumes.append(float(data[pos+5]))
    n=pos+9+natoms[-1]
    a=data[pos+9:n]
    atoms.append(map(lambda x:map(float,x.split()),a))
    msdNeighb.append(float(data[n+1]))
    msdIdeal.append(float(data[n+3]))
    data=data[n+5:]

optimE={38:-173.252,76:-402.384,104:-582.038}
nState=len(msdNeighb)
na=natoms[-1]

#Figure 1: MSD & Energy vs step #
pl.subplot(211)
pl.title("%d Atoms"%na)
pl.ylabel("Energy")
pl.plot([0,nState],[optimE[na]/na,optimE[na]/na],color="black",ls="--",label="Optimal Energy/atom")
pl.plot(range(nState),[i/na for i in energies],label="Energy/atom")
pl.subplot(212)
pl.plot(range(nState),msdIdeal,label="MSD-ideal")
pl.legend(loc=0)
pl.xlabel("Step")


#Figure 2: Visualize Lowest Energy Configuration
fig2=pl.figure()
#a2=fig2.add_subplot(111,project='3d')

iatoms=atoms[energies.index(min(energies))]
plotsimulation(basis,iatoms,[na],fig2)
pl.show()
