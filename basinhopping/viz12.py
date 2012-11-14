#!/usr/bin/python

import pylab as pl
import sys
#mine
from poscarsPlot import plotsimulation
from struct_tools import neighbors

def usage():
    print "%s <natom>"%sys.argv[0].split("/")[-1]

def loadData(fname):
    data=open(fname,"r").readlines()
    if data[0][0]=="#":
        data.pop(0)
    natoms=int(data[1])
    es=list()
    vs=list()
    ats=list()
    mN=list()
    mI=list()
    pos=14+natoms
    for i,line in enumerate(data):
        j=i%pos
        if(j==3):
            es.append(float(line))
        if(j==5):
            vs.append(float(line))
        if(j==9):
            ats.append(list())
        if(j>=9 and j<9+natoms):
            ats[-1].append(map(float,line.split()))
        if(j==natoms+10):
            mN.append(float(line))
        if(j==natoms+12):
            mI.append(float(line))
    return es,vs,ats,mN,mI

def printMins(energies,msdI,ti):
    dexE=map(lambda x:x.index(min(x)),energies)
    dexM=map(lambda x:x.index(min(x)),msdI)

    print "Mins by Energy: "+ti
    print dexE
    print [energies[i][j] for i,j in enumerate(dexE)]
    print [msdI[i][j] for i,j in enumerate(dexE)]

    print "Mins by MSD(I): "+ti
    print [energies[i][j] for i,j in enumerate(dexM)]
    print [msdI[i][j] for i,j in enumerate(dexM)]


basis=[[11.0/2,0,0],[0,11.0/2,0],[0,0,11.0/2]]
optimE={38:-173.928,76:-402.894,104:-582.086}

nState=5000

nAtom=int(sys.argv[1])
logAnames=["data/finalstate_N"+str(nAtom)+"_A_"+str(i)+".dat" for i in range(5)]
logBnames=["data/finalstate_N"+str(nAtom)+"_B_"+str(i)+".dat" for i in range(5)]

energiesA,volumesA,atomsA,msdNA,msdIA=zip(*map(loadData,logAnames))
printMins(energiesA,msdIA,"A")
energiesB,volumesB,atomsB,msdNB,msdIB=zip(*map(loadData,logBnames))
printMins(energiesB,msdIB,"B")

"""
pl.subplot(211)
for x in energiesB:
    pl.plot(x)
pl.subplot(212)
for x in msdIB:
    pl.plot(x)
pl.show()
"""
"""
nState=len(msdNeighb)

#Figure 1: MSD & Energy vs step #
pl.subplot(211)
pl.title("%d Atoms"%natoms)
pl.ylabel("Energy")
pl.plot([0,nState],[optimE[natoms]/natoms,optimE[natoms]/natoms],color="black",ls="--",label="Optimal Energy/atom")
pl.plot(range(nState),[i/natoms for i in energies],label="Energy/atom")
pl.subplot(212)
pl.plot(range(nState),msdIdeal,label="MSD-ideal")
pl.legend(loc=0)
pl.xlabel("Step")


#Figure 2: Visualize Lowest Energy Configuration
fig2=pl.figure()
#a2=fig2.add_subplot(111,project='3d')
mc=energies.index(min(energies))
iatoms=atoms[mc]
ienergy=energies[mc]
print ienergy
aa=plotsimulation(basis,iatoms,[natoms],fig2)
bnd=[[0.0,11.0],[0.0,11.0],[0.0,11.0]]
halfNeigh=neighbors(iatoms,bnd,1.5,"half")
for a in range(natoms):
    for b in halfNeigh[a]:
        x,y,z=zip(*[iatoms[a],iatoms[b]])
        aa.plot(x,y,z,c="black",lw=1.5)
pl.show()
"""
