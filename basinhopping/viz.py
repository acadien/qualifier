#!/usr/bin/python

import pylab as pl
import sys

data=open(sys.argv[1],"r").readlines()
data.pop(0) #drop the header
atoms=list()
msds=list()
energies=list()
volumes=list()
natoms=list()

bounds=[[11.0,0,0],[0,11.0,0],[0,0,11.0]]
pos=0
while len(data)>1:
    natoms.append(int(data[pos+1]))
    energies.append(float(data[pos+3]))
    volumes.append(float(data[pos+5]))
    n=pos+9+natoms[-1]
    a=data[pos+9:n]
    atoms.append(map(lambda x:map(float,x.split()),a))
    msds.append(float(data[n+1]))
    data=data[n+2:]

nState=len(msds)
nAtom=natoms[-1]
pl.plot(range(nState),msds,label="MSD")
pl.plot(range(nState),[i/nAtom for i in energies],label="Energy")
pl.show()
