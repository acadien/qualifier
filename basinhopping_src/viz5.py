#!/usr/bin/python

import pylab as pl
import sys
import matplotlib.gridspec as gridspec
#mine
from poscarsPlot import plotsimulation
from struct_tools import neighbors

def usage():
    print "%s"%sys.argv[0].split("/")[-1]

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
#logAnames={{}}
lognames={}
for k in [38,76,104]:
    lognames[(k,"A")]=["../basinhopping_data/finalstate_N"+str(k)+"_A_"+str(i)+".dat" for i in range(5)]
    lognames[(k,"B")]=["../basinhopping_data/finalstate_N"+str(k)+"_B_"+str(i)+".dat" for i in range(5)]
gs=[gridspec.GridSpec(1, 1),gridspec.GridSpec(1, 1)]
gs[0].update(left=0.1, right=0.52)
gs[1].update(left=0.54, right=0.96)
a=[]
a.append(pl.subplot(gs[0][:,:]))
a.append(pl.subplot(gs[1][:,:]))
width=0.01
for k in [38,76,104]:
    for index,h in enumerate(["A","B"]):
        pl.sca(a[index])
        energies,volumes,atoms,msdN,msdI=zip(*map(loadData,lognames[(k,h)]))
        #mean consecutive local minima distance
        meanCLMD=[sum(i)/len(i) for i in msdN] 

        #minimum Root-Mean-Square distance to Ideal case
        minRMSI=map(min,msdI)

        avgMCLMD=sum(meanCLMD)/len(meanCLMD)

        minMRMSI=min(minRMSI)
        maxMRMSI=max(minRMSI)
        avgMRMSI=sum(minRMSI)/len(minRMSI)
        
        print avgMCLMD,avgMRMSI
        
        pl.bar(avgMCLMD-width/2., maxMRMSI-minMRMSI, width=width, bottom=minMRMSI,lw=2,color="none")
        pl.text(avgMCLMD-width/2.,maxMRMSI*1.01,"N="+str(k))
        pl.scatter(avgMCLMD,avgMRMSI,marker="x",c="black",s=50,linewidths=2,zorder=10)
        pl.plot([avgMCLMD-width/2,avgMCLMD+width/2],[avgMRMSI,avgMRMSI],c="black",lw=2)


for i,h in enumerate(["A","B"]):
    pl.sca(a[i])
    pl.axis([0,0.15,0,3.5])    
    pl.xlabel("$\overline{\Delta RMS}_{Neighbor}$",size=14)
    if i==0:
        pl.xticks(pl.xticks()[0][:-1])
        pl.ylabel("$min(RMS_{ideal})$",size=14)
    if i==1:
        pl.yticks(pl.yticks()[0],[""]*len(pl.yticks()[0]))
    pl.title("Method %s"%h)

pl.show()
