#usr/bin/env/python

## This program is for a 2D Monte Carlo simulation of two species (peptides
## and PEG) on a surface, with various interaction energies.

import time,sys, math, gc
from math import log, sqrt, pow,exp
import random, numpy
from random import shuffle
traj=open("traj.gro",'w')
NHBout=open("NHB.xvg",'w')
NPEGout=open("NPEG.xvg",'w')
prep =open("percentrepl.xvg",'w')

dim = 80 # nm
area = dim*dim
apl = 0.64  # N^2/pep
Nlipid = area/apl  ###
#Rpep = 1.05  ##nm
#Rwig = 0.5
#Rhard = 2*Rpep-Rwig
#########       INPUT  		#############

x1 = 0.05  ### number of peptides per pepid
x2 = 0.05   ### number of PEGs per pepid
			
Npep = int(Nlipid*x1)
NPeg = int(Nlipid*x2)
Nbet = Npep
N = Npep + Nbet  + NPeg 
NU = Npep + NPeg
type = dict([x,1] for x in range(Npep))
type.update(dict([x,2] for x in range(Npep,Npep+NPeg)))
type.update(dict([x,3] for x in range(Npep+NPeg,N)))
pi = math.acos(-1)
#######     INITIAL COORDINATES   ########3
Arad = 0.7 ### 1.5 for low pH, 1.2 for high pH
RepFac = -0.4 ## 0.5 : repulsion at 0.5nm (low pH), 0.25: repulsion at 0.75nm (high pH)
### 0: fusion peptide
### 1-4: LE amino acid pair

def checkdim(z):
	if (z>dim): zz = z - dim
	elif (z<0): zz = z + dim
	else: zz = z
	return zz
def checkd(z):
	if (z>dim/2): zz = dim - z
	else: zz = z
	return zz

x = [0]*N
y = [0]*N 
for i in range(NU):
	x[i] = random.random()*dim
	y[i] = random.random()*dim
for j in range(NU,N):
	ang = random.random()*2*pi
	rad = random.random()*Arad
	x[j] = checkdim(x[j-NU] + math.cos(ang)*rad)
	y[j] = checkdim(y[j-NU] + math.sin(ang)*rad)


########    INTERACTION POTENTIALS   ########

####### PEG-PEG

Prange = 5
A1 = 192.1
A2 = 3.906
B1 = 1.05
B2 = 3.62
eps = 1
sigll = 0.4

######### pep-pep

Erange = 3.5 
pm2 = 6.26
pm1 = -33.53
p0 = 74.41
p1 = -64.92
p2 = 20.76
p3 = -2.18
Fpp = 0.8

######## pep-PEG

Pm2 = 1.47
Pm1 = -6.15
P0 = 15.21
P1 = -7.67
P2 = 1.085

FPp = 0.6 * Fpp  
Brange = 1.12
VBrange = 1.2
#HBstrength = 6.72
HBstrength = 4.57
def Bpot(r):
	energy = 300*(r-1)*(r-1) - HBstrength
	return energy

def Vpot(r):
	energy = 1/pow(r+RepFac,12) + Fpp*(pm2/(r*r) + pm1/r + p0 + p1*r + p2*r*r + p3*r*r*r)
	return energy

def VBpot(r):
	energy = 1/pow(r,12)
	return energy

def Ppot(r):
	energy = A1*exp(-r/(B1*sigll))+A2*exp(-r/(B2*sigll))+4*eps*pow((sigll/r),12)
	return energy

def pPpot(r):
	energy = 1/pow(r+RepFac,12) + FPp*(pm2/(r*r) + pm1/r + p0 + p1*r + p2*r*r + p3*r*r*r)
	#energy = 0
	return energy

######   DEFINE PEPTIDE CLUSTERS  #######

Crange = 5

def cluster(x,y):
	tc=time.time()
	ppc = [[0 for k in range(Npep/2)] for i in range(Npep)]
	z = [1 for k in range(Npep)]
	for pep in range(Npep):
		ppc[pep][0]=pep
		for pep2 in range(Npep):
			dx = abs(x[pep]-x[pep2])
			dy = abs(y[pep]-y[pep2])
			if dx>dim/2: dx=dim-dx
			if dy>dim/2: dy=dim-dy
			r = math.sqrt(dx*dx+dy*dy)
			if r<Crange and pep!=pep2:
				gc.disable()
				ppc[pep][z[pep]]=pep2
				z[pep]+=1
				gc.enable()
	ppc = [ppc[pep][:z[pep]] for pep in range(Npep)] 
	i = 0
	g = []
	while i<len(ppc):
		g.append(ppc[i]) 
		j=i+1
		while j<len(ppc):
			if list(set(ppc[i]) & set(ppc[j]))!=[]:
				gc.disable()
				g[i] = list(set(g[i] + ppc[j]))
				del(ppc[j])
				gc.enable()
				j-=1
			j+=1
		i+=1
	
	return g


g = cluster(x,y)
#print g
######   CREATE NEIGHBOR LISTS   ######

Nrange = 30
Nupdate = 15 #Nupdate < Nrange - Erange

def neighbors(x,y):
	nlist = [[] for i in range(NU)]
	for i in range(NU):
		for j in range(NU):
			dx = abs(x[i]-x[j])
			dy = abs(y[i]-y[j])
			if dx>dim/2: dx=dim-dx
			if dy>dim/2: dy=dim-dy
			r = math.sqrt(dx*dx+dy*dy)
			if r<Nrange and i!=j:
				gc.disable()
				nlist[i].append(j)
				gc.enable()
	#t2=time.time()
	#print "neighbor time",t2-tn1
	return nlist

nlist = neighbors(x,y)

def PEGincluster(x,y):

	nP = 0
	for i in range(Npep,NPEG+Npep):
		n = 0
		for j in nlist[i]:
			if j < Npep:
				dx = checkd(abs(x[i]-x[j])
				dy = checkd(abs(y[i]-y[j])
				if math.sqrt(dx*dx+dy*dy)<2.5:
					n += 1
		if n > 4:
		nP += 1
	
	nfrac = float(nP)/float(NPeg)
	return nfrac


def tenergy(i,j,ilist):	

	tx = x[i]; ty = y[i]
	if i>=Npep:
		dx = checkd(abs(x[j+NU]-x[j]))
		dy = checkd(abs(y[j+NU]-y[j]))
		xb = checkdim(x[i]+dx)
		yb = checkdim(y[i]+dy)
	else: 
		xb = x[i+NU]
		yb = y[i+NU]

	Ep = 0
	EP = 0
	for j in ilist:
		if i!=j:
			dxold = checkd(abs(tx-x[j]))  #### pep-pep, pep-PEG
			dyold = checkd(abs(ty-y[j]))
			rold = math.sqrt(dxold*dxold + dyold*dyold)

			if rold < Erange:
				if j<Npep:

					Ep += Vpot(rold)      #### pep-pep
					dx = checkd(abs(xb-x[j+NU])) #### beta-beta
					dy = checkd(abs(yb-y[j+NU]))
					dB = math.sqrt(dx*dx+dy*dy)
					if dB < Brange: Ep += Bpot(dB)

					dx = checkd(abs(xb-x[j]))  ###pep-beta
					dy = checkd(abs(yb-y[j]))
					dVB = math.sqrt(dx*dx+dy*dy)
					if dVB < VBrange: Ep += VBpot(dVB)
				
					dx = checkd(abs(tx-x[j+NU]))  ###pep-beta
					dy = checkd(abs(ty-y[j+NU]))
					dVB = math.sqrt(dx*dx+dy*dy)
					if dVB < VBrange: Ep += VBpot(dVB)
				
					EP += pPpot(rold)    ####  PEG -pep
					dx = checkd(abs(tx-x[j+NU])) ### PEG-beta
					dy = checkd(abs(ty-y[j+NU]))
					dVB = math.sqrt(dx*dx+dy*dy)
					if dVB < VBrange: EP += VBpot(dVB)
						
				else:
					Ep += pPpot(rold) ##### pep-PEG
					EP += Ppot(rold)

					dx = checkd(abs(xb-x[j]))
					dy = checkd(abs(yb-y[j]))
					dVB = math.sqrt(dx*dx+dy*dy)
					if dVB < VBrange: Ep += VBpot(dVB)

					

	return Ep, EP

#####   MONTE CARLO   ######

steps=500000
delta = 2
NHB = [0]*Npep

for s in range(steps):
	if s%Nupdate==0:
		nlist = neighbors(x,y)
		#g = cluster(x,y)

	# REPLACEMENT MOVES

	#i = int(random.random()*Npep)
	#j = int(random.random()*NPeg)
	#cEp, nEP = tenergy(i,j,nlist[i])
	#nEp, cEP = tenergy(j,i,nlist[j])

	#if (nEp + nEP) < (cEp + cEP) or exp(-((nEp + nEP)-(cEp + cEP))) > random.random():
					#print (nEp+nEP),(cEp+cEP)
					
	#	dxb = checkdim(x[i+NU] - x[i])
	#	dyb = checkdim(y[i+NU] - y[i])
	#	x[i+NU] = checkdim(x[j] + dxb)
	#	y[i+NU] = checkdim(y[j] + dyb)
	#	tx = x[i]
	#	ty = y[i]
	#	x[i] = x[j]
	#	y[i] = y[j]
	#	x[j] = tx
	#	y[j] = ty
#

    ####### Single Moves
	Rlist = [int(z) for z in range(NU)]
	random.shuffle(Rlist)
	NHBtotal = 0
	for i in Rlist:

		Et = 0
		nx = 0
		ny = 0
		nxb = 0
		nyb = 0
		r = random.random()*delta

		theta = random.random()*2*pi
		dx = math.cos(theta)*r
		dy = math.sin(theta)*r

		nx = checkdim(x[i] + dx)
		ny = checkdim(y[i] + dy)

		if i < Npep: 
			rad = random.random()*Arad
			theta = random.random()*2*pi
			dx = math.cos(theta)*rad
			dy = math.sin(theta)*rad

			nxb = checkdim(nx + dx)
			nyb = checkdim(ny + dy)
	
		nE = 0
		cE = 0
		NHBold = [0,0]
		NHBnew = [0,0]
		tlist = nlist[i]
		random.shuffle(tlist)
		for j in tlist:

			dxold = checkd(abs(x[i]-x[j]))  #### pep-pep, pep-PEG
			dyold = checkd(abs(y[i]-y[j]))
			dxnew = checkd(abs(nx-x[j]))
			dynew = checkd(abs(ny-y[j]))
			rold = math.sqrt(dxold*dxold + dyold*dyold)
			rnew = math.sqrt(dxnew*dxnew + dynew*dynew)

			if i < Npep:
				if j < Npep:
					if rold < Erange:
						cE += Vpot(rold)
						dx = checkd(abs(x[i+NU]-x[j+NU])) #### beta-beta
						dy = checkd(abs(y[i+NU]-y[j+NU]))
						dB = math.sqrt(dx*dx+dy*dy)
						if dB < Brange: 
							tE = Bpot(dB)
							if tE < max(NHBold):
								NHBold[NHBold.index(max(NHBold))] = tE
								cE += tE

						dx = checkd(abs(x[i+NU]-x[j]))
						dy = checkd(abs(y[i+NU]-y[j]))
						dVB = math.sqrt(dx*dx+dy*dy)
						if dVB < VBrange: cE += VBpot(dVB)

						dx = checkd(abs(x[i]-x[j+NU]))
						dy = checkd(abs(y[i]-y[j+NU]))
						dVB = math.sqrt(dx*dx+dy*dy)
						if dVB < VBrange: cE += VBpot(dVB)

					if rnew < Erange:
						nE += Vpot(rnew)
						dx = checkd(abs(nxb-x[j+NU])) #### beta-beta
						dy = checkd(abs(nyb-y[j+NU]))
						dB = math.sqrt(dx*dx+dy*dy)
						if dB < Brange: 
							tE = Bpot(dB)
							if tE < max(NHBold):
								NHBnew[NHBnew.index(max(NHBnew))] = tE
								cE += tE

						dx = checkd(abs(nxb-x[j]))
						dy = checkd(abs(nyb-y[j]))
						dVB = math.sqrt(dx*dx+dy*dy)
						if dVB < VBrange: nE += VBpot(dVB)

						dx = checkd(abs(nx-x[j+NU]))
						dy = checkd(abs(ny-y[j+NU]))
						dVB = math.sqrt(dx*dx+dy*dy)
						if dVB < VBrange: nE += VBpot(dVB)


				if j > Npep:
					if rold < Erange:
						cE += pPpot(rold)
						dx = checkd(abs(x[i+NU] - x[j]))
						dy = checkd(abs(y[i+NU] - y[j]))
						dVB = math.sqrt(dx*dx+dy*dy)
						if dVB < VBrange: cE += VBpot(dVB)

					if rnew < Erange:
						nE += pPpot(rnew)
						dx = checkd(abs(nxb - x[j]))
						dy = checkd(abs(nyb - y[j]))
						dVB = math.sqrt(dx*dx+dy*dy)
						if dVB < VBrange: nE += VBpot(dVB)

			elif i >= Npep:
				if j < Npep:
					if rold < Erange:
						cE += pPpot(rold)
						dx = checkd(abs(x[i] - x[j+NU]))
						dy = checkd(abs(y[i] - y[j+NU]))
						dVB = math.sqrt(dx*dx+dy*dy)
						if dVB < VBrange: cE += VBpot(dVB)

					if rnew < Erange:
						nE += pPpot(rnew)
						dx = checkd(abs(nx - x[j+NU]))
						dy = checkd(abs(nx - y[j+NU]))
						dVB = math.sqrt(dx*dx+dy*dy)
						if dVB < VBrange: nE += VBpot(dVB)


				elif j >= Npep:
					if rold < Prange:
						cE += Ppot(rold)
					if rnew < Prange:
						nE += Ppot(rnew)


		if nE < cE or exp(-(nE-cE))>random.random():
			x[i] = nx
			y[i] = ny
			if i < Npep :
				x[i+NU] = nxb
				y[i+NU] = nyb
			L = len([k for k in NHBnew if k!=0])
			NHBtotal += L
		else: 
			L = len([k for k in NHBold if k!=0])
			NHBtotal += L



##### Print to files

	if (s%5==0):
		print>>NHBout, float(NHBtotal)/float(Npep)
		PEGfrac = PEGincluster(x,y)
		print>>NPEGout,  PEGfrac

 	if (s%5==0):				
		#print "s",s
#		print "single move time",t4-t3
#		print "grp move time",t5-t4
#	a	print "cluster time",t2-t1
		print>>traj,"conf"
		print>>traj,N
		for z in range(Npep):
			print >> traj,"   1XXX   H   ",z*2+1,"     ",x[z],"    ",y[z]+dim+10,"   0"
			print >> traj,"   1XXX   O   ",z*2+2,"     ",x[z+NU],"    ",y[z+NU]+dim+10,"   0"
		for z in range(Npep,Npep+NPeg):
			print >> traj,"   1XXX   C   ",z+Npep+1,"     ",x[z],"    ",y[z]+dim+10,"   0"
							
		print>>traj,dim,"  ",dim,"  ",1

		#NgLo=open("gL.xvg",'a')
		#print>>NgLo, log(T) , len(g)
		
		#print >> sizeo,"conf",Nmin
		#print >> sizeo,int(N/8)-1
		#print >> sizeo,"   1XXX   O   ",1,"     ",1,"    ",N1-sum(NP),"   0"
		#for i in range(2,int(N/8)):
	#		print >> sizeo,"   1XXX   O   ",i,"     ",i,"    ",NP.count(i),"   0"
	#	print>>sizeo,dim,"  ",dim,"  ",1

	#	Navg = float(N1)/len(g)
	#	Navgo=open("Navg.xvg",'a')
	#	print>>Navgo, log(T) , Navg
		t6 =time.time()
		#print "print time",t6-t5


