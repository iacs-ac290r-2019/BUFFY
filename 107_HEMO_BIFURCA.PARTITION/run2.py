#!/usr/bin/env python

from MagicUniverse import *
import sys

# preprocess_parallel_mesh(4); sys.exit(1)

NU = 0.1
RHO = 1
NSTEP = 1000
NDIAG = 100
NVTKFREQ = 100
U_BAR = 0.01
Re = 5.0
Pe = 10
C0 = 0.01
C1 = 1.0
UNFREEZE_TIME = 500

R = Re * NU / 2 / U_BAR
DIFFUSIVITY = U_BAR * R / Pe

MagicBegins()

# define universe

u = Universe()
s = Scale()
m = Mesh()
f = Fluid()
c = Fluid()
t = Tracker()

u.addItems([s,m,f, c,t])

u.setTitle('Test_Re=5_Steps=100')
u.setNumberOfSteps(NSTEP)
u.setTemperature(0.0)
u.setStateRestart(False)
u.setStateDumpFrequency(-1)

u.create()

# set params

s.set(name='MonoScale', mesh=m,  actors=[f,c,t])

m.setRegularMesh(True)

m.setPeriodicity('111')
m.setDomainDecomposition(3)


t.setDiagnosticFrequency(NDIAG)
# t.setDataShow(velocity=True)
# t.setMapDirections('zx')

t.setVtkDump(True, meshtype='unstructured', frequency=NVTKFREQ)

f.setName('BloodFlow')
f.setCollisionType('BGK')
f.setDensityUniform(RHO)
f.setViscosity(NU)

# NEWLY ADDED
c.setName('Bolus')
c.setCollisionType('BGK')
c.setDiffusivity(DIFFUSIVITY)
c.setAdvector(f)
c.setADR(True)

f.setInletOutletMethod('zouhe')
c.setInletOutletMethod('zouhe')
# f.setInletOutletMethod('equilibrium')
# f.setInletOutletMethod('closed')

f.setStabilizeLB(True)

u.decorate()
f.setIOValue('inlet', 1, U_BAR)
f.setIOValue('outlet', 2, U_BAR)
c.setIOValue('inlet', 1, C0)
c.setIOValue('outlet', 2, C0)

nx, ny, nz = m.getBox()
nx = int(nx); ny = int(ny); nz = int(nz)
profile2 = c.getArray(nx*ny*nz)

for k in xrange(1,nz+1):
	for j in xrange(1,ny+1):
		for i in xrange(1,nx+1):
			ifl = m.getLocator(i,j,k)
			
			if k > nz/8 - 3 and k < nz/8. + 3: # create BOLUS
				profile2[ifl] = C1
			else:
				profile2[ifl] = C0

c.setDensityProfile(profile2)
c.setFreeze(True)

for itime in u.cycle():
	if itime == UNFREEZE_TIME:
		c.setFreeze(False)			

	u.animate()

MagicEnds()
