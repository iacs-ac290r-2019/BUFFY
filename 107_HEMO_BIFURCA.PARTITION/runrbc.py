#!/usr/bin/env python

from MagicUniverse import *
import sys

# preprocess_parallel_mesh(4); sys.exit(1)

NU = 0.1
RHO = 1
NSTEP = 5000
NDIAG = 500
NVTKFREQ = 500
U_BAR = 0.01
Re = 5.0
Pe = 10
C0 = 0.01
C1 = 1.0
UNFREEZE_TIME = 3000
GROWTIME = 1000 #1000
RESTART = False

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

m.setPeriodicity('100')
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
myid = get_myproc()

if not RESTART:
    profile2 = c.getArray(nx*ny*nz)
    for k in range(1,nz+1):
        for j in range(1,ny+1):
            for i in range(1,nx+1):

                ifl = m.getLocator(i,j,k)

                if i > nx/8 - 3 and i < nx/8. + 3: # create BOLUS
                    profile2[ifl] = C1
                else:
                    profile2[ifl] = C0

    c.setDensityProfile(profile2)

    # freeze the fluid and drug until released
    f.setFreeze(True)
    c.setFreeze(True)

    # cap RBC forces to a large roof to avoid instabilities
    a.setCapForces(True, forcecap=1.e4, torquecap=1.e4, velcap=0.4, angvelcap=0.6)
    # set a robust friction coefficient for initial equilibration of RBC
    a.setGamma(gammaT=0.1, gammaR=0.1)
    a.setZeroVelocity()

for itime in u.cycle():

  if not RESTART:

    # gently increase excluded volume
    if itime <= GROWTIME:
        tscale = 0.5 * (1. + itime/GROWTIME) 
        if myid==0 and itime%10==0: 
            print 'Rescaling interaction to', tscale
            sys.stdout.flush()
        a.scaleVdwParameters(tscale)

    #  increase plasma-RBC coupling
    elif itime == int(1.5*GROWTIME):
        f.setFreeze(False)
        a.setGamma(gammaT=0.001, gammaR=0.001)
        a.setZeroVelocity()

    # RBC free to move with right coupling
    elif itime == int(2.*GROWTIME):
        a.setGamma(gammaT=0.01, gammaR=0.01)
        a.setZeroVelocity()

    # allow the drug to move 
    elif itime == int(2.5*GROWTIME):
        c.setFreeze(False)

  # BC: stripe of given drug density
  c.setDensityStripe('x', nx, C0) 

  u.animate()

MagicEnds()