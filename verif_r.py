#!/usr/bin/env python
from math import log

Lz = 0.6 # m
#Lz = 0.2 # m
#r_max=0.048 # m
#r_min=0.012 # m
r_max=0.20 # m
r_min=0.08 # m
mass = 8e-8 # kg

r0=0.18 # m
#r0=0.04 # m
dt= 2e-12
t_end=2e-6
n_steps=int(t_end/dt)

V0 = 6e5    # in V
I0 = 2e5    # initial current in Amp
I20 = 2e5   #
R1 = 0.0    # Ohm
C1 = 1e-6   # Farad
L1 = 7e-7   # Henry
mu0=1.256637e-6
mu0_4pi=1e-7
mu0_2pi=2e-7

print "# (1)-t (2)-z in cm (3)-I (4)-Vd (5)-Vc (6)-vz in cm/s"

# initialize

t=0
vr=0
r=r0
I=I0
q=V0*C1

for i in range(n_steps):
    dqdt=-I
    drdt=vr
    I2=I20*log(r0/r_min)/log(r/r_min)
    dvrdt=mu0_4pi/mass*(-I*I+I2*I2)*Lz/r
    didt=(q/C1-I*R1+I*mu0_2pi*Lz*vr/r)/(L1+mu0_2pi*Lz*log(r_max/r))
    Vd=q/C1-I*R1-didt*L1
    if(i%1000==0): print "%.9e %.9e %.9e %.9e %.9e %.9e"%(t,r*100,I,Vd,q/C1,vr*100)
    q+=dt*dqdt
    I+=dt*didt
    r+=dt*drdt
    vr+=dt*dvrdt
    t+=dt
    

    