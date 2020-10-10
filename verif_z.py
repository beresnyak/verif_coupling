#!/usr/bin/env python

Lz = 0.779 # m
#Lz = 0.769 #m
Lxy = 1 #Lx/Ly
#Lz= 76.9 #cm
# mass = 8e-5 g
mass = 8e-8 #kg

#z0=1.0 cm
z0=0.01 #m
#z0= 0.0
dt= 1e-10
t_end=2e-6
n_steps=int(t_end/dt)

V0 = 6e5    # in V
I0 = 2e5    # initial current in Amp
I20 = 2e5   #
R1 = 0.0    # Ohm
C1 = 1e-6   # Farad
L1 = 7e-7   # Henry
mu0=1.256637e-6

print "# (1)-t (2)-z in cm (3)-I (4)-Vd (5)-Vc (6)-vz in cm/s"

# initialize

t=0
vz=0
z=z0
I=I0
q=V0*C1

for i in range(n_steps):
    dqdt=-I
    dzdt=vz
    I2=I20*(Lz-z0)/(Lz-z)
    dvzdt=mu0/2/mass*(I*I-I2*I2)*Lxy
    didt=(q/C1-I*R1-I*mu0*Lxy*vz)/(L1+mu0*z*Lxy)
    Vd=q/C1-I*R1-didt*L1
    if(i%100==0): print "%.9e %.9e %.9e %.9e %.9e %.9e"%(t,z*100,I,Vd,q/C1,vz*100)
    q+=dt*dqdt
    I+=dt*didt
    z+=dt*dzdt
    vz+=dt*dvzdt
    t+=dt
    

    