//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
/*============================================================================*/
/*! \file verif_z.cpp
 *  \brief cartesian flux compressor for verification of coupling */
/*  ideal, athena++  */
/* compile with */
/*  ./configure.py -b --flux hlle --prob verif_z --cxx icc  */
/*============================================================================*/

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include <stdlib.h>

void        z_drive(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

//initial params of the fluid
Real mass,z0;
Real rho_min,E_min,dt_min;
FILE* circ_output;

Real MyTimeStep(MeshBlock *pmb){ return dt_min; }


//circuit params
Real I0,I20;
Real V0;// V
Real R1;//Ohm
Real C1;//F
Real L1;//H

// dynamical variables of the circuit
Real tot_curr;// in Amp, not Gauss
Real curr_rate=0;// 1/I dI/dt to estimate dVz/dz
Real q; // charge on capacitor, coulomb
Real voltage;//on the device

//globals
Real gamma1;

//#define PRESSURE  

const double speed_of_light=3e10;
const double m2cm=1e2;
const double m2cm12=10;
const double J2erg=1e7;
const double curr_s2g=speed_of_light/m2cm12;//3e9
const double Ath2gauss=3.544907702;//sqrt(4*M_PI)
const double volt_s2g=m2cm12*J2erg/speed_of_light;//3.33e-3
const double eV2J=1.6e-19; // charge e in SI

Real x_max,x_min,y_max,y_min,z_max;

void Mesh::InitUserMeshData(ParameterInput *pin)
{
 mass=pin->GetReal("problem","mass");
 rho_min=pin->GetReal("problem","rho_min");
 E_min=pin->GetReal("problem","E_min");
 z0=pin->GetReal("problem","z0");
 dt_min=pin->GetReal("problem","dt_min");

 I0=pin->GetReal("problem","I0");//initial current in Amp
 I20=pin->GetReal("problem","I20");//inside current in Amp
 V0=pin->GetReal("problem","V0");// in Volts
 R1=pin->GetReal("problem","R1"); //in Ohm
 C1=pin->GetReal("problem","C1"); // in Farad
 L1=pin->GetReal("problem","L1"); // in Henry

 x_max = mesh_size.x1max;
 x_min = mesh_size.x1min;
 y_max = mesh_size.x2max;
 y_min = mesh_size.x2min;
 z_max = mesh_size.x3max;

 EnrollUserBoundaryFunction(INNER_X3, z_drive);
 EnrollUserTimeStepFunction(MyTimeStep);

 circ_output=fopen("circuit.txt","w");
 fprintf(circ_output,"# 1) time  2) V_C1  3) I  4) V_device 5) x 6) vx\n");

}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
   gamma1 = peos->GetGamma() - 1.0;

q=V0*C1; //initialize circuit
tot_curr=I0;
Real Ld= 4*M_PI*1e-7*z0*(x_max-x_min)/(y_max-y_min)/m2cm;
printf("Ld=%e\n",Ld);
curr_rate=(tot_curr<1e3)? 0:(q/C1-tot_curr*R1)/(L1+Ld)/tot_curr; // 1/I dI/dt at t=0
voltage=curr_rate*tot_curr*Ld; 
printf("voltage=%e\n",voltage);

  return;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

//initialize v and B 
Real a_coeff=-sqrt(4*M_PI/m2cm)/(2*M_PI);
Real dz=pcoord->dx3v(ks); //uniform grid only!
int k_mass=ks+int(z0/dz+0.5);
Real rho_mass=mass/dz/(x_max-x_min)/(y_max-y_min);
//initial velocity of mass=0;
Real B_phi1=sqrt(4*M_PI/m2cm)*tot_curr/(y_max-y_min); 
Real B_phi2=sqrt(4*M_PI/m2cm)*I20/(y_max-y_min); 

for(int k=ks; k<=ke+1; k++)
for(int j=js; j<=je+1; j++)
for(int i=is; i<=ie+1; i++)
   {
     Real rho=rho_min;
     Real vz=0;
     Real Bphi=0;
     Real z=pcoord->x3v(k);
     Real P=0;
     //if(k>=k_mass-1 && k<=k_mass+1) rho=rho_mass;
     if(k==k_mass) rho=rho_mass;
     if(k<k_mass) { vz=curr_rate*(z0-z); Bphi=B_phi1;}
     if (k>k_mass){
#ifndef PRESSURE         
         Bphi=B_phi2;
#else
         P=0.5*SQR(B_phi2);
#endif         
         }
      pfield->b.x1f(k,j,i) = 0;
      pfield->b.x2f(k,j,i) = Bphi;
      pfield->b.x3f(k,j,i) = 0;
      phydro->u(IDN,k,j,i)=rho;
      phydro->u(IM1,k,j,i)=0;
      phydro->u(IM2,k,j,i)=0;
      phydro->u(IM3,k,j,i)=rho*vz;
      phydro->u(IEN,k,j,i)=0.5*SQR(Bphi)+0.5*rho*SQR(vz)+E_min+P/gamma1;
   }
  return;
}

//calculate voltage on the device, flux, energy, etc
//update the circuit
void MeshBlock::UserWorkInLoop(void)
{
Real volt=0.0;//inflow electric potential
int jm=(js+je)/2;
for(int i=is; i<=ie; i++)
 {
   Real Vz1=0.5*(phydro->u(IM3,ks-1,jm-1,i)/phydro->u(IDN,ks-1,jm-1,i)+
                 phydro->u(IM3,ks-1,jm  ,i)/phydro->u(IDN,ks-1,jm  ,i));
   Real Vz2=0.5*(phydro->u(IM3,ks,jm-1,i)/phydro->u(IDN,ks,jm-1,i)+
                 phydro->u(IM3,ks,jm  ,i)/phydro->u(IDN,ks,jm  ,i));
   Real Bphi1=pfield->b.x2f(ks-1,jm,i);
   Real Bphi2=pfield->b.x2f(ks  ,jm,i);
   Real Ey=0.5*(Vz1*Bphi1+Vz2*Bphi2);
   volt+=Ey*pcoord->dx1v(i);
 }

voltage=volt*sqrt(4*M_PI)/(m2cm12*J2erg); //from Athena units to SI

Real position=0.0;
Real velocity=0.0;
Real mass=0.0;

for(int k=ks; k<=ke; k++)
{
 mass+=phydro->u(IDN,k,js,is);
 position+=pcoord->x3v(k)*phydro->u(IDN,k,js,is);
 velocity+=phydro->u(IM3,k,js,is);          
}

position/=mass;
velocity/=mass;

//update the circuit

q-=tot_curr*pmy_mesh->dt;
curr_rate=(tot_curr<1e3)? 0:(q/C1-tot_curr*R1-voltage)/L1/tot_curr; // 1/I dI/dt
tot_curr+=(q/C1-tot_curr*R1-voltage)*(pmy_mesh->dt)/L1;

fflush(circ_output);
// 1) time  2) V_C1  3) I  4) V_device 5) x 6) vx
fprintf(circ_output,"%.10e %.10e %.10e %.10e %.10e %.10e\n",
        pmy_mesh->time,q/C1,tot_curr,voltage,position,velocity);
  return;
}


void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
if(Globals::my_rank==0){
  fclose(circ_output);
}
  return;
}

void z_drive(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
   const Real a_coeff=-sqrt(4*M_PI/m2cm)/(2*M_PI);
   const int kl = ks - ngh;

//collect vz for all inflow zones
int nzones=0;
Real vz_total=0.0;

for(int i=is; i<=ie; i++)
for(int j=js; j<=je; j++)
          {
           vz_total+=prim(IVZ,ks,j,i);
           nzones++;
          };
if(nzones==0) throw std::runtime_error("error, total number of inflow zones ==0");
Real vz_avr=vz_total/nzones;

Real B_phi1=sqrt(4*M_PI/m2cm)*tot_curr/(y_max-y_min); 

for(int i=is; i<=ie; i++)  //energy from generator, "almost inflow"
for(int k=kl; k<ks;  k++)
for(int j=js; j<=je; j++)
   {
           prim(IDN,k,j,i) = rho_min;
           prim(IVX,k,j,i) = 0.0;
           prim(IVY,k,j,i) = 0.0;
           prim(IVZ,k,j,i) = vz_avr+curr_rate*(ks-k)*pcoord->dx3v(k);
           prim(IPR,k,j,i) = E_min*gamma1;
           b.x3f(k,j,i) = 0.0;
   }

for(int i=is; i<=ie+1; i++)  
for(int k=kl; k<ks;  k++)
for(int j=js; j<=je; j++)
   {
           b.x1f(k,j,i) = 0.0;
   }

for(int i=is; i<=ie; i++)  
for(int k=kl; k<ks;  k++)
for(int j=js; j<=je+1; j++)
   {
           b.x2f(k,j,i) = B_phi1;
   }

}
