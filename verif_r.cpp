//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
/*============================================================================*/
/*! \file verif_r
 *  \brief cylindical flux compressor for verification of coupling */
/*  ideal, athena++  */
/* compile with */
/*  ./configure.py -b --flux hlle --coord cylindrical --prob verif_r --cxx icc  */
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

//prototypes
void rout_bc(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

//float rng() { return ((float)rand()/(float)(RAND_MAX));}

//initial params of the fluid
Real mass,r0;
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

const double speed_of_light=3e10;
const double m2cm=1e2;
const double m2cm12=10;
const double J2erg=1e7;
const double curr_s2g=speed_of_light/m2cm12;//3e9
const double Ath2gauss=3.544907702;//sqrt(4*M_PI)
const double volt_s2g=m2cm12*J2erg/speed_of_light;//3.33e-3
const double eV2J=1.6e-19; // charge e in SI

Real r_max,r_min,z_max;


void Mesh::InitUserMeshData(ParameterInput *pin)
{
 mass=pin->GetReal("problem","mass");
 rho_min=pin->GetReal("problem","rho_min");
 E_min=pin->GetReal("problem","E_min");
 r0=pin->GetReal("problem","r0");
 dt_min=pin->GetReal("problem","dt_min");

 I0=pin->GetReal("problem","I0");//initial current in Amp
 I20=pin->GetReal("problem","I20");//inside current in Amp
 V0=pin->GetReal("problem","V0");// in Volts
 R1=pin->GetReal("problem","R1"); //in Ohm
 C1=pin->GetReal("problem","C1"); // in Farad
 L1=pin->GetReal("problem","L1"); // in Henry

 r_max = mesh_size.x1max;
 r_min = mesh_size.x1min;
 z_max = mesh_size.x3max;

     EnrollUserBoundaryFunction(OUTER_X1, rout_bc);
     EnrollUserTimeStepFunction(MyTimeStep);

 circ_output=fopen("circuit.txt","w");
 fprintf(circ_output,"# 1) time  2) V_C1  3) I  4) V_device 5) x 6) vx\n");
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
   gamma1 = peos->GetGamma() - 1.0;

q=V0*C1; //initialize circuit
tot_curr=I0;
Real Ld= 2e-7*z_max*log(r_max/r0)/m2cm;
printf("Ld=%e\n",Ld);
curr_rate=(tot_curr<1e3)? 0:(q/C1-tot_curr*R1)/(L1+Ld)/tot_curr; // 1/I dI/dt at t=0
voltage=curr_rate*tot_curr*Ld; 
printf("voltage=%e\n",voltage);

  return;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

//Real rho_fill=mass/(z_max*M_PI*r0*r0);
const Real a_coeff=-sqrt(4*M_PI/m2cm)/(2*M_PI);
Real dr=pcoord->dx1v(is); //uniform grid only!
Real dr_mass=(r_max-r_min)/800; // beres somewhat arbitrary check later
int i_mass2=ie-int((r_max-r0-dr_mass/2)/dr+0.5);
int i_mass1=ie-int((r_max-r0+dr_mass/2)/dr+0.5);
printf("n_mass=%d",i_mass2-i_mass1+1);
Real Vmass=z_max*M_PI*(SQR(pcoord->x1f(i_mass2+1))-SQR(pcoord->x1f(i_mass1)));
Real rho_mass=mass/Vmass;

for (int k=ks; k<=ke; k++)
for (int j=js; j<=je+1; j++)
for (int i=is; i<=ie+1; i++) 
  {
             Real r=pcoord->x1v(i);

             Real my_rho=rho_min;
             Real my_vr=0;
             Real my_bphi=0;
             Real E=E_min;

             if(i>i_mass2) 
              {
                my_vr   = -curr_rate*r*log(r/r0);  
                my_bphi =  tot_curr*a_coeff/r;
              }
             else if(i<i_mass1)
                   {
                     my_bphi =  I20*a_coeff/r;
                   }
                  else
                   {
                     my_rho=rho_mass;
                     E=0.25*(SQR(tot_curr)+SQR(I20))*SQR(a_coeff/r)/gamma1;
                   }       
             phydro->u(IDN,k,j,i) = my_rho;
             phydro->u(IM1,k,j,i) = my_rho*my_vr;
             phydro->u(IM2,k,j,i) = 0.0;
             phydro->u(IM3,k,j,i) = 0.0;
             phydro->u(IEN,k,j,i) = 0.5*my_rho*SQR(my_vr)+0.5*SQR(my_bphi)+E;
             pfield->b.x1f(k,j,i) = 0.0;
             pfield->b.x2f(k,j,i) = my_bphi;
             pfield->b.x3f(k,j,i) = 0.0;
            
  };

  return;
}

//calculate voltage on the device
//update the circuit
void MeshBlock::UserWorkInLoop(void)
{

//calculate voltage
Real volt=0.0;//inflow electric potential
int jm=(js+je)/2;
for(int k=ks; k<=ke; k++)
 {
   Real Vr1=phydro->u(IM1,k,jm,ie)/phydro->u(IDN,k,jm,ie);
   Real Vr2=phydro->u(IM1,k,jm,ie+1)/phydro->u(IDN,k,jm,ie+1);
   Real Bphi1=pfield->b.x2f(k,jm,ie);
   Real Bphi2=pfield->b.x2f(k,jm,ie+1);
   Real Ez=0.5*(Vr1*Bphi1+Vr2*Bphi2);
   volt+=Ez*pcoord->dx3v(k);
 }

voltage=volt*sqrt(4*M_PI)/(m2cm12*J2erg); //from Athena units to SI

Real position=0.0;
Real velocity=0.0;
Real mass=0.0;

for (int i=is; i<=ie; i++) 
{
 mass+=phydro->u(IDN,ks,js,i);
 position+=pcoord->x1v(i)*phydro->u(IDN,ks,js,i);
 velocity+=phydro->u(IM1,ks,js,i);          
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


void rout_bc(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

const Real a_coeff=-sqrt(4*M_PI/m2cm)/(2*M_PI);

//collect vr for all inflow zones
int nzones=0;
Real vr_total=0.0;
for(int k=ks; k<=ke; k++)
for(int j=js; j<=je; j++)
          {
           vr_total+=prim(IVX,k,j,ie);
           nzones++;
          };
if(nzones==0) throw std::runtime_error("error, total number of inflow zones ==0");
Real vr_avr=vr_total/nzones;

for(int k=ks; k<=ke; k++) //energy from generator, "almost inflow"
      {
          for(int j=js; j<=je; j++)
          for(int i=ie+1; i<=ie+ngh;  i++)
          {
           prim(IDN,k,j,i) = rho_min;
           double delta_r=(i-ie)*pcoord->dx1v(ie);
           double r=pcoord->x1v(i);
           prim(IVX,k,j,i) = vr_avr+(vr_avr/r-curr_rate)*delta_r;
           prim(IVY,k,j,i) = 0.0;
           prim(IVZ,k,j,i) = 0.0;
           prim(IPR,k,j,i) = E_min*gamma1;
          }

          for(int j=js; j<=je; j++)
          for(int i=ie+2; i<=ie+ngh+1; i++) b.x1f(k,j,i) = 0.0;
          for(int j=js; j<=je+1; j++)
          for(int i=ie+1; i<=ie+ngh;  i++) b.x2f(k,j,i) = tot_curr*a_coeff/pcoord->x1v(i);
          for(int j=js; j<=je; j++)
          for(int i=ie+1; i<=ie+ngh;  i++) b.x3f(k,j,i) = 0.0;

      }
}

