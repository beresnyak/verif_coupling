//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
/*============================================================================*/
/*! \file verif_magnoh
 *  \brief magnoh for verification of coupling */
/*  ideal, athena++  */
/* compile with */
/*  ./configure.py -b --flux hlle --coord cylindrical --prob verif_magnoh --cxx icc  */
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
#include <stdio.h>

//prototypes
void rout_bc(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

//initial params of the fluid
Real rho_min,E_min,dt_min;
static Real alpha, beta, rho0, P0, pcoeff;
static Real bphi0, bz;

//globals
Real gamma1;

Real r_max,r_min,z_max;
Real rt;
std::string sss_file;
AthenaArray<Real> sss_table;
int sss_n_entries;
Real sss_rmin,sss_rmax;
Real sss_logrmax,sss_logstep;

const double sss_time=3e-8;
const double vshock=1e7;

void Mesh::InitUserMeshData(ParameterInput *pin)
{
// for coupling 
 rho_min=pin->GetReal("problem","rho_min");
 E_min=pin->GetReal("problem","E_min");
// dt_min=pin->GetReal("problem","dt_min");
// for magnoh 
 alpha  = pin->GetReal("problem", "alpha");
 beta  = pin->GetReal("problem", "beta");
 pcoeff= pin->GetReal("problem", "pcoeff");
 rho0  = pin->GetReal("problem", "d");
 bphi0 = pin->GetReal("problem", "bphi")/sqrt(4*M_PI);//conversion to Athena units
 bz = pin->GetReal("problem", "bz")/sqrt(4*M_PI);//conversion to Athena units
 P0 = 4*M_PI*pcoeff*(bphi0*bphi0+bz*bz);
 sss_file= pin->GetString("problem", "file");

 r_max = mesh_size.x1max;
 r_min = mesh_size.x1min;
 z_max = mesh_size.x3max;

 rt=r_max;
 
     EnrollUserBoundaryFunction(OUTER_X1, rout_bc);
}


double sss_value(double r, double t, int var)
{
    if(var<1||var>3){
     fprintf(stderr,"invalid variable in sss_value\n");
     exit(1);   
    }
    if(t==0||sss_time*r/t>=sss_rmax) 
    { 
      switch(var)
      {  
         case 1: return 0;
         case 2: return bphi0*pow(r,beta);
         case 3: return bz*pow(r,beta);
      }      
    }    
    else{
        double scale_factor=1.0;
        if(var>1) scale_factor=pow(t/sss_time,beta)/sqrt(4*M_PI); // note conversion to Athena units
        double ind=(log(sss_time*r/t)-sss_logrmax)/sss_logstep; //always positive
        int i=(int)ind; double di=ind-i;
        if(i>=sss_n_entries-1) return 0.0;
        else return scale_factor*(sss_table(i,var)+di*(sss_table(i+1,var)-sss_table(i,var))); //log linear interp
    }   
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  gamma1 = peos->GetGamma() - 1.0;
  return;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  AthenaArray<Real> az;
  az.NewAthenaArray(nx2,nx1);
  for (int j=js; j<=je+1; ++j) 
    for (int i=is; i<=ie+1; ++i)       
      az(j,i) = (bphi0/(beta+1))*pow(pcoord->x1f(i),beta+1);
    
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
    pfield->b.x1f(k,j,i) = (az(j+1,i) - az(j,i))/pcoord->dx2f(j)/pcoord->x1f(i);
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x2f(k,j,i) = -(az(j,i) - az(j,i+1))/pcoord->dx1f(i);
  }}}
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x3f(k,j,i) = bz*pow(pcoord->x1v(i),beta);
  }}}
  az.DeleteAthenaArray();

  // initialize conserved variables
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
             // Volume centered coordinates
             Real  rho=rho0*pow(pcoord->x1v(i),alpha);
             Real  P  =  P0*pow(pcoord->x1v(i),2*beta);
             phydro->u(IDN,k,j,i) = rho;
             phydro->u(IM1,k,j,i) = 0.0;
             phydro->u(IM2,k,j,i) = 0.0;
             phydro->u(IM3,k,j,i) = 0.0;
             phydro->u(IEN,k,j,i) = P/(peos->GetGamma()-1.0);
          if (NON_BAROTROPIC_EOS && MAGNETIC_FIELDS_ENABLED) {
             phydro->u(IEN,k,j,i)+=
                0.5*0.25*(SQR(pfield->b.x1f(k,j,i)+pfield->b.x1f(k,j,i+1))
                        + SQR(pfield->b.x2f(k,j,i)+pfield->b.x2f(k,j+1,i))
                        + SQR(pfield->b.x3f(k,j,i)+pfield->b.x3f(k+1,j,i)));
          } // NON_BAROTROPIC_EOS && MAGNETIC_FIELDS_ENABLED
  }}}; //main loop
  
//read sss_file
 FILE *s_file;
 char line[200];
 if ((s_file=fopen (sss_file.c_str(),"r"))==NULL)
 { 
  fprintf(stderr,"Unable to open self-similar solution: %s\n",sss_file.c_str());
  exit(1);
 }
 int num_lines=0;
 int dummy1; double dummy2,dummy3,dummy4,dummy5,dummy6,dummy7,dummy8;
 while (fgets(line, 200, s_file)) {
 if (*line == '#') continue;
 if (sscanf(line,"%d%lf%lf%lf%lf%lf%lf%lf",
     &dummy1,&dummy2,&dummy3,&dummy4,&dummy5,&dummy6,&dummy7,&dummy8)!= 8) {
       fprintf(stderr,"invalid format in %s\n",sss_file.c_str()); 
       exit(1);
     } else {
         num_lines++;
     }
 }
 sss_table.NewAthenaArray(num_lines,4);// r,vr,bphi,bz
 fseek(s_file, 0, SEEK_SET );
 int counter=0;
 while (fgets(line, 200, s_file) && counter<num_lines) {
 if (*line == '#') continue;
 if (sscanf(line,"%d%lf%lf%lf%lf%lf%lf%lf",
     &dummy1,&sss_table(counter,0),&dummy3,&sss_table(counter,1),
     &sss_table(counter,2),&sss_table(counter,3),&dummy7,&dummy8)!= 8) {
       fprintf(stderr,"invalid format in %s\n",sss_file.c_str()); 
       exit(1);
     } else {
         counter++;
     }
 }
 sss_n_entries=counter;
 fprintf(stderr,"Success in reading self similar solution table from %s with %d entries\n",
                 sss_file.c_str(),sss_n_entries);
 sss_rmin=sss_table(sss_n_entries-1,0);
 sss_rmax=sss_table(0,0);
 sss_logrmax=log(sss_rmax);
 sss_logstep=(log(sss_rmin)-sss_logrmax)/(sss_n_entries-1); // negative
 fclose(s_file);

 fprintf(stderr,"Testing self similar solution table\n");
 fprintf(stderr,"t=%.4e r=%.4e v=%.4e bphi=%.4e bz=%.4e\n",
                 0.0,0.1,sss_value(0.1,0,1),sss_value(0.1,0,2),sss_value(0.1,0,3));
 fprintf(stderr,"t=%.4e r=%.4e v=%.4e bphi=%.4e bz=%.4e\n",
                 3e-8,2.5e4,sss_value(2.5e4,3e-8,1),sss_value(2.5e4,3e-8,2),sss_value(2.5e4,3e-8,3));
 fprintf(stderr,"t=%.4e r=%.4e v=%.4e bphi=%.4e bz=%.4e\n",
                 3e-8,3.01e-1,sss_value(3.01e-1,3e-8,1),sss_value(3.01e-1,3e-8,2),sss_value(3.01e-1,3e-8,3));
  return;
}


void MeshBlock::UserWorkInLoop(void)
{
rt+=sss_value(rt,pmy_mesh->time,1)*pmy_mesh->dt;
if(pmy_mesh->ncycle%50==0) printf("rt=%e\n",rt);
  return;
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  return;
}


void rout_bc(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{

//self-similar values of Bphi and Bz;
    
Real Bphi_rt=sss_value(rt,pmb->pmy_mesh->time,2);
Real Bz_rt=sss_value(rt,pmb->pmy_mesh->time,3);
    
for(int k=ks; k<=ke; k++) // inflow with vr and rho_min
      {
          for(int j=js; j<=je; j++)
          for(int i=ie+1; i<=ie+ngh;  i++)
          {
           prim(IDN,k,j,i) = rho_min;
           prim(IVX,k,j,i) = prim(IVX,k,j,ie);
           prim(IVY,k,j,i) = 0.0;
           prim(IVZ,k,j,i) = 0.0;
           prim(IPR,k,j,i) = E_min*gamma1;
          }

          for(int j=js; j<=je; j++)
          for(int i=ie+2; i<=ie+ngh+1; i++) b.x1f(k,j,i) = 0.0;
          for(int j=js; j<=je+1; j++)
          for(int i=ie+1; i<=ie+ngh;  i++) b.x2f(k,j,i) = Bphi_rt*rt/pmb->pcoord->x1v(i); // vacuum field
      }      
      
          for(int k=ks; k<=ke+1; k++)
          for(int j=js; j<=je; j++)
          for(int i=ie+1; i<=ie+ngh;  i++) b.x3f(k,j,i) = Bz_rt;

      
}

