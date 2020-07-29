/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "bond_harmonic_fibers.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondHarmonicFibers::BondHarmonicFibers(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondHarmonicFibers::~BondHarmonicFibers()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
  }
}

/* ---------------------------------------------------------------------- */

void BondHarmonicFibers::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr,rk,rer,xFac;
  double E1,E2,R1,R2,m1,m2,n1,n2,eps1,eps2,ebreak;
  //Parameters in micro units 
  E1=3.5e3;
  //E2=0.77e9;
  eps1=0.92;
  //eps2=0.08;
  ebreak=0.97;
  R1=3.20e5; //valore della forza a eps1 dello stress-strain
  //R2=0.3311e9; //valore della forza a eps2 dello stress-strain
  //m1=E1-2*R1/(eps1-eps2);
  //m2=E2-2*R2/(eps2-eps1);
  //n1=R1-(E1-2*R1/(eps1-eps2))*eps1;
  //n2=R2-(E2-2*R2/(eps2-eps1))*eps2;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    //read each bond and compute its length
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx,dely,delz);

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    //r0 of the type of bond: normal or crosslink
    dr = r - r0[type];
    //define the engeneering strain
    rer = dr/r0[type];
    //Form factor derived from a sort of fermi-dirac behavior of the fiber to smooth the breaking point 
    //coefficient 500 choose upon what criteria?
    xFac=1/(exp((rer-ebreak)*500)+1);

    if (rer<=eps1) {
	//linear elastic response of the force in the first part of the curve
    fbond=E1*rer;}
    /*if (rer>eps1 && rer<=eps2) {
	//response of the force fitted from the exerimental data with opportune parameters 
    fbond=(m1*rer+n1)*((rer-eps2)/(eps1-eps2))*((rer-eps2)/(eps1-eps2))+(m2*rer+n2)*((rer-eps1)/(eps2-eps1))*((rer-eps1)/(eps2-eps1));}*/
    if (rer>eps1) { //fbond=0;
	//linear elastic response of the force in the last part of the curve
    fbond=R1+E1*(rer-eps1);}
    //k of the type of bond normal or crosslink as crossection area. 
    //here is the right k lammps code usaul took 2k as parameter.
    fbond=fbond*k[type]*xFac;

    // computing the bond force between bonded atoms
 
    fbond = -fbond/r;

    //computing energy between bonded atoms 
    if (eflag) ebond = 0.5*fbond*dr;
 
    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fbond;
      f[i1][1] += dely*fbond;
      f[i1][2] += delz*fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fbond;
      f[i2][1] -= dely*fbond;
      f[i2][2] -= delz*fbond;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* -----------------------------------------------------------------------------
 
	here took r0[type] and k[type] from input parameters in lammps command 

   -----------------------------------------------------------------------------*/

void BondHarmonicFibers::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(r0,n+1,"bond:r0");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondHarmonicFibers::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double r0_one = force->numeric(FLERR,arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length 
------------------------------------------------------------------------- */

double BondHarmonicFibers::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void BondHarmonicFibers::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void BondHarmonicFibers::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondHarmonicFibers::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],r0[i]);
}

/* ---------------------------------------------------------------------- */


double BondHarmonicFibers::single(int type, double rsq, int i, int j, double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double rk = k[type] * dr;
  fforce = 0;
  if (r > 0.0) fforce = -2.0*rk/r;
  return rk*dr;
}
/* ----------------------------------------------------------------------
    Return ptr to internal members upon request.
------------------------------------------------------------------------
void *BondHarmonicFibers::extract( char *str, int &dim )
{
  dim = 1;
  if( strcmp(str,"kappa")==0) return (void*) k;
  if( strcmp(str,"r0")==0) return (void*) r0;
  return NULL;
}*/
