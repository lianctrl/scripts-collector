
#include "header.h"

typedef char NAME[6];    /* Chr: a name; Lth: 50           */
#define NMAXRES   1000
#define NMAXFRAME 10000
//#define RCUT 2.4
#define DIM 3

void readtoendofline(FILE *);
void Usage(int argc, char **argv);

/*******************************************************************************

                              MAIN ROUTINE

*******************************************************************************/

int main (int argc, char **argv)
{
/*============================================================================*/
/*                LOCAL VARIABLES DECLARATION                                 */
/*============================================================================*/

    real       **x, x_ini[4], x_old[4], x_temp[4], t, dt, *msd;
    real       RCUT, xref[4];
    char       **FileName = argv, *File, *num;
    FILE       *fp, *fp2, *fp3;
    int        nframe_tot, nframe, nskip, num_traj, nstart;
    int        n, i, ii, d, step;
    real       hx, hy, hz, xorig, yorig, zorig, *pot, *pot2, Ez;
    real       **xelpot, xg[DIM+1], R, dR[DIM+1];
    int        nx, ny, nz, nxgrid, nygrid, nzgrid, ntotgrid, nelpot;

/*============================================================================*/

    if (argc != 4) Usage(argc, argv);

    if ((fp = fopen(FileName[1],"r"))==NULL) {
      fprintf(stdout,"Cannot open %s file!\n", FileName[1]);
      exit(1);
    }

    if ((fp2 = fopen(FileName[2],"w"))==NULL) {
      fprintf(stdout,"Cannot open %s file!\n", FileName[2]);
      exit(1);
    }

    if ((fp3 = fopen(FileName[3],"w"))==NULL) {
      fprintf(stdout,"Cannot open %s file!\n", FileName[3]);
      exit(1);
    }

    RCUT       = atof(argv[4]); //A°

/*
    Efield     = atof(argv[5]);

*/

/*==================================================================*/
/*  II)       Read input data  
   
    data from M373E grid space 0.5 A°
    object 1 class gridpositions counts 288 288 176
    origin 0.20929 0.20929 0.116493
    delta 0.463768 0 0
    delta 0 0.463768 0
    delta 0 0 0.482554

    grid space 0.25 A°
    object 1 class gridpositions counts 540 540 360
    origin 0.101082 0.101082 -0.00682449
    delta 0.247343 0 0
    delta 0 0.247343 0
    delta 0 0 0.235915

    ------------------------------------------------
    data from WT grid space 0.5 A°
    object 1 class gridpositions counts 288 288 168
    origin 0.25119 0.221191 0.141869
    delta 0.467588 0 0
    delta 0 0.467588 0
    delta 0 0 0.498181

    grid space 0.25 A°
    object 1 class gridpositions counts 540 540 336
    origin 0.142082 0.112083 0.0173264
    delta 0.24938 0 0
    delta 0 0.24938 0
    delta 0 0 0.249091

*/

    nxgrid = 540;
    nygrid = 540;
    nzgrid = 360;
    ntotgrid = nxgrid * nygrid * nzgrid;

    xorig = 0.101082;
    yorig = 0.101082;
    zorig = -0.0068245;

    hx = 0.247343;
    hy = 0.247343;
    hz = 0.235915;

    pot  = AllocVecR(ntotgrid);
    pot2 = AllocVecR(ntotgrid);

//  Electric Field unit in NAMD is kcal/mol/˚A/e
//  one unit is equivalent to 0.0434 Volts/˚A
//  Electric Potential unit in DX fiormat is kt/e
//  one unit is equivalent to 0.0258 Volts

//  Electric field from NAMD input (from command line)
/*  Ef for WT    = 0.2718

    Ef for M373E = 0.2721
*/


    Ez = 0.2721;           // Electric Field in kcal/mol/˚A/e

    fprintf(stdout,"\n");
    fprintf(stdout,"  Box x,y,z                  :  %12.8lf   %12.8lf  %12.8lf\n", hx * nxgrid, hy * nygrid, hz * nzgrid);
    fprintf(stdout,"  Box center x,y,z           :  %12.8lf   %12.8lf  %12.8lf\n", hx * nxgrid / 2.0 + xorig, hy * nygrid / 2.0 + yorig, hz * nzgrid / 2.0 + zorig);
    fprintf(stdout,"  Box orig (down left) x,y,z :  %12.8lf   %12.8lf  %12.8lf\n", xorig, yorig, zorig);


    fprintf(stdout,"\n");
    fprintf(stdout,"  Electric Field applied     :  %12.8lf kcal/mol/˚A/e  %12.8lf Volts/˚A\n", Ez, 0.0434 * Ez);
    fprintf(stdout,"  Electric Potential applied :  %12.8lf Volts\n", 0.0434 * Ez * hz * nzgrid);
    fflush(stdout);

    Ez = Ez * 0.0434;    // Electric Field in Volts/˚A
    Ez = Ez / 0.0258;    // Electric Field in kt/e/˚A


    /* Add external electric potential to electric potential given by the system only as provided by VMD .dx */
    for (nx = 0; nx < nxgrid; nx++) {
       for (ny = 0; ny < nygrid; ny++) {
          for (nz = 0; nz < nzgrid; nz++) {

             n = nx * nygrid * nzgrid + ny * nzgrid + nz + 1;

             fscanf(fp,"%lf", &pot[n]);

             pot2[n] = pot[n] - nz * hz * Ez;

          }
       }
    }

//  Put the center of Calpha of the residues of the selectivity filter
//  oc_cac_dim.py script has been used

/*  Wild Type select filter center
    ('X:', 66.52740)
    ('Y:', 67.90098)
    ('Z:', 50.62985)

    M373E select filter center
    ('X:', 67.70304)
    ('Y:', 66.54276)
    ('Z:', 50.37272)
*/ 

    xref[1] = 67.70304;
    xref[2] = 66.54276;
    xref[3] = 50.37272;

    for (nx = 0; nx < nxgrid; nx++) {
       for (ny = 0; ny < nygrid; ny++) {
          for (nz = 0; nz < nzgrid; nz++) {

             n = nx * nygrid * nzgrid + ny * nzgrid + nz + 1;

             xg[1] = nx * hx + xorig;
             xg[2] = ny * hy + yorig;
             xg[3] = nz * hz + zorig;

             for (d = 1; d <= DIM; d++) {
                 dR[d] = xref[d] - xg[d];
             }
             R = sqrt(dR[1]*dR[1] + dR[2]*dR[2]);


             /* Write coordinate and electric field in Volt/Ang of the old and corrected one*/
             if (R < RCUT) {
               fprintf(fp2," %12.5lf   %12.5lf   %12.5lf       %12.5lf  %12.5lf\n", xg[1], xg[2], xg[3], pot[n] * 0.0258, pot2[n] * 0.0258);
               fflush(fp2);
             }

          }
       }
    }


    /* Write the new electric potential map that now includes the external potential */
    for (n = 1; n <= ntotgrid; n=n+3) {

       fprintf(fp3,"%.5lf %.5lf %.5lf\n", pot2[n], pot2[n+1], pot2[n+2]);
    }

 exit(1);


}
/*******************************************************************************

                              ROUTINE: Usage

*******************************************************************************/

void Usage (int argc, char **argv)
{
     fprintf(stdout, "\n     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     fprintf(stdout,   "     @                                                @\n");
     fprintf(stdout,   "     @               JJ_ELECPOTENTIAL                 @\n");
     fprintf(stdout,   "     @                                                @\n");
     fprintf(stdout,   "     @                                                @\n");
     fprintf(stdout,   "     @               by  G. Brancato                  @\n");
     fprintf(stdout,   "     @                                                @\n");
     fprintf(stdout,   "     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n");
     fprintf(stdout,   "  Usage:                                               \n");
     fprintf(stdout,   "  jj_elecpot_3D_to_2D old.dx pot_coord.xyz new.dx RCUT \n\n");
     fflush(stdout);
     exit(1);
}

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* readtoendofline: Function to read to end of line in read_coord files     */
/*==========================================================================*/
void readtoendofline(FILE *fp){
  int eol,ch;
  eol = (int )'\n';
  ch = eol+1;
  while(ch!=eol&&ch!=EOF){ch=fgetc(fp);}
  if(ch==EOF){
      printf("\n@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("ERROR: End of file reached                     \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/
}/* end routine */
/*==========================================================================*/
