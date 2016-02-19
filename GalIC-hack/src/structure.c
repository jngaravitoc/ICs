/*******************************************************************************
 * This file is part of the GALIC code developed by D. Yurin and V. Springel.
 *
 * Copyright (c) 2014
 * Denis Yurin (denis.yurin@h-its.org) 
 * Volker Springel (volker.springel@h-its.org)
 *******************************************************************************/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

#include "allvars.h"
#include "proto.h"



static double fc(double c)
{
  return c * (0.5 - 0.5 / pow(1 + c, 2) - log(1 + c) / (1 + c)) / pow(log(1 + c) - c / (1 + c), 2);
}

static double jdisk_int(double x, void *param)
{
  double vc2, Sigma0, vc, y;

  if(x > 1.0e-10 * All.Halo_A)
    vc2 = All.G * (halo_get_mass_inside_radius(x) + bulge_get_mass_inside_radius(x)) / x;
  else
    vc2 = 0;

  if(vc2 < 0)
    terminate("vc2 < 0");

  Sigma0 = All.Disk_Mass / (2 * M_PI * All.Disk_H * All.Disk_H);
  y = x / (2 * All.Disk_H);

  if(y > 1e-4)
    vc2 +=
      x * 2 * M_PI * All.G * Sigma0 * y * (gsl_sf_bessel_I0(y) * gsl_sf_bessel_K0(y) -
					   gsl_sf_bessel_I1(y) * gsl_sf_bessel_K1(y));

  vc = sqrt(vc2);

  return pow(x / All.Disk_H, 2) * vc * exp(-x / All.Disk_H);
}


static double gc_int(double x, void *param)
{
  return pow(log(1 + x) - x / (1 + x), 0.5) * pow(x, 1.5) / pow(1 + x, 2);
}




void structure_determination(void)
{
  double jhalo, jdisk, jd;
  double hnew, dh, f_c, mhmvir;


  /* total galaxy mass */
  All.Mvir = pow(All.Vvir, 3) / (6.9713700231733506 * All.G * All.Hubble);
  /* virial radius of galaxy */
  All.Rvir = All.Vvir / (6.9713700231733506 * All.Hubble);

  //fprintf(stdout, "Mvir I %f, Rvir I %f and the constants G %f and H %f\n", All.Mvir, All.Rvir, All.G, All.Hubble);
  /* halo scale radius */
  All.Halo_Rs = All.Rvir / All.Halo_C;
	  
  /* set the scale factor of the hernquist halo */
  f_c = log(1.0 + All.Halo_C) - All.Halo_C / (1.0 + All.Halo_C); 
  All.Halo_A = All.Halo_Rs * 1.0 / (1.0/sqrt(2 * f_c) -  1.0 / All.Halo_C);
	  
  /* set the Hernquist equivalent mass */
  mhmvir = pow(All.Halo_A / All.Halo_Rs, 2) / (2 * f_c);
  All.Mvir = All.Mvir * mhmvir;

 
  All.LowerDispLimit = pow(0.01 * All.Vvir, 2);


  /* determine the masses of all components */
  All.Disk_Mass = All.MD * All.Mvir;
  All.Bulge_Mass = All.MB * All.Mvir;
  
  /* Halo mass */
  All.Halo_Mass = All.Mvir;
  mpi_printf("Halo Mass %f and A parameter %f, disk mass %f and Bulge Mass %f  \n", All.Halo_Mass, All.Halo_A, All.Disk_Mass, All.Bulge_Mass);


  All.BH_Mass = All.MBH;
  if(All.MBH > 0)
    All.BH_N = 1;
  else
    All.BH_N = 0;
  
  
  jhalo = All.Lambda * sqrt(All.G) * pow(All.Mvir, 1.5) * sqrt(2 * All.Rvir / fc(All.Halo_C));
  jdisk = All.JD * jhalo;

  double halo_spinfactor =
    1.5 * All.Lambda * sqrt(2 * All.Halo_C / fc(All.Halo_C)) * pow(log(1 + All.Halo_C) -
								   All.Halo_C / (1 + All.Halo_C),
								   1.5) / structure_gc(All.Halo_C);

  mpi_printf("\nStructural parameters:\n");
  mpi_printf("Rvir            = %g\n", All.Rvir);
  mpi_printf("Mvir            = %g  (this is the total mass of the equivalent Hernquist)\n", All.Mvir);
  mpi_printf("A (halo)        = %g\n", All.Halo_A);
  mpi_printf("halo_spinfactor = %g\n", halo_spinfactor);

  /* first guess for disk scale length */
  All.Disk_H = sqrt(2.0) / 2.0 * All.Lambda / fc(All.Halo_C) * All.Rvir; 
  All.Disk_Z0 = All.DiskHeight * All.Disk_H;	/* sets disk thickness */ 

  All.Bulge_A = All.BulgeSize;	/* this will be used if no disk is present */

  MType[1] = All.Halo_Mass;
  MType[2] = All.Disk_Mass;
  MType[3] = All.Bulge_Mass;

  NType[1] = All.Halo_N;
  NType[2] = All.Disk_N;
  NType[3] = All.Bulge_N;


  if(All.Disk_Mass > 0)
    {
      do
	{
	  jd = structure_disk_angmomentum();	/* computes disk momentum */

	  hnew = jdisk / jd * All.Disk_H;

	  dh = hnew - All.Disk_H;

	  if(fabs(dh) > 0.5 * All.Disk_H)
	    dh = 0.5 * All.Disk_H * dh / fabs(dh);
	  else
	    dh = dh * 0.1;

	  All.Disk_H = All.Disk_H + dh;

	  /* mpi_printf("Jd/J=%g   hnew: %g  \n", jd / jhalo, All.Disk_H);
	   */

	  All.Disk_Z0 = All.DiskHeight * All.Disk_H;	/* sets disk thickness */
	}
      while(fabs(dh) / All.Disk_H > 1e-5);
    }

  mpi_printf("H  (disk)       = %g\n", All.Disk_H);
  mpi_printf("Z0 (disk)       = %g\n", All.Disk_Z0);
  mpi_printf("A (bulge)       = %g\n", All.Bulge_A);
}


double structure_disk_angmomentum(void)
{
  gsl_function F;
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &jdisk_int;

  double result, abserr;

  gsl_integration_qag(&F, 0, dmin(30 * All.Disk_H, All.Rvir),
		      0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  result *= All.Disk_Mass;

  gsl_integration_workspace_free(workspace);

  return result;
}


double structure_gc(double c)
{
  gsl_function F;
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &gc_int;

  double result, abserr;

  gsl_integration_qag(&F, 0, c, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return result;
}
