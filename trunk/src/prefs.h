/*
*				preflist.h
*
* Include file for preflist.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 1999-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	Stuff is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	Stuff is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with Stuff. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		03/03/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _GALAXIES_H_
#include "galaxies.h"
#endif
#ifndef _STARS_H_
#include "stars.h"
#endif
#ifndef _SED_H_
#include "sed.h"
#endif

typedef enum {DETECT_PHOTONS, DETECT_ENERGY}	detenum;

/*------------------------------- preferences -------------------------------*/
typedef struct
  {
  char		*(cat_name[SED_MAXNPB]);/* Filename(s) of output catalog */
  int		ncat_name;		/* Number of catalogs */
  int		width[SED_MAXNPB];	/* Image width */
  int		nwidth;			/* Number of width entries */
  int		height[SED_MAXNPB];	/* Image height */
  int		nheight;	       	/* Number of height entries */
  double	pixscale[SED_MAXNPB];	/* Pixel scales (arcsec) */
  int		npixscale;	       	/* Number of pixel size entries */
  double	maglim[2];	/* Magnitude limits of the simulation */
/* Cosmology */
  double	h0;		/* Hubble constant (km.Mpc-1) */
  double	omegam;		/* Omega "matiere" */
  double	omegal;		/* Omega "lambda" */
  double	integ_zstep;	/* Integration step along z (h-1.Mpc) */
/* Passbands */
  char		refpb_name[MAXCHAR];	/* Reference passband specs */
  detenum	refdetect_type;		/* How is the ref. flux integrated */
  char		*(pb_name[SED_MAXNPB]);		/* Observed passband specs */
  int		npb_name;		/* Number of bandpasses */
  detenum	obsdetect_type[SED_MAXNPB];
					/* How is the obs. flux integrated */
  int		nobsdetect_type;	/* Number of obsdetect_types */
  double	gain[SED_MAXNPB];	/* Conversion factor in e-/ADU */
  int		ngain;			/* Number of gains */
  double	area[SED_MAXNPB];	/* Collecting area for computing ZPs */
  int		narea;			/* Number of areas */
/* SEDs */
  char		refcalibsed_name[MAXCHAR];	/* Ref. system SED filename */
  char		*(pbcalibsed_name[SED_MAXNPB]);	/* Obs. system SED filenames */
  int		npbcalibsed_name;		/* Number of bandpasses */
  char		backsed_name[MAXCHAR];		/* Sky backgr. SED filename */
  char		*(gal_sedname[GAL_MAXNSED]);	/* Galaxy SED filenames */
  int		ngal_sedname;			/* Number of galaxy SEDs */
  char		*(star_sedname[STAR_MAXNSED]);	/* Star SED filenames */
  int		nstar_sedname;			/* Number of star SEDs */
  int		gal_bsedno[GAL_MAXNTYPE];	/* SED-indices for bulges */
  int		ngal_bsedno;			/* Number of indices */
  int		gal_dsedno[GAL_MAXNTYPE];	/* SED-indices for disks */
  int		ngal_dsedno;			/* Number of indices */
/* Extinction */
  char		extinct_name[MAXCHAR];		/* Extinction Al/Av filename */
  double	gal_iextinc[GAL_MAXNTYPE];	/* External extinc. */
  int		ngal_iextinc;			/* Number of alpha's */
/* Luminosity functions */
  double	lf_phistar[GAL_MAXNTYPE]; /* Schechter's phi* (h-3.Mpc-3) */
  int		nlf_phistar;		  /* Number of phi*'s */
  double	lf_phistarevol[GAL_MAXNTYPE];	  /* P = dln(phi*) / dln(1+z)*/
  int		nlf_phistarevol;	  /* Number of P's */
  double	lf_mstar[GAL_MAXNTYPE];	  /* Schechter's M* (mag+5.log10(h)) */
  int		nlf_mstar;		  /* Number of M*'s */
  double	lf_mstarevol[GAL_MAXNTYPE];	  /* Q = dM* / dln(1+z) */
  int		nlf_mstarevol;		  /* Number of Q's */
  double	lf_alpha[GAL_MAXNTYPE];	  /* Schechter's alpha */
  int		nlf_alpha;		  /* Number of alpha's */
  double	lf_mlim[2];	/* Bounds to the LF (mag+5.log10(h) */
/* Morphology */
  double	gal_hubtype[GAL_MAXNTYPE];/* Hubble type */
  int		ngal_hubtype;		  /* Number of Hubble types */
  double	gal_bt[GAL_MAXNTYPE];	  /* B/T ratios */
  int		ngal_bt;		  /* Number of B/T's */
  int		gal_ntype;		  /* Total number of galaxy types */
  double	gal_dbeta;		  /* Tully-Fisher's beta */
  double	gal_drstar;		  /* Disk effective radius at M* */
  double	gal_drevol;		  /* Disk radius evolution parameter */
  double	gal_dsigmalambda;	  /* Dispersion of disk radii */
  double	gal_bmknee;		  /* Binggeli's abs. magnitude at knee*/
  double	gal_brknee;		  /* Binggeli's eff. radius at knee */
  double	gal_brevol;		  /* Bulge radius evolution parameter */
/* Galaxy clusters */
  char		clusterlist_name[MAXCHAR];/* Cluster list filename (input) */
  char		clusterlistout_name[MAXCHAR];/*Cluster list filename (output)*/
/* Lensing */
  double	lens_kappa;		  /* Convergence lensing parameter */
  double	lens_gamma[2];		  /* Shear lensing parameters */
  int		nlens_gamma;		  /* Number of shear parameters */
/* Stars */
  int		starflag;	/* Include stars? */
  double	galcoord[2];	/* Galactic coordinates */
  int		starseed;	/* Seed for the random generator of stars */
  int		galseed;	/* Seed for the random generator of galaxies */
/* Multithreading */
  int		nthreads;			/* Number of active threads */
/* miscellaneous */
  char		datadir_name[MAXCHAR];	/* Global path for SEDs, filters,... */
  enum {QUIET, NORM, LOG, FULL} verbose_type;	/* display type */
  char		sdate_start[12];		/* SCAMP start date */
  char		stime_start[12];		/* SCAMP start time */
  char		sdate_end[12];			/* SCAMP end date */
  char		stime_end[12];			/* SCAMP end time */
  double	time_diff;			/* Execution time */
  }	prefstruct;

prefstruct	prefs;

/*----------------------------- Internal constants --------------------------*/

#define		MAXCHARL	16384   /* max. nb of chars in a string list */
#define		MAXLIST		64	/* max. nb of list members */
#define		MAXLISTSIZE	2000000	/* max size of list */

/*-------------------------------- protos -----------------------------------*/
extern char	*list_to_str(char *listname);

extern int	cistrcmp(char *cs, char *ct, int mode);

extern void	dumpprefs(int state),
		endprefs(void),
		readprefs(char *filename,char **argkey,char **argval,int narg),
		useprefs(void);

