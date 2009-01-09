/*
 				lf.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	Stuff
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for lf.c
*
*	Last modify:	11/09/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef	_LF_H_
#define	_LF_H_

/*--------------------------------- constants -------------------------------*/

#define	MABS_MIN	(-20.0)	/* Minimum absolute magnitude - M* */
#define	MABS_MAX	20.0	/* Maximum absolute magnitude - M* */
#define	MABS_NSTEP	1000	/* Number of integration steps */

/*-------------------------- structure definitions --------------------------*/

typedef struct
  {
  double	phistar;	/* Schechter's phi* parameter (m-3) */
  double	mstar;		/* Schechter's M* parameter (mag) */
  double	dmstar;		/* Change in M* due to z or other factor */
  double	alpha;		/* Schechter's alpha parameter */
  double	mabsmin,mabsmax;/* Bounds of the LF */
  double	mstarevol;	/* M* evolution factor dm* / ln(1+z) */
  double	phistarevol;	/* phi* evolution factor dln(phi*) / dln(1+z)*/
  } lfstruct;


/*-------------------------------- protos -----------------------------------*/

lfstruct	*lf_evol(lfstruct *lf, double z, lfstruct *lfevol);

double	lf_intlumschechter(lfstruct *lf, double mmin, double mmax),
	lf_intschechter(lfstruct *lf, double mmin, double mmax),
	lf_rndschechter(lfstruct *lf, double mmin, double mmax),
	lf_schechter(lfstruct *lf, double m);

#endif
