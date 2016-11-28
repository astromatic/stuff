/*
*				lf.h
*
* Include file for lf.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 1999-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		26/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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
  double	dmstar_faceon;	/* Change in M* due to extinction */
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
