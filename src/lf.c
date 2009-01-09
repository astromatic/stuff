/*
                                  lf.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        Stuff
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       Functions dealing with luminosity functions.
*
*       Last modify:    11/09/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "globals.h"
#include "cosmo.h"
#include "lf.h"
#include "random.h"


/****************************** lf_intschechter ******************************/
/*
The Schechter (1976) galaxy luminosity function integral from mmin to mmax.
*/
double	lf_intschechter(lfstruct *lf, double mmin, double mmax)
  {
   static double	phi[MABS_NSTEP], alphas;
   static int		flag;
   double		dm, di, phimin,phimax, moffset;
   int			i, over;

  dm = (double)(MABS_MAX-MABS_MIN)/MABS_NSTEP;

  if (!flag || lf->alpha!=alphas)
    {
/*-- Never passed through here before! Let's initialize arrays */
    moffset = MABS_MIN + lf->mstar + 0.5*dm;
    for (i=0; i<MABS_NSTEP; i++)
      phi[i] = lf_schechter(lf, i*dm + moffset) + (i?phi[i-1]:0.0);
    flag = 1;
    alphas = lf->alpha;
    }

/* Now make the interpolation */
/* Lower bound */
  mmin -= lf->mstar;
  if (mmin < MABS_MIN)
    mmin = MABS_MIN;
  di = (mmin-MABS_MIN)/dm-0.4999;
  i = (int)(di+0.5);
  di -= (double)i;
  if (i>=MABS_NSTEP-1)
    {
    over = i+2-MABS_NSTEP;
    i -= over;
    di += (double)over;
    }
  phimin = phi[i]+(phi[i+1]-phi[i])*di;

/* Upper bound */
  mmax -= lf->mstar;
  if (mmax < MABS_MIN)
    mmax = MABS_MIN;
  di = (mmax-MABS_MIN)/dm-0.4999;
  i = (int)(di+0.5);
  di -= (double)i;
  if (i>=MABS_NSTEP-1)
    {
    over = i+2-MABS_NSTEP;
    i -= over;
    di += (double)over;
    }
  phimax = phi[i]+(phi[i+1]-phi[i])*di;
  return 0.921034*lf->phistar*(phimax-phimin)*dm;
  }


/**************************** lf_intlumschechter *****************************/
/*
The luminosity density from Schechter (1976) galaxy luminosity function,
integrated from mmin to mmax.
*/
double	lf_intlumschechter(lfstruct *lf, double mmin, double mmax)
  {
   static double	phi[MABS_NSTEP], alphas;
   static int		flag;
   double		dm, di, phimin,phimax, moffset, mag;
   int			i, over;

  dm = (double)(MABS_MAX-MABS_MIN)/MABS_NSTEP;

  if (!flag || lf->alpha!=alphas)
    {
/*-- Never passed through here before! Let's initialize arrays */
    moffset = MABS_MIN + lf->mstar + 0.5*dm;
    for (i=0; i<MABS_NSTEP; i++)
      {
      mag = i*dm + moffset;
      phi[i] = lf_schechter(lf, mag)*exp(-0.921034*mag) + (i?phi[i-1]:0.0);
      }
    flag = 1;
    alphas = lf->alpha;
    }

/* Now make the interpolation */
/* Lower bound */
  mmin -= lf->mstar;
  if (mmin < MABS_MIN)
    mmin = MABS_MIN;
  di = (mmin-MABS_MIN)/dm-0.4999;
  i = (int)(di+0.5);
  di -= (double)i;
  if (i>=MABS_NSTEP-1)
    {
    over = i+2-MABS_NSTEP;
    i -= over;
    di += (double)over;
    }
  phimin = phi[i]+(phi[i+1]-phi[i])*di;

/* Upper bound */
  mmax -= lf->mstar;
  if (mmax < MABS_MIN)
    mmax = MABS_MIN;
  di = (mmax-MABS_MIN)/dm-0.4999;
  i = (int)(di+0.5);
  di -= (double)i;
  if (i>=MABS_NSTEP-1)
    {
    over = i+2-MABS_NSTEP;
    i -= over;
    di += (double)over;
    }
  phimax = phi[i]+(phi[i+1]-phi[i])*di;
  return 0.921034*lf->phistar*(phimax-phimin)*dm;
  }


/******************************* lf_schechter ********************************/
/*
The Schechter (1976) dphi/dm galaxy luminosity function.
*/
double	lf_schechter(lfstruct *lf, double m)
  {
  return exp(0.921034*(lf->alpha+1.0)*(lf->mstar-m)
	- exp(0.921034*(lf->mstar-m)));
  }


/****************************** lf_rndschechter ******************************/
/*
Return a random absolute magnitude following the Schechter (1976)
distribution.
*/
double  lf_rndschechter(lfstruct *lf, double mmin, double mmax)
  {
   double	mrnd, prnd;

  prnd  = lf_schechter(lf, mmax);
  if (1.0>prnd && mmax>lf->mstar)
    prnd = 1.0;
  do
    {
    mrnd = mmin + (mmax-mmin)*random_double();
/*-- Enveloppe */
    } while (prnd*random_double() > lf_schechter(lf, mrnd));

  return mrnd;
  }


/********************************** lf_evol **********************************/
/*
Luminosity function evolution as a function of redshift .
*/
lfstruct 	*lf_evol(lfstruct *lf, double z, lfstruct *lfevol)
  {
   double	lzp1;

  *lfevol = *lf;
  lzp1 = log(z + 1.0);
  lfevol->mstar += (lfevol->dmstar = lf->mstarevol*lzp1);
  lfevol->phistar *= exp(lf->phistarevol*lzp1);

  return lfevol;
  }


