/*
*				galaxies.c
*
* Manage galaxies.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 1999-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		28/05/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "define.h"
#include "globals.h"
#include "cosmo.h"
#include "lf.h"
#include "galaxies.h"
#include "prefs.h"
#include "random.h"
#include "sed.h"


/******************************** galtype_init *******************************/
/*
Create type-dependent galaxy structures.
*/
galtypestruct	*galtype_init(double hubtype, double bt, double extinct,
			lfstruct *lf,
			sedstruct *bsed, sedstruct *dsed,
			sedstruct *tau_i, sedstruct *refpb)
  {
   galtypestruct	*galtype;
   sedstruct		*dised;
   double		sum1,sum2;

/* Allocate memory for the galaxy type structure */
  QCALLOC(galtype, galtypestruct, 1);

/* Copy properties */
  galtype->hubtype = hubtype;
  galtype->bt = bt;
/*-- Pre-compute magnitude-offsets for bulge and disk */
  galtype->bdm = bt>0.0? -2.5*log10(bt) : 100.0;
  galtype->ddm = bt<1.0? -2.5*log10(1.0-bt) : 100.0;
/*-- Initialize the luminosity function */
  QMALLOC(galtype->lf, lfstruct, 1);
  *(galtype->lf) = *lf;
/* Correct M* for disk internal extinction */
  galtype->lf->mstar -= 2.5*log10(1.0+(1.0-bt)
			*DEXP(0.12041*extinct));
/* Subtract from M* the average disk extinction to get the face-on value */
  galtype->beta = prefs.gal_dbeta;
/* The 1.67835 is here to convert the effective radius to a scale length */
  galtype->disk_rstar = prefs.gal_drstar*KPC/H/1.67835;
  galtype->disk_sigmalambda = prefs.gal_dsigmalambda;
  galtype->bulge_rstar = prefs.gal_brknee*KPC/H;
  galtype->bsed = sed_dup(bsed);
  galtype->dsed = sed_dup(dsed);
  galtype->tau_i = sed_dup(tau_i);
/* Normalization of extinction */
  sum1 = sed_mul(dsed, 1.0, refpb, 1.0, &dised);
  sum2 = sed_mul(dised, 1.0, tau_i, 1.0, NULL);
  sed_end(dised);
  galtype->disk_taufact = sum2>0.0? -0.4*extinct*sum1/sum2
				    : 0.0;
  return galtype;
  }


/******************************* galtype_end *********************************/
/*
Terminate a galaxy type structure.
*/
void	galtype_end(galtypestruct *galtype)
  {

  if (galtype->lf)
    free(galtype->lf);
  if (galtype->bsed)
    sed_end(galtype->bsed);
  if (galtype->dsed)
    sed_end(galtype->dsed);
  if (galtype->tau_i)
    sed_end(galtype->tau_i);

  free(galtype);

  return;
  }


/********************************* gal_end ***********************************/
/*
Terminate a galaxy structure.
*/
void	gal_end(galstruct *gal)
  {
  free(gal->bt);
  free(gal->mag);

  free(gal);

  return;
  }


/******************************* gal_bulgesize *******************************/
/*
Returns the average bulge size (in m) as a function of the absolute magnitude.
From Fig. 7 in Binggeli et al. (1984).
*/
double gal_bulgesize(galtypestruct *galtype, double mabs)
  {
/* Subtract the log10(re) = f(Mabs) knee from the absolute magnitude */
  mabs -= prefs.gal_bmknee + deltaMH;

  return galtype->bulge_rstar*DEXP(-(mabs<0.0? 0.3:0.1)*mabs)/H;
  }


/******************************* gal_bulgeflat *******************************/
/*
Returns the (random) apparent bulge flattening as a function of
cos(inclination) (Sandage et al. 1970). Warning: i is here in quadrature with
the original Sandage definition.
*/
double gal_bulgeflat(double cosi)
  {
   double	q;

/* Distribution of intrinsic flattening */
  q = GAL_BULGEQMEAN + random_gauss(GAL_BULGEQSIG);
  if (q>1.0)
    q = 1.0;

  return sqrt((q*q-1.0)*(1.0-cosi*cosi) + 1.0);
  }


/******************************* gal_disksize ********************************/
/*
Returns the (random) disk scale-size (in m) as a function of the absolute
magnitude (de Jong & Lacey 2000).
*/
double gal_disksize(galtypestruct *galtype, double mabs)
  {
/* A log-normal distribution */
  return galtype->disk_rstar*exp(0.921034*galtype->beta
		*(mabs-galtype->lf->mstar)
	+ random_gauss(galtype->disk_sigmalambda));
  }


/********************************* gal_cosi **********************************/
/*
Returns the (random) cos(inclination).
*/
double gal_cosi(void)
  {
/* Flat angular distribution */
  return cos(random_double()*GAL_MAXINCLIN*DEG);
  }


/********************************* gal_init **********************************/
/*
"Render" a galaxy at a given cosmological redshift with a given luminosity
function.
*/
galstruct *gal_init(galtypestruct *galtype, double z, double mabsmax,
			sedstruct **pb, int npb)
  {
   galstruct	*gal;
   sedstruct	*dsed;
   lfstruct	lfevol;
   double	mabs,mabsb,mabsd, da, kb,kd,kt, bt, dm;
   int		p;

  QCALLOC(gal, galstruct, 1);
  QMALLOC(gal->bt, double, npb);
  QMALLOC(gal->mag, double, npb);

/* Evolve the luminosity function */
  lf_evol(galtype->lf, z, &lfevol);

/* Absolute magnitude in the reference band */
  gal->mabs = mabs = lf_rndschechter(&lfevol, lfevol.mabsmin, mabsmax);
  mabsb = mabs + galtype->bdm;	  /* Bulge */
  mabsd = mabs + galtype->ddm;	  /* Disk */
/* Shape */
  gal->bposang = gal->dposang = 360.0*random_double() - 180.0;
  gal->dflat = gal_cosi();
  gal->bflat = gal_bulgeflat(gal->dflat);
/* Disk SED (absorption) */
  dsed = sed_extinc(galtype->dsed, galtype->tau_i,
		galtype->disk_taufact*log(gal->dflat));
/* Angular distance element (in m/arcsec)*/
  da = cosmo_dlum(z)*ARCSEC/((1.0+z)*(1.0+z));
/* Compensate for increase of surface brightness due to luminosity evolution */
  dm = lfevol.dmstar;
  gal->bsize = gal_bulgesize(galtype, mabsb-dm)/da*pow(1.0+z, prefs.gal_brevol);
  gal->dsize = gal_disksize(galtype, mabsd-dm)/da*pow(1.0+z, prefs.gal_drevol);
  bt = galtype->bt;
/* Now in each observed passband */
  for (p=0; p<npb; p++)
    {
/*-- Linear k-corrections for bulge, disk and total */
    kb = bt>0.0? sed_kcor(galtype->bsed, pb[p], z) : 1.0;
    kd = bt<1.0? sed_kcor(dsed, pb[p], z) : 1.0;
    if ((kt=bt*(kb-kd)+kd) < 1/BIG)
      kt = 1/BIG;
/*-- Apparent B/T in this passband */
    gal->bt[p] = bt*kb/kt;
/*-- Apparent magnitude */
    gal->mag[p] = mabs + cosmo_dmodul(z) - 2.5*log10(kt);
    }

  sed_end(dsed);

  return gal;
  }


/********************************* gal_shear ********************************/
/*
Apply a shear term to a galaxy structure.
*/
void	gal_shear(galstruct *gal, double kappa, double *gamma, int npb)
  {
   double	a,b,d, ct,st,mx2,my2,mxy, s11,s12,s21,s22, smx2,smy2,smxy, dm;
   int		p;

  s11 = 1.0 - kappa - gamma[0];
  s12 = s21 = -gamma[1];
  s22 = 1.0 - kappa + gamma[0];

/* disk */
  if (gal->dflat != 0.0 && gal->dsize!= 0.0)
    {
    a = gal->dsize/sqrt(gal->dflat);
    b = a*gal->dflat;
    ct = cos(0.0174533*gal->dposang);
    st = sin(0.0174533*gal->dposang);
    mx2 = (a*a*ct*ct+b*b*st*st);
    my2 = (a*a*st*st+b*b*ct*ct);
    mxy = (a*a-b*b)*ct*st;
    smx2 = mx2*s11*s11+my2*s12*s12+mxy*2.0*s11*s12;
    smy2 = mx2*s21*s21+my2*s22*s22+mxy*2.0*s21*s22;
    smxy = mx2*s11*s21+my2*s12*s22+mxy*(s11*s22+s12*s21);
    d = sqrt(0.25*(smx2-smy2)*(smx2-smy2)+smxy*smxy);
    a = sqrt(0.5*(smx2+smy2)+d);
    b = 0.5*(smx2+smy2)-d;
    b = b>=0? sqrt(b) : 0.0;
    gal->dposang = 28.64789*atan2(2.0*smxy, smx2-smy2);
    gal->dsize = sqrt(a*b);
    gal->dflat = b/a;
    }
/* bulge */
  if (gal->bflat != 0.0 && gal->bsize!= 0.0)
    {
    a = gal->bsize/sqrt(gal->bflat);
    b = a*gal->bflat;
    ct = cos(0.0174533*gal->bposang);
    st = sin(0.0174533*gal->bposang);
    mx2 = (a*a*ct*ct+b*b*st*st);
    my2 = a*a*st*st+b*b*ct*ct;
    mxy = (a*a-b*b)*ct*st;
    smx2 = mx2*s11*s11+my2*s12*s12;
    smy2 = mx2*s21*s21+my2*s22*s22;
    smxy = mx2*s11*s21+my2*s21*s22+mxy*(s11*s22+s12*s22);
    d = sqrt(0.25*(smx2-smy2)*(smx2-smy2)+smxy*smxy);
    a = sqrt(0.5*(smx2+smy2)+d);
    b = 0.5*(smx2+smy2)-d;
    b = b>=0? sqrt(b) : 0.0;
    gal->bposang = 28.64789*atan2(2*smxy, smx2-smy2);
    gal->bsize = sqrt(a*b);
    gal->bflat = b/a;
    }

  dm = 2.5*log10(fabs(s11*s22-s21*s12));
  for (p=0; p<npb; p++)
    gal->mag[p] += dm;

  return;
  }


