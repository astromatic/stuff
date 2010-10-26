/*
*				cosmo.c
*
* Manage cosmology.
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

/******************************** cosmo_dconf ********************************/
/*
Conformal distance. Numerical integration along geodesic, allowing for a
cosmological constant.
*/
double	cosmo_dconf(double z0)
  {
   static double	chi[Z_NSTEP], omegas, lambdas, dlz;
   static int		flag;
   double		dz,z, di;
   int			i, over;

  if (!flag || OmegaM!=omegas || OmegaL!=lambdas)
    {
/*-- Never passed through here before! Let's initialize arrays */
    dlz = log(Z_MAX/Z_MIN)/Z_NSTEP;
    for (i=0; i<Z_NSTEP; i++)
      {
      z = Z_MIN*exp((i+0.5)*dlz);
      dz = z*dlz;
      chi[i] = i?
	 chi[i-1] + 1/sqrt((1+z)*(1+z)*(1+OmegaM*z)-z*(2+z)*OmegaL)*dz
	:z;	/* Initialize with the Euclidean(flat) approximation */
      }
    flag = 1;
    omegas = OmegaM;
    lambdas = OmegaL;
    }

/* Now make the interpolation */
  if (z0>=Z_MIN)
    {
    di = (log(z0/Z_MIN)/dlz-0.5);
    i = (int)(di+0.5);
    di -= (double)i;
    if (i>=Z_NSTEP-1)
      {
      over = i+2-Z_NSTEP;
      i -= over;
      di += (double)over;
      }
    return chi[i]+(chi[i+1]-chi[i])*di;
    }
  else
    return chi[0]*z0/Z_MIN;
  }


/******************************** cosmo_dlum *********************************/
/*
Luminosity distance, allowing for a cosmological constant.
See Carroll (1992).
*/
double	cosmo_dlum(double z)
  {
   double	chi, omegak;

/* Compute conformal distance */
  chi = cosmo_dconf(z);

/* Compute the curvature */
  omegak = 1.0 - OmegaM - OmegaL;
  if (omegak>0.0)
/*-- Open Universe */
    return C/H0*(1+z)*sinh(sqrt(omegak)*chi)/sqrt(omegak);
  else if (omegak<0.0)
/*-- Closed Universe */
    return C/H0*(1+z)*sin(sqrt(-omegak)*chi)/sqrt(-omegak);
  else
/*-- Flat Universe */
    return C/H0*(1+z)*chi;
  }


/****************************** cosmo_dvol ***********************************/
/*
Return the comoving volume element.
see Peebles (1993).
*/
double	cosmo_dvol(double z)
  {
   double	chi, omegak, sk;

/* Compute conformal distance */
  chi = cosmo_dconf(z);

/* Compute the curvature */
  omegak = 1.0 - OmegaM - OmegaL;
  if (omegak>0.0)
/*-- Open Universe */
    {
    sk = sinh(sqrt(omegak)*chi);
    return C*C*C*sk*sk/(H0*H0*H0*omegak
	*sqrt((1+z)*(1+z)*(1+OmegaM*z)-z*(2+z)*OmegaL));
    }
  else if (omegak<0.0)
/*-- Closed Universe */
    {
    sk = sin(sqrt(-omegak)*chi);
    return -C*C*C*sk*sk/(H0*H0*H0*omegak
	*sqrt((1+z)*(1+z)*(1+OmegaM*z)-z*(2+z)*OmegaL));
    }
  else
/*-- Flat Universe */
    return C*C*C*chi*chi/(H0*H0*H0*sqrt(OmegaM*(1+z)*(1+z)*(1+z)+OmegaL));
  }


/****************************** cosmo_dlc ***********************************/
/*
Return the comoving (physical) line element:
dlc/dz = dV/dOmega/dz*((1+z)/D_lum)**2
*/
double	cosmo_dlc(double z)
  {
  return C/(H0*sqrt((1+z)*(1+z)*(1+OmegaM*z)-z*(2+z)*OmegaL));
  }


/******************************* cosmo_dmodul ********************************/
/*
Compute the distance modulus at a given z (excluding k- and e-corrections).
*/
double	cosmo_dmodul(double z)
  {
/* The first term is -2.5*log10(10.0/4*PI) (removed) */
  return /* 0.24802466 + */ 5*log10(cosmo_dlum(z)/(10*PC));
  }


/******************************* cosmo_dmoduldif *****************************/
/*
Modulus difference (for computing z_max).
*/
double	cosmo_dmoduldif(double z, lfstruct *lf, double mappmax)
  {
   lfstruct	lfevol;
   double	dm;

  dm = mappmax - lf_evol(lf,z, &lfevol)->mabsmin
	- cosmo_dmodul(z)		/* Distance modulus */
    /*+ 2.5*log10(1+z)*/;			/* Basic K-correction */
  return dm;
  }


/********************************* cosmo_zmax ********************************/
/*
Perform minimisation through line-search. Based on algorithm described in
Numerical Recipes.
*/

#define SHFT(a,b,c,d)   {(a)=(b);(b)=(c);(c)=(d);}      /* For line-search */
#define SIGN(a,b)       ((b)>0.0? fabs(a) : -fabs(a))   /* For line-search */
#define CGOLD		0.3819660	/* Complement to the golden section */
#define TINY		1e-20		/* Almost nothing */
#define ITMAX		100		/* Max. nb of iter. in line-search */
#define TOL		1e-3		/* Fract. tolerance in line-search */

double	cosmo_zmax(lfstruct *lf, double mappmax)
  {
   double	a,b,d,e,etemp, fu,fv,fw,fx, p, q,r, tol1,tol2, u,v,w,x,xm,
		ax,bx,cx, dm;
   int		iter;

  ax = 0.0;
  cx = Z_MAX;
  bx = (ax+cx)/2.0;

/* Brent's algorithm for finding the minimum */
  d= e = 0.0;
  a = (ax < cx) ? ax : cx;
  b = (ax > cx) ? ax : cx;
  x = w = v = bx;
  dm = cosmo_dmoduldif(x, lf, mappmax);
  fw = fv = fx = dm*dm;
  for (iter=ITMAX; iter--;)
    {
    xm = 0.5*(a+b);
    tol2 = 2 * (tol1=TOL*fabs(x)+TINY);
    if (fabs(x-xm) <= (tol2-0.5*(b-a)))
      return x;
    if (fabs(e) > tol1)
      {
      r = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v)*q - (x-w)*r;
      q = 2*(q-r);
      if (q > 0.0)
        p = -p;
      q = fabs(q);
      etemp = e;
      e = d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
        d = CGOLD*(e=(x >= xm ? a-x : b-x));
      else
        {
        d = p/q;
        u = x+d;
        if (u-a < tol2 || b-u < tol2)
          d = SIGN(tol1,xm-x);
        }
      }
    else
      d = CGOLD*(e=(x >= xm ? a-x : b-x));
    u = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    dm = cosmo_dmoduldif(u, lf, mappmax);
    if ((fu=dm*dm) <= fx)
      {
      if (u >= x)
        a = x;
      else
        b = x;
      SHFT(v, w, x, u);
      SHFT(fv, fw, fx, fu);
      }
    else
      {
      if (u < x)
        a = u;
      else
        b = u;
    if (fu <= fw || w == x)
      {
      v = w;
      w = u;
      fv = fw;
      fw = fu;
      }
    else if (fu <= fv || v == x || v == w)
      {
      v = u;
      fv = fu;
      }
    }
  }

  warning("Too many iterations in ", "zmax()");
  return  x;
  }

#undef	SHFT
#undef	SIGN

#undef CGOLD
#undef TINY
#undef ITMAX
#undef TOL

