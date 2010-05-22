/*
                                  random.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        Stuff
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       functions returning random numbers.
*
*       Last modify:    23/05/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "define.h"
#include "globals.h"
#include "random.h"

#ifdef  HAVE_LGAMMA
#define LOGGAMMA	lgamma
#else
#define LOGGAMMA	gammln
static double		gammln();
#endif

static int		seed_time(void);

/****** random_gauss ********************************************************
PROTO   double random_gauss(double sigma)
PURPOSE Generate a random number with a (centered) Gaussian pdf.
INPUT   Standard deviation.
OUTPUT  Gauss-distributed random number.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 29/10/97
*/
double	random_gauss(double sigma)
  {
   double	x,y,z, r;

  while((z=pow(x=random_double()-0.5,2.0) + pow(y=random_double()-0.5,2.0))
	> 0.25);
  while ((r=random_double()) <= 0.0);

  return sigma*sqrt(-2.0*log(r)/z)*x;
  }


/****** random_int **********************************************************
PROTO   int random_int(void)
PURPOSE Generate a random integer over the range [0,2^31].
INPUT   -.
OUTPUT  Random integer number with uniform distribution.
NOTES   The actual upper bound of the range is implementation-dependent.
AUTHOR  E. Bertin (IAP)
VERSION 23/05/2010
*/
int	random_int(void)
  {
#if defined	HAVE_LRAND48
  return (int)lrand48();
#else
  return (int)rand();
#endif
  }


/****** random_double ********************************************************
PROTO   double random_double(void)
PURPOSE Generate a random number with uniform distribution over [0.0,1.0[
INPUT   -.
OUTPUT  Random double with uniform distribution.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 22/05/2010
*/
double	random_double(void)
  {

#if defined	HAVE_DRAND48
  return (double)drand48();
#else
  return (double)rand() / RAND_MAX;
#endif
  }


/****i* seed_time **********************************************************
PROTO   int seed_time(void)
PURPOSE Return a random generator seed based on current time.
INPUT   -.
OUTPUT  Integer based on current number of microseconds.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 23/05/2010
*/
int	seed_time(void)
  {
   struct timeval	tp;
   struct timezone	tzp;
   int			dummy;

  dummy = gettimeofday(&tp,&tzp);
  return (int)tp.tv_usec;
  }


/****** init_random **********************************************************
PROTO   void init_random(int seed)
PURPOSE Initialize the random number generator.
INPUT   Seed.
OUTPUT  -.
NOTES   The seed is used to initialize the random sequence at a particular
        position, which is implementation-dependent. If seed = 0, then the
        actual seed is taken from the time() function (which varies each
        second).
AUTHOR  E. Bertin (IAP)
VERSION 23/05/2010
*/
void	init_random(int seed)
  {
#if defined	HAVE_SRAND48
  if (seed)
    srand48((long)seed);
  else
    srand48((long)seed_time());
#else
  if (seed)
    srand((unsigned int)seed);
  else
    srand((unsigned int)seed_time());
#endif

  return;
  }


/****i* gammln ***************************************************************
PROTO   double gammln(double xx)
PURPOSE Returns the log of the Gamma function (from Num. Recipes in C, p.168).
INPUT   A double.
OUTPUT  Log of the Gamma function.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 29/10/97
*/
double	gammln(double xx)

  {
   double		x,tmp,ser;
   static double	cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
   int			j;

  tmp=(x=xx-1.0)+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<6;j++)
    ser += cof[j]/(x+=1.0);

  return log(2.50662827465*ser)-tmp;
  }


/****** random_poisson *******************************************************
PROTO   double random_double(double xm)
PURPOSE Returns a random number with Poisson deviate (from Num. Recipes in C.,
        p.222) centered on xm.
INPUT   Mean of the Poisson distribution.
OUTPUT  A double containing the integer (!) variable with Poisson deviate.
NOTES   I am still searching for a faster algorithm!!
AUTHOR  E. Bertin (IAP)
VERSION 23/05/2010
*/
double	random_poisson(double xm)
  {
   static double	sq,alxm,g,oldm=(-1.0);
   double		em,t,y;
   double		gammln();

  if (xm < 12.0)
    {
    if (xm != oldm)
      {
      oldm=xm;
      g=exp(-xm);
      }
    em = -1.0;
    t=1.0;
    do
      {
      em += 1.0;
      t *= random_double();
      } while (t > g);
    }
  else
    {
    if (xm != oldm)
      {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-LOGGAMMA(xm+1.0);
      }
    do
      {
      do
        {
        y=tan(PI*random_double());
        em=sq*y+xm;
        } while (em < 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-LOGGAMMA(em+1.0)-g);
      } while (random_double() > t);
    }

  return em;
  }

