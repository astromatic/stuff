/*
*				igm.c
*
* Reproduce the effects of the inter-galactic medium.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		27/05/2013
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
#include "igm.h"
#include "sed.h"

/* Wavelengths and absorption coefficients for Lyman hydrogen lines */
/* from xspec code by Martin Still (NASA Ames) and Frank Marshall (NASA/GSFC):*/
/* http://heasarc.gsfc.nasa.gov/xanadu/xspec/models/zigm.html */

const double lyman_wave[] = {1215.67, 1025.72, 972.537, 949.743, 937.803,
			     930.748, 926.226, 923.150, 920.963, 919.352,
			     918.129, 917.181, 916.429, 915.824, 915.329,
			     914.919, 914.576, 911.75};
const double lyman_coeff[] = {0.0036, 0.0017, 0.0011846, 0.0009410, 0.0007960,
			      0.0006967,0.0006236,0.0005665,0.0005200,0.0004817,
			      0.0004487,0.0004200,0.0003947,0.000372 ,0.0003520,
 			      0.0003334,0.00031644};

/****** sed_igmmadauextinct **************************************************
PROTO   sedstruct *sed_igmmadauextinct(double z)
PURPOSE Compute the extinction curve (in the observer's restframe) produced by
	the IGM on a source at a given redshift
INPUT   Source redshift.
OUTPUT  Extinction curve.
NOTES   Based on Madau et al. 1995
	http://adsabs.harvard.edu/abs/1995ApJ...441...18M
	with input from Damien Leborgne <leborgne@iap.fr>
AUTHOR  E. Bertin (IAP)
VERSION 27/05/2013
*/
sedstruct	*sed_igmmadauextinct(double z)
  {
   sedstruct	*taused;
   double	*tauwave,*taudata,
		w, dtau,tau, xem,xem018,xem046,xem168,xem132,
		xc,xc300,xc150,xc018,xc046,xc168,xc132;
   int		i,l, nlambda,nlambdamax;

  xem = 1.0+z;
  xem018 = pow(xem,0.18);
  xem046 = pow(xem,0.46);
  xem168 = pow(xem,1.68);
  xem132 = xem168/(xem*xem*xem);
  nlambdamax = (IGM_WAVE_END - IGM_WAVE_START) / IGM_WAVE_STEP + 1;
  taused = sed_new("IGM extinction", nlambdamax);
  tauwave = taused->wave;
  taudata = taused->data;
  w = IGM_WAVE_START;
  for (l=0; l<nlambdamax; l++)
    {
    tau = IGM_METAL_COEFF *pow(w/lyman_wave[0], IGM_METAL_POWER);
    for (i=0; (xc=w/lyman_wave[i]) < xem && i<IGM_NLYMAN; i++)
      tau += lyman_coeff[i] * pow(xc, IGM_LYMAN_POWER);
    if (!i)
      break;
    if (i==IGM_NLYMAN)
      {
      xc300 = xc*xc*xc;
      xc150 = sqrt(xc300);
      xc018 = pow(xc,0.18);
      xc046 = pow(xc,0.46);
      xc168 = xc150*xc018;
      xc132 = xc018/xc150;
      if ((dtau = 0.25*xc300*(xem046-xc046) + 9.4*xc150*(xem018-xc018)
		-0.7*xc300*(xc132-xem132) - 0.023*(xem168-xc168)) > 0.0)
        tau += dtau;
      }
    *(tauwave++) = w*ANGSTROEM;
    *(taudata++) = tau;
    w += IGM_WAVE_STEP;
    }

  nlambda = l;
  if (nlambda<2)
    error(EXIT_FAILURE, "*Internal Error*:",
	"source SED out of IGM extinction curve limits in sed_avgigmextinct()");
  if (nlambda<nlambdamax)
    {
    *(tauwave++) = w*ANGSTROEM;
    *(taudata++) = 0.0;
    nlambda++;
    }
  if (nlambda<nlambdamax)
    {
    *tauwave = IGM_WAVE_END*ANGSTROEM;
    *taudata = 0.0;
    nlambda++;
    }

  taused->wavemin = *taused->wave;
  taused->wavemax = *tauwave;

  taused->ndata = nlambda;
  QREALLOC(taused->wave, double, taused->ndata);
  QREALLOC(taused->data, double, taused->ndata);

  return taused;
  }


