/*
*				sed.c
*
* Manage Spectral Energy Distributions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 1999-2018 IAP/CNRS/UPMC
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
*	Last modified:		09/04/2018
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
#include "prefs.h"
#include "cosmo.h"
#include "lf.h"
#include "sed.h"


/****** sed_dup **************************************************************
PROTO	sedstruct *sed_dup(sedstruct *sed)
PURPOSE	Duplicate a SED or passband.
INPUT	SED structure pointer.
OUTPUT	Pointer to a copy of the input sed.
NOTES 	-.
AUTHOR	E. Bertin (IAP)
VERSION	27/05/2013
*/
sedstruct	*sed_dup(sedstruct *sed)
  {
   sedstruct	*sedout;

  QMALLOC(sedout, sedstruct, 1);
  *sedout = *sed;
  if (sedout->wave)
    {
    QMEMCPY(sed->wave, sedout->wave, double, sedout->ndata);
    }
  if (sedout->data)
    {
    QMEMCPY(sed->data, sedout->data, double, sedout->ndata);
    }

  return sedout;
  }


/****** sed_load *************************************************************
PROTO	sedstruct *sed_load(char *seddir_name, char *sed_name)
PURPOSE	Load, reformat, and combine SEDs or passbands.
INPUT	SED structure pointer.
OUTPUT	Pointer to the loaded sed.
NOTES	Negative lambda indicate that the data are stored in f_nu instead of
	f_lambda units.
AUTHOR	E. Bertin (IAP)
VERSION	25/09/2016
*/
sedstruct	*sed_load(char *seddir_name, char *sed_name)
  {
   sedstruct		*sedcomp[SED_MAXNCOMP];
   const char		notokstr[] = {" \t=,;\n\r\""};
   char			fullname[MAXCHAR], shortname[MAXCHAR],
			*tokptr1, *tokptr2,
			*cstr, *cstr2;
   sedstruct		*sed, *sed2;
   FILE			*infile;
   double		*wavet, *datat,
			waveval, dataval;
   int			j,k,last, dataflag, fnuflag, ncomp, ndatat;

  cstr = strtok_r(sed_name, notokstr, &tokptr1);
/* Decompose passband components */
  for (j=0; (cstr2 = strtok_r(j?NULL:cstr, "*", &tokptr1)); j++)
    {
/*--  Include directory path */
    strcpy(shortname, cstr2);
    if (*cstr2 != '/')
      {
      strcpy(fullname, seddir_name);
      strcat(fullname, cstr2);
      cstr2 = fullname;
      }
    if ((infile = fopen(cstr2,"r")) == NULL)
      {
      sprintf(gstr, "%s.pb", cstr2);
      if ((infile = fopen(gstr,"r")) == NULL)
        {
        sprintf(gstr, "%s.sed", cstr2);
        if ((infile = fopen(gstr,"r")) == NULL)
          {
          sprintf(gstr, "%s.ext", cstr2);
          if ((infile = fopen(gstr,"r")) == NULL)
            error(EXIT_FAILURE,"*ERROR*: can't read ", cstr2);
          }
        }
      }
    dataflag = 0;
    ndatat = SED_NDATA;
    sedcomp[j] = sed_new(strtok_r(shortname, notokstr, &tokptr2), ndatat);
    wavet = sedcomp[j]->wave;
    datat = sedcomp[j]->data;
    fnuflag = SED_FLAMBDADATA;
    for (k=1; fgets(gstr, MAXCHAR, infile);)
      {
      if (k>ndatat)
        {
        ndatat += SED_NDATA;
        QREALLOC(sedcomp[j]->wave, double, ndatat);
        wavet = sedcomp[j]->wave+k-1;
        QREALLOC(sedcomp[j]->data, double, ndatat);
        datat = sedcomp[j]->data+k-1;
        }
      if ((*gstr!=0) && (*gstr!=(char)'#')
	&& sscanf(gstr, "%lf %lf", &waveval, &dataval) == 2)
        {
        if (waveval < 0.0)
          {
          waveval = -waveval;
          fnuflag = SED_FNUDATA;
          }
        if (waveval < 1/BIG)
          waveval = 1/BIG;
        if (dataval>0.0)
          dataflag = 1;
        waveval*=ANGSTROEM;
        if (k>1 && waveval<*(wavet-1))
          error(EXIT_FAILURE, "*ERROR*: decreasing wavelength in ", cstr2);
        *(wavet++) = waveval;
        *(datat++) = dataval;
        k++;
        }
      }
    fclose(infile);

/*-- Save some unused memory */
    last = (sedcomp[j]->ndata = k-1) - 1;
    if (k<=ndatat)
      {
      QREALLOC(sedcomp[j]->wave, double, sedcomp[j]->ndata);
      QREALLOC(sedcomp[j]->data, double, sedcomp[j]->ndata);
      }

/*-- Characterize the sampling */
    if (!dataflag)
      warning("Empty SED or pass-band in ", cstr2);
    sedcomp[j]->wavemin = sedcomp[j]->wave[0];
    sedcomp[j]->wavemax = sedcomp[j]->wave[last];
    sedcomp[j]->fnu_flag = fnuflag;
    }

/* Combine passband components */
  ncomp = j;
  sed = sedcomp[0];
  if (ncomp>1)
    {
/*-- Multiply all component SEDs */
    for (j=1; j<ncomp; j++)
      {
      sed_mul(sed, 1.0, sedcomp[j], 1.0, &sed2);
      sed_end(sed);
      sed_end(sedcomp[j]);
      sed = sed2;
      }
    }

  return sed;
  }


/****** sed_new **************************************************************
PROTO	sedstruct *sed_new(char *name, int ndata)
PURPOSE	Create an empty SED or passband, filled with zeroes;
INPUT	SED name,
	number of SED points.
OUTPUT	Pointer to the new sed.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	27/05/2013
*/
sedstruct	*sed_new(char *name, int ndata)
  {
   sedstruct *sed;

  QMALLOC(sed, sedstruct, 1);

  sed->ndata = ndata;
  QCALLOC(sed->wave, double, ndata);
  QCALLOC(sed->data, double, ndata);
  strcpy(sed->name, name);
  sed->wavemin = sed->wavemax = 0.0;
  sed->fnu_flag = SED_FLAMBDADATA;

  return sed;
  }


/****** sed_end **************************************************************
PROTO	void sed_end(sedstruct *sed)
PURPOSE	Terminate SED or passband.
INPUT	Pointer to the sed.
OUTPUT	-.
NOTES	Memory is freed upon exit.
AUTHOR	E. Bertin (IAP)
VERSION	09/04/2018
*/
void	sed_end(sedstruct *sed)
  {
  if (sed) {
    free(sed->wave);
    free(sed->data);
    free(sed);
  }

  return;
  }


/****** sed_mul **************************************************************
PROTO	double sed_mul(sedstruct *sed1, double wf1,
			sedstruct *sed2, double wf2, sedstruct **psedo)
PURPOSE	Multiply and integrate 2 SEDs or passbands without degrading resolution.
INPUT	Pointer to SED #1,
	multiplicative factor for wavelengths of SED #1,
	pointer to SED #2,
	multiplicative factor for wavelengths of SED #2,
	pointer to the SED product.
OUTPUT	Sum of the SED product.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/05/2013
*/
double	sed_mul(sedstruct *sed1, double wf1, sedstruct *sed2, double wf2,
		sedstruct **psedo)
  {
   sedstruct	*sedo;
   double	*wsedo,*dsedo, *wsed1,*dsed1, *wsed2,*dsed2,
		wavemin,wavemax, wsed1a,dsed1a, wsed2a,dsed2a,
		ssed1a,ssed1b,ssed2a,ssed2b, integ, wa,wai,wb,w, cwf1,cwf2;
   int		nsedo,nsed1,nsed2, ndatamax, fnu1,fnu2;

/* Find the smallest range */
  wavemin = sed1->wavemin*wf1;
  if ((w=sed2->wavemin*wf2) > wavemin)
    wavemin = w;
  wavemax = sed1->wavemax*wf1;
  if ((w=sed2->wavemax*wf2) < wavemax)
    wavemax = w;
  ndatamax = sed1->ndata+sed2->ndata;

/* Allocate memory for destination if requested */
  if (psedo)
    {
    QMALLOC(sedo, sedstruct, 1);
    *psedo = sedo;
    QMALLOC(sedo->wave, double, ndatamax);
    QMALLOC(sedo->data, double, ndatamax);
    wsedo = sedo->wave;
    dsedo = sedo->data;
    nsedo = 0;
    }
  else
    {
    sedo = NULL;
    wsedo = dsedo = NULL;	/* Avoid gcc -Wall warnings */
    nsedo = 0;
    }

/* Now go through each available wavelength bin */
  integ = 0.0;
  fnu1 = sed1->fnu_flag;
  fnu2 = sed2->fnu_flag;
  cwf1 = C*wf1*wf1;
  cwf2 = C*wf2*wf2;
  dsed1a = dsed2a = 0.0;
  wsed1a = wsed2a = 0.0;
  ssed1a = ssed2a = 0.0;
  nsed1 = sed1->ndata;
  nsed2 = sed2->ndata;
  wsed1 = sed1->wave;
  dsed1 = sed1->data;
  wsed2 = sed2->wave;
  dsed2 = sed2->data;

  for (wa=wavemin*(1-2*SED_DLAMBDA); wa<=wavemax; wa=wb)
    {
    wai = wa/wf1;
    while (nsed1 && *wsed1/wai < 1.0+SED_DLAMBDA)
        {
        wsed1a = *(wsed1++);
        dsed1a = *(dsed1++);
        nsed1--;
        }
    wai = wa/wf2;
    while (nsed2 && *wsed2/wai < 1.0+SED_DLAMBDA)
        {
        wsed2a = *(wsed2++);
        dsed2a = *(dsed2++);
        nsed2--;
        }

/*-- Exit if no more data point available */
    if (!(nsed1 && nsed2))
      break;

/*-- Select the lowest wavelength bin as the next to come */
      wb = *wsed1*wf1;
      if ((w=*wsed2*wf2) < wb)
        wb = w;

/*-- Interpolate data points */
    ssed1b =  dsed1a + (*dsed1-dsed1a)*(wb/wf1-wsed1a)/(*wsed1-wsed1a);
    if (fnu1)
      ssed1b *= cwf1/(wb*wb);
    ssed2b =  dsed2a + (*dsed2-dsed2a)*(wb/wf2-wsed2a)/(*wsed2-wsed2a);
    if (fnu2)
      ssed2b *= cwf2/(wb*wb);
/*-- Add the contribution of the trapezium product to the integral */
    integ += (wb-wa)*((ssed1a*ssed2b+ssed1b*ssed2a)/2.0
		+ (ssed1b-ssed1a)*(ssed2b-ssed2a)/3.0);
    if (psedo)
      {
/*---- Multiply data */
      *(wsedo++) = wb;
      *(dsedo++) = ssed1b*ssed2b;
      nsedo++;
      }
/*-- The current SED values become the previous ones! */
    ssed1a = ssed1b;
    ssed2a = ssed2b;
    }

  if (psedo)
    {
/*-- Keep only what is necessary from what was allocated */
    sedo->ndata = nsedo;
    if (nsedo)
      {
      QREALLOC(sedo->wave, double, nsedo);
      QREALLOC(sedo->data, double, nsedo);
      sedo->wavemin = sedo->wave[0];
      sedo->wavemax = sedo->wave[nsedo-1];
      }
    else
      {
      free(sedo->wave);
      free(sedo->data);
      sedo->wave = sedo->data = NULL;
      sedo->wavemin = sedo->wavemax = 1.0;
      }
    sprintf(sedo->name, "%s * %s",  sed1->name, sed2->name);
    sedo->fnu_flag = fnu1|fnu2;
    }

  return integ;
  }


/****** sed_kcor *************************************************************
PROTO	double sed_kcor(sedstruct *sed, sedstruct *pb, double z)
PURPOSE	Compute k-correction (+band correction) for SED sed through passband pb.
INPUT	Pointer to the source SED,
	pointer to the passband "SED",
	source redshift.
OUTPUT	k-correction (in linear units).
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	27/05/2013
*/
double	sed_kcor(sedstruct *sed, sedstruct *pb, double z)
  {
   double	kcor;

  kcor = sed_mul(sed, 1.0+z, pb, 1.0, NULL)/(1.0+z);
/* Prevent infinite K-corrections */
  if (kcor<1/BIG)
    kcor = 1/BIG;
  return kcor;
  }


/****** sed_calib ************************************************************
PROTO	double sed_calib(sedstruct *sed, sedstruct *pb)
PURPOSE	Calibrate an SED: normalize it through a specific passband.
INPUT	Pointer to the source SED,
	pointer to the passband "SED".
OUTPUT	Computed normalization factor.
NOTES	The input SED is normalized with the computed normalization factor.
AUTHOR	E. Bertin (IAP)
VERSION	27/05/2013
*/
double	sed_calib(sedstruct *sed, sedstruct *pb)
  {
   double	*data, integ;
   int		i;

  integ = sed_mul(sed, 1.0, pb, 1.0, NULL);
  if (integ == 0)
    {
    sprintf(gstr, "Rest-frame SED %s has no overlap\n         with reference"
	" passband (%s)", sed->name, pb->name);
    error(EXIT_FAILURE, "*Error*: ", gstr);
    }
  data = sed->data;
  for (i=sed->ndata; i--;)
    *(data++) /= integ;

  return integ;
  }


/****** pb_calib *************************************************************
PROTO	void pb_calib(sedstruct **pb, sedstruct **pbcalibsed, int npb, 
		sedstruct *refpb, sedstruct *refcalibsed,
		sedstruct *backsed)
PURPOSE	Calibrate passbands: normalize them to a magnitude system.
INPUT	Pointer to the array of input passbands,
	pointer to the array if calibration passbands (e.g., Vega, AB,..),
	number of passbands,
	pointer to the reference passband,
	pointer to the reference calibration passband (e.g., Vega, AB,..),
	pointer to the background SED.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	27/05/2013
*/
void	pb_calib(sedstruct **pb, sedstruct **pbcalibsed, int npb, 
		sedstruct *refpb, sedstruct *refcalibsed,
		sedstruct *backsed)
  {
   sedstruct		*normsed, *pbt;
   double		*data, *wave, magzp, backinteg;
   int			i, j;

/* Normalize reference spectrum */
  normsed = sed_new("5556 A normalization", 3);
  normsed->wavemin = normsed->wave[0] = 5555.0*ANGSTROEM;
  normsed->wave[1] = 5556*ANGSTROEM;
  normsed->wavemax = normsed->wave[2] = 5557.0*ANGSTROEM;

/* The "line" FWHM is indeed 1 A, and we want the result to be REF_PHOTRATE */
  normsed->data[1] = 1.0/(REF_ENERGY*ANGSTROEM);
  sed_calib(refcalibsed, normsed);
  for (j=0; j<npb; j++)
    sed_calib(pbcalibsed[j], normsed);
  sed_end(normsed);

/* Normalize background spectrum */
  if (backsed)
    {
    data = backsed->data;
/*-- We divide by hnu */
    for (i=backsed->ndata; i--;)
      *(data++) /= ANGSTROEM;
    }


/* Convert fluxes to number of photons if DETECTION_TYPE is PHOTONS */
  if (prefs.refdetect_type == DETECT_PHOTONS)
    {
    data = refpb->data;
    wave = refpb->wave;
/*-- We divide by hnu */
    for (i=refpb->ndata; i--;)
      *(data++) *= *(wave++)/REF_PHOTENERGY/REF_WAVELENGTH;
    }

  for (j=0; j<npb; j++)
    if (prefs.obsdetect_type[j] == DETECT_PHOTONS)
      {
      pbt = pb[j];
      data = pbt->data;
      wave = pbt->wave;
/*---- We divide by hnu */
      for (i=pbt->ndata; i--;)
        *(data++) *= *(wave++)/REF_PHOTENERGY/REF_WAVELENGTH;
      }

/* Compute dot-product of the reference passband and the reference spectrum */
  sed_calib(refpb, refcalibsed);

  NPRINTF(OUTPUT, "------------- Passband name                 zero-point"
		" background (mag.arcsec-2)\n");

  for (j=0; j<npb; j++)
    {
    magzp = 2.5*log10(sed_calib(pb[j], pbcalibsed[j])*prefs.area[j]/prefs.gain[j]);
    backinteg = (backsed? sed_mul(backsed, 1.0, pb[j], 1.0, NULL): 1.0);
    NPRINTF(OUTPUT, strlen(pb[j]->name)>30?
		  "Passband #%-3d %-27.27s...  %+6.3f  %+6.3f\n"
		: "Passband #%-3d %-30.30s  %+6.3f  %+6.3f\n",
		j+1, pb[j]->name, magzp,-2.5*log10(backinteg/4.254e10));
    }

  return;
  }


/****** sed_extinc ***********************************************************
PROTO	double sed_extinc(sedstruct *sed, sedstruct *tau, double taufact
				sedstruct **psedo)
PURPOSE	Apply an extinction law to an SED.
INPUT	Pointer to the source SED,
	pointer to the extinction "SED" (tau),
	multiplicative factor for extinction "SED".
	pointer to the extinguished SED.
OUTPUT	Sum of the extincted SED.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	21/11/2016
*/
double	sed_extinc(sedstruct *sed, sedstruct *tau, double taufact,
				sedstruct **psedo)
  {
   sedstruct	*extinc;
   double	*extdata, *taudata,
		prod;
   int		i;

  extinc = sed_dup(tau);

/* Compute exp(-tau) */
  extdata = extinc->data;
  taudata = tau->data;
  for (i=extinc->ndata; i--;)
    *(extdata++) = exp(-taufact**(taudata++));

/* Compute the new spectrum */
  prod = sed_mul(sed, 1.0, extinc, 1.0, psedo);

  sed_end(extinc);

  return prod;
  }



