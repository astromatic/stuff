/*
*				makeit.c
*
* Main loop.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 1999-2016 IAP/CNRS/UPMC
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
*	Last modified:		23/11/2016
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifdef USE_THREADS
#include <pthread.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "define.h"
#include "globals.h"
#include "prefs.h"
#include "celsys.h"
#include "clusters.h"
#include "cosmo.h"
#include "galaxies.h"
#include "igm.h"
#include "lf.h"
#include "random.h"
#include "sed.h"

/********************************** makeit ***********************************/
void	makeit(void)
  {
   celsysstruct		*celsys;
   clusterstruct	**clusters;
   galtypestruct	**galtype;
   galstruct		*gal;
   sedstruct		*pb[SED_MAXNPB],*pbcorr[SED_MAXNPB],
			*pbcalibsed[SED_MAXNPB],
			*galsed[GAL_MAXNSED], *starsed[STAR_MAXNSED],
			*refcalibsed, *refpb, *tau_i, *tau_igm,
			*bsed,*dsed, *backsed, *sed;
   lfstruct		*lf, *lfevol;
   time_t		thetime, thetime2;
   struct tm		*tm;
   FILE			*(catfile[SED_MAXNPB]);
   static char		seddir_name[MAXCHAR], filterdir_name[MAXCHAR],
			extdir_name[MAXCHAR];
   double		pos[2], mabsmin, mabsmax, dn, bt, kb,kd,kt,
			width,height, xscale,yscale, xoffset,yoffset,
			x,y, radius, omcradius, domega, mag, fluxratio,
			z, dz,zlow,zhigh, zmax, dlc, galdens, dtime;
   long			nsource;
   int			c,i,g,n,p,s, npb, ngaltype, ngalsed, nstarsed, ncluster,
			clustindex, shearflag, angcoordflag;

/* Processing start date and time */
  thetime = time(NULL);
  dtime = counter_seconds();
  tm = localtime(&thetime);
  sprintf(prefs.sdate_start,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_start,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);

  NFPRINTF(OUTPUT, "");
  QPRINTF(OUTPUT,
	"----- %s %s started on %s at %s with %d thread%s\n\n",
		BANNER,
		MYVERSION,
		prefs.sdate_start,
		prefs.stime_start,
		prefs.nthreads,
		prefs.nthreads>1? "s":"");

/* Set frequently used variables */
  NFPRINTF(OUTPUT, "Initializing constants...")
/* Cosmological parameters */
  H0 = prefs.h0*KM/MPC;
  H = prefs.h0/100.0;
  OmegaM = prefs.omegam;
  OmegaL = prefs.omegal;
  deltaMH = 5*log10(H);	/* Action of H0 on absolute magnitudes */

  npb = prefs.npb_name;

/* Constant shear */
  shearflag = fabs(prefs.lens_gamma[0])>1e-6 || fabs(prefs.lens_gamma[1])>1e-6
		|| fabs(prefs.lens_kappa)>1e-6 ;

/* Coordinates */
  angcoordflag = (prefs.coord_type != COORD_PIXEL);
  if (angcoordflag) {
/*-- Initialize celestial coordinates */
    celsys = celsys_init(prefs.coord_center);
    radius = prefs.field_size[0] / 2.0 * DEG;
    omcradius = 1.0 - cos(radius);
    domega = 2.0 * PI * omcradius;
  } else {
/*-- Cartesian position offsets and scales */
    width = prefs.pixscale * prefs.field_size[0] * 1.1;
    height = prefs.pixscale * prefs.field_size[1] * 1.1;
    xscale = width / prefs.pixscale;
    yscale = height / prefs.pixscale;
    xoffset = (xscale - prefs.field_size[0]) / 2.0;
    yoffset = (yscale - prefs.field_size[1]) / 2.0;
    domega = width * height * ARCSEC * ARCSEC;
  }

/* Set default directories */
  sprintf(seddir_name, "%s/seds/", prefs.datadir_name);
  NFPRINTF(OUTPUT,"");
  NPRINTF(OUTPUT, "SED directory:     %s\n", seddir_name);
  sprintf(filterdir_name, "%s/filters/", prefs.datadir_name);
  NFPRINTF(OUTPUT,"");
  NPRINTF(OUTPUT, "Filter directory:  %s\n", filterdir_name);
  sprintf(extdir_name, "%s/extinct/", prefs.datadir_name);
  NFPRINTF(OUTPUT,"");
  NPRINTF(OUTPUT, "Extinct directory: %s\n", extdir_name);

/* Get the observed passbands */
  NFPRINTF(OUTPUT, "Loading passbands...")
  for (p=0; p<npb; p++)
    pb[p] = sed_load(filterdir_name, prefs.pb_name[p]);

/* Get the reference passband */
  refpb = sed_load(filterdir_name, prefs.refpb_name);

/* Get the calibration SED (which defines the magnitude system) */
  NFPRINTF(OUTPUT, "Loading the calibration SEDs...")
  refcalibsed = sed_load(seddir_name, prefs.refcalibsed_name);
  for (p=0; p<npb; p++)
    pbcalibsed[p] = sed_load(seddir_name, prefs.pbcalibsed_name[p]);

/* Get the sky background SED */
  if (*prefs.backsed_name != '*')
    {
    NFPRINTF(OUTPUT, "Loading the background SED...")
    backsed = sed_load(seddir_name, prefs.backsed_name);
    }
  else
    backsed = NULL;

/* Calibrate the passbands in the chosen magnitude system */
  NFPRINTF(OUTPUT, "Calibrating passbands...")
  pb_calib(pb, pbcalibsed, npb, refpb, refcalibsed, backsed);

/* Forget the calibration SED to gain memory */
  sed_end(refcalibsed);
  for (p=0; p<npb; p++)
    sed_end(pbcalibsed[p]);

/* Get the star SEDs */
  NFPRINTF(OUTPUT, "Loading star SEDs...")
  nstarsed = prefs.nstar_sedname;
  for (p=0; p<nstarsed; p++)
    starsed[p] = sed_load(seddir_name, prefs.star_sedname[p]);

/* Calibrate the star SEDs */
  NFPRINTF(OUTPUT, "Calibrating star SEDs...")
  for (i=0; i<nstarsed; i++)
    sed_calib(starsed[i], refpb);

/* Get the galaxy SEDs */
  NFPRINTF(OUTPUT, "Loading galaxy SEDs...")
  ngalsed = prefs.ngal_sedname;
  for (p=0; p<ngalsed; p++)
    galsed[p] = sed_load(seddir_name, prefs.gal_sedname[p]);

/* Calibrate the galaxy SEDs */
  NFPRINTF(OUTPUT, "Calibrating galaxy SEDs...")
  for (i=0; i<ngalsed; i++)
    sed_calib(galsed[i], refpb);

/* Get the internal extinction law (in tau units) */
  NFPRINTF(OUTPUT, "Loading the internal extinction law...")
  tau_i = sed_load(extdir_name, prefs.extinct_name);

/* Initialize field galaxy types (luminosity functions, SEDs) */
  NFPRINTF(OUTPUT, "Initializing galaxy types...")
  ngaltype = prefs.gal_ntype;
  QMALLOC(galtype, galtypestruct *, ngaltype);
  QMALLOC(lf, lfstruct, 1);	/* Used as a template only */
  lf->mabsmax = prefs.lf_mlim[1]+deltaMH;	/* This is indep. of gal.type*/
  lf->mabsmin = prefs.lf_mlim[0]+deltaMH;	/* This is indep. of gal.type*/
  for (g=0; g<ngaltype; g++)
    {
    lf->phistar = prefs.lf_phistar[g]*H*H*H/(MPC*MPC*MPC);
/*-- Subtract from M* the average disk extinction to get the face-on value */
    lf->mstar = prefs.lf_mstar[g] + deltaMH;
    lf->alpha = prefs.lf_alpha[g];
    lf->phistarevol = prefs.lf_phistarevol[g];
    lf->mstarevol = prefs.lf_mstarevol[g];
/*-- Spectral Energy Distributions (SEDs) */
/*-- Bulge component */
    if ((n=prefs.gal_bsedno[g]-1) >= ngalsed)
      error(EXIT_FAILURE, "*Error*: no such Bulge SED-index", "");
    bsed = galsed[n];
/*-- Disk component */
    if ((n=prefs.gal_dsedno[g]-1) >= ngalsed)
      error(EXIT_FAILURE, "*Error*: no such Disk SED-index", "");
    dsed = galsed[n];
    galtype[g] = galtype_init(prefs.gal_hubtype[g], prefs.gal_bt[g],
			prefs.gal_iextinc[g], lf, bsed, dsed, tau_i, refpb);
    }

  free(lf);
  QMALLOC(lfevol, lfstruct, 1);	/* Used as a template only */

/* Use an appropriate step (constant in comoving radial coordinate) */
  dlc = prefs.integ_zstep*MPC/H;

/* Initialize the random generator */
  NFPRINTF(OUTPUT, "Initializing galaxy random generator...")
  init_random(prefs.galseed);

/* Check out if a list of clusters is provided */
  ncluster = 0;
  if (strcmp(prefs.clusterlist_name, "NONE"))
    {
    NFPRINTF(OUTPUT, "Reading cluster list...")
    cluster_read(prefs.clusterlist_name, &clusters, &ncluster);
    NFPRINTF(OUTPUT,"");
    NPRINTF(OUTPUT, "%d cluster%s found\n", ncluster, ncluster>1?"s":"");
    if (ncluster)
      cluster_sortz(clusters, ncluster);
    }

/* Open the output catalogs */
  NFPRINTF(OUTPUT, "Creating catalogs...")
  for (p=0; p<npb; p++)
    if ((catfile[p] = fopen(prefs.cat_name[p],"w")) == NULL)
      error(EXIT_FAILURE,"*ERROR*: can't create ", prefs.cat_name[p]);

/* Integrate from "far away" to 0 redshift */
  NFPRINTF(OUTPUT, "Generating galaxies...")
  zmax = Z_MAX;
  nsource = 0;
  clustindex = 0;
  for (zhigh=zmax; zhigh>0.0; zhigh=zlow)
    {
/*-- Use an appropriate step (constant in comoving radial coordinate) */
    dz = dlc/cosmo_dlc(zhigh);
    zlow = zhigh - dz;
    if (zlow<0.0)
      {
      zlow = 0.0;
      dz = zhigh - zlow;
      }
    z = (zlow+zhigh)/2.0;

    sprintf(gstr, "Generating galaxies... z =%6.2f", z);
    NFPRINTF(OUTPUT, gstr)

/*-- Compute corrected passbands because of IGM opacity */
    if (prefs.igm_type != IGM_NONE)
      {
      tau_igm = sed_igmmadauextinct(zlow);
      for (p=0; p<npb; p++)
        sed_extinc(pb[p], tau_igm, 1.0, &pbcorr[p]);
      sed_end(tau_igm);
      }
    else
      for (p=0; p<npb; p++)
        pbcorr[p] = pb[p];

/*-- Normalize cluster densities */
    for (c=clustindex ;c<ncluster && clusters[c]->z > zlow; c++)
      cluster_normdens(clusters[c], galtype, ngaltype);

/*-- Loop over galaxy types */
    for (g=0; g<ngaltype; g++)
      {
      lf = galtype[g]->lf;
      bt = galtype[g]->bt;
      if (zlow>0.0)
        {
/*------ Compute a lower limit to the brigthness of detectable galaxies */
        mabsmax = prefs.maglim[1] - gal_mag(galtype[g], 0.0, zlow, refpb,
					NULL, NULL);
        if (lf->mabsmax < mabsmax)
          mabsmax = lf->mabsmax;
        }
      else
        mabsmax = lf->mabsmax;

/*---- Make the luminosity function evolve with redshift */
      lf_evol(lf, z, lfevol);

/*---- Field galaxies */
/*---- Compute the volume element and the local density -> dN = n*dV */
      galdens = lf_intschechter(lfevol, lfevol->mabsmin, mabsmax);
      dn = cosmo_dvol(z)*dz*domega*galdens;
      n = random_poisson(dn);
      nsource += n;
/*---- Generate field galaxies */
      for (i=n; i--;)
        {
/*------ Redshift */
        z = random_double()*(zhigh-zlow)+zlow;
        gal = gal_init(galtype[g], z, mabsmax, refpb, pbcorr, npb);
        if (gal->refmag < prefs.maglim[0] || gal->refmag > prefs.maglim[1])
          {
          gal_end(gal);
          continue;
          }
/*------ Add constant shear */
        if (shearflag)
          gal_shear(gal, prefs.lens_kappa, prefs.lens_gamma, npb);

/*------ Randomize coordinates */
        x = random_double();
        y = random_double();
        if (angcoordflag) {
/*-------- Position in celestial (angular) coordinates */
          pos[0] = x * 360.0;
          pos[1] = asin(1.0 - y * omcradius) / DEG;
/*-------- Convert to equatorial */
          celsys_to_eq(celsys, pos);
        } else {
/*-------- Position in Cartesian (pixel) coordinates */
          pos[0] = x * xscale - xoffset;
          pos[1] = y * yscale - yoffset;
        }
/*------ Now in each observed passband */
        for (p=0; p<npb; p++)
          {
          fprintf(catfile[p],
		angcoordflag ? "200 %11.7f %+11.7f %8.4f %5.3f "
			"%9.3f %5.3f %+7.2f %9.3f %5.3f %+7.2f %8.5f %+5.1f\n"
			  : "200 %11.4f %11.4f %8.4f %5.3f "
			"%9.3f %5.3f %+7.2f %9.3f %5.3f %+7.2f %8.5f %+5.1f\n",

		pos[0], pos[1],
		gal->mag[p], gal->bt[p],
		gal->bsize, gal->bflat, gal->bposang,
		gal->dsize, gal->dflat, gal->dposang,
		z, galtype[g]->hubtype);
          }
        gal_end(gal);
        }
/*---- Cluster galaxies */
      for (c=clustindex ;c<ncluster && clusters[c]->z > zlow; c++)
        {
/*------ Compute the number of galaxies to simulate */
        n = random_poisson(clusters[c]->lf_eqvol*galdens);
        nsource += n;
        clusters[c]->real_ngal += n;
/*------ Generate cluster galaxies */
        for (i=n; i--;)
          {
          cluster_rndgalpos(clusters[c], &x, &y, &z);
          if (z<=0.0)
            continue;
          gal = gal_init(galtype[g], z, mabsmax, refpb, pbcorr, npb);
          if (gal->refmag < prefs.maglim[0] || gal->refmag > prefs.maglim[1])
            {
            gal_end(gal);
            continue;
            }
          clusters[c]->real_mabs += exp(-0.921034*gal->mabs);
/*-------- Add constant shear */
          if (shearflag)
            gal_shear(gal, prefs.lens_kappa, prefs.lens_gamma, npb);

/*-------- Position in pixel coordinates*/
          pos[0] = clusters[c]->x + x / prefs.pixscale;
          pos[1] = clusters[c]->y + y / prefs.pixscale;
/*-------- Now in each observed passband */
          for (p=0; p<npb; p++)
            {
            fprintf(catfile[p],
		"200 %11.4f %11.4f %8.4f %5.3f "
		"%9.3f %5.3f %+7.2f %9.3f %5.3f %+7.2f %8.5f %+5.1f\n",
		pos[0], pos[1],
		gal->mag[p], gal->bt[p],
		gal->bsize, gal->bflat, gal->bposang,
		gal->dsize, gal->dflat, gal->dposang, z, galtype[g]->hubtype);
            }
          gal_end(gal);
          }
        }
      }
/*-- Compute cluster luminosity stats */
    for (c=clustindex ;c<ncluster && clusters[c]->z > zlow; c++)
      if (clusters[c]->real_ngal)
        clusters[c]->real_mabs = -2.5*log10(clusters[c]->real_mabs);
    clustindex = c;
    for (p=0; p<npb; p++)
      if (pbcorr[p] != pb[p])
        sed_end(pbcorr[p]);
    }

  NFPRINTF(OUTPUT,"");
  NPRINTF(OUTPUT, "%ld galax%s generated\n", nsource, nsource>1?"ies":"y");

/* Generate field stars */
  if ((prefs.starflag)) {
//-- Initialize the random generator
    NFPRINTF(OUTPUT, "Initializing star random generator...")
    init_random(prefs.starseed);

    dn = domega*10000000000.0;
    n = (int)random_poisson(dn);
    nsource += n;
    sprintf(gstr, "Adding %d stars...", n);
    NFPRINTF(OUTPUT, gstr);

    for (i=n; i--;) {
//---- Randomize coordinates */
      x = random_double();
      y = random_double();
      if (angcoordflag) {
//------ Position in celestial (angular) coordinates */
        pos[0] = x * 360.0;
        pos[1] = asin(1.0 - y * omcradius) / DEG;
//------ Convert to equatorial */
        celsys_to_eq(celsys, pos);
      } else {
//------ Position in Cartesian (pixel) coordinates */
        pos[0] = x * xscale - xoffset;
        pos[1] = y * yscale - yoffset;
      }
      mag = prefs.maglim[1] + log10(random_double()+1e-8) / 0.4;

//---- Randomize the stellar SED */
      while ((s=(int)(random_double()*nstarsed)) >= nstarsed);
      sed = starsed[s];

//---- Now in each observed passband */
      for (p=0; p<npb; p++) {
        fluxratio = sed_mul(sed, 1.0, pb[p], 1.0, NULL);
        fprintf(catfile[p],
		angcoordflag ? "100 %11.7f %+11.7f %8.4f\n"
			  : "100 %11.4f %11.4f %8.4f\n",
		pos[0], pos[1],
		mag + (fluxratio > 1.0/BIG? -2.5 * log10(fluxratio) : 99.0));
      }
    }

    NFPRINTF(OUTPUT,"");
    NPRINTF(OUTPUT, "%ld star%s generated\n", n, n>1?"s":"");
  }

/* Close catalogs */
  NFPRINTF(OUTPUT, "Closing catalogs...")
  for (p=0; p<npb; p++)
    if (fclose(catfile[p]))
      error(EXIT_FAILURE,"*ERROR*: can't close ", prefs.cat_name[p]);

/* Write summary output cluster file */
  if (ncluster)
    cluster_write(prefs.clusterlistout_name, clusters, ncluster);

/* Free memory */
  NFPRINTF(OUTPUT, "Freeing memory...")
  if (celsys)
    celsys_end(celsys);
  for (c=0; c<ncluster; c++)
    free(clusters[c]);
  if (ncluster)
    free(clusters);
  for (g=0; g<ngaltype; g++)
    galtype_end(galtype[g]);
  free(galtype);
  for (p=0; p<ngalsed; p++)
    sed_end(galsed[p]);
  for (p=0; p<nstarsed; p++)
    sed_end(starsed[p]);
  for (p=0; p<npb; p++)
    sed_end(pb[p]);
  sed_end(refpb);
  sed_end(tau_i);
  if (backsed)
    sed_end(backsed);
  free(lfevol);

/* Processing end date and time */
  thetime2 = time(NULL);
  tm = localtime(&thetime2);
  sprintf(prefs.sdate_end,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_end,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);
  prefs.time_diff = counter_seconds() - dtime;

  return;
  }

/****** counter_seconds *******************************************************
PROTO	double counter_seconds(void)
PURPOSE	Count the number of seconds (with an arbitrary offset).
INPUT	-.
OUTPUT	Returns a number of seconds.
NOTES	Results are meaningful only for tasks that take one microsec or more.
AUTHOR	E. Bertin (IAP)
VERSION	24/09/2009
 ***/
double	counter_seconds(void)
  {
   struct timeval	tp;
   struct timezone	tzp;
   int			dummy;

  dummy = gettimeofday(&tp,&tzp);
  return (double) tp.tv_sec + (double) tp.tv_usec * 1.0e-6;
  }


