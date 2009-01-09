/*
                                  makeit.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        Stuff
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       Main loop
*
*       Last modify:    09/02/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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

#include "define.h"
#include "globals.h"
#include "prefs.h"
#include "clusters.h"
#include "cosmo.h"
#include "galaxies.h"
#include "lf.h"
#include "random.h"
#include "sed.h"

/********************************** makeit ***********************************/
void	makeit(void)
  {
   clusterstruct	**clusters;
   galtypestruct	**galtype;
   galstruct		*gal;
   sedstruct		*pb[SED_MAXNPB],*galsed[GAL_MAXNSED],
			*starsed[STAR_MAXNSED],
			*calibsed, *refpb, *tau_i, *bsed,*dsed, *backsed;
   lfstruct		*lf, *lfevol;
   time_t		thetime, thetime2;
   struct tm		*tm;
   FILE			*(catfile[SED_MAXNPB]);
   static char		seddir_name[MAXCHAR], filterdir_name[MAXCHAR],
			extdir_name[MAXCHAR];
   static double	xscale[SED_MAXNPB], xoffset[SED_MAXNPB],
			yscale[SED_MAXNPB], yoffset[SED_MAXNPB];
   double		mabsmax, dn, bt, kb,kd,kt,
			width,height, widthmax,heightmax, x,y, xp,yp, domega,
			z, dz,zlow,zhigh, zmax, dlc, galdens;
   long			nsource;
   int			c,i,g,n,p, npb, ngaltype, ngalsed, nstarsed, ncluster,
			clustindex;

/* Processing start date and time */
  thetime = time(NULL);
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

/* First find the largest spread in x and y */
  widthmax = heightmax = -BIG;
  for (p=0; p<npb; p++)
    {
    if ((width = prefs.pixscale[p]*prefs.width[p])>widthmax)
      widthmax = width;
    if ((height = prefs.pixscale[p]*prefs.height[p])>heightmax)
      heightmax = height;
    }
  widthmax *= 1.1;
  heightmax *= 1.1;
/* Position offsets and scales for each passband */
  for (p=0; p<npb; p++)
    {
    xscale[p] = widthmax/prefs.pixscale[p];
    yscale[p] = heightmax/prefs.pixscale[p];
    xoffset[p] = (xscale[p]-prefs.width[p])/2.0;
    yoffset[p] = (yscale[p]-prefs.height[p])/2.0;
    }
  domega = widthmax*heightmax*ARCSEC*ARCSEC;

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
  NFPRINTF(OUTPUT, "Loading the calibration SED...")
  calibsed = sed_load(seddir_name, prefs.calibsed_name);

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
  pb_calib(pb, npb, refpb, calibsed, backsed);

/* Forget the calibration SED to gain memory */
  sed_end(calibsed);

/* Get the galaxy SEDs */
  NFPRINTF(OUTPUT, "Loading galaxy SEDs...")
  ngalsed = prefs.ngal_sedname;
  for (p=0; p<ngalsed; p++)
    galsed[p] = sed_load(seddir_name, prefs.gal_sedname[p]);

/* Get the star SEDs */
  NFPRINTF(OUTPUT, "Loading star SEDs...")
  nstarsed = prefs.nstar_sedname;
  for (p=0; p<nstarsed; p++)
    starsed[p] = sed_load(seddir_name, prefs.star_sedname[p]);

/* Calibrate the SEDs */
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
/*------ Compute the minimum absolute magnitude of detectable galaxies */
/*------ Linear k-corrections for bulge and disk */
        kb = bt>0.0? sed_kcor(galtype[g]->bsed, refpb, zlow) : 1.0;
        kd = bt<1.0? sed_kcor(galtype[g]->dsed, refpb, zlow) : 1.0;
        if ((kt=bt*(kb-kd)+kd) < 1/BIG)
          kt = 1/BIG;
        mabsmax = prefs.maglim[1]-cosmo_dmodul(zlow)+2.5*log10(kt);
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
        gal = gal_init(galtype[g], z, mabsmax, pb, npb);

/*------ Position */
        x = random_double();
        y = random_double();
/*------ Now in each observed passband */
        for (p=0; p<npb; p++)
          {
/*-------- Position in pixel coordinates*/
          xp = x*xscale[p]-xoffset[p];
          yp = y*yscale[p]-yoffset[p];
          fprintf(catfile[p],
		"200 %10.3f %10.3f %8.4f %5.3f "
		"%9.3f %5.3f %+7.2f %9.3f %5.3f %+7.2f %8.5f %+5.1f\n",
		xp,yp,
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
          gal = gal_init(galtype[g], z, mabsmax, pb, npb);
          clusters[c]->real_mabs += exp(-0.921034*gal->mabs);
/*-------- Now in each observed passband */
          for (p=0; p<npb; p++)
            {
/*-------- Position in pixel coordinates*/
            xp = clusters[c]->x + x/prefs.pixscale[p];
            yp = clusters[c]->y + y/prefs.pixscale[p];
            fprintf(catfile[p],
		"200 %10.3f %10.3f %8.4f %5.3f "
		"%9.3f %5.3f %+7.2f %9.3f %5.3f %+7.2f %8.5f\n",
		xp,yp,
		gal->mag[p], gal->bt[p],
		gal->bsize, gal->bflat, gal->bposang,
		gal->dsize, gal->dflat, gal->dposang, z);
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
    }

  NFPRINTF(OUTPUT,"");
  NPRINTF(OUTPUT, "%ld galax%s generated\n", nsource, nsource>1?"ies":"y");

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
  prefs.time_diff = difftime(thetime2, thetime);

  return;
  }

