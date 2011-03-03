/*
*				preflist.h
*
* Configuration keyword definitions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	Stuff
*
*	Copyright:		(C) 1999-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		03/03/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "key.h"
#ifndef _GALAXIES_H_
#include "galaxies.h"
#endif
#ifndef _STARS_H_
#include "stars.h"
#endif

#ifdef	USE_THREADS
#define	THREADS_PREFMAX THREADS_NMAX
#else
#define	THREADS_PREFMAX 65535
#endif

int idummy;

pkeystruct key[] =
 {
  {"BULGE_FRACTION", P_FLOATLIST, prefs.gal_bt, 0,0, 0.0, 1.0,
   {""}, 1, GAL_MAXNTYPE, &prefs.ngal_bt},
  {"CALIBSED_REF", P_STRING, prefs.refcalibsed_name},
  {"CALIBSED_OBS", P_STRINGLIST, prefs.pbcalibsed_name, 0,0, 0.0,0.0,
   {""}, 1, SED_MAXNPB, &prefs.npbcalibsed_name},
  {"CATALOG_NAME", P_STRINGLIST, prefs.cat_name, 0,0, 0.0,0.0,
   {""}, 1, SED_MAXNPB, &prefs.ncat_name},
  {"CLUSTER_LIST", P_STRING, prefs.clusterlist_name},
  {"CLUSTER_LISTOUT", P_STRING, prefs.clusterlistout_name},
  {"COLLECT_AREA", P_FLOATLIST, prefs.area, 0,0, 0.0, BIG,
   {""}, 1, SED_MAXNPB, &prefs.narea},
  {"DATA_DIRECTORY", P_STRING, prefs.datadir_name},
  {"DISK_BETA", P_FLOAT, &prefs.gal_dbeta, 0,0, -100,-0.001},
  {"DISK_EXTINCT", P_FLOATLIST, prefs.gal_iextinc, 0,0, 0.0, BIG,
   {""}, 1, GAL_MAXNTYPE, &prefs.ngal_iextinc},
  {"DISK_REVOL", P_FLOAT, &prefs.gal_drevol, 0,0, -10.0,10.0},
  {"DISK_RSTAR", P_FLOAT, &prefs.gal_drstar, 0,0, 0.0, 1e12},
  {"DISK_SIGMAL", P_FLOAT, &prefs.gal_dsigmalambda, 0,0, 0.0, 1e12},
  {"DISTANCE_STEP", P_FLOAT, &prefs.integ_zstep, 0,0, 0.0, 10000.0},
  {"EXTINCT_NAME", P_STRING, prefs.extinct_name},
  {"GAIN", P_FLOATLIST, prefs.gain,  0,0, 0.0, BIG,
   {""}, 1, SED_MAXNPB, &prefs.ngain},
  {"GALACTIC_COORDS", P_FLOAT, prefs.galcoord, 0,0, 0.0, 360.0,
   {""}, 2,2, &idummy},
  {"H0", P_FLOAT, &prefs.h0, 0,0, 1.0, 1000.0},
  {"HUBBLE_TYPE", P_FLOATLIST, prefs.gal_hubtype, 0,0, -6.0, 10.0,
   {""}, 1, GAL_MAXNTYPE, &prefs.ngal_hubtype},
  {"IMAGE_HEIGHT", P_INTLIST, prefs.height, 1, 1000000000, 0.0,0.0,
   {""}, 1, SED_MAXNPB, &prefs.nheight},
  {"IMAGE_WIDTH", P_INTLIST, prefs.width, 1, 1000000000, 0.0,0.0,
   {""}, 1, SED_MAXNPB, &prefs.nwidth},
  {"INCLUDE_STARS", P_BOOL, &prefs.starflag},
  {"LENS_KAPPA", P_FLOAT, &prefs.lens_kappa, 0,0, -1.0, 1.0},
  {"LENS_GAMMA", P_FLOATLIST, prefs.lens_gamma, 0,0, -1.0, 1.0,
   {""}, 2, 2, &prefs.nlens_gamma},
  {"LF_ALPHA", P_FLOATLIST, prefs.lf_alpha, 0,0, -100.0, 10.0,
   {""}, 1, GAL_MAXNTYPE, &prefs.nlf_alpha},
  {"LF_MAGLIMITS", P_FLOATLIST, prefs.lf_mlim, 0,0, -100.0, 100.0,
   {""}, 2,2, &idummy},
  {"LF_MSTAR", P_FLOATLIST, prefs.lf_mstar, 0,0, -100.0, 100.0,
   {""}, 1, GAL_MAXNTYPE, &prefs.nlf_mstar},
  {"LF_MSTAREVOL", P_FLOATLIST, prefs.lf_mstarevol, 0,0, -10.0, 10.0,
   {""}, 1, GAL_MAXNTYPE, &prefs.nlf_mstarevol},
  {"LF_PHISTAR", P_FLOATLIST, prefs.lf_phistar, 0,0, 0.0, 1e12,
   {""}, 1, GAL_MAXNTYPE, &prefs.nlf_phistar},
  {"LF_PHISTAREVOL", P_FLOATLIST, prefs.lf_phistarevol, 0,0, -10.0, 10.0,
   {""}, 1, GAL_MAXNTYPE, &prefs.nlf_phistarevol},
  {"MAG_LIMITS", P_FLOATLIST, prefs.maglim, 0,0, -30.0, 50.0,
   {""}, 2,2, &idummy},
  {"OBSDETECT_TYPE", P_KEYLIST, prefs.obsdetect_type, 0,0, 0.0,0.0,
   {"PHOTONS", "ENERGY", ""}, 1, SED_MAXNPB, &prefs.nobsdetect_type},
  {"OMEGA_LAMBDA", P_FLOAT, &prefs.omegal, 0,0, 0.0, 10.0},
  {"OMEGA_M", P_FLOAT, &prefs.omegam, 0,0, 0.0, 10.0},
  {"NTHREADS", P_INT, &prefs.nthreads, 0, THREADS_PREFMAX},
  {"PASSBAND_OBS", P_STRINGLIST, prefs.pb_name, 0,0, 0.0,0.0,
   {""}, 1, SED_MAXNPB, &prefs.npb_name},
  {"PASSBAND_REF", P_STRING, prefs.refpb_name},
  {"PIXEL_SIZE", P_FLOATLIST, prefs.pixscale, 0,0, 1e-3,1e3,
   {""}, 1, SED_MAXNPB, &prefs.npixscale},
  {"REFDETECT_TYPE", P_KEY, &prefs.refdetect_type, 0,0, 0.0,0.0,
   {"PHOTONS", "ENERGY", ""}},
  {"SED_BACKGROUND", P_STRING, prefs.backsed_name},
  {"SED_GALAXIES", P_STRINGLIST, prefs.gal_sedname, 0,0, 0.0,0.0,
   {""}, 1, GAL_MAXNSED, &prefs.ngal_sedname},
  {"SED_STARS", P_STRINGLIST, prefs.star_sedname, 0,0, 0.0,0.0,
   {""}, 1, STAR_MAXNSED, &prefs.nstar_sedname},
  {"SEDINDEX_BULGE", P_INTLIST, prefs.gal_bsedno, 1,GAL_MAXNSED+1, 0.0, 0.0,
   {""}, 1, GAL_MAXNTYPE, &prefs.ngal_bsedno},
  {"SEDINDEX_DISK", P_INTLIST, prefs.gal_dsedno, 1,GAL_MAXNSED+1, 0.0, 0.0,
   {""}, 1, GAL_MAXNTYPE, &prefs.ngal_dsedno},
  {"SEED_GALAXIES", P_INT, &prefs.galseed, 0, 0x7fffffffL},
  {"SEED_STARS", P_INT, &prefs.starseed, 0, 0x7fffffffL},
  {"SPHEROID_REVOL", P_FLOAT, &prefs.gal_brevol, 0,0, -10.0,10.0},
  {"SPHEROID_RKNEE", P_FLOAT, &prefs.gal_brknee, 0,0, 0.0, 1e12},
  {"SPHEROID_MKNEE", P_FLOAT, &prefs.gal_bmknee, 0,0, -100.0,100.0},
  {"VERBOSE_TYPE", P_KEY, &prefs.verbose_type, 0,0, 0.0,0.0,
   {"QUIET","NORMAL","FULL",""}},

  {""}
 };

char			keylist[sizeof(key)/sizeof(pkeystruct)][32];
static const char	notokstr[] = {" \t=,;\n\r\""};

char *default_prefs[] =
 {
"# Default configuration file for " BANNER " " MYVERSION,
"# EB " DATE,
"#",
" ",
"#--------------------------------- Image -------------------------------------",
" ",
"CATALOG_NAME    g.list,r.list,i.list   # output catalog file name(s)",
"IMAGE_WIDTH     2048            # width of the simulated frame",
"IMAGE_HEIGHT    2048            # height of the simulated frame",
"PIXEL_SIZE      0.2             # pixel size (arcsec)",
"MAG_LIMITS      16.0,28.0       # allowed range of apparent magnitudes",
" ",
"#------------------------ Zero-points and Background -------------------------",
" ",
"COLLECT_AREA    1.0             # Effective collecting area(s) (in m2)",
"GAIN            1.0             # Detector conversion factor(s) (in e-/ADU)",
"SED_BACKGROUND  zodiacal        # Background SED",
" ",
"#-------------------------------- Passbands ----------------------------------",
" ",
"PASSBAND_REF    couch/Bj        # Reference passband",
"CALIBSED_REF    Vega            # SED for ref. passband mag.system (AB or Vega)",
"REFDETECT_TYPE  ENERGY          # Ref. detector type: \"PHOTONS\" or \"ENERGY\"",
"PASSBAND_OBS    megaprime/g,megaprime/r,megaprime/i     # Observed passband(s)",
"CALIBSED_OBS    AB              # SED(s) for obs. mag.system(s) (AB or Vega)",
"OBSDETECT_TYPE  PHOTONS         # Obs. detector type: \"PHOTONS\" or \"ENERGY\"",
" ",
"#------------------------------- Cosmology -----------------------------------",
" ",
"H0              70.0            # Hubble constant (km.s-1.Mpc-1)",
"OMEGA_M         0.3             # Matter density in units of critical density",
"OMEGA_LAMBDA    0.7             # Cosmol constant in units of critical density",
"*DISTANCE_STEP   5.0             # Integration step along z (h-1.Mpc)",
" ",
"#----------------------- Spectral Energy Distributions -----------------------",
" ",
"SED_GALAXIES    E,Sbc,Scd,Im    # SEDs for galaxy components",
" ",
"SEDINDEX_BULGE  1,1,1,1,1,1     # bulge SED indices in SED_GALAXIES",
"SEDINDEX_DISK   1,1,2,2,3,4     # disk SED indices in SED_GALAXIES",
" ",
"#----------------------------- Galaxy types ----------------------------------",
" ",
"HUBBLE_TYPE     -6.0,-2.0,2.0,4.0,6.0,10.0      # Simien & deVaucouleurs 86",
"BULGE_FRACTION  1.0,0.57,0.32,0.16,0.049,0.0    # Simien & deVaucouleurs 86",
"                                # B/T in ref. band for each galaxy component",
" ",
"#------------------------ Galaxy luminosity functions ------------------------",
" ",
"LF_PHISTAR      4.95e-3,4.95e-3,7.2e-3,5.0e-3,1.2e-3,1.2e-3",
"                                # Schechter's phi* density parameter (h3.Mpc-3)",
"LF_MSTAR        -19.58,-19.58,-19.53,-19.17,-19.15,-19.15",
"                                # Schechter's M* absolute magnitude",
"LF_ALPHA        -0.54,-0.54,-0.99,-1.24,-1.50,-1.50",
"                                # Schechter's alpha parameter",
"LF_MAGLIMITS    -27.0,-13.0     # bounds to the luminosity function",
" ",
"LF_PHISTAREVOL  -1.7,-1.7,-1.2,-1.2,-1.2,1.9",
"                                # Density evolution factor dln(phi*)/dln(1+z)",
"LF_MSTAREVOL    -1.0,-1.0,-1.0,-1.0,-1.0,-1.5",
"                                # Luminosity evolution factor dM*/dln(1+z)",
" ",
"#------------------------------- Extinction ----------------------------------",
" ",
"EXTINCT_NAME    LMC.ext         # Internal extinction law",
"DISK_EXTINCT    0.0,0.75,1.23,1.47,1.47,1.23",
"                                # de Vaucouleurs' alpha (extinction parameter)",
" ",
"*#------------------------------- Disk sizes ---------------------------------",
"*",
"*DISK_BETA       -0.214          # Beta parameter of the disk size distribution",
"*DISK_RSTAR      3.85            # re* parameter of the disk size distribution",
"*                                # (h-1.kpc) (see de Jong & Lacey 2000)",
"*DISK_REVOL      -0.40           # re* evolution factor dln(re*)/dln(1+z)",
"*                                # (see e.g. Trujillo et al. 2006)",
"*DISK_SIGMAL     0.36            # sigma_lambda parameter of the disk size",
"*                                # distribution (see de Jong & Lacey 2000)",
"*",
"*#----------------------------- Spheroid sizes -------------------------------",
"*",
"*SPHEROID_MKNEE  -20.0           # Abs. mag. at knee (Binggeli et al. 1984)",
"*SPHEROID_RKNEE  1.58            # re* parameter of spheroid size distribution",
"*                                # (h-1.kpc) (see Binggeli et al. 1984)",
"*SPHEROID_REVOL  -0.50           # re* evolution factor dln(re*)/dln(1+z)",
"*                                # (see e.g. Trujillo et al. 2006)",
"*",
"*#--------------------------- Galaxy clusters ---------------------------------",
"*",
"*CLUSTER_LIST    NONE            # ASCII File containing a list of clusters:",
"*                                # x(pixels) y(pixels) z M beta rc(h.^-1 Mpc)...",
"*                                # ... rmax(h^-1 Mpc) sig_v(h^-1 km/s)",
"*CLUSTER_LISTOUT cluster_out.list # output file containing the list of clusters:",
"*                                # x(pixels) y(pixels) z M beta rc(h.^-1 Mpc)...",
"*                                # ... rmax(h^-1 Mpc) sig_v(h^-1 km/s) M_real...",
"*                                # ... Ngal rc('') rmax('') ...",
"*",
"*#------------------------------- Lensing -------------------------------------",
"*LENS_KAPPA      0.0             # weak lensing (constant) convergence parameter",
"*LENS_GAMMA      0.0,0.0         # weak lensing (constant) shear parameters ",
"*",
"#------------------------------ Stellar field --------------------------------",
" ",
"INCLUDE_STARS   N               # allow addition of a stellar field?",
"GALACTIC_COORDS 270,30.0        # galactic coordinates",
"*SED_STARS       Vega            # stellar SEDs",
" ",
"#----------------------------- Random Seeds ----------------------------------",
" ",
"SEED_STARS      0               # random seed for the stellar field (0=time)",
"SEED_GALAXIES   0               # random seed for the galaxy field (0=time)",
" ",
"#----------------------------- Miscellaneous ---------------------------------",
" ",
"DATA_DIRECTORY  " DATA_DIR,
"                                # Top of directory tree containing Stuff data",
"VERBOSE_TYPE    NORMAL          # \"QUIET\",\"NORMAL\", \"LOG\" or \"FULL\"",
#ifdef USE_THREADS
"NTHREADS        0               # Number of simultaneous threads for the SMP",
"                                # version of " BANNER " (0 = automatic)",
#else
"NTHREADS        1               # 1 single thread",
#endif
""};

