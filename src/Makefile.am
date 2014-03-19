#
#				Makefile.am
#
# src Makefile. Process this file with automake to generate a Makefile
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	Stuff
#
#	Copyright:		(C) 2005-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
#
#	License:		GNU General Public License
#
#	Stuff is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	Stuff is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with Stuff. If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		27/05/2013
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if USE_PLPLOT
#CPLOTSOURCE	= cplot.c
#endif

bin_PROGRAMS	= stuff
stuff_SOURCES	= clusters.c cosmo.c $(CPLOTSOURCE) galaxies.c igm.c lf.c \
		  main.c makeit.c misc.c prefs.c random.c sed.c stars.c \
		  clusters.h cosmo.h define.h galaxies.h globals.h igm.h key.h \
		  lf.h preflist.h prefs.h random.h sed.h stars.h threads.h \
		  types.h
DATE=`date +"%Y-%m-%d"`
