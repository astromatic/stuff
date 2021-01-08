#
#				configure.ac
#
# Process this file with autoconf to produce a configure script.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	Stuff
#
#	Copyright:		(C) 2005-2021 IAP/CNRS/SorbonneU
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
#	along with Stuff.  If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		08/01/2021
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# First, disable the annoying config.cache
define([AC_CACHE_LOAD],)
define([AC_CACHE_SAVE],)

# This is your standard Bertin source code...
AC_INIT(Stuff, 2.0.1, [https://github.com/astromatic/stuff],
	stuff, [http://astromatic.net/software/stuff])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR(src/makeit.c)
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE
date=`date +%Y-%m-%d`
date1=`date -R`
date2=`date +"%a %b %d %Y"`
date3=`date +"%B %Y"`
AC_DEFINE_UNQUOTED(DATE, "$date", [Archive creation date])
AC_SUBST(PACKAGER, "AstrOmatic")
AC_SUBST(DATE1, "$date1")
AC_SUBST(DATE2, "$date2")
AC_SUBST(DATE3, "$date3")

# Include macros
sinclude(acx_prog_cc_optim.m4)
sinclude(acx_pthread.m4)
sinclude(acx_urbi_resolve_dir.m4)

# Provide a special option for setting the package release number
AC_ARG_WITH(release,
	[AS_HELP_STRING([--with-release=<release number>],
	[set the package release number (default = 1)])],
	[],
	[with_release="no"])
if test "$with_release" = "no"; then
  with_release=1
fi
AC_SUBST(PACKAGE_RELEASE, "$with_release")

# Display pakage and version number
AC_MSG_RESULT([******** Configuring:  $PACKAGE_NAME $PACKAGE_VERSION - $PACKAGE_RELEASE ($date) ********])

# Provide special option for choosing automatically the compilation flags.
AC_MSG_CHECKING([if compilation flags are set automatically])
AC_ARG_ENABLE(auto-flags,
	[AS_HELP_STRING([--enable-auto-flags],
	[Let the configure script choose the compilation flags (default = no)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Provide special option for the Linux Intel C compiler
AC_MSG_CHECKING([whether the INTEL compiler is enabled])
AC_ARG_ENABLE(icc,
	[AS_HELP_STRING([--enable-icc],
	[Use the Intel compiler (default = no)])],
        CC="icc"
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Checks for programs.
AC_LANG(C)

AC_SEARCH_LIBS([strerror],[cposix])
if test "$enable_auto_flags" = "yes"; then
  CFLAGS=""
  LDFLAGS=""
  ACX_PROG_CC_OPTIM
fi
AC_PROG_RANLIB
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB(m, sin)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for INTEL math header files.
if test "$enable_icc" = "yes"; then
  AC_CHECK_HEADERS(mathimf.h)
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([getenv floor pow sqrt strchr strcspn strtol lgamma \
	drand48 lrand48 srand48 gettimeofday])

# Set the data directory to a true absolute path
datadir2=$(URBI_RESOLVE_DIR([$datadir]))

# Provide a special option for the default XSLT URL
with_stuff_datadir="$datadir2/$PACKAGE_TARNAME/data"
AC_ARG_WITH(stuff-datadir,
        [AS_HELP_STRING([--with-stuff-datadir=<default path for Stuff data>],
        [Provide an alternative path for Stuff data])])

AC_DEFINE_UNQUOTED([DATA_DIR], "$with_stuff_datadir",
	[Default path for Stuff data])

# Set flags for multithreading
n_pthreads=1024
AC_ARG_ENABLE(threads,
	[AS_HELP_STRING([--enable-threads@<:@=<max_number_of_threads>@:>@],
	[Enable multhreading (on with up to 1024 threads by default)])],
    if test "$enableval" = "no"; then
      use_pthreads="no"
    else
      use_pthreads="yes"
      if test "$enableval" != "yes"; then
        n_pthreads=$enableval
      fi
    fi,
    use_pthreads=yes
    )

# Provide special option for profiling
AC_MSG_CHECKING([for profiler mode])
AC_ARG_ENABLE(profiling,
	[AS_HELP_STRING([--enable-profiling],
	[Enable special mode for profiling (default = no)])]
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Enable linking options for making the executable as portable as possible.
AC_MSG_CHECKING([best linking option])
AC_ARG_ENABLE(best-link,
	[AS_HELP_STRING([--enable-best-link],
	[Choose the right combination of static and dynamic linking to make \
the executable as portable as possible (default = no)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

################# Actions to complete in case of multhreading ################
AC_DEFINE_UNQUOTED(THREADS_NMAX, $n_pthreads,[Maximum number of POSIX threads])
if test "$use_pthreads" = "yes"; then
  AC_MSG_CHECKING([for multithreading])
  AC_MSG_RESULT([maximum of $n_pthreads thread(s)])
  AC_DEFINE(USE_THREADS, 1, [Triggers multhreading])
# CC, CFLAGS and LIBS are system and compiler-dependent
  ACX_PTHREAD
  CC="$PTHREAD_CC"
  [AM_CFLAGS="$AM_CFLAGS $PTHREAD_CFLAGS -D_REENTRANT"]
  LIBS="$LIBS $PTHREAD_LIBS"
fi
AM_CONDITIONAL(USE_THREADS, [test $use_pthreads = "yes"])

# Compile with profiling option
if test "$enable_profiling" = "yes"; then
  if test "$enable_icc" = "yes"; then
    AM_CFLAGS="$AM_CFLAGS -pq"
  else
    AM_CFLAGS="$AM_CFLAGS -pg"
  fi
  enable_best_link="no"
fi

# "Best" linking option
if test "$enable_best_link" = "yes"; then
  if test "$enable_icc" = "yes"; then
    AM_LDFLAGS="-static-intel -qopenmp-link static -shared-libgcc \
	-static-libtool-libs -no-intel-extensions -avoid-version $AM_LDFLAGS"
  else
    AM_LDFLAGS="-shared-libgcc -static-libtool-libs -avoid-version $AM_LDFLAGS"
  fi
fi

AC_SUBST(AM_CFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)

# Display compiler and linker flags
AC_MSG_RESULT([***************************************************************])
AC_MSG_RESULT([Compile cmdline:  $CC $AM_CPPFLAGS $CPPFLAGS $AM_CFLAGS $CFLAGS])
AC_MSG_RESULT([Link    cmdline:  $CC $AM_LDFLAGS $LDFLAGS $LIBS])
AC_MSG_RESULT([Default XSLT URL: $xsl_url])
AC_MSG_RESULT([***************************************************************])

AC_CONFIG_FILES([
Makefile
stuff.spec
debian/Makefile
debian/changelog
doc/Makefile
doc/src/conf.py
man/Makefile
man/stuff.1
src/Makefile])
AC_OUTPUT
