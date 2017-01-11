#
#				Makefile.am
#
# Main Makefile. Process this file with automake to generate a Makefile
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	Stuff
#
#	Copyright:		(C) 2005-2017 IAP/CNRS/UPMC
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
#	Last modified:		11/01/2017
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ACLOCAL_AMFLAGS		= -I m4
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= src man
EXTRA_DIST		= data doc autogen.sh \
			  AUTHORS BUGS ChangeLog COPYRIGHT \
			  HISTORY INSTALL LICENSE README.md THANKS
RPM_ROOTDIR		= `rpmbuild -E %_topdir`
RPM_SRCDIR		= $(RPM_ROOTDIR)/SOURCES
dist-hook:
	rm -rf `find $(distdir) -name .svn`

rpm:	dist
	cp -f $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean --nodeps $(PACKAGE_TARNAME).spec

rpm-icc:	dist
	cp -f $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	USE_ICC="1" rpmbuild -ba --clean --nodeps $(PACKAGE_TARNAME).spec

rpm-best:	dist
	cp -f $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	USE_BEST="1" rpmbuild -ba --clean --nodeps $(PACKAGE_TARNAME).spec

install-data-local:
	if test -e $(DESTDIR)$(pkgdatadir); then \
	  $(RM) -rf $(DESTDIR)$(pkgdatadir); \
	fi;
	if ! test -d $(DESTDIR)$(pkgdatadir); then \
	  $(mkinstalldirs) $(DESTDIR)$(pkgdatadir); \
	fi;
	cp -rf data $(DESTDIR)$(pkgdatadir)/data;
	$(RM) -rf `find $(DESTDIR)$(pkgdatadir) -name .svn`

