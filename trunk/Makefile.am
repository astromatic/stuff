# Main Makefile for Stuff
# Copyright (C) 2005-2009 Emmanuel Bertin.
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= man src
EXTRA_DIST		= data doc \
			  AUTHORS BUGS ChangeLog COPYRIGHT \
			  HISTORY INSTALL README THANKS \
			  acx_pthread.m4 acx_prog_cc_optim.m4
RPM_ROOTDIR		= `rpmbuild --nobuild -E %_topdir`
RPM_SRCDIR		= $(RPM_ROOTDIR)/SOURCES
dist-hook:
	rm -rf `find $(distdir) -name .svn`

rpm:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean --nodeps $(PACKAGE_NAME).spec

rpm-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	USE_ICC="1" rpmbuild -ba --clean --nodeps $(PACKAGE_NAME).spec

install-data-local:
	if test -e $(DESTDIR)$(pkgdatadir); then \
	  $(RM) -rf $(DESTDIR)$(pkgdatadir); \
	fi;
	if ! test -d $(DESTDIR)$(pkgdatadir); then \
	  $(mkinstalldirs) $(DESTDIR)$(pkgdatadir); \
	fi;
	cp -rf data $(DESTDIR)$(pkgdatadir)/data;
	$(RM) -rf `find $(DESTDIR)$(pkgdatadir) -name .svn`

