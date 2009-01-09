# Main Makefile for Stuff
# Copyright (C) 2005-2008 Emmanuel Bertin.
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= man src
EXTRA_DIST		= data doc \
			  AUTHORS BUGS ChangeLog COPYING \
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

rpm:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec

rpm-opteron:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall -march=opteron" rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec                                                               

rpm-athlon:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall -m32 -march=i686 -msse -mfpmath=sse -mtune=athlon" rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec

rpm-opteron-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc" TPXFLAGS="-O3 -axWP -ip -unroll" rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec                                                               

rpm-athlon-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc32" TPXFLAGS="-O3 -axWP -ip -unroll" rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec

icc:
	$(MAKE) CC="icc" CFLAGS="-O3 -axWP -ip -unroll -pthread" LIBS="-lm -lpthread"

debug:
	$(MAKE) CFLAGS="-O3 -funroll-loops -fomit-frame-pointer -Wall -m32 -g" LIBS="-L/usr/lib/ -lpthread -lm"

