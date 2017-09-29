.. File Installing.rst

Installing the software
=======================

Hardware requirements
---------------------

|Stuff| runs in (ANSI) text-mode from a shell. A window system is
not necessary for basic operation.

The amount of memory required by |Stuff| is very low, a few tens of megabytes
at most. The current version of |Stuff| does not take advantage of multiple CPU
cores.

Obtaining |Stuff|
-----------------

For Linux users, the simplest way to have |Stuff| up and running is to install
the standard binary package the comes with your Linux distribution. Run, e.g.,
``apt-get stuff`` (on Debian) or ``dnf stuff`` (Fedora) and
|Stuff|, as well as all its dependencies, will automatically be installed. If
you decided to install the package this way you may skip the following and move
straight to the :ref:`next section <using_stuff>`.

However if |Stuff| is not available in your distribution, or to obtain the most
recent version, the |Stuff| source package can be downloaded from `the official
GitHub repository <https://github.com/astromatic/stuff>`_ . One may choose
`one of the stable releases <https://github.com/astromatic/stuff/releases>`_,
or for the fearless, `a copy of the current master development branch
<https://github.com/astromatic/stuff/archive/master.zip>`_.

Software requirements
---------------------

|Stuff| has been developed on `GNU/Linux <http://en.wikipedia.org/wiki/Linux>`_
machines and should compile on any
`POSIX <http://en.wikipedia.org/wiki/POSIX>`_-compliant system (this includes
|OSX|_ and `Cygwin <http://www.cygwin.com>`_ on |Windows|_, at the price of
some difficulties with the configuration).

Installation
------------

To install from the |GitHub| source package, you must first uncompress the
archive:

.. code-block:: console

  $ unzip stuff-<version>.zip

A new directory called :file:`stuff-<version>` should now appear at the current
location on your disk. Enter the directory and generate the files required by
the `autotools <http://en.wikipedia.org/wiki/GNU_Build_System>`_, which the
package relies on:

.. code-block:: console

  $ cd stuff-<version>
  $ sh autogen.sh

A :program:`configure` script is created. This script has many options, which
may be listed with the ``--help`` option:

.. code-block:: console

  $ ./configure --help

No options are required for compiling with the default GNU C compiler
(:program:`gcc`) if all the required libraries are installed at their default
locations:

.. code-block:: console

  $ ./configure

Compared to :program:`gcc` and the librairies above, using the
|Intel| compiler (:program:`icc`) can give the |Stuff|
executable a strong boost in performance, thanks to better vectorized code.
If :program:`icc` is installed on your system [#geticc]_ , you can take
advantage of it using

.. code-block:: console

  $ ./configure --enable-icc

Additionally, if the |Stuff| binary is to be run on a different machine
that does not have :program:`icc`  installed (e.g., a cluster
computing node), you must configure a partially statically linked executable
using

.. code-block:: console

  $ ./configure --enable-icc --enable-auto-flags --enable-best-link

In all cases, |Stuff| can now be compiled with

.. code-block:: console

  $ make -j

An :file:`src/stuff` executable is created. For system-wide installation, run
the usual

.. code-block:: console

  $ sudo make install

You may now check that the software is properly installed by simply
typing in your shell:

.. code-block:: console

  $ stuff

which will return the version number and other basic information (note that
some shells require the :program:`rehash` command to be run before making a
freshly installed executable accessible in the execution path).

.. [#mac_install] Mac OS X |.dmg|_ packages should be available soon.
.. [#geticc] The Linux versions of the |Intel| compiler and |MKL| are
   `available for free to academic researchers, students, educators and open
   source contributors <http://software.intel.com/qualify-for-free-software>`_.

.. include:: keys.rst

