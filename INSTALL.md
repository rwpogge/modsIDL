# Installation Instructions

## System Requirements

The modsIDL pipeline was developed on a Linux machine with 8 Intel
Core i7 3.4GHz processors and 8Gb of RAM.  It has also been tested 
and run on a MacBook Air with 2GHz Intel Core i7 and 8Gb of RAM.  
We recommend a system with at least a 3GHz processor and 8Gb of RAM.
The modsIDL pipeline has also been run successfully on machines with 
slower processors and less memory than this recommended limit.  However, 
given the large size of processed MODS files (~300Mb after sky 
subtraction), a minimum of 4Gb of RAM is recommended.

You need your own IDL license.  We developed and tested
modsIDL with IDL version 8.1. 

For Python we recommend the latest Anaconda distribution of Python, which we
have found works best for most Linux distributions and Mac OS/X computers. Currently,
however, the modsIDL code is only compatible with Python 2.7 (a Python 3 version
will be hopefully released later in 2019):

 > [Anaconda Python Distribution](https://www.anaconda.com/).

The modsIDL programs are incorporated into the XIDL software developed by
Jason X. Prochaska and Joe Hennawi and others:

 > [XIDL GitHub Repository](https://github.com/profxj/xidl)

This makes use of the numerous tools that have already been tested by numerous 
projects.  This version of modsIDL was developed using the SVN installation 
of XIDL from the fall of 2012(!). We have not been tracking the on-going
development of XIDL, so therefore provide a frozen copy of
the entire XIDL package that was used in developing modsIDL.

Advanced XIDL users may wish to install the MODS scripts into their personal versions
of XIDL.  We only recommended this path for those users who are familiar with the
structure of XIDL and are aware that changes may have occurred in XIDL that
may not work with the current pipeline. For users who wish to employ this 
option, the MODS files are located in the appropriate directory under the 
Longslit branch of XIDL.

## Unpacking 

Unpack the tarball (e.g., modsIDL_vX.Y.tgz) in the usual way

> tar -xvzpf modsIDL_vX.Y.tgz

This will create the modsIDL/ directory and fill it with the requisite
programs and support files. 

The modsIDL directory has five (5) subdirectories:

<pre>
   coyote/    -  a version of the Coyote IDL Library
   idlspec2d/ -  IDL package built by  David J. Schlegel & Scott Burles
   idlutils/  -  (djs) IDL package built and distributed by David J. Schlegel
   xidl/      -  the XIDL package built and distributed by Jason X. Prochaska
   Docs/      -  copy of this manual (PDF) and any supplementary documentation
</pre>

## Installation
	
To install the modsIDL package, the modsIDL directory must be in your IDL
path.  In general, most IDL codes are kept together in a top level directory
such as ~/myidl.  All code is then placed under this directory.  The user 
will also need to define the necessary environment variables and update their
IDL_PATH as necessary.  The following is recommended to be added to your .cshrc 
file:
<pre>
   if ( -d ~/myidl/modsIDL ) then
      setenv IDLUTILS_DIR ~/myidl/modsIDL/idlutils
      setenv IDLSPEC2D_DIR ~/myidl/modsIDL/idlspec2d
      setenv XIDL_DIR ~/myidl/modsIDL/xidl
      setenv LONGSLIT_DIR $XIDL_DIR/Spec/Longslit
      setenv IDL_PATH +$IDL_DIR\/lib:+$IDL_DIR\/examples:+~/myidl:+pro
   endif
</pre>
Or alternatively, add the following to your .bash file:
<pre>
   if [ -d ~/myidl/modsIDL ] then
      export IDLUTILS_DIR=~/myidl/modsIDL/idlutils
      export IDLSPEC2D_DIR=~/myidl/modsIDL/idlspec2d
      export XIDL_DIR=~/myidl/modsIDL/xidl
      export LONGSLIT_DIR=$XIDL_DIR/Spec/Longslit
      export IDL_PATH=+$IDL_DIR\/lib:+$IDL_DIR\/examples:+~/myidl:+pro
   endif
</pre>
To verify that modsIDL has been included in your path you can type:
<pre>
   env | grep –i modsIDL
</pre>
At this point you are ready to use the modsIDL pipeline.

## Updates

At present updates must be obtained by checking the modsIDL github repository
and updating by hand. We are discussing how to (or whether to) incorporate later 
versions of modsIDL into the GitHub distribution of XIDL, but have not yet
committed to a particular course of action.
