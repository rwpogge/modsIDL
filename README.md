# modsIDL - MODS Spectral Data Reduction Pipeline

IDL Data Reduction Pipeline for the Multi-Object Double Spectrographs at the Large Binocular Telescope Observatory

<B>This is work in progress, not yet ready for release or download!!!</B>

## Overview

modsIDL is a suite of programs built on Jason Prochaska's XIDL package for the reduction of MODS long-slit and multi-slit spectroscopy. Long-slit spectra are considered to be a subset of multi-slit spectra. Both grating and prism mode spectra can be reduced using modsIDL, but for the latter please contact us directly for advice, as prism mode is not as widely used and those parts of the pipeline are still in an advanced stage of development.

Sky subtraction, critical for faint-object spectroscopy, is accomplished using the 2D B-Spline algorithm developed by Dan Kelson [2003, PASP, 115, 688](http://adsabs.harvard.edu/abs/2003PASP..115..688K) . The modsIDL pipeline maintains the original detector pixel grid throughout the steps leading up to 1D science spectrum extraction without geometric rectification. This allows us to robustly estimate and propagate errors throughout the reduction process, resulting in internally-consistent and robust 2D and 1D error arrays for each science target spectrum.

The ultimate output of modsIDL is a Multi-Extension FITS format file containing row-stacked 1D wavelength- and flux-calibrated sky-subtracted spectra with associated error arrays, sky spectra, and other related data. 

## User's Manual

Please download and read the user's manual on this repository for instructions on how to use and install the modsIDL package. This manual describes essential details on how the programs work to take you from raw MODS CCD images to wavelength and flux calibrated sky-subtracted 1D spectra. and the 2D reduction steps required. Be sure to read it carefully before using these programs to process MODS raw images. 

## Installation

See [INSTALL.md](INSTALL.md) for installation instructions.

## System Requirements

 * Linux (CentOS 5++, FC++, RHE9++), but note you need different 32- and 64-bit distributions unless you are prepared to recompile. 
 * Mac OS X (10.7 or later) 
 * Computer should have at least a 3GHz Processor and 8Gb of RAM. More cores and more RAM is always better. 
 * IDL v8.1 or later (requires a license) installed. 
 * [modsCCDRed](https://github.com/rwpogge/modsCCDRed) Basic 2D reduction programs (python) 
 * Python 2.7, preferrably the free Anaconda distribution from Continuum Analytics. You will at least need numpy and astropy which come as part of Anaconda. It will not work with Python 3 yet. We cannot guarantee proper operation if you do not use Anaconda. 

## Example Data

If you are running the worked example in chapter 5 of the manual, this is where you download the raw \
and processed data files:

> [modsIDL Example Data](http://www.astronomy.ohio-state.edu/MODS/Software/modsIDL/Data/index.html)

Please read the webpage at this link before downloading these files, as many of the packages are quite 
large (few Gb) and you need to know what you need before committing. 

## Terms & Conditions of Support

The modsIDL package is provided as-is, with no warranty or offer of individual user support. We are willing to fix bugs, answer basic questions, and provide updates, but individual instruction or hand-holding is beyond our limited human resources to provide.

We expect that users have a basic working knowledge of IDL at the level of knowing the basic command syntax and workflow of IDL pre-written procedures (a knowledge of IDL programming is helpful but not required for basic use). Similarly, we expect that users are well-acquainted with the basic principles of CCD spectroscopy in general and MODS spectra in particular. 

## Acknowledging modsIDL 

modsIDL was developed for reducing data obtained with the MODS1 and MODS2 instruments at the Large Binocular Telescope
Observatory.  The MODS were built with with major support provided by grants from the U.S. National Science Foundation's
Division of Astronomical Sciences  Advanced Technologies and Instrumentation (AST-9987045), the NSF/NOAO TSIP Program,
and matching funds provided by the Ohio State University Office of Research and the Ohio Board of Regents. Additional 
support was provided by NSF Grant AST-1108693.

If you used modsIDL in your research work, we ask that you follow emerging
[software citation principles](https://doi.org/10.7717/peerj-cs.86) being adopted by the astronomical community
to ensure the proper citation of software in scientific publications. 

The modsIDL pipeline was developed independently of the MODS instrument project as part of NSF Grant AST-1108693 which
supported the Chemical Abundances of Spirals (CHAOS) project (PIs: Evan Skillman, UM, and Rick Pogge, OSU). It grew out of the
need of the CHAOS project for an automated multi-object reduction and analysis pipeline for MODS data, and was intended from
the start to be made available to the rest of the LBTO community as one of the "deliverables" of this NSF project. The modsIDL
development project was led by Dr. Kevin Croxall, an OSU postdoc from 2012 until 2016.

If you publish MODS data reduced with the modsIDL pipeline, we ask that you add this line to the acknowledgments section of
your paper:

    This paper made use of the modsIDL spectral data reduction reduction pipeline developed in part with 
    funds provided by NSF Grant AST-1108693 and a generous gift from OSU Astronomy alumnus David G.
    Price through the Price Fellowship in Astronomical Instrumentation. 

A software DOI has been registered with zenodo.org:

 > [![DOI](https://zenodo.org/badge/167826611.svg)](https://zenodo.org/badge/latestdoi/167826611)

We and our funding agencies thank you for properly acknowledging their support in your papers. 

## Author & History

Dr. Kevin Croxall, the lead on the modsIDL project, was a postdoctoral research fellow at The Ohio State University
Department of Astronomy from 2012 until 2016, supported by funds provided by NSF Grant AST-1108693 and by the David G. Price 
Postdoctoral Fellowship in Astronomical Instrumentation. Kevin left OSU and Astronomy in May 2016 and is making a career in
data science. 

Support for modsIDL will continue on an informal basis as we can find resources. In the meantime, MODS instrument PI 
Prof. Richard Pogge (pogge.1@osu.edu) will address any bugs or questions as time permits. 

## Revision History

### Current Version: v1.0.1 - 2019 Feb 10 - First binocular release
