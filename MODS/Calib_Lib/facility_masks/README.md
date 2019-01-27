# MODS Facility Slit Masks

Last Update: 2019 Jan 24 [housekeeping]

This directory contains the Gerber files for cutting copies of the
MODS facility slit masks.  For the facility long-slit masks we also
include a pseudo-MMS file describing the masks for use with multi-slit
reduction programs (this treats the segmented long-slit masks as a
kind of multi-slit mask, which they are).

The Gerber files are all in RS-274X (aka "Extended") Gerber Format as
described by the Ucamco Reference Document US-GXXX-TB-01-EN-E (Dec
2010).

File Extensions:
<pre>
   .gbr = Gerber RS-274X File
   .mms = MMS File (MODS Mask Specification File)
</pre>
See the MODS Manual for details of the facility long-slit specifications.

## Manifest:
<pre>
   ls5x60x0.3.xyz  --  0.3-arcsec segmented long-slit mask
   ls5x60x0.6.xyz  --  0.6-arcsec segmented long-slit mask
   ls5x60x0.8.xyz  --  0.8-arcsec segmented long-slit mask
   ls5x60x1.0.xyz  --  1.0-arcsec segmented long-slit mask
   ls5x60x1.2.xyz  --  1.2-arcsec segmented long-slit mask
   ls5x60x2.4.xyz  --  2.4-arcsec segmented long-slit mask
   ls60x5.xyz      --  5-arcsec x 60-arcsec spectrophotometric wide slit mask
   sudokuMask.xyz  --  9x9 Sudoko pinhole mask (MOS x,y,lambda Calibration)
   pinholeSlit.xyz --  Linear pinhole slit
   sieveMask.xyz   --  Pinhole Sieve Mask (imaging calibration)
</pre>
There are no mms versions of the pinhole and sieve mask gerber files.
