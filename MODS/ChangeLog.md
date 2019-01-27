# modsIDL Change Log

Directory: xidl/Spec/Longslit/pro/LBT/MODS

Make note of all changes to scripts in this directory in reverse
chronological order (recent changes push down old).  Note the
script(s), briefly describe the change, and remember to date and
"sign" the changes.

## 2019 Jan 26: Version 1.0 - First Binocular Release [rwp/osu]

First release with a verified set of Sudoku transforms derived for
the MODS2 grating spectroscopy modes.

Version was tested with MODS2 2-channel grating spectra which
uncovered a number of bugs related to MODS2 reduction that are
addressed in this release.

NOTE: Not covered in this release:
 * MODS2 Prism Mode transforms
 * MODS2 quicklook reduction (requires substantial rewrites)

### Specific Changes:
<pre>    
      Transforms/blue_map.pro
      Transforms/red_map.pro
          Added MODS2 Blue & Red Grating sudoku transform coefficients

      mods_mask.pro
      mods_wstruct.pro
         Changed tests of the TELESCOP header keyword that previous
         tested for '*LBT-SX*' to test for '*LBT-* instead to properly
         detect the MODS2 case.

      mods_plan.pro
         Added a missing test for case INSTRUME=MOD2* that resulted in not
         building plan files correctly for MODS2 reductions.

      mods_extract_singlechan.pro 
      mods_fluxstand_singlechan.pro
      mods_skyfit2d_singlechan.pro
         Fixed bug for the MODS2R case where it was setting otherwise undefined
         variable red_channel_code instead of channel_code.

      modsadcalc.pro
         Fixed bug in case where the target declination was greater than
         the site latitude, which resulted in the wrong sign for the 
         atmospheric refraction deflection along the slit axis, resulting
         in the wrong spectral trace shape.

      mods_extract_dualchan.pro
         Fixed bugs in MODS2 channel detection.

      Calib_Lib/facility_masks/
         Added mms and gerber files for the LS5x60x2.4 facility slit mask:
         ls5x60x2.4.mms/.gbr 
</pre>