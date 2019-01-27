;+
; NAME:
;   MKMMS_WS
;
; VERSION:
;   1.0 (Feb, 2014)
;
; PURPOSE:
;   Make an MMS file with the appropriate header info so that the
;   slits all have RA and DEC information to pass along to the 
;   mods-pipeline extraction software.
;
; REFERENCE:
;   MODS_MANUAL?
;
; CALLING SEQUENCE:
;   mkmms_ws,infile
;
; OUTPUT:
;   A file named 5arcsecSpectrophotometryWideSlit.mms is saved to the 
;   current directory.  This file contains all the information 
;   necessary for the MODS pipline to construct a wode-slit
;   structure.  The header information from the input file is
;   inserted to allow for the passing of coordinates.  The adopted 
;   name is what is looked for by the pipe-line for a WS 
;   observation.
;
; EXAMPLE:
;   mkmms_ws,mods1b.20130213.0016_otf.fits
;
; HISTORY:
;   2014-02-17 (KV Croxall): Written.
;-
;##############################################################################

;=========================================================================
;   stringad  -- convert sexigesimal RA and Dec to degrees
;       09-AUG-90 Version 1 written by Kerry McQuade
;       20-AUG-90 Put code to account for '-0' back in after it was
;               removed by someone.  E. Deutsch
;       17-JUL-95 Added support for coordinates separated by colons, e.g.
;               17:00:45.2 25:4:32.4, which IRAF uses.  E. Deutsch
;=========================================================================
pro stringad,coords,ra,dec
  On_error,2

  arg = N_params()
  if ( arg LT 1 ) then begin
    print,'Call: IDL> STRINGAD,coord_string,ra,dec'
    print,"e.g.: IDL> STRINGAD,'17 00 45.2  25 4 32.4',ra,dec"
    print," or : IDL> STRINGAD,'17:00:45.2  25:4:32.4',ra,dec"
    return
  endif

; Remove any gaps between '-' or '+' and numeral  
  I = strpos(coords,'+ ')
  if ( I GE 0 ) then strput,coords,'  +', I-1
  J = strpos(coords,'- ')
  if ( J GE 0 ) then strput,coords,'  -',J-1

; Replace colons with spaces - Added by Deutsch 7/17/95
  i=0
  while (i ne -1) do begin
    i=strpos(coords,':')
    if (i ne -1) then strput,coords,' ',i
    endwhile

  radec = getopt(coords,'F')
  if ( N_elements(radec) LT 6 ) then message, $
        'Coordinate format should be HR MIN SEC DEG DMIN DSEC'

  ra = ten(radec(0:2)*15.)      ;Convert to decimal degrees
  dec = ten(radec(3:5))

; Some formats write this: '12 34 15.33 -0 12 45.3'  Make this convert properly
  if ((strpos(coords,'-') NE -1) and (dec gt 0)) then dec = -dec
  if (arg LT 2) then print,'Decimal coords:   ',ra,dec

  return
  end

;=========================================================================
;   mkmms_ws - write the MMS file for a given MODS LONGSLIT spectrum
;=========================================================================
pro mkmms_ws,infile

  outfile = '5arcsecSpectrophotometryWideSlit.mms'
  openw,42,outfile,WIDTH=950

;read the input header
  in = mrdfits(infile,0,hdr)
  targalpha = strcompress(sxpar(hdr[*, 0], 'TELRA'), /rem)
  targdelta = strcompress(sxpar(hdr[*, 0], 'TELDEC'), /rem)
  coords = targalpha + ' ' + targdelta
  stringad,coords,objra,objdec
  remchar,targalpha,':'
  remchar,targdelta,':'

  targpa = strcompress(sxpar(hdr[*, 0], 'POSANGLE'), /rem) 
  targeq = strcompress(sxpar(hdr[*, 0], 'OBJEPOCH'), /rem)
  targorig = strcompress(sxpar(hdr[*, 0], 'OBSERVAT'), /rem)
  fitsinst = strcompress(sxpar(hdr[*, 0], 'INSTRUME'), /rem)
  fitsobj = strcompress(sxpar(hdr[*, 0], 'OBJNAME'), /rem)
  currenttime = systime(/utc)

;calculate the slit centers
;   Convenience function to evaluate the standard Cartesian 2D
;   coordinate system rotation.  Note that for applying this to
;   astronomical standard coordinates (xi,eta), the helicity of rotAng
;   has the opposite sign (e.g., compute xi,eta when rotating by
;   celestial position angle posAng, use rotAng=-posAng).
  xs = [0,0,0,0,0]
  ys = [0,-63,63,-126,126]
  rotAng = -targpa
  sinPA = sin(rotAng * !dpi/180)
  cosPA = cos(rotAng * !dpi/180)
  xr = xs*cosPA - ys*sinPA
  yr = xs*sinPA + ys*cosPA

;   Given standard coordinates (xi,eta) of an object and the celestial
;   coordinates of the field center (ra0,dec0), compute the celestial
;   coordinates (ra,dec) of the object on the celestial sphere.  This
;   geometric problem is described in W.W. Smart, Textbook on Spherical
;   Astronomy, Chapter XII, sections 160 and 161.
  xi_r = xr/3600.0 * !dpi/180
  eta_r = yr/3600.0 * !dpi/180
  ra0_r = objra * !dpi/180

  cosd0 = cos(objdec * !dpi/180)
  tand0 = tan(objdec * !dpi/180)
  denom = 1.0 - eta_r*tand0
  dra = atan(xi_r / (cosd0*denom))
  raslits = (ra0_r + dra) *180/!dpi
  decslits = 180/!dpi * (atan((cos(dra)*(eta_r+tand0))/denom))

;create the header
  printf,42,'# created by MKMMS_LS.pro -- in the MODS-XIDL package'
  printf,42,'# orig filename  : 1.0arcsecSegmentedLongSlit.mms'
  printf,42,'# user-date      : ',systime(/utc)
  printf,42,'# fits file      : ',infile
  printf,42,'# fits instrument: ',fitsinst
  printf,42,'# fits origin    : ',targorig
  printf,42,'# fits object    : ',fitsobj
  printf,42,'# target setup for MODS'
  printf,42,'# --------------------------------------------'
  printf,42,''
  printf,42,''
  printf,42,'TEL.TARG.ALPHA    ',targalpha
  printf,42,'TEL.TARG.DELTA    ',targdelta
  printf,42,'TEL.ROT.OFFANGLE  ',targpa
  printf,42,'TEL.TARG.EQUINOX  ',targeq
  printf,42,''
  printf,42,'INS.SLIT.NUMBER    5'
  printf,42,''
  printf,42,'INS.TARG101.NAME               CENTER'
  printf,42,'INS.TARG101.SHAPE             STRAIGHT'
  printf,42,'INS.TARG101.WID    5.0'
  printf,42,'INS.TARG101.LEN   60.0'
  printf,42,'INS.TARG101.ROT    0.0'
  printf,42,'INS.TARG101.ALPHA  ',targalpha
  printf,42,'INS.TARG101.DELTA  ',targdelta
  printf,42,'INS.TARG101.WIDMM    3.600'
  printf,42,'INS.TARG101.LENMM   36.000'
  printf,42,'INS.TARG101.XMM      0.000'
  printf,42,'INS.TARG101.YMM      0.000'
  printf,42,''
  printf,42,'SEQ.PREP.MASK_NA M010_LS60x5'
  printf,42,'SEQ.PREP.MASK_ID     500001'

close,42
end
