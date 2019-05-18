;+
; NAME:
;   mods_plan
;
; PURPOSE:
;   Create plan file(s) for running the mods pipeline.  This code
;   parses headers, does image stats, etc.
;
; CALLING SEQUENCE:
;   mods_plan, [ fileexpr, indir, planfile= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   fileexpr   - File names in the input directory; default to '*.fits*'
;   indir      - Input directory(s) for reading files;
;                default to current directory
;   planfile   - Output plan file; default to 'plan.par'.
;                This file is put in the same directory as the raw data files.
;
; OUTPUT:
;
; COMMENTS:
;   One plan file is made for each input directory.
;
;   The following flavors of images are listed:
;     bias
;     domeflat
;     iflat (internal flat)
;     twiflat
;     arc
;     science
;
; EXAMPLES:
; mods_plan,'*.fits','/b/martell/data_arx/09072005/Raw/',planfile='plan-master.par'
; BUGS:
;
; PROCEDURES CALLED:
;   fileandpath()
;   headfits()
;   idlutils_version() - removed 2019 May 18 [rwp]
;   splog
;   sxpar()
;   yanny_write
;
; INTERNAL SUPPORT ROUTINES:
;   mods_plan_struct()
;
;  To Do: for twilight images taken with a grism in but no mask,
;  assume they are slitless, and treat them as pixel flat field images 
;  not twiflats. 
;
;
;
; REVISION HISTORY:
;   13-Mar-2005  Written by David Schlegel, LBL.
;-
;------------------------------------------------------------------------------
function mods_plan_struct, nfile
   planstr = create_struct(name='lexp', $
    'FILENAME'    , '', $
    'FLAVOR'      , '', $
    'TARGET'      , '', $
    'EXPTIME'     , 0., $
    'INSTRUMENT'  , '', $
    'GRATING'     , '', $
    'WAVE'        , '', $                       
    'MASKNAME'    , ''  )
   return, replicate(planstr, nfile)
end
;------------------------------------------------------------------------------
pro mods_plan, fileexpr, indir, planfile=planfile

   if not keyword_set(verbose) then !Quiet=1
   COMMON SITE, lat, lng, tzone
   DRADEG = 180.d0/!dpi
   
   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(fileexpr)) then fileexpr = '*.fits*'
   if (NOT keyword_set(planfile)) then planfile = 'plan.par'

   if not keyword_set(INDIR) then begin
       spawn, '\ls -d '+indir, dirlist
       if (NOT keyword_set(dirlist)) then begin
           splog, 'No input directories found'
           return
       endif
       ndir = n_elements(dirlist)
   endif else begin
       dirlist = [indir]
       ndir = n_elements(dirlist)
   endelse
   
   ;----------
   ; Loop over each input directory

   for idir = 0L, ndir-1L do begin
       splog, 'Working on directory ', dirlist[idir]
       cd, dirlist[idir], current = olddir
       if (idir EQ 0) then origdir = olddir
       filenames = findfile(fileexpr, count = nfile)
       splog, 'Number of FITS files found: ', nfile
       
       if (nfile GT 0) then begin
           planstr = mods_plan_struct(nfile)
           planstr.filename = filenames
           for i = 0L, nfile-1L do begin
               hdr = xheadfits(filenames[i])
               if (size(hdr, /tname) EQ 'STRING') then begin
                  if (strmatch(sxpar(hdr, 'INSTRUME'), 'LRIS*')) then begin
                      ;;---------------------------
                      ;; This is the LRIS on Keck I
                      ;;---------------------------
                       exptime = sxpar(hdr, 'ELAPTIME')
                       trapdoor = strtrim(sxpar(hdr, 'TRAPDOOR'))
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       ;; Figure out if telescope is pointing
                       ;; throught the dome slit (to within 4
                       ;; degrees), or at the
                       ;; inside of the dome.  Note that domeflats are
                       ;; generally taken with a telescope
                       ;; position 90 degrees from the slit,
                       ;; but we won't force this case to be true.
                       slitposn = sxpar(hdr, 'AZ')
                       domeposn = sxpar(hdr, 'DOMEPOSN')
                       slitposn = (slitposn + 360.) mod 360.
                       domeposn = (domeposn + 360.) mod 360.
                       if (abs(domeposn - slitposn) lt 4.) then $
                          throughslit = 1 else throughslit = 0
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       observatory, 'keck', obs_struct
                       tzone = obs_struct.tz
                       lng = 360.d0 - obs_struct.longitude
                       lat = obs_struct.latitude
                       jd = 2400000.5D + sxpar(hdr, 'MJD-OBS')
                       sunpos, jd, ra1, dec1 ; returns degrees
                       zenpos, jd, ra2, dec2 ; returns radians
                       ra2 = ra2 * DRADEG
                       dec2 = dec2 * DRADEG
                       sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                       lamps = sxpar(hdr, 'LAMPS')
                       if (lamps EQ '0,0,0,0,0,1') then $
                         planstr[i].flavor = 'iflat' $
                       else if (strmid(lamps, 0, 9) NE '0,0,0,0,0') then $
                         planstr[i].flavor = 'arc' $
                       else if (trapdoor EQ 'closed') then $
                         planstr[i].flavor = 'bias' $
                       else if (throughslit eq 0) then $
                         planstr[i].flavor = 'domeflat' $
                       else if (sun_angle GT -12.) then $
                         planstr[i].flavor = 'twiflat' $
                       else if exptime LT 120.0 THEN $
                         planstr[i].flavor = 'std' $
                       else planstr[i].flavor = 'science'
                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'TARGNAME'))
                       planstr[i].maskname = strtrim(sxpar(hdr, 'SLITNAME'))
                  ; Note that GRISNAME will only be defined for blue CCDs,
                  ; and GRANAME will only be defined for red CCDs. 
                  ; MSWAVE is only defined for red CCDs
                       instrument = strtrim(sxpar(hdr, 'INSTRUME'))
                       planstr[i].instrument = instrument
                       planstr[i].grating = $
                         + (instrument EQ 'LRISBLUE' ? sxpar(hdr, 'GRISNAME') : '') $
                         + (instrument EQ 'LRIS' ? sxpar(hdr, 'GRANAME') : '')
                       planstr[i].wave = $
                         + (instrument EQ 'LRISBLUE' ? '' : '') $
                         + (instrument EQ 'LRIS' ? $
                            strcompress(sxpar(hdr, 'MSWAVE'), /rem): '')
                       
                    endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*CA-3.5*') $
                    then begin
                       ;;---------------------------
                       ;; This is the CAHA 3.5m
                       ;;---------------------------
                       exptime = sxpar(hdr, 'EXPTIME')
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       observatory, 'ca', obs_struct
                       tzone = obs_struct.tz
                       lng = 360.d0 - obs_struct.longitude
                       lat = obs_struct.latitude
                       jd = sxpar(hdr, 'JUL-DATE')
                       sunpos, jd, ra1, dec1 ; returns degrees
                       zenpos, jd, ra2, dec2 ; returns radians
                       ra2 = ra2 * DRADEG
                       dec2 = dec2 * DRADEG
                       sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                       caha_hh = caha_hdr(hdr)
                       imagetyp = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                       IF imagetyp EQ 'dark' OR imagetyp EQ 'bias' $
                       THEN planstr[i].flavor = 'bias' $
                       ELSE IF imagetyp EQ 'simulation' OR imagetyp EQ 'arc'THEN $
                          planstr[i].flavor = 'arc' $
                       ELSE IF imagetyp EQ 'flat' THEN $
                          planstr[i].flavor = 'domeflat' $
                       ELSE IF imagetyp EQ 'science' THEN BEGIN 
                          IF sun_angle GT -12. THEN $
                             planstr[i].flavor = 'twiflat' $
                          ELSE IF exptime LT 350.0 THEN $
                             planstr[i].flavor = 'std' $
                          ELSE planstr[i].flavor = 'science'
                       ENDIF ELSE message,  'Unrecognized obstype'
                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'OBJECT'))
                       planstr[i].maskname = strtrim(sxpar(caha_hh, 'SLIT_WID'))
                       instr1 = strcompress(sxpar(hdr, 'INSTRUME'), /rem)
                       path = strcompress(sxpar(caha_hh, 'PATH'), /rem)
                       instrument = instr1 + '-' + path
                       planstr[i].instrument = instrument
                       planstr[i].grating = $
                          + (path EQ 'BLUE' ? sxpar(caha_hh, 'GRAT1_NA') : '') $
                          + (path EQ 'RED'  ? sxpar(caha_hh, 'GRAT2_NA') : '')
                       planstr[i].wave = $
                          + (path EQ 'BLUE' ? sxpar(caha_hh, 'GRAT1_WL') : '') $
                          + (path EQ 'RED'  ? sxpar(caha_hh, 'GRAT2_WL') : '')
                    endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*CA-2.2*') $
                    then begin
                       ;;---------------------------
                       ;; This is the CAHA 2.2m
                       ;;---------------------------
                       exptime = sxpar(hdr, 'EXPTIME')
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       observatory, 'ca', obs_struct
                       tzone = obs_struct.tz
                       lng = 360.d0 - obs_struct.longitude
                       lat = obs_struct.latitude
                       jd = sxpar(hdr, 'JUL-DATE')
                       sunpos, jd, ra1, dec1 ; returns degrees
                       zenpos, jd, ra2, dec2 ; returns radians
                       ra2 = ra2 * DRADEG
                       dec2 = dec2 * DRADEG
                       sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                       ;; This may be unnecessary
                       caha_hh = hdr
                       imagetyp = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                       IF imagetyp EQ 'dark' OR imagetyp EQ 'bias' $
                       THEN planstr[i].flavor = 'bias' $
                       ELSE IF imagetyp EQ 'simulation' OR imagetyp EQ 'arc'THEN $
                          planstr[i].flavor = 'arc' $
                       ELSE IF imagetyp EQ 'flat' THEN $
                          planstr[i].flavor = 'domeflat' $
                       ELSE IF imagetyp EQ 'science' THEN BEGIN 
                          IF sun_angle GT -12. THEN $
                             planstr[i].flavor = 'twiflat' $
                          ELSE IF exptime LT 350.0 THEN $
                             planstr[i].flavor = 'std' $
                          ELSE planstr[i].flavor = 'science'
                       ENDIF ELSE message,  'Unrecognized obstype'
                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'OBJECT'))
                       planstr[i].maskname = strtrim(sxpar(caha_hh, 'INSAPDY')) ;; Microns not arcsec
                       instr1 = strcompress(sxpar(hdr, 'INSTRUME'), /rem)
                       path = strcompress(sxpar(caha_hh, 'PATH'), /rem)
                       instrument = instr1 + '-' + path
                       planstr[i].instrument = instrument
                       planstr[i].grating = strtrim(sxpar(caha_hh,'INSGRNAM'),2)
                       planstr[i].wave = float(sxpar(caha_hh,'INSGRWL0'))
                    endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*ESO-NTT*') $
                    then begin
                       ;;---------------------------
                       ;; This is EFOSC on the NTT 3.6m
                       ;;---------------------------
                       exptime = sxpar(hdr, 'EXPTIME')
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       observatory, 'eso', obs_struct
                       tzone = obs_struct.tz
                       lng = 360.d0 - obs_struct.longitude
                       lat = obs_struct.latitude
                       jd = 2400000.5D + sxpar(hdr, 'MJD-OBS') 
                       sunpos, jd, ra1, dec1 ; returns degrees
                       zenpos, jd, ra2, dec2 ; returns radians
                       ra2 = ra2 * DRADEG
                       dec2 = dec2 * DRADEG
                       sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                       ind_catg = $
                          WHERE(stregex(hdr, 'HIERARCH ESO DPR CATG', /bool))
                       catg = strcompress(repstr(strmid(hdr[ind_catg], 30, 14), "'", ''), /rem)
                       ind_type = $
                          WHERE(stregex(hdr, 'HIERARCH ESO DPR TYPE', /bool))
                       type = strcompress(repstr(strmid(hdr[ind_type], 30, 14), "'", ''), /rem)
                       IF catg EQ 'CALIB' THEN BEGIN
                          IF type EQ 'DARK' OR type EQ 'BIAS' $
                          THEN planstr[i].flavor = 'bias' $
                          ELSE IF type EQ 'WAVE' THEN $
                             planstr[i].flavor = 'arc' $
                          ELSE IF type EQ 'FLAT' THEN $
                             planstr[i].flavor = 'domeflat' 
                       ENDIF ELSE IF catg EQ 'SCIENCE' THEN BEGIN
                          IF sun_angle GT -12. THEN $
                             planstr[i].flavor = 'twiflat' $
                          ELSE IF exptime LT 350.0 OR type EQ 'STD' THEN $
                             planstr[i].flavor = 'std' $
                          ELSE planstr[i].flavor = 'science'
                       ENDIF ELSE message,  'Unrecognized obstype'
                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'OBJECT'))
                       ind_slit = $
                          WHERE(stregex(hdr, 'HIERARCH ESO INS SLIT1 NAME' $
                                        , /bool))
                       slit = strcompress(repstr(strmid(hdr[ind_slit], 30, 14), "'", ''), /rem)
                       slit = repstr(slit, '#', '')
                       planstr[i].maskname = slit
                       planstr[i].instrument = $
                          strcompress(sxpar(hdr, 'INSTRUME'), /rem)
                       ind_gris = $
                          WHERE(stregex(hdr, 'HIERARCH ESO INS GRIS1 NAME' $
                                        , /bool))
                       gris = strcompress(repstr(strmid(hdr[ind_gris], 30, 14), "'", ''), /rem)
                       gris = repstr(gris, '#', '')
                       planstr[i].grating = gris
                       planstr[i].wave = '0.0'
           endif else if  (strmatch(sxpar(hdr, 'INSTRUME'), 'GMOS-N*')) OR $
                     (strmatch(sxpar(hdr, 'INSTRUME'), 'GMOS-S*')) $
                     then begin
                    ;---------------------------
                    ; This is GMOS on Gemini-N/S
                    ;---------------------------
                       exptime = strcompress(sxpar(hdr, 'EXPOSURE'), /rem)
                       obstype = strcompress(sxpar(hdr, 'OBSTYPE'), /rem)
                       object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                       instrument = strcompress(sxpar(hdr, 'INSTRUME'), /rem)
                       
                       if obstype EQ 'FLAT' then planstr[i].flavor = 'domeflat' $
                       else if obstype EQ 'ARC' then planstr[i].flavor = 'arc' $
                       else if obstype EQ 'BIAS'then planstr[i].flavor = 'bias' $
                       else if (obstype EQ 'OBJECT' AND object EQ 'Twilight') $
                         then planstr[i].flavor = 'twiflat' $
                       else if obstype EQ 'MASK'then planstr[i].flavor $
                          = 'mask' $
                       else if obstype EQ 'DARK'then planstr[i].flavor $
                          = 'dark' $
                       else if (obstype EQ 'OBJECT' AND object NE 'Twilight') $
                         then BEGIN
                           if exptime LT 120.0 THEN planstr[i].flavor = 'std' $
                           ELSE  planstr[i].flavor = 'science'
                       ENDIF ELSE message,  'Unrecognized obstype'
                       planstr[i].exptime = exptime
                       planstr[i].target = object
                       planstr[i].maskname = strtrim(sxpar(hdr, 'MASKNAME'))
                       instrument = strtrim(sxpar(hdr, 'INSTRUME'))
                  planstr[i].instrument = instrument
                  ;; XGMOS ?
                  if sxpar(hdr,'XGMOS') then $
                    planstr[i].instrument = planstr[i].instrument+'X'
                  planstr[i].grating =  strcompress(sxpar(hdr, 'GRATING') $
                                                    , /rem)
                  planstr[i].wave = strcompress(sxpar(hdr, 'GRWLEN'), /rem)
              endif else if (strmatch(sxpar(hdr, 'INSTRUME'), 'DEIMOS*')) then begin
                  ;---------------------------
                  ; This is the DEIMOS on Keck I
                  ;---------------------------
                       exptime = sxpar(hdr, 'ELAPTIME')
                       trapdoor = strtrim(sxpar(hdr, 'TRAPDOOR'))
                  ; Decide if this is twilight: sun_angle is the
                  ; angle of the sun *above* the horizon, so must
                  ; be less than -12 if darker than 12-degree twi.
                       observatory, 'keck', obs_struct
                       tzone = obs_struct.tz
                       lng = 360.d0 - obs_struct.longitude
                       lat = obs_struct.latitude
                       jd = 2400000.5D + sxpar(hdr, 'MJD-OBS')
                       sunpos, jd, ra1, dec1 ; returns degrees
                       zenpos, jd, ra2, dec2 ; returns radians
                       ra2 = ra2 * DRADEG
                       dec2 = dec2 * DRADEG
                       sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)

                       lamps = sxpar(hdr, 'LAMPS')
                       if strmatch(lamps,'Qz*') then $
                         planstr[i].flavor = 'domeflat' $
                       else if (strmid(lamps, 0, 3) NE 'Off') then $
                         planstr[i].flavor = 'arc' $
                       else if strmatch(trapdoor,'closed*') then $
                         planstr[i].flavor = 'bias' $
                       else if (sun_angle GT -12.) then $
                         planstr[i].flavor = 'twiflat' $
                       else begin
                           if exptime LT 120.0 THEN $
                             planstr[i].flavor = 'std' $
                           else planstr[i].flavor = 'science'
                       endelse 
                       planstr[i].exptime = exptime
                       if not strmatch(sxpar(hdr,'MOSMODE'),'Spectral') then $
                         planstr[i].flavor = 'unknown'
                       
                       planstr[i].target = strtrim(sxpar(hdr, 'TARGNAME'))
                       planstr[i].maskname = strtrim(sxpar(hdr, 'SLMSKNAM'))

                       instrument = strtrim(sxpar(hdr, 'INSTRUME'))
                       planstr[i].instrument = instrument
                       planstr[i].grating = strtrim(sxpar(hdr, 'GRATENAM'),2)
                       gpos = strtrim(sxpar(hdr, 'GRATEPOS'),2)

                       planstr[i].wave = sxpar(hdr,'G'+strtrim(gpos,2)+'TLTWAV')
              endif else if strmatch(sxpar(hdr, 'TELESCOP'), 'mmt*') THEN BEGIN
                    ;---------------------------
                    ; This is Blue/Red Channel on MMT 
                    ;---------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  obstype = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  ccd    =  strcompress(sxpar(hdr, 'DETECTOR'), /rem)
                  if ccd EQ 'ccd35' or ccd EQ 'mmtbluechan' then begin
                     instrument = 'MMT Blue Channel' 
                     dispers = sxpar(hdr, 'DISPERSE')
                     IF KEYWORD_SET(DISPERS) THEN $
                        planstr[i].grating = strcompress(dispers, /rem) $
                     ELSE planstr[i].grating =  '831 2nd order'
                     planstr[i].wave = ''
                  ENDIF else if (ccd EQ 'ccd34') OR strmatch(ccd,'mmtredchan') then BEGIN
                     instrument = 'MMT Red Channel' 
                     ;; These next lines might not work with the old
                     ;; chip
                     dispers = sxpar(hdr, 'DISPERSE')
                     prs = strsplit(dispers,'-',/extract)
                     planstr[i].grating =  prs[0]
                     planstr[i].wave = prs[1]
                  Endif ELSE message, 'Unrecognized MMT instrument'
                  if obstype EQ 'flat' then planstr[i].flavor = 'domeflat' $
                  else if obstype EQ 'comp' then planstr[i].flavor = 'arc' $
                  else if obstype EQ 'zero'then planstr[i].flavor = 'bias' $
                  else if obstype EQ 'object' then begin
                      if exptime LT 120.0 THEN planstr[i].flavor = 'std' $
                      else planstr[i].flavor = 'science' 
                  endif else message,  'Unrecognized obstype'
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  planstr[i].maskname = object
                  planstr[i].instrument = instrument
               endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*kp4m*') $
               THEN BEGIN
                  ;;---------------------------
                  ;; This is RC Spectrograph on the KPNO 4m
                  ;;--------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  obstype = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  ccd    =  strcompress(sxpar(hdr, 'DETECTOR'), /rem)
                  instrument = 'RC Spectrograph'
                  dispers = sxpar(hdr, 'DISPERSE')
                  planstr[i].grating = strcompress(dispers, /rem) 
                  planstr[i].wave = strcompress(sxpar(hdr, 'TILTPOS'), /rem)
                  if obstype EQ 'flat' then planstr[i].flavor = 'domeflat' $
                  else if obstype EQ 'comp' then planstr[i].flavor = 'arc' $
                  else if obstype EQ 'zero'then planstr[i].flavor = 'bias' $
                  else if obstype EQ 'object' then begin
                     if exptime LT 200.0 THEN planstr[i].flavor = 'std' $
                     else planstr[i].flavor = 'science' 
                  endif else if obstype EQ 'focus' THEN $
                     planstr[i].flavor = 'focus' $
                  else message,  'Unrecognized obstype'
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  planstr[i].maskname =  $
                     strcompress(sxpar(hdr, 'APERTURE'), /rem)
                  planstr[i].instrument = instrument
               endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*Baade*') $
               THEN BEGIN
                    ;---------------------------
                    ; This is IMACS on Magellan 
                    ;---------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  obstype = strcompress(sxpar(hdr, 'EXPTYPE'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  instrument = 'IMACS'
                  maskname = strcompress(sxpar(hdr, 'SLITMASK'), /rem)
                  dispers = sxpar(hdr, 'DISPERSE')
                  planstr[i].grating = strcompress(sxpar(hdr, 'DISPERSR'), /rem)
                  planstr[i].wave = ''
                  IF exptime GT 10.0 THEN planstr[i].flavor = 'science' $
                  ELSE  planstr[i].flavor = '????'
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  planstr[i].maskname = maskname
                  planstr[i].instrument = instrument
              endif else if  strmatch(sxpar(hdr, 'DETECTOR'), 'mars*') $
                then begin
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  obstype = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  instrument = 'KPNO Mars'
                  if obstype EQ 'flat' then planstr[i].flavor = 'domeflat' $
                  else if obstype EQ 'comp' then planstr[i].flavor = 'arc' $
                  else if obstype EQ 'zero'then planstr[i].flavor = 'bias' $
                  else if obstype EQ 'object' then $
                  planstr[i].flavor = 'science' $
                  else message,  'Unrecognized obstype'
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  planstr[i].maskname = object
                  planstr[i].instrument = instrument
                  planstr[i].grating =  '8050'
                  planstr[i].wave = ''
;              endif else if (strmatch(strtrim(sxpar(hdr, 'INSTRUME'),2), 'Kast') OR strmatch (strtrim(sxpar(hdr, 'INSTRUME'),2), 'Kast blue arm') OR  strmatch(strtrim(sxpar(hdr, 'INSTRUME'),2), 'KAST blue') OR strmatch(strtrim(sxpar(hdr, 'INSTRUME'),2), 'KAST'))  then begin

              endif else if (stregex(sxpar(hdr,'INSTRUME'),'.*kast.*',$
                                    /boolean,/fold_case) eq 1) or $
                (stregex(sxpar(hdr,'VERSION'),'kast*', /bool, /fold_case) EQ 1) then begin
                  ;;---------------------------
                  ;; This is Kast on Shane 3m (Mt. Hamilton)
                  ;;---------------------------
                  data = xmrdfits(filenames[i],/fscale)
                  exptime = float(sxpar(hdr, 'EXPOSURE')) > $
                            float(sxpar(hdr, 'EXPTIME'))
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  planstr[i].target = object
                  planstr[i].exptime = exptime
                  ;; Blue or Red?
                  if strtrim(sxpar(hdr,'SPSIDE'),2) EQ 'blue' OR $
                    strmid(sxpar(hdr,'VERSION'),0,5) EQ 'kastb' then begin
                      planstr[i].instrument = 'KAST-B' 
                      planstr[i].grating = strtrim(sxpar(hdr,'GRISM_N'),2)
                      planstr[i].wave = sxpar(hdr,'GRISM_P')
                  endif else begin
                      planstr[i].instrument = 'KAST-R' 
                      planstr[i].grating = strtrim(sxpar(hdr,'GRATNG_N'),2)
                      planstr[i].wave = sxpar(hdr,'GRTILT_P')
                  endelse
                  ;; Type
                  if strmatch(planstr[i].grating,'open') then begin
                      print, 'Imaging ', filenames[i], '.  Skipping..'
                      planstr[i].flavor = 'IMG'
                      continue
                  endif
                  if strtrim(sxpar(hdr,'SHUTTER'),2) NE 'OPEN' AND $
                    strtrim(sxpar(hdr,'OBSTYPE'),2) NE 'OBJECT' then begin
                      if exptime LE 2 then planstr[i].flavor = 'bias' $
                      else planstr[i].flavor = 'dark'
                  endif else begin
                      if exptime LT 120. then begin
                          if median(data) GT 500. then planstr[i].flavor = 'domeflat' $
                          else planstr[i].flavor = 'arc'  
                      endif else planstr[i].flavor = 'science'  
                  endelse
              endif else if strmatch(strtrim(sxpar(hdr, 'INSTRUME'), 2) $
                                     , 'DIS') then begin
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJNAME'), /rem)
                  IF strmatch(object, '*TestSlew*') THEN $
                    planstr[i].target = 'calib' $
                  ELSE planstr[i].target = object
                  planstr[i].exptime = exptime
                  ;; Blue or Red?
                  if strmatch(sxpar(hdr, 'DETECTOR'), '*blue*') then $
                    planstr[i].instrument = 'DIS-B' $
                  else if strmatch(sxpar(hdr, 'DETECTOR'), '*red*') then $
                    planstr[i].instrument = 'DIS-R' $
                  else planstr[i].instrument = ' '
                  planstr[i].grating = strtrim(sxpar(hdr, 'GRATING'), 2)
                  planstr[i].wave = sxpar(hdr, 'DISPWC')
                  planstr[i].maskname = strtrim(sxpar(hdr, 'SLITMASK'))
                  if strmatch(sxpar(hdr, 'IMAGETYP'), '*zero*')  $
                    then planstr[i].flavor = 'bias' $
                  else if  strmatch(sxpar(hdr, 'IMAGETYP'), '*flat*') $
                    then planstr[i].flavor = 'iflat' $
                  else if strmatch(sxpar(hdr, 'IMAGETYP'), '*object*') AND $
                    exptime GT 120.0 THEN planstr[i].flavor = 'science' $
                  else if strmatch(sxpar(hdr, 'IMAGETYP'), '*object*') AND $
                    exptime LE 120.0 THEN planstr[i].flavor = 'arc' $
                  else planstr[i].flavor = 'unknown'  
               endif else if strmatch(strtrim(sxpar(hdr, 'TELID'), 2), '200') $
                  THEN BEGIN
                    ;---------------------------
                    ; This is DBSP on Palomar 200"
                    ;---------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  planstr[i].exptime = exptime
                  ;; Blue or Red?
                  IF strtrim(sxpar(hdr, 'FPA'), 2) EQ 'DBSP_BLUE' THEN $
                     planstr[i].instrument = 'P200-B' $
                  ELSE IF strtrim(sxpar(hdr, 'FPA'), 2) EQ 'DBSP_RED' THEN $
                     planstr[i].instrument = 'P200-R' $
                  ELSE planstr[i].instrument = ' '
                  planstr[i].grating = strtrim(sxpar(hdr, 'GRATING'), 2)
                  planstr[i].wave =  strtrim(sxpar(hdr, 'ANGLE'), 2)
                  planstr[i].maskname = 'slit' + strtrim(sxpar(hdr, 'APERTURE'))
                  lamps  = strtrim(sxpar(hdr, 'LAMPS'), 2)
                  lamp_vec = long(strmid(lamps, lindgen(7), 1))
                  turret = strtrim(sxpar(hdr, 'TURRET'), 2)
                  imgtype=strtrim(sxpar(hdr,'IMGTYPE'))
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  IF strmatch(imgtype, 'flat') THEN BEGIN
                     planstr[i].flavor = 'iflat' 
                     planstr[i].target = object
                  ENDIF ELSE IF strmatch(imgtype, 'cal') THEN BEGIN
                     planstr[i].flavor = 'arc'
                     lamp_str = ''
                     lamp_names = ['D', 'FeAr', 'Hg', 'Ar', 'Ne', 'He'$
                                   , 'InCand']
                     FOR j = 0L, 6 DO BEGIN
                        IF lamp_vec[j] GT 0 THEN BEGIN
                           IF NOT KEYWORD_SET(lamp_str) THEN $
                              lamp_str = lamp_names[j] $
                           ELSE lamp_str = lamp_str + '-' + lamp_names[j]
                        ENDIF
                     ENDFOR
                     planstr[i].target = lamp_str
                  ENDIF ELSE IF exptime EQ 0 THEN BEGIN
                     planstr[i].flavor = 'bias' 
                     planstr[i].target = 'none'
                  ENDIF ELSE IF strmatch(imgtype, 'object') THEN BEGIN
                     if exptime LT 120.0 THEN planstr[i].flavor = 'std' $
                     else planstr[i].flavor = 'science' 
                     planstr[i].target = object
                  ENDIF ELSE planstr[i].flavor = 'unknown'  

               endif else if (strmatch(sxpar(hdr, 'INSTRUME'), 'MODS1*')) then begin
                      ;;---------------------------
                      ;; This is MODS1 on LBT
                      ;;---------------------------
                       exptime = sxpar(hdr, 'EXPTIME')
                       imtype = strtrim(sxpar(hdr, 'IMAGETYP'))
                       
                       if imtype eq 'BIAS' then planstr[i].flavor = 'bias' $
                       else if imtype eq 'COMP' then planstr[i].flavor = 'arc' $
                       else if imtype eq 'FLAT' then planstr[i].flavor = 'domeflat' $
                       else if imtype eq 'SKY' then planstr[i].flavor = 'twiflat' $
                       else if imtype eq 'STD' then planstr[i].flavor = 'std' $
                       else if imtype eq 'OBJECT' then planstr[i].flavor = 'science'

                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'OBJECT'))
                       planstr[i].maskname = strtrim(sxpar(hdr, 'MASKNAME'))
                       instrument = strtrim(sxpar(hdr, 'INSTRUME'))
                       planstr[i].instrument = instrument
                       planstr[i].grating = strtrim(sxpar(hdr, 'GRATNAME'))
                       
               endif else if (strmatch(sxpar(hdr, 'INSTRUME'), 'MODS2*')) then begin
                      ;;---------------------------
                      ;; This is MODS2 on LBT
                      ;;---------------------------
                       exptime = sxpar(hdr, 'EXPTIME')
                       imtype = strtrim(sxpar(hdr, 'IMAGETYP'))
                       
                       if imtype eq 'BIAS' then planstr[i].flavor = 'bias' $
                       else if imtype eq 'COMP' then planstr[i].flavor = 'arc' $
                       else if imtype eq 'FLAT' then planstr[i].flavor = 'domeflat' $
                       else if imtype eq 'SKY' then planstr[i].flavor = 'twiflat' $
                       else if imtype eq 'STD' then planstr[i].flavor = 'std' $
                       else if imtype eq 'OBJECT' then planstr[i].flavor = 'science'

                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'OBJECT'))
                       planstr[i].maskname = strtrim(sxpar(hdr, 'MASKNAME'))
                       instrument = strtrim(sxpar(hdr, 'INSTRUME'))
                       planstr[i].instrument = instrument
                       planstr[i].grating = strtrim(sxpar(hdr, 'GRATNAME'))
                       
               endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*bok*') $
                    THEN BEGIN
                  ;;---------------------------
                  ;; This is the B&C spectrograph at the Bok 2.3-meter
                  ;; telescope; added by Moustakas on 2011-Jun-08
                  ;;---------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  obstype = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  ccd    =  strcompress(sxpar(hdr, 'DETECTOR'), /rem)
                  instrument = 'Bok B&C'
                  dispers = sxpar(hdr, 'DISPERSE')
                  planstr[i].grating = strcompress(dispers, /rem) 
                  planstr[i].wave = strcompress(sxpar(hdr, 'TILTPOS'), /rem)
                  if obstype EQ 'flat' then planstr[i].flavor = 'domeflat' $
                  else if obstype EQ 'comp' then planstr[i].flavor = 'arc' $
                  else if obstype EQ 'zero'then planstr[i].flavor = 'bias' $
                  else if obstype EQ 'object' then begin
                     planstr[i].flavor = 'science' 
;                    if exptime LT 200.0 THEN planstr[i].flavor = 'std' $
;                    else planstr[i].flavor = 'science' 
                  endif else if obstype EQ 'focus' THEN $
                     planstr[i].flavor = 'focus' $
                  else message,  'Unrecognized obstype'
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  planstr[i].maskname =  $
                     strcompress(sxpar(hdr, 'APERTURE'), /rem)
                  planstr[i].instrument = instrument
               ENDIF ELSE $
                  splog, 'WARNING: Unknown instrument for ', filenames[i]
            ENDIF
      ENDFOR              ; End loop over files
           
      logfile = repstr(planfile, '.par', '') + '.log'
      plotfile = repstr(planfile, '.par', '') + '.ps'

         hdr = ''
;         hdr = [hdr, '# Biases grouped by INSTRUMENT']
;         hdr = [hdr, '# Pixel flats grouped by INSTRUMENT+GRATING+MASKNAME+WAVE']
;         hdr = [hdr, '# Wavelength solutions from arcs grouped by INSTRUMENT+GRATING+MASKNAME+WAVE']
;         hdr = [hdr, '# Illumination flats grouped by INSTRUMENT+GRATING+MASKNAME+WAVE']
         hdr = [hdr, ' ']
         hdr = [hdr, "logfile '" + logfile + "'   # Log file"]
         hdr = [hdr, "plotfile '" + plotfile + "'   # Plot file"]
         hdr = [hdr, "indir '" + dirlist[idir] + "'   # Raw data directory"]
         hdr = [hdr, "tempdir Temp     # Temporary (working) directory"]
         hdr = [hdr, "scidir  Science  # Science output directory"]
;         if n_elements(maxobj) GT 0 then $
;            hdr = [hdr, "maxobj  "+string(maxobj,format='(i)') + $
;                        " # Maximum number of Objects per slit"] $
;         else 
         hdr = [hdr, "maxobj     5  # Maximum number of Objects per slit"] 
         hdr = [hdr, "minslit 20  # Minimum slit width"]
;         hdr = [hdr, "slitthresh 0.02 # Sets threshold for slit identification. 0.1 works best for Longslit spectra"]
         hdr = [hdr, "reduxthresh 0.01 # Sets the fraction of the brightest objects on each slit that is reduced"]
         ; hdr = [hdr, "idlutilsVersion '" + idlutils_version() + "'  # Version of idlutils when building plan file"]
         hdr = [hdr, "LongslitVersion '" + longslit_version() $
          + "'  # Version of Longslit when building plan file"]

         ; Write this plan file
         cd, olddir
         yanny_write, planfile, ptr_new(planstr), hdr=hdr, /align

      endif
   endfor ; End loop over directories

   cd, origdir

   return
end
;------------------------------------------------------------------------------
