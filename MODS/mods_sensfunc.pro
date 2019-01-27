;+
; NAME:
;   long_sensfunc
;
; PURPOSE:
;   Use one or more standard star spectra to determine the
;   spectroscopic response function. 
;
; CALLING SEQUENCE:
;
; INPUTS:
;   scifile         - file containing object structure which has spectrum 
;                     for standard star
;
;   standard_name   - name of the standard star
;   
;   sensfuncfile    - File to write the sensitivity function out to
;
;
; OPTIONAL INPUTS:
;   OBJID           - object id of standar star in the object structure. Default
;                     is to the first object. 
;   nresln          - Break point spacing in resolution elemnets
;                     (default=20)
;   /MSK_BALM       - Mask Balmer lines (recommended but not default)
;   /NOGREY         - do not grey-shift (applicable when multiple
;     standards are passed)
;
; OUTPUTS:
;   mag_set        - structure containing b-spline info for sensitivity function
;
; OPTIONAL OUTPUTS:
;  sensfunc         - sensitivity function evaluated
;
; COMMENTS:
;                   See README file in /apps2/iraf211/iraf/noao/lib/onedstds/
;                   for list of standard stars and the names of the
;                   associated files
;
; EXAMPLES:
;
; BUGS:
;                   Does not take into account atmospheric extinction!!!
;                   Leaves out first and last wavelength bins of
;                   sensitivity function
;
; PROCEDURES CALLED:
;   traceset2xy (idlutils)
;   xy2traceset (idlutils)
;   splog       (idlutils)
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   01-Oct-2005  Written by J. Hennawi UC Berkeley
;   08-Jun-2011  J. Moustakas, UCSD - optionally deal with a list of
;     standards, and generate a QAplot
;------------------------------------------------------------------------------
FUNCTION mods_sensfunc, stdstar,sensfuncfile, SCHHDR = SCIHDR $
		      , std_name = std_name , STDFILE = STDFILE $
		      , MSK_BALM = msk_balm $
                      , LINE_BALM = LINE_BALM $
                      , BALM_MASK_WID = BALM_MASK_WID $
		      , AIRMASS = AIRMASS $
		      , EXPTIME = EXPTIME $
                      , sensfunc = sensfunc, sensfit = sensfit $
                      , wave = wave, flux=flux, nresln = nresln $
                      , chk = chk $
                      , sciind = sciind1, std_flux = flux_std_int $
                      ,  INMASK = INMASK, FINALMASK = FINALMASK $
                      , IVAR=ivar, resln=resln, nofit=nofit, nogrey=nogrey $
			, no_tell_mask=no_tell_mask,verbose=verbose,alt=alt

  if not keyword_set(verbose) then !Quiet=1

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'sfunc = long_sensfunc(scifil, outfil, STD_NAME=) [v1.0]'
      return, -1
  endif 

; stay five pixels away from edge
  IF NOT KEYWORD_SET(NRESLN) THEN NRESLN = 10  ;; JXP -- Avoid 'narrow' features
  IF NOT KEYWORD_SET(LO_BUFFER) THEN LO_BUFFER = 8L
  IF NOT KEYWORD_SET(HI_BUFFER) THEN HI_BUFFER = 5L
  IF NOT KEYWORD_SET(BALM_MASK_WID) THEN BALM_MASK_WID = 5.0D

  maxf = max(stdstar.flux,sciind)

  longslit_dir = getenv('LONGSLIT_DIR')
  if strmid(std_name,2,3,/reverse_offset) eq '10a' then begin
	std_file = longslit_dir+ '/pro/LBT/MODS/Calib_Lib/' + std_name + '.dat'
	; read in standar star spectrum ;ABmag0 = -48.6 c=2.99792458d18A/s
	readcol,std_file,wave_std,abmag_std,step
	flux_std = 10^((abmag_std+48.6)/(-2.5))*(2.99792458d18/(wave_std^2))*1.0d17
	plot,wave_std,flux_std,psym=4

	;alter parameters for 10A spacing...
	NRESLN=1.D
	MSK_BALM = 1
	BALM_MASK_WID = 3
	line_balm = [4102.8 $                        ;; (H9,CaII H, Hf,H-delta)
                   , 4341.6 $
                   , 4862.7 $ ;; (H-gamma,H-beta,H-alpha,Fe)
                   , 5407.0, 6564.6, 8224.8, 8239.2]
  endif else begin
	std_file = longslit_dir+ '/calib/standards/calspec/' + std_name + '.fits.gz'
	; read in standar star spectrum
	std = xmrdfits(std_file, 1,/silent)
	wave_std = std.WAVELENGTH
	flux_std = 1.0d17*std.FLUX   ; fluxes are in units of 1.0e-17 erg/cm^2/s/A
;	if (std_name eq 'g191b2b_stisnic_002') $
;	  or (std_name eq 'g191b2b_mod_005') $
;	  or (std_name eq 'feige34_stis_001') $
;          or (std_name eq 'bd_28d4211_stis_001') $
;	  then begin 
;		MSK_BALM = 1
;	endif
	if (std_name eq 'bd_33d2642_fos_003') then begin
		NRESLN=5
		MSK_BALM = 1
		BALM_MASK_WID=1.
	endif
        if (std_name eq 'feige66_002') then begin
                NRESLN=5
                MSK_BALM = 1
                BALM_MASK_WID=1.
        endif

  endelse 

;addin far red pieces
  if (std_name eq 'bd_33d2642_fos_003') then begin
        std_file2 = longslit_dir+ '/pro/LBT/MODS/Calib_Lib/bd33d2642_10a.dat'
        readcol,std_file2,wave_std2,abmag_std,step
        flux_std2 = 10^((abmag_std+48.6)/(-2.5))*(2.99792458d18/(wave_std2^2))*1.0d17

	plot,wave_std,flux_std,xrange=[6000,10000]
	oplot,wave_std2,flux_std2,psym=4	
	addin = where(wave_std2 gt 9100)
	oplot,wave_std2[addin],flux_std2[addin],psym=4,color=cgcolor('red')

	orig = where(wave_std lt 9100)
	print,'Supplementing FOS spectrum in the far red'
	wave_std = [wave_std[orig],wave_std2[addin]]
	flux_std = [flux_std[orig],flux_std2[addin]]
  endif
  if (std_name eq 'feige66_002') then begin
        std_file2 = longslit_dir+ '/pro/LBT/MODS/Calib_Lib/feige66_10a.dat'
        readcol,std_file2,wave_std2,abmag_std,step
        flux_std2 = 10^((abmag_std+48.6)/(-2.5))*(2.99792458d18/(wave_std2^2))*1.0d17

        plot,wave_std,flux_std,xrange=[6000,10000]
        oplot,wave_std2,flux_std2,psym=4
        addin = where(wave_std2 gt 9100)
        oplot,wave_std2[addin],flux_std2[addin],psym=4,color=cgcolor('red')

        orig = where(wave_std lt 9100)
        print,'Supplementing FOS spectrum in the far red'
        wave_std = [wave_std[orig],wave_std2[addin]]
        flux_std = [flux_std[orig],flux_std2[addin]]
   endif




; uncalibrated observed spectrum
wave = stdstar.wave
flux = stdstar.flux
ivar = stdstar.ivar
alt_flux = stdstar.flux2

IF KEYWORD_SET(INMASK) THEN BEGIN
   IF n_elements(INMASK) NE n_elements(ivar) THEN $
      message, 'Your mask must be aligned with the std spectrum'
   badpix = WHERE(INMASK EQ 0, nbad)
   IF nbad GT 0 THEN ivar[badpix] = 0.0
ENDIF

; parse headers and read in extinction file
extfile = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/mods_extiction.dat'
readcol,extfile,ext_wave,ext_val,/silent
ext_coef = bspline_iterfit(ext_wave,ext_val,bkspace=3,nord=3)
ext = bspline_valu(wave,ext_coef)
ext = 10.0D^(0.4D*ext*airmass)/exptime  ;; TAKE OUT EXPTIME TOO
;ext = long_extinct(wave, scihdr, AIRMASS = AIRMASS, EXPTIME = EXPTIME)

; extinction correct data and divide by exposure time
flux = flux*ext
ivar = ivar/ext^2
alt_flux = alt_flux*ext

; find the min and max of the calibration spectrum
wave_min_std = min(wave_std)
wave_max_std = max(wave_std)
; find the min and max of the standard
ind_sort = sort(wave)
wave_min_obs = min(wave)        ;wave[ind_sort[LO_BUFFER-1]]
wave_max_obs = max(wave)        ;wave[ind_sort[nwave-1L-HI_BUFFER]]

wave_min = wave_min_std > wave_min_obs
wave_max = wave_max_std < wave_max_obs 

calib_inds = WHERE(wave GE wave_min AND wave LE wave_max)
wave = wave[calib_inds]
flux = flux[calib_inds]
ivar = ivar[calib_inds]
alt_flux = alt_flux[calib_inds]

sort_ind = sort(wave)
wave = wave[sort_ind]
flux = flux[sort_ind]
ivar = ivar[sort_ind]
alt_flux = alt_flux[sort_ind]

nwave = n_elements(wave)
;; Don't use the edges in the fit
;IF KEYWORD_SET(LO_BUFFER) THEN ivar[0:LO_BUFFER] = 0
;IF KEYWORD_SET(HI_BUFFER) THEN ivar[nwave-1L-HI_BUFFER:nwave-1L] = 0

;interpolate calbiration spectrum onto observed wavelengths
flux_std_int = interpol(flux_std, wave_std, wave)

; Compute an effective resolution for the standard. This could be improved
; to setup an array of breakpoints based on the resolution. At the 
; moment we are using only one number
std_res = 2.0*djs_median(abs(wave_std - shift(wave_std, 1)))
resln = std_res

if keyword_set(MSK_BALM) then begin
   IF NOT KEYWORD_SET(LINE_BALM) THEN $
      line_balm = [3836.4 $
                   , 3890.1 $
                   , 3969.6 $
                   , 4102.8 $                        ;; (H9,CaII H, Hf,H-delta)
                   , 4341.6 $
;;           , 4687.3 $
                   , 4862.7 $ ;; (H-gamma,H-beta,H-alpha,Fe)
                   , 5407.0 $
                   , 6564.6 $
                   , 8224.8 $
                   , 8239.2] 
   nbalm = n_elements(line_balm)
   for qq = 0L, nbalm-1 do begin
      mskwv = where(abs(wave-line_balm[qq]) LE BALM_MASK_WID*resln, nbad)
      if nbad NE 0 then ivar[mskwv] = 0.
   endfor
endif


if std_name eq 'bd_33d2642_fos_003' then $
        line_balm = (wave GE 3714.0D AND wave LE 3723.0D) OR $
               (wave GE 3650.0D AND wave LE 3820.0D) OR $
               (wave GE 3725.0D AND wave LE 3735.0D) OR $
               (wave GE 3740.0D AND wave LE 3755.0D) OR $
               (wave GE 3760.0D AND wave LE 3775.0D) OR $
               (wave GE 3785.0D AND wave LE 3806.0D) OR $
               (wave GE 3810.0D AND wave LE 3820.0D) OR $
               (wave GE 3824.0D AND wave LE 3841.0D) OR $
	       (wave GE 3880.0D AND wave LE 3895.0D) OR $
               (wave GE 3957.0D AND wave LE 3979.0D) OR $
               (wave GE 4000.0D AND wave LE 4030.0D) OR $
               (wave GE 4087.0D AND wave LE 4120.0D) OR $
               (wave GE 4135.0D AND wave LE 4145.0D) OR $
               (wave GE 4328.0D AND wave LE 4355.0D) OR $
               (wave GE 4677.0D AND wave LE 4692.0D) OR $
               (wave GE 4830.0D AND wave LE 4931.0D) OR $
               (wave GE 5402.0D AND wave LE 5417.0D) OR $
	       (wave GE 6535.0D AND wave LE 6590.0D) $
  else line_balm = (wave GE 3718.0D AND wave LE 3727.0D) OR $
               (wave GE 3731.0D AND wave LE 3740.0D) OR $
               (wave GE 3745.0D AND wave LE 3757.0D) OR $
               (wave GE 3764.0D AND wave LE 3779.0D) OR $
               (wave GE 3790.0D AND wave LE 3808.0D) OR $
               (wave GE 3816.0D AND wave LE 3824.0D) OR $
               (wave GE 3824.0D AND wave LE 3851.0D) OR $
               (wave GE 3959.0D AND wave LE 3979.0D) OR $
               (wave GE 4087.0D AND wave LE 4118.0D) OR $
               (wave GE 4328.0D AND wave LE 4355.0D) OR $
               (wave GE 4677.0D AND wave LE 4692.0D) OR $
               (wave GE 4830.0D AND wave LE 4931.0D) OR $
               (wave GE 5402.0D AND wave LE 5417.0D) OR $
               (wave GE 6535.0D AND wave LE 6590.0D)


        balm_ind = where(line_balm,nbalm)
        IF nbalm GT 0 THEN ivar[balm_ind] = 0.0

;; Mask telluric absorption
if not keyword_set(no_tell_mask) then begin
	tell = (wave GE 7580.0D AND wave LE 7720.0D) OR $
	       (wave GE 7160.0D AND wave LE 7340.0D) OR $
	       (wave GE 6850.0  AND wave LE 6930.0D)
	tell_ind = where(tell, ntell)
	IF ntell GT 0 THEN ivar[tell_ind] = 0.0
endif
finalmask = where(ivar GT 0)

;plot,wave,flux,/ylog,yrange=[100,2e4],/ystyle
;oplot,wave[finalmask],flux[finalmask],color=cgcolor('dodger blue')
;oplot,wave[finalmask],alt_flux[finalmask],color=cgcolor('red')
;stop


if (keyword_set(nofit) eq 0) then begin ; jm09dec19ucsd
   mag_set = bspline_magfit(wave, flux, ivar, flux_std_int $
     , bkspace = resln*nresln $
     , maxiter = 10, maxrej = 5, upper = 3, lower = 3 $
     , sensfit = sensfit, sensfunc = sensfunc $
     , wave_min = wave_min, wave_max = wave_max $
     , outmask = outmask) 

   mag_set2 = bspline_magfit(wave, alt_flux, ivar, flux_std_int $
     , bkspace = 1*resln*nresln $
     , maxiter = 10, maxrej = 5, upper = 3, lower = 3 $
     , sensfit = sensfit2, sensfunc = sensfunc2 $
     , wave_min = wave_min, wave_max = wave_max $
     , outmask = outmask)

;;;;; Test Plots ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plot,wave,flux
; plot,wave,sensfit,color=cgcolor('dodger blue')   
; oplot,wave,flux_std_int/sensfit,color=cgcolor('purple')
; stop  

; std_file
;file3= '/usr/custom/xidl/xidl/Spec/Longslit/calib/standards/calspec/g191b2b_stisnic_002.fits.gz'
;file3= '/usr/custom/xidl/xidl/Spec/Longslit/calib/standards/calspec/feige34_stis_001.fits.gz'
;file3= '/usr/custom/xidl/xidl/Spec/Longslit/calib/standards/calspec/bd_33d2642_fos_003.fits.gz'
;file3= '/usr/custom/xidl/xidl/Spec/Longslit/calib/standards/calspec/feige66_002.fits.gz'
;file3= '/usr/custom/xidl/xidl/Spec/Longslit/calib/standards/calspec/hz43_mod_005.fits.gz'
;s3 = mrdfits(file3,1,hdr)
;std = s3.flux*1.0d17
;w1 = indgen(5601)*0.5+3000.D
;w1 = indgen(9001)*0.5 + 5500.D
;flux_std_hst = interpol(std, s3.wavelength, w1)

; plot,wave,flux_std_int,xrange=[3200,5900],xstyle=1,xtitle='Wavelength [A]',ytitle='Flux*1E17',charsize=1.5
;plot,wave,flux_std_int,xrange=[4300,4900],xstyle=1,xtitle='Wavelength [A]',ytitle='Flux*1E17',charsize=1.5
;plot,wave,flux_std_int,xrange=[5500,10100],xstyle=1,xtitle='Wavelength [A]',ytitle='Flux*1E17',charsize=1.5
;plot,wave,flux_std_int,xrange=[6500,8000],xstyle=1,xtitle='Wavelength [A]',ytitle='Flux*1E17',charsize=1.5
;plot,wave,flux_std_int,xrange=[3300,4200],xstyle=1,xtitle='Wavelength [A]',ytitle='Flux*1E17',charsize=1.5,/ystyle

; oplot,wave,flux*sensfit,color=cgcolor('red')
; oplot,wave,alt_flux*sensfit2,color=cgcolor('orange')
; oplot,[5755,5755],[0,10000],linestyle=1
; oplot,w1,flux_std_hst,color=cgcolor('gold')
;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if keyword_set(CHK) then x_splot, wave, sensfit, /blo 
   IF KEYWORD_SET(sensfuncfile) THEN mwrfits, mag_set, sensfuncfile, /create
   IF KEYWORD_SET(sensfuncfile) THEN ssflen = STRLEN(sensfuncfile)
   IF KEYWORD_SET(sensfuncfile) THEN sensfuncfile2 = STRMID(sensfuncfile, 0, ssflen-8) + 'ssf2.fits'
   IF KEYWORD_SET(sensfuncfile) THEN mwrfits, mag_set2, sensfuncfile2, /create
endif else return, 1

alt = mag_set2
RETURN, mag_set
END
