;+
; modsADCalc - compute the differential atmospheric refraction for an observation
;
; Computes the differential atmospheric refraction (aka atmospheric dispersion)
; given an observing site (latitude), target declination, and midpoint of hour angle.
;
; For a given slit position angle and width, it computes the predicted xy
; deflections. 
;
; Inputs from the header:
;     siteID    Site name text (e.g., "LBT")
;     siteLat   Site latitude in decimal degrees
;     P         Atmospheric Pressure in hPa
;     T         Air Temperature in Celsius
;     RH        Relative Humidity in %
;     objDec    Object declination, may be decimal or sexagesimal
;     slitPA    Slit position angle in decimal degrees, North-to-East
;     slitWidth Slit width in arcseconds
;     haStart   Starting Hour Angle in decimal degrees
;     lamGuide  Guide reference wavelength in nm
; 
; Models:
;   The atmospheric refractive index model is from Owens, J.C. 1967,
;   Ap.Opt. 6, 51 and uses the Buck 1996 formula for the water vapor
;   saturation pressure (Buck, A.L. 1981 J.App.Met 20, 1527 the 1996
;   updated formula).  Inputs required are the air pressure in hPa,
;   temperature in Celsius, and relative humidity.
;
;   For most applications, mean values derived from historical
;   meteorology data or at least pressure and temperature of the US
;   Standard Atmosphere model (1976 Model) are appropriate.  The
;   sensitivity to humidity is small if the site is reasonably dry
;   (e.g. 10-30% RH).
;
; EXAMPLE:
;    modsadcalc,'Science/2Dsub_blue_NGC628_F1.fits',5007,blue_wav1d,xdelt,ydelt
;
; Author:
;    K. Croxall, OSU Astronomy Dept.
;    croxall@astronomy.ohio-state.edu
;    2013 May 20
;    based on modsADCalc.pl by R. Pogge 
;##############################################################################

; Atmospheric Refractive Index Model
;   Sources:  Owens, J.C. 1967, Ap.Opt., 6, 51
;             Buck, A.L. 1981, J.App.Met. 20, 1527; 1996 Update Formula

;=========================================================================
; satPres() - compute the saturation pressure of water given T in C
;             uses the Buck 1996 formula valid for T=-50..+80 C
;=========================================================================

function satPres,Tc
    Psat = 6.1121*exp((18.679-(Tc/234.5))*(Tc/(257.4+Tc)))
    return,Psat
end

;=========================================================================
; dryDens() - compute the dry-air density factor given T in C, and Pdry in hPa
;=========================================================================

function dryDens,degC,Pdry
    degK = degC + 273.16;   convert C to K
    Ds = (pDry/degK)*(1.0 + pDry*((5.790E-7) - (9.3250E-4/degK) + (0.25844/(degK^2))));
    return,Ds
end

;=========================================================================
; wetDens() - compute the water-vapor density factor given T in C and Pwet in hPa
;=========================================================================

function wetDens,degC,Pwet
    degK = degC + 273.16;   convert C to K
    Dw = (pWet/degK)*(1.0 + pWet*(1.0 + 3.7E-4*pWet)*(-2.37321E-3 + (2.23366/degK)$
                - (710.792/(degK^2)) + (77514.1/(degK^3))))
    return,Dw
end

;=========================================================================
; refrac() - compute the index of refraction given lambda in nanometers,
;            temperature in C, pressure in hPa, and relative humidity in %
;=========================================================================

function refrac,lam,degC,Pair,RH
    Psat = satPres(degC)
    Pwet = Psat*(RH/100.0)
    Pdry = Pair-Pwet
    Ddry = dryDens(degC,Pdry)
    Dwet = wetDens(degC,Pwet)
    sig2 = (1000/lam)^2;  # sig has units of 1/micron
    nDry = Ddry*(2371.34 + (683939.7/(130-sig2)) + (4547.3/(38.9-sig2)))
    nWet = Dwet*(6487.31 + 58.058*sig2 - 0.7115*(sig2^2) + 0.08851*(sig2^3));
    n = 1.0e-8*(nDry + nWet)
    return,n
end

;=========================================================================
; stringad - Converts a string of sexigesimal coordinates into decimal degrees.
;	09-AUG-90 Version 1 written by Kerry McQuade
;	20-AUG-90 Put code to account for '-0' back in after it was
;		removed by someone.  E. Deutsch
;	17-JUL-95 Added support for coordinates separated by colons, e.g.
;		17:00:45.2 25:4:32.4, which IRAF uses.  E. Deutsch
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

;                        Remove any gaps between '-' or '+' and numeral  

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
; parAngle(lat,dec,HA) - compute the parallactic angle from lat, dec, & HA
;                        lat = site latitude in decimal degrees
;                        dec = target declination in decimal degrees
;                         HA = hour angle in decimal hours
;
; -parAngle is the celestial position angle along parallactic
;=========================================================================

function parAngle,siteLat,dec,ha
    latRad = siteLat * !dpi/180
    decRad = dec * !dpi/180
    haRad  = (15.0*ha) * !dpi/180 ; convert h->deg->rad
    top = sin(haRad)
    bot = tan(latRad)*cos(decRad) - sin(decRad)*cos(haRad)
    parRad = atan(top,bot)
    parAng = parRad * 180/!dpi
    if (parAng lt -90.0) then parAng = -(parAng + 180.0)
    if (parAng gt 90.0) then parAng = 180.0 - parAng
   return,parAng
end

;=========================================================================
; airmass(lat,dec,HA) - compute airmass (secZ) from lat, dec, and HA
;                       lat = site latitude in decimal degrees
;                       dec = target declination in decimal degrees
;                        HA = hour angle in decimal hours
;=========================================================================

function airmass,siteLat,dec,ha
    latRad = siteLat * !dpi/180
    decRad = dec * !dpi/180
    haRad  = (15.0*ha) * !dpi/180
    cosZD = sin(latRad)*sin(decRad) + cos(latRad)*cos(decRad)*cos(haRad)
    if (cosZD gt 0) then am = 1.0/cosZD
    if (am gt 99) then am = 99.99
    return,am
end

;=========================================================================
; modsadcalc - get the dispersion offsets
;=========================================================================

pro modsadcalc,file,lamGuide,lamarr,xdelt,ydelt,verbose=verbose
  if not keyword_set(verbose) then !Quiet=1

  version = 'modsADCalc IDL v1.3.1';
  verDate = '2013-05-17';
  image = mrdfits(file,0,hdr,/silent)

	; These the entries we care about

  siteID = strcompress(sxpar(hdr[*, 0], 'TELESCOP'), /rem)
  siteLat = strcompress(sxpar(hdr[*, 0], 'LATITUDE'), /rem)
  P = strcompress(sxpar(hdr[*, 0], 'LBTPRES'), /rem)
  T = strcompress(sxpar(hdr[*, 0], 'LBTTEMP'), /rem)
  RH = strcompress(sxpar(hdr[*, 0], 'LBTHUM'), /rem)
  objRA = strcompress(sxpar(hdr[*, 0], 'TELRA'), /rem)
  objDec = strcompress(sxpar(hdr[*, 0], 'TELDEC'), /rem)
  slitPA = strcompress(sxpar(hdr[*, 0], 'POSANGLE'), /rem)
  slitWidth = 1.0 						; make input
  haStart = strcompress(sxpar(hdr[*, 0], 'HA'), /rem)
  exptime = strcompress(sxpar(hdr[*, 0], 'EXPTIME'), /rem)
  amhead = strcompress(sxpar(hdr[*, 0], 'AIRMASS'), /rem)
; need number of exposures possibly to caluclate a middle dispersion
; need an input 1d wavelength vector to give out a correction vector for each slit
; need to pass in the ref/centering wavelength

; Additional input processing
  coords = objRA + ' ' + objDec
  stringad,coords,ra,dec
  decStr = objDec
  decRad = dec * !dpi/180.
  haStart = strsplit(haStart,':',/extract)
  haSign = strmid(haStart[0],0,1)
   if haSign eq '-' then haStart = haStart[0] - haStart[1]/60.D - haStart[2]/3600.D $
   else haStart = haStart[0] + haStart[1]/60.D + haStart[2]/3600.D
  hastop = haStart + exptime / 3600.D
  haMid = haStart + exptime / 3600.D / 2.D

; Validate the hour angle range and step
  if (haStart lt -12) or (haStart gt 12) then begin
    print,'ERROR: Invalid Starting Hour Angle' ;
    print,'Hour angles must be between -12 and +12 hours.';
    stop
  endif

; Validate the Slit PA
  if slitPA lt -180.0 then slitPA = slitPA + 360.
  if slitPA gt 180.0 then slitPA = slitPA -360.
  paRad = slitPA * !dpi/180.

; Validate the atmospheric properties

  if (T lt -50.0) or (T gt 80.0) then begin
    print,'ERROR: Invalid Air Temperature' 
    print,'The Buck (1996) formula for the water vapor saturation pressure is'
    print,'only valid between -50 and +80 C.'
    stop
    endif

  if (P lt 400.0) or (P gt 1200.0) then begin
    print,'ERROR: Invalid Air Pressure'
    print,'The air pressure should be between 400 and 1200 hPa'
    stop
    endif

  if (RH lt 0.0) or (RH gt 100.0) then begin
    print,'ERROR: Invalid Relative Humidity' 
    print,'The relative humidity must be between 0 and 100%'
    stop
    endif

; Compute the refractive index at the guiding wavelength
  lamGuide /= 10.
  nG = refrac(lamGuide,T,P,RH)

; Set the wavelength arrays
  lam = lamarr/10. ; convert wavelengths to nm
  numLam = n_elements(lam)

; calculate the delta array
  x = fltarr(numLam)
  y = fltarr(numLam)
  for j=0,numLam-1 do begin ; loop over wavelengths
    n = refrac(lam[j],T,P,RH)
    yLam=0
    parAng = parAngle(siteLat,dec,haMid)
    paUp = parAng
	if (dec gt siteLat) then begin
	    paUp = -parAng - 180
	    if (paUp gt 0) then  paDown = paUp - 180 $
	    else paDown = paUp + 180
	endif else begin
            paUp = parAng
	    if (paUp gt 0) then paDown = paUp - 180 $
	    else paDown = paUp + 180
	endelse
	am = airmass(siteLat,dec,haMid)
	if (am gt 0)  && (am lt 99) then begin
	    tanZD = tan(acos(1/am))
	    R = 206264.8*tanZD*(n-nG)
	    x[j] =  R*sin((paUp - slitPA) * !dpi/180.)
	    y[j] = -R*cos((paUp - slitPA) * !dpi/180.)
	endif
	if y[j] ne y[j] then y[j] =0
  endfor
  xdelt=x
  ydelt=y
  return
end

