;+
; NAME:
;   mods_slits2mask
;
; PURPOSE:
;   convert the slit-structure into a slit-mask image.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   
; REVISION HISTORY:
;   2013         Re-Written for MODS data by K. Croxall
;-  
;------------------------------------------------------------------------------

function mods_slits2mask, tset_slits1, xshift=xshift, nslit=nslit $
                          , VERBOSE = VERBOSE, ONLYSLITS = ONLYSLITS

   if (size(tset_slits1,/tname) NE 'STRUCT' $
    OR n_elements(tset_slits1) NE 2) then $
    message, 'TSET_SLITS must be a 2-element structure'

   tset_slits = tset_slits1
   if (keyword_set(xshift)) then $
    tset_slits.coeff[0,*] = tset_slits.coeff[0,*] + xshift

   dims = tset_slits[0].dims
   nx = dims[0]
   ny = dims[1]

   xx1 = tset_slits[0].xx1
   xx2 = tset_slits[1].xx2

   if (size(xx1,/n_dimen) EQ 1) then nslit = 1 $
    else nslit = (size(xx1,/dimens))[1]

;fill out curve for prism mode....
   IF (size(xx1[*,0],/dim) lt 8192) then begin
	xx1_prime = dblarr(8192,nslit)
        xx2_prime = dblarr(8192,nslit)
        yy1_prime = dblarr(8192,nslit)
        yy2_prime = dblarr(8192,nslit)

	initial_pix = yy1[0,*]
	for incslit = 0,nslit-1 do begin
	  ip = initial_pix[incslit]
	  for pix=ip,ip+489 do xx1_prime[pix,incslit] = xx1[pix-ip,incslit]
          for pix=ip,ip+489 do xx2_prime[pix,incslit] = xx2[pix-ip,incslit]
          for pix=ip,ip+489 do yy1_prime[pix,incslit] = yy1[pix-ip,incslit]
          for pix=ip,ip+489 do yy2_prime[pix,incslit] = yy2[pix-ip,incslit]
	endfor
	xx1 = xx1_prime
        xx2 = xx2_prime
        yy1 = yy1_prime
        yy2 = yy2_prime
   ENDIF
;;;;

   IF KEYWORD_SET(ONLYSLITS) THEN BEGIN
      nloop = n_elements(onlyslits)
      doslits = onlyslits - 1L
   ENDIF ELSE BEGIN
      nloop = nslit
      doslits = lindgen(nslit)
   ENDELSE

   ; Generate the mask image
   slitmask = intarr(dims)
   for ii=0L, nloop-1L do begin
      islit = doslits[ii]
      for iy=0L, ny-1L do begin
         x1 = round(xx1[iy,islit])
         x2 = round(xx2[iy,islit])
         if (x1 GE x2) AND KEYWORD_SET(VERBOSE) THEN $
          splog, 'WARNING: Slit start and end positions appear to cross!'
         x1 = x1 > 0
         x2 = x2 < (nx-1)
         if (x1 LE x2) then begin
             if (total(slitmask[x1:x2, iy]) GT 0) AND KEYWORD_SET(VERBOSE) then $
               splog, 'WARNING: Slits appear to overlap!'
             slitmask[x1:x2, iy] = islit+1
         endif
      endfor
   endfor

   return, slitmask
end
;------------------------------------------------------------------------------
  
