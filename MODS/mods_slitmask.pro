;+
; NAME:
;   mods_slitmask
;
; PURPOSE:
;   Generate slitmask structure for MODS multislit observations
;
; CALLING SEQUENCE:
;   mods_slitmask, mms,nx,ny
;
; INPUTS:
;   mms - name of mms file from the observation
;   nx  - number of pixels in x
;   ny  - number of pixels in y
;   CHANNEL - Red or blue channel observations
;   GRATNAME - name of the grating or Prism used
;   FILENAME - name of the filename
;   INTERACTIVE_SLITS - interactively ajust the locations of the slits
;   AUTOTUNE_SLITS - autotune the slit placement
;   TRANSPOSE_IM - XIDL uses transposed MODS images.  If the input image 
;	has not yet been transposed use this.
;
; OUTPUTS:
;   tset_slits - 2-element array of trace sets, where the first defines
;                the starting slit positions, and the second one defines
;                the ending slit positions
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   This will break if the binning is not 1x1
;
; PROCEDURES CALLED:

; REVISION HISTORY:
;   Dec-2012  Written by K. Croxall
;   Jun-2013  Revised to include prism mode and improve fits
;              Also added vertical shifts to center the mask/slits
;-
;------------------------------------------------------------------------------

;=========================================================================
;  slit_twiddle -- automatically tighten slit placement 
;=========================================================================
FUNCTION slit_twiddle,lower_trace,upper_trace,image
 topcurve = fltarr(101)
 botcurve = fltarr(101)
 offset_find_u = fltarr(101,24)
 offset_find_1 = fltarr(101,24)

;upper edge
for i=0,23 do offset_find_u[*,i] = image[upper_trace[*,0],upper_trace[*,1]+i-15]
delta_u = fltarr(101,24)
for i=1,23 do delta_u[*,i-1] = abs(offset_find_u[*,i] - offset_find_u[*,i-1])
median_delta_u=fltarr(24)
for i=0,23 do median_delta_u[i]=median(delta_u[*,i])
offset_u = where(median_delta_u eq max(median_delta_u))+1
offset_u = offset_u[0] - 15.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;de-bug plots
;tmp=indgen(24)-15.
;plot,tmp,offset_find_u[50,*]
;oplot,[offset_u,offset_u],[10,1e6]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;lower edge
for i=0,23 do offset_find_1[*,i] = image[lower_trace[*,0],lower_trace[*,1]+i-15]
delta_l = fltarr(101,24)
for i=1,23 do delta_l[*,i-1] = abs(offset_find_1[*,i] - offset_find_1[*,i-1])
median_delta_l=fltarr(24)
for i=0,23 do median_delta_l[i]=median(delta_l[*,i])
offset_l = where(median_delta_l eq max(median_delta_l))+1
offset_l = offset_l[0] - 15.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;de-bug plots
;tmp=indgen(24)-15.
;plot,tmp,offset_find_1[50,*]
;oplot,[offset_l,offset_l],[10,1e6]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

offset = (offset_u + offset_l)/2.
oldoffset = offset

 tol = 2.e29
 while tol gt 0 do begin
	for i=0,100 do topcurve[i] = image[upper_trace[i,0],upper_trace[i,1]+offset]
        for i=0,100 do botcurve[i] = image[lower_trace[i,0],lower_trace[i,1]+offset]
         if abs(total(topcurve-botcurve)) lt tol then begin
         	tol = abs(total(topcurve-botcurve))
                oldoffset = offset
		if tol GT 1e6 then step = 3. else step = 1.
                if total(topcurve-botcurve) lt 0 then offset-= step
                if total(topcurve-botcurve) gt 0 then offset+= step
         endif else begin
         	offset = oldoffset
                tol=-1
         endelse
 endwhile
; if abs(offset) gt 3 then begin
;	splog,'Found a large auto-offset of '+string(offset)
;	splog,'This may be caused by light from a nearby slit. Adopting no shift'
;	offset = 0
; endif

 splog,'Auto-offset of '+string(offset)+' found'
 RETURN,offset
END

;=========================================================================
;  man_twiddle -- manually check and adjust slit placement  
;=========================================================================
FUNCTION man_twiddle,xarray,lower,upper,image,slitno,beg_x,end_x
 offset=0
 ; plot slit and traces
 ncolors = !D.Table_Size
 device,DECOMPOSED=0
 LoadCT,39
 min_val = min(lower[beg_x:end_x])-10. > 0.
 max_val = max(upper[beg_x:end_x])+10. <3087
 slit_image = (BytScl(image[beg_x:end_x,min_val:max_val]));min=-5,max=60,/nan))

 Ok=0
 LeSigh: ;print,'Applying manual offset of ',offset,' to slit #',slitno

 TVImage,slit_image,Position=[0.05,0.08,0.99,0.99],/erase
 plot,xarray,lower-min_val,/noerase,Position=[0.05,0.08,0.99,0.99],$
     	xstyle=1,xrange=[beg_x,end_x],ystyle=1,yrange=[0,max_val-min_val],$
	xtitle='Pixels Along Dispersion Axis',ytitle='Pixels Relative To Spatial Axis'
 plot,xarray,upper-min_val,/noerase,Position=[0.05,0.08,0.99,0.99],$
        xstyle=1,xrange=[beg_x,end_x],ystyle=1,yrange=[0,max_val-min_val]

 plot,xarray,lower-min_val+offset,/noerase,Position=[0.05,0.08,0.99,0.99],$
        xstyle=1,xrange=[beg_x,end_x],ystyle=1,yrange=[0,max_val-min_val],linestyle=2
 plot,xarray,upper-min_val+offset,/noerase,Position=[0.05,0.08,0.99,0.99],$
        xstyle=1,xrange=[beg_x,end_x],ystyle=1,yrange=[0,max_val-min_val],linestyle=2

 while Ok ne 1 do begin
 	print,'Accept Slit? (y/n/abort)'
 	redo = ''
 	read,': ',redo
 	if redo eq 'n' then begin
 		Print,'Vertically Offset slit by (pixels):'
 	        read,': ',offset
 	        GOTO,LeSigh
 	endif else if redo eq 'y' then begin
 	        print,'Slit Accepted with a manual offset of ',offset
	        Ok = 1
	endif else if redo eq 'abort' then begin
                stop
        endif else begin
                print,'I do not understand ',redo
        endelse
 endwhile
 RETURN,offset
END

;=========================================================================
;  mods_slitmask -- Create the MODS-XIDL slit structure 
;=========================================================================
FUNCTION MODS_SLITMASK, mms, nx, ny, $
	CHANNEL= CHANNEL,GRATNAME=GRATNAME,$
	FILENAME=FILENAME,INTERACTIVE_SLITS=INTERACTIVE_SLITS, $
	AUTOTUNE_SLITS=AUTOTUNE_SLITS,TRANSPOSE_IM=TRANSPOSE_IM


IF keyword_set(INTERACTIVE_SLITS) then $
	if not WindowAvailable(1) then window,1,retain=2,$
	title='Interactive Slit Alignment'

image = mrdfits(filename,0,maskhd,/silent)
IF keyword_set(TRANSPOSE_IM) then image = transpose(image)

;if strmid(gratname,0,1) eq 'P' then readcol,mms,input,val,FORMAT='A,A',/silent,skipline=4 $
;	else readcol,mms,input,val,FORMAT='A,A',/silent,skipline=13
readcol,mms,input,val,FORMAT='A,A',/silent,skipline=13

global = 0
global_offset = 0
rot_b=0.
rot_r=0.
yshift_b=0.
yshift_r=0.

col1 = input
col2 = input
col3 = input
dim=size(input, /dimensions)
for n=0,dim[0]-1 do begin
        parts=str_sep(input[n],'.')
        col1[n]=parts[0]
        col2[n]=parts[1]
        col3[n]=parts[2]
endfor

tarRA_sex = 'null'
tarDEC_sex = 'null'
tarRA_deg = 999
tarDEC_deg = 999
pangle = 999

; extract the rotation and target center
for n=0,dim[0]-1  do begin
    if col2[n] eq 'ROT' then pangle = strtrim(string(val[n],Format='(F20.1)'),2)
    if col2[n] eq 'TARG' then begin
        if col3[n] eq 'ALPHA' then begin
            ;taralpha = string(val[n],Format='(F20.3)')
            ;tardelta = string(val[n+1],Format='(F20.3)')
            taralpha = val[n]
            tardelta = val[n+1]
		radot = strpos(taralpha,'.')
		ra1=strmid(taralpha,radot-6,2)
		ra2=strmid(taralpha,radot-4,2)
		ra3=strmid(taralpha,radot-2)
		decdot = strpos(tardelta,'.')
		dec1=strmid(tardelta,decdot-7,3)
		dec2=strmid(tardelta,decdot-4,2)
		dec3=strmid(tardelta,decdot-2)
		tarRA_deg = 15 * (double(ra1) + double(ra2)/60 + double(ra3)/3600)
		tarDEC_deg = double(dec1) + double(dec2)/60 + double(dec3)/3600
		if strmid(dec1,0,1) eq '-' then tarDEC_deg = $
			double(dec1) - double(dec2)/60 - double(dec3)/3600
		tarRA_sex = ra1 + ':' + ra2 + ':' + ra3
		tarDEC_sex = dec1 + ':' + dec2 + ':' + dec3
        endif
    endif
endfor

; count the number of science slits
slitno=0
for n=0,dim[0]-1  do begin
        if col3[n] eq 'WIDMM' then begin
                WIDMM=val[n]
                LENMM=val[n+1]
                if WIDMM gt 0.6 or LENMM gt 0.6 and WIDMM ne LENMM then slitno+=1
        endif
endfor
splog,'Found '+string(slitno)+' slits in the mms file.'

; set up the slit structure
tset_proto = $
  { func    :    'legendre', $		;used in fitting coeffs.... outdated?
    xmin    :    0.0, $			
    xmax    :    float(ny-1), $
;    coeff   :    dblarr(2, slitno), $	;used in fitting coeffs.... outdated?
    dims    :    long([nx, ny]), $
    tarRA   :    tarRA_sex, $
    tarDEC  :    tarDEC_sex, $
    apRA    :    strarr(slitno), $
    apDEC   :    strarr(slitno), $
    tarRA2  :    tarRA_deg, $
    tarDEC2 :    tarDEC_deg, $
    apRA2   :    dblarr(slitno), $
    apDEC2  :    dblarr(slitno), $
    posang  :    pangle, $ 
    length  :    0.0D, $
    width   :    0.0D, $
    pix_beg :    dblarr(slitno), $
    pix_end :    dblarr(slitno),  $
    xx1     :    dblarr(8192,slitno), $
    xx2     :    dblarr(8192,slitno), $
    XMM_arr :    dblarr(slitno), $
    YMM_arr :    dblarr(slitno), $
    LENMM_arr:   dblarr(slitno), $
    WIDMM_arr:   dblarr(slitno) $
  }
tset_slits = replicate(tset_proto, 2)

; fill the structure with the slit traces and positions
slitno=0
for n=0,dim[0]-1  do begin
    if col3[n] eq 'WIDMM' then begin
	; get slit info from MMS
	if val[n-6] eq 'refslit' then continue
        WIDMM=val[n]
        LENMM=val[n+1]
        XMM=val[n+2]
        YMM=val[n+3]
	WID = strtrim(string(val[n-5],Format='(F20.1)'),2)
	LEN = strtrim(string(val[n-4],Format='(F20.1)'),2)
	ALPHA=val[n-2]
	DELTA=val[n-1]
	;select slits of interest
        if WIDMM gt 0.6 or LENMM gt 0.6 and WIDMM ne LENMM then begin
            if ((CHANNEL eq 'MODS1B') or (CHANNEL eq 'MODS2B')) then begin
                blue_lam=fltarr(101)
                coord_pix1=fltarr(101,2)
                coord_pix0=fltarr(101,2)
		if GRATNAME eq 'G450L' then begin
		  splog,'No G450L grating exists.'
		  splog,'there may have been an issue in the header (this is known to have been mis-entered once in the past).'
		  splog,'please check your data to be sure you are using the MODS Blue grating (G400L)'
		  print,'To adopt G400L, you may .continue'
		  stop
		  GRATNAME = 'G400L'
		endif
                if GRATNAME eq 'G400L' then begin
		  tset_slits[0].pix_beg[slitno] = 0.0
		  tset_slits[0].pix_end[slitno] = 8191
                  tset_slits[1].pix_beg[slitno] = 0.0
                  tset_slits[1].pix_end[slitno] = 8191
		  ;calculate the slit trace
                  for i=0,100 do begin
                      blue_lam[i]=3000+28*i
                      coord_pix1[i,*]=blue_map(CHANNEL,XMM,-(YMM-LENMM/2.+XMM*tan(rot_b)),blue_lam[i])
                      coord_pix0[i,*]=blue_map(CHANNEL,XMM,-(YMM+LENMM/2.+XMM*tan(rot_b)),blue_lam[i])
                  endfor
                  ;catch slits going off the detector
                  if max(coord_pix1[*,0]) gt 8191 then begin
                        print,'The trace for Slit number ',slitno,' goes off the detector'
			print,'The trace goes off the red end of the chip.'
			print,'PROCEDE WITH CAUTION - suggestion: plot,coord_pix1[*,0]'
                        stop
                  endif
                  if max(coord_pix1[*,1]) gt 3088 then begin
                        print,'The trace for Slit number ',slitno,' goes off the detector'
			print,'The top of the slit leaves the top of the chip.'
                        print,'PROCEDE WITHsuggestion: suggestion: plot,coord_pix1[*,1]'
                        stop
                  endif
                  if min(coord_pix0[*,0]) lt 0 then begin
                        print,'The trace for Slit number ',slitno,' goes off the detector'
			print,'The traces goes of the blue end of the chip.'
                        print,'PROCEDE WITH CAUTION - suggestion: plot,coord_pix0[*,0]'
                        print,'possible fix if only one pixel is off: coord_pix0[0,0] = 0'
                        stop
                  endif
                  if min(coord_pix0[*,1]) lt 0 then begin
                        print,'The trace for Slit number ',slitno,' goes off the detector'
			print,'The trace goes off the bottom of the chip.'
                        print,'PROCEDE WITH CAUTION - suggestion: plot,coord_pix0[*,1]'
                        stop
                  endif

                  ; account for y-shift of mask
		  IF keyword_set(AUTOTUNE_SLITS) then begin
		  splog,'Autotune location of slit'+string(slitno + 1)
		  offset = slit_twiddle(coord_pix0,coord_pix1,image)
                  	coord_pix1[*,1] += offset
                  	coord_pix0[*,1] += offset
		  endif

                  ;interpolate curve edges
		  xarr = indgen(8192)
		  interp_xx1 = interpol(coord_pix0[*,1],coord_pix0[*,0],xarr) + global_offset
                  interp_xx2 = interpol(coord_pix1[*,1],coord_pix1[*,0],xarr) + global_offset

		  IF keyword_set(INTERACTIVE_SLITS) then begin
			offset_int = man_twiddle(xarr,interp_xx1,interp_xx2,image,slitno,0,8191)
			Ok=0
			while Ok ne 1 do begin 
				print,'Apply Global Offset? (y/n/abort)'
				global = ''
				read,': ',global
				if global eq 'y' then begin
					INTERACTIVE_SLITS = 0
					AUTOTUNE_SLITS = 0
					global_offset = offset + offset_int
					ok=1
				endif else if global eq 'n' then begin
					ok = 1
				endif else if redo eq 'abort' then begin
					stop
				endif else begin
					print,'I do not understand ',redo 
				endelse
			endwhile
        	        interp_xx1 += offset_int
	                interp_xx2 += offset_int
		  endif
		  tset_slits[0].xx1[*, slitno] = interp_xx1
		  tset_slits[1].xx2[*, slitno] = interp_xx2

                  ;set the fit coeffs for slit tracing
                  coeffs1=svdfit(coord_pix1[*,0],coord_pix1[*,1],2,func='flegendre')
;                  tset_slits[1].coeff[*, slitno] = coeffs1
                  coeffs0=svdfit(coord_pix0[*,0],coord_pix0[*,1],2,func='flegendre')
;                  tset_slits[0].coeff[*, slitno] = coeffs0

                endif else if GRATNAME eq 'P450L' then begin
		  ; find the prism trace for the slit
                  for i=0,100 do begin
                      blue_lam[i]=3200+35*i
                      coord_pix1[i,*]=blue_map_prism(CHANNEL,XMM,-(YMM-LENMM/2.+XMM*tan(rot_b)),blue_lam[i])
                      coord_pix0[i,*]=blue_map_prism(CHANNEL,XMM,-(YMM+LENMM/2.+XMM*tan(rot_b)),blue_lam[i])
                  endfor
                  ;catch slits going off the detector
                  if max(coord_pix1[*,0]) gt 8191 then begin
                        print,'The trace for Slit number ',slitno+1,' goes off the detector'
                        stop
                  endif
                  if max(coord_pix1[*,1]) gt 3088 then begin
                        print,'The trace for Slit number ',slitno+1,' goes off the detector'
                        stop
                  endif
                  if min(coord_pix0[*,0]) lt 0 then begin
                        print,'The trace for Slit number ',slitno+1,' goes off the detector'
                        stop
                  endif
                  if min(coord_pix0[*,1]) lt 0 then begin
                        print,'The trace for Slit number ',slitno+1,' goes off the detector'
                        stop
                  endif
		  ; account for prism mode not using the full chip
                  coord_pix0[*,0] += 2048
                  coord_pix1[*,0] += 2048
		  ; account for y-shift of mask
                  IF keyword_set(AUTOTUNE_SLITS) then begin
                  splog,'Autotune location of slit'+string(slitno + 1)
                  offset = slit_twiddle(coord_pix0,coord_pix1,image)
                        coord_pix1[*,1] += offset
                        coord_pix0[*,1] += offset
                  endif
		  ;set the begining and end of the prism slits
                  tset_slits[0].pix_beg[slitno] = coord_pix0[0,0]
                  tset_slits[0].pix_end[slitno] = coord_pix0[100,0]
                  tset_slits[1].pix_beg[slitno] = coord_pix0[0,0]
                  tset_slits[1].pix_end[slitno] = coord_pix0[100,0]

		  ;interpolate curve edges
                  xarr = indgen(8192)
                  interp_xx1 = interpol(coord_pix0[*,1],coord_pix0[*,0],xarr) + global_offset
                  interp_xx2 = interpol(coord_pix1[*,1],coord_pix1[*,0],xarr) + global_offset
	          ;blank the begining and end of the trace (otherwise prism slits overlap)
                  IF keyword_set(INTERACTIVE_SLITS) then begin
                        offset_int = man_twiddle(xarr,interp_xx1,interp_xx2,image,slitno,coord_pix0[0,0],coord_pix0[100,0])
                        Ok=0
                        while Ok ne 1 do begin
                                print,'Apply Global Offset? (y/n/abort)'
                                global = ''
                                read,': ',global
                                if global eq 'y' then begin
                                        INTERACTIVE_SLITS = 0
                                        AUTOTUNE_SLITS = 0
                                        global_offset = offset + offset_int
                                        ok=1
                                endif else if global eq 'n' then begin
                                        ok = 1
                                endif else if redo eq 'abort' then begin
                                        stop
                                endif else begin
                                        print,'I do not understand ',redo
                                endelse
                        endwhile
                        interp_xx1 += offset_int
                        interp_xx2 += offset_int
                  endif
		  for i = 0,coord_pix0[0,0]-1 do interp_xx1[i] = 0
		  for i = coord_pix0[100,0]+1,8191 do interp_xx1[i] = 0
                  for i = 0,coord_pix0[0,0]-1 do interp_xx2[i] = 0
                  for i = coord_pix0[100,0]+1,8191 do interp_xx2[i] = 0
                  tset_slits[0].xx1[*, slitno] = interp_xx1
                  tset_slits[1].xx2[*, slitno] = interp_xx2

                  ;set the fit coeffs for slit tracing
;                  coeffs1=poly_fit(coord_pix1[*,0],coord_pix1[*,1],1)
;                  tset_slits[1].coeff[*, slitno] = coeffs1
;                  coeffs0=poly_fit(coord_pix0[*,0],coord_pix0[*,1],1)
;                  tset_slits[0].coeff[*, slitno] = coeffs0

                endif else begin
                      print,'Unknown Grating/Prism for the Blue Channel'
                      stop
                endelse
		
;;;;;;;;;;;PLOT DIAG;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;		ypos = polyleg(coord_pix0[*,0],tset_slits[0].coeff[*,slitno])
;               yposb = polyleg(coord_pix1[*,0],tset_slits[1].coeff[*,slitno]) 
;
;   print,coord_pix0[0,0],ypos(0),yposb(0)
;   print,coord_pix0[100,0],ypos(100),yposb(100)
;   plot,coord_pix0[*,0],coord_pix0[*,1],yrange=[min(coord_pix0[*,1])-10,max(coord_pix0[*,1])+10],xrange=[0,8192]
;   oplot,coord_pix0[*,0],ypos,color=cgColor('red')
;   oplot,coord_pix1[*,0],coord_pix1[*,1]
;   oplot,coord_pix1[*,0],yposb,color=cgColor('red')
;   xarr = indgen(8192)
;   result = interpol(coord_pix0[*,1],coord_pix0[*,0],xarr)
;   oplot,xarr,result,color=cgColor('green'),linestyle=2
;   oplot,xarr,interp_xx1,color=cgColor('blue'),linestyle=3
;   stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	    ; RED CHANNEL BRANCH
	    endif else if ((CHANNEL eq 'MODS1R') or (CHANNEL eq 'MODS2R')) then begin
	        red_lam=fltarr(101)
		coord_pix1=fltarr(101,2)
		coord_pix0=fltarr(101,2)
		if GRATNAME eq 'G670L' then begin
                  tset_slits[0].pix_beg[slitno] = 0.0
                  tset_slits[0].pix_end[slitno] = 8191
                  tset_slits[1].pix_beg[slitno] = 0.0
                  tset_slits[1].pix_end[slitno] = 8191
 		  for i=0,100 do begin
		      red_lam[i]=5500+40*i
                      coord_pix1[i,*]=red_map(CHANNEL,XMM,-(YMM-LENMM/2.+XMM*tan(rot_r)),red_lam[i])
		      coord_pix0[i,*]=red_map(CHANNEL,XMM,-(YMM+LENMM/2.+XMM*tan(rot_r)),red_lam[i])
		  endfor
		  ;catch slits going off the detector
		  if max(coord_pix1[*,0]) gt 8191 then begin
                        print,'The trace for Slit number ',slitno+1,' goes off the detector'
                        print,'The trace goes off the red end of the chip.'
                        print,'PROCEDE WITH CAUTION - suggestion: plot,coord_pix1[*,0]'
			print,'Official recommendation is to remove this slit from the MMS file.'
                        stop
                  endif
                  if max(coord_pix1[*,1]) gt 3088 then begin
                        print,'The trace for Slit number ',slitno+1,' goes off the detector'
                        print,'The top of the slit leaves the top of the chip.'
                        print,'PROCEDE WITHsuggestion: suggestion: plot,coord_pix1[*,1]'
                        print,'Official recommendation is to remove this slit from the MMS file.'
                        stop
                  endif
                  if min(coord_pix0[*,0]) lt 0 then begin
                        print,'The trace for Slit number ',slitno+1,' goes off the detector'
                        print,'The traces goes of the blue end of the chip.'
                        print,'PROCEDE WITH CAUTION - suggestion: plot,coord_pix0[*,0]'
                        print,'possible fix if only one pixel is off: coord_pix0[0,0] = 0'
                        print,'Official recommendation is to remove this slit from the MMS file.'
                        stop
                  endif
                  if min(coord_pix0[*,1]) lt 0 then begin
			print,'The trace for Slit number ',slitno+1,' goes off the detector'
                        print,'The trace goes off the bottom of the chip.'
                        print,'PROCEDE WITH CAUTION - suggestion: plot,coord_pix0[*,1]'
                        print,'Official recommendation is to remove this slit from the MMS file.'
			stop
		  endif
                  IF keyword_set(AUTOTUNE_SLITS) then begin
                  splog,'Autotune location of slit'+string(slitno + 1)
                  offset = slit_twiddle(coord_pix0,coord_pix1,image)
                        coord_pix1[*,1] += offset
                        coord_pix0[*,1] += offset
                  endif
                  xarr = indgen(8192)
                  interp_xx1 = interpol(coord_pix0[*,1],coord_pix0[*,0],xarr) + global_offset
                  interp_xx2 = interpol(coord_pix1[*,1],coord_pix1[*,0],xarr) + global_offset
                  IF keyword_set(INTERACTIVE_SLITS) then begin
                        offset_int = man_twiddle(xarr,interp_xx1,interp_xx2,image,slitno,0,8191)
                        Ok=0
                        while Ok ne 1 do begin
                                print,'Apply Global Offset? (y/n/abort)'
                                global = ''
                                read,': ',global
                                if global eq 'y' then begin
                                        INTERACTIVE_SLITS = 0
                                        AUTOTUNE_SLITS = 0
                                        global_offset = offset + offset_int
                                        ok=1
                                endif else if global eq 'n' then begin
                                        ok = 1
                                endif else if redo eq 'abort' then begin
                                        stop
                                endif else begin
                                        print,'I do not understand ',redo
                                endelse
                        endwhile
                        interp_xx1 += offset_int
                        interp_xx2 += offset_int
                  endif
                  tset_slits[0].xx1[*, slitno] = interp_xx1
                  tset_slits[1].xx2[*, slitno] = interp_xx2

                endif else if GRATNAME eq 'P700L' then begin
                  for i=0,100 do begin
                      red_lam[i]=5600+40*i
                      coord_pix1[i,*]=red_map_prism(CHANNEL,XMM,-(YMM-LENMM/2.+XMM*tan(rot_r)),red_lam[i])
                      coord_pix0[i,*]=red_map_prism(CHANNEL,XMM,-(YMM+LENMM/2.+XMM*tan(rot_r)),red_lam[i])
                  endfor
                  ;catch slits going off the detector
                  if max(coord_pix1[*,0]) gt 8191 then begin
                        print,'The trace for Slit number ',slitno,' goes off the detector'
                        stop
                  endif
                  if max(coord_pix1[*,1]) gt 3088 then begin
                        print,'The trace for Slit number ',slitno,' goes off the detector'
                        stop
                  endif
                  if min(coord_pix0[*,0]) lt 0 then begin
                        print,'The trace for Slit number ',slitno,' goes off the detector'
                        stop
                  endif
                  if min(coord_pix0[*,1]) lt 0 then begin
                        print,'The trace for Slit number ',slitno,' goes off the detector'
                        stop
                  endif
                  ; account for prism mode not using the full chip
                  coord_pix0[*,0] += 2048
                  coord_pix1[*,0] += 2048
                  ; account for y-shift of mask
                  IF keyword_set(AUTOTUNE_SLITS) then begin
                  splog,'Autotune location of slit'+string(slitno + 1)
                  offset = slit_twiddle(coord_pix0,coord_pix1,image)
                        coord_pix1[*,1] += offset
                        coord_pix0[*,1] += offset
                  endif
                  tset_slits[0].pix_beg[slitno] = coord_pix0[0,0]
                  tset_slits[0].pix_end[slitno] = coord_pix0[100,0]
                  tset_slits[1].pix_beg[slitno] = coord_pix0[0,0]
                  tset_slits[1].pix_end[slitno] = coord_pix0[100,0]

                  ;interpolate curve edges
                  xarr = indgen(8192)
                  interp_xx1 = interpol(coord_pix0[*,1],coord_pix0[*,0],xarr) + global_offset
                  interp_xx2 = interpol(coord_pix1[*,1],coord_pix1[*,0],xarr) + global_offset
                  IF keyword_set(INTERACTIVE_SLITS) then begin
                        offset_int = man_twiddle(xarr,interp_xx1,interp_xx2,image,slitno,coord_pix0[0,0],coord_pix0[100,0])
                        Ok=0
                        while Ok ne 1 do begin
                                print,'Apply Global Offset? (y/n/abort)'
                                global = ''
                                read,': ',global
                                if global eq 'y' then begin
                                        INTERACTIVE_SLITS = 0
                                        AUTOTUNE_SLITS = 0
                                        global_offset = offset + offset_int
                                        ok=1
                                endif else if global eq 'n' then begin
                                        ok = 1
                                endif else if redo eq 'abort' then begin
                                        stop
                                endif else begin
                                        print,'I do not understand ',redo
                                endelse
                        endwhile
                        interp_xx1 += offset_int
                        interp_xx2 += offset_int
                  endif
                  ;blank the begining and end of the trace (otherwise prism slits overlap)
                  for i = 0,coord_pix0[0,0]-1 do interp_xx1[i] = 0
                  for i = coord_pix0[100,0]+1,8191 do interp_xx1[i] = 0
                  for i = 0,coord_pix0[0,0]-1 do interp_xx2[i] = 0
                  for i = coord_pix0[100,0]+1,8191 do interp_xx2[i] = 0
                  tset_slits[0].xx1[*, slitno] = interp_xx1
                  tset_slits[1].xx2[*, slitno] = interp_xx2
                  
                  ;set the fit coeffs for slit tracing
                  coeffs1=poly_fit(coord_pix1[*,0],coord_pix1[*,1],1)
                  ;tset_slits[1].coeff[*, slitno] = coeffs1
                  coeffs0=poly_fit(coord_pix0[*,0],coord_pix0[*,1],1)
                  ;tset_slits[0].coeff[*, slitno] = coeffs0
                endif else begin
                      print,'Unkown Grating/Prism for the Red Channel'
                      stop
                endelse
	    endif

	    ;propogate the slit coordinates from the MMS file
	    radot = strpos(ALPHA,'.')
            ra1=strmid(ALPHA,radot-6,2)
            ra2=strmid(ALPHA,radot-4,2)
            ra3=strmid(ALPHA,radot-2)
            decdot = strpos(DELTA,'.')
            dec1=strmid(DELTA,decdot-7,3)
            dec2=strmid(DELTA,decdot-4,2)
            dec3=strmid(DELTA,decdot-2)
            apRA_deg = 15.D * (double(ra1) + double(ra2)/60.D + double(ra3)/3600.D)
            apDEC_deg = double(dec1) + double(dec2)/60.D + double(dec3)/3600.D
            apRA_sex = ra1 + ':' + ra2 + ':' + ra3
            apDEC_sex = dec1 + ':' + dec2 + ':' + dec3
            tset_slits[0].tarRA = tarRA_sex
            tset_slits[0].tarDEC = tarDEC_sex
            tset_slits[0].apRA[slitno] = apRA_sex
            tset_slits[0].apDEC[slitno] = apDEC_sex
            tset_slits[0].tarRA2 = tarRA_deg
            tset_slits[0].tarDEC2 = tarDEC_deg
            tset_slits[0].apRA2[slitno] = apRA_deg
            tset_slits[0].apDEC2[slitno] = apDEC_deg
	    tset_slits[0].width = WID
            tset_slits[0].length = LEN
            tset_slits[0].XMM_arr[slitno] = XMM
            tset_slits[0].YMM_arr[slitno] = YMM
            tset_slits[0].LENMM_arr[slitno] = LENMM
            tset_slits[0].WIDMM_arr[slitno] = WIDMM
	    slitno+=1
        endif
    endif
endfor


if global eq 'y' then begin
	INTERACTIVE_SLITS = 1
	AUTOTUNE_SLITS = 1
endif 

RETURN, tset_slits
END
