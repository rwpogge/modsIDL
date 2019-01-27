;=========================================================================
;=========================================================================
pro readnum,value
        satisfactory = 0
        input = ''
        repeat begin
                read,': ',input
                if not valid_num(input) then print,'I need a number' else begin
                        value = float(input)
                        satisfactory = 1
                endelse
        endrep until satisfactory
return
end

;=========================================================================
; center - Compute the center of the desired line in the slit
;		after centering on the strongest line in the
;		specified wavelength regions
;=========================================================================
pro center,file,wave,slit_num,wav_min,wav_max,result_x,$
	result_y,goodwav,wav1d,stack,stack2,xx1,xx2,$
	trim1,trim2

        ;load the 2D reticfied spectrum for centering
        xx1_ap = xx1[*,slit_num-1]+trim1
        xx2_ap = xx2[*,slit_num-1]-trim2

; working to trim for prism-----------------------
	trim_disp_dir = where(xx1_ap gt 15)
;        ap_im = x_ordrectify(file,xx1_ap[*],xx2_ap[*],/nocorrect)
;        ap_wave = x_ordrectify(wave,xx1_ap[*],xx2_ap[*],/nocorrect)
        ap_im = x_ordrectify(file[*,trim_disp_dir],xx1_ap[trim_disp_dir],xx2_ap[trim_disp_dir],/nocorrect)
        ap_wave = x_ordrectify(wave[*,trim_disp_dir],xx1_ap[trim_disp_dir],xx2_ap[trim_disp_dir],/nocorrect)
;---------------------------------------------
        ap_im = transpose(ap_im)
        ap_wave = transpose(ap_wave)
        sl_dims=size(ap_im,/dim)
        sl_nx=sl_dims[0]
        sl_ny=sl_dims[1]
        ;find the line peak to center on
        wav1d=ap_wave[*,sl_ny/2]
        goodwav = where(wav1d lt wav_max and wav1d gt wav_min)
        tmp = ap_im[goodwav,*]
        stack = total(tmp,2)
        result=where(stack eq max(stack))
        result_x = goodwav[result]
        ;find which part of the slit is at the median... mininimum pulled out the slit edge
        tmp = ap_im[result_x-4:result_x+5,*]
        stack2 = total(tmp,1)
	result_y=where(stack2 eq max(stack2))
return
end

;=========================================================================
; extractionplot - plot the variosu cuts and extractions
;			for the selected aperture
;=========================================================================
pro extractionplot,bluew,redw,blueg,redg, $
	bluest,redst,bluex,redx,bluest2,redst2,bluey,redy, $
	bluesp,redsp,lower,upper

common exclude1, blank_lowb, blank_lowr, blank_hib, blank_hir

;Plot up the centering and extractions
;bottom panel
xb0 = bluew[blueg[0]]-20
xr0 = redw[redg[0]]-20
xb1 = max(bluew[blueg])+20
xr1 = max(redw[redg])+20
plot,bluew[blueg],bluest,Position=[0.09,0.05,0.45,0.20],XRANGE=[xb0,xb1],xstyle=1
line_x = [bluew[bluex],bluew[bluex]]
line_x1 = [bluew[bluex-9],bluew[bluex-9]]
line_x2 = [bluew[bluex+10],bluew[bluex+10]]
minval = !Y.CRANGE[0]
maxval = !Y.CRANGE[1]
line_y = [minval,maxval]
oplot,line_x,line_y,Color=cgColor('dodger blue')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')
plot,redw[redg],redst,Position=[0.59,0.05,0.95,0.20],$
        /noerase, XRANGE=[xr0,xr1],xstyle=1
line_x = [redw[redx],redw[redx]]
line_x1 = [redw[redx-9],redw[redx-9]]
line_x2 = [redw[redx+10],redw[redx+10]]
minval = !Y.CRANGE[0]
maxval = !Y.CRANGE[1]
line_y = [minval,maxval]
oplot,line_x,line_y,Color=cgColor('dodger blue')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')

;middle pannel
stacksz = size(bluest2,/dim)
plot,indgen(stacksz),bluest2,XRANGE=[0,stacksz],$
	position=[0.09,0.27,0.45,0.60],$
        /noerase,yrange=[min(stacksz),max(stacksz)],xstyle=1
line_x = [bluey,bluey]
minval = !Y.CRANGE[0]
maxval = !Y.CRANGE[1]
line_y = [minval,maxval]
line_x1 = [bluey+upper/0.12,bluey+upper/0.12]
line_x2 = [bluey-lower/0.12,bluey-lower/0.12]
niter = size(blank_lowb,/dim)
niter = niter[0]-1
for i=0,niter do polyfill,[blank_lowb[i],blank_lowb[i],blank_hib[i],blank_hib[i]], $
	[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data
plot,indgen(stacksz),bluest2,XRANGE=[0,stacksz],$
        position=[0.09,0.27,0.45,0.60],$
        /noerase,yrange=[min(stacksz),max(stacksz)],xstyle=1,xtitle='Pixels'
oplot,line_x,line_y,Color=cgColor('dodger blue')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')
axis,xaxis=1,XRANGE=[0,stacksz*0.12],xtitle='ARCSEC'
stacksz = size(redst2,/dim)
plot,indgen(stacksz),redst2,XRANGE=[0,stacksz],$
	position=[0.59,0.27,0.95,0.60],$
        /noerase,yrange=[min(stacksz),max(stacksz)],xstyle=1
line_x = [redy,redy]
minval = !Y.CRANGE[0]
maxval = !Y.CRANGE[1]
line_y = [minval,maxval]
line_x1 = [redy+upper/0.123,redy+upper/0.123]
line_x2 = [redy-lower/0.123,redy-lower/0.123]
for i=0,niter do polyfill,[blank_lowr[i],blank_lowr[i],blank_hir[i],blank_hir[i]],$
	[minval,maxval,maxval,minval],color=cgColor('firebrick')
oplot,line_x,line_y,Color=cgColor('dodger blue')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')
plot,indgen(stacksz),redst2,XRANGE=[0,stacksz],$
        position=[0.59,0.27,0.95,0.60],$
        /noerase,yrange=[min(stacksz),max(stacksz)],xstyle=1,xtitle='Pixels'
axis,xaxis=1,XRANGE=[0,stacksz*0.123],xtitle='ARCSEC'

;top panel
plot,bluew,bluesp,position=[0.09,0.72,0.45,0.99],$
        /noerase,XRANGE=[xb0-100,xb1+50];,YRANGE=[200,900]
plot,redw,redsp,position=[0.59,0.72,0.95,0.99],$
	/noerase,XRANGE=[xr0-100,xr1+100];,YRANGE=[200,900]

return
end

;=========================================================================
; mods_extract1d - extract and calibrate 1D dual channel MODS spectra
;			from 2D spectra
;=========================================================================
pro mods_extract_prism,red_scifile,blue_scifile, red_cal,blue_cal, z = z, $
	apertures=apertures,outname=outname,wave_blue=wave_blue,wave_red=wave_red, $
	centerLine=centerLine, trim_rt = trim_rt, trim_rb = trim_rb, trim_bb = trim_bb, trim_bt = trim_bt,$
	noPlotAp=noPlotAp, red_slits = red_slits, blue_slits = blue_slits,$
	man_scale = man_scale

;deal with parameters that were not entered
if n_elements(outname) EQ '' then begin
	redstub = STRSPLIT(red_scifile,'_',/extract)
	redout = 'Science/1Dsub_' + redstub[1] + '.fits'
	redoutcal = 'Science/1Dsub_' + redstub[1] + '_cal.fits'
	bluestub = STRSPLIT(blue_scifile,'_',/extract)
	blueout = 'Science/1Dsub_' + bluestub[1] + '.fits'
	blueoutcal = 'Science/1Dsub_' + bluestub[1] + '_cal.fits'
endif else begin
	redout = 'Science/1Dsub_red_' + outname + '.fits'
	redoutcal = 'Science/1Dsub_red_' + outname + '_cal.fits'
	blueout = 'Science/1Dsub_blue_' + outname + '.fits'
	blueoutcal = 'Science/1Dsub_blue_' + outname + '_cal.fits'
endelse


if n_elements(z) EQ '' then begin
   print,''
   print, 'You gave me no redshift.  I will assume (1+z)=1'
   print,''
   wait,0.5
   z=0
endif
z = z+1 ;later code uses z+1 so increment up

if n_elements(trim_rt) EQ 0 then trim_rt = 1
if n_elements(trim_bt) EQ 0 then trim_bt = 1

if n_elements(trim_rb) EQ 0 then trim_rb = 1
if n_elements(trim_bb) EQ 0 then trim_bb = 1

if n_elements(centerLine) EQ '' then begin
	rlow_cent = 5650*z
	rhi_cent = 7800*z
	blow_cent = 3890*z
	bhi_cent =  5170*z
endif else begin
	rlow_cent = centerLine[2]*z
        rhi_cent = centerLine[3]*z
        blow_cent = centerLine[0]*z
        bhi_cent =  centerLine[1]*z
endelse

common exclude1, blank_lowb, blank_lowr, blank_hib, blank_hir
blank_lowb = [0]
blank_lowr = [0]
blank_hib = [0]
blank_hir = [0]

; Load the images
; Science image
verify_red = file_test(red_scifile)
verify_blue = file_test(blue_scifile)
if verify_red then begin
	red_sciimg = mrdfits(red_scifile,0,red_hdr)
        red_skyimg = mrdfits(red_scifile,2,red_hdrsky)
        red_varimg = mrdfits(red_scifile,1,rheadvar)
        red_varimg = red_skyimg * red_varimg + red_sciimg * red_varimg
 endif else begin
                print,'------------------------------------------------------'
                print,'You need a valid RED science image'
                print,red_scifile,' does not exist. Check your Science folder'
                stop
 endelse
if verify_blue then begin
        blue_sciimg = mrdfits(blue_scifile,0,blue_hdr)
        blue_skyimg = mrdfits(blue_scifile,2,blue_hdrsky)
        blue_varimg = mrdfits(blue_scifile,1,bheadvar)
        blue_varimg = blue_skyimg * blue_varimg + blue_sciimg * blue_varimg
 endif else begin
                print,'------------------------------------------------------'
                print,'You need a valid BLUE science image'
                print,blue_scifile,' does not exist. Check your Science folder'
                stop
 endelse

; Get Mask Info
red_instrument = strcompress(sxpar(red_hdr[*,0], 'INSTRUME'), /rem)
blue_instrument = strcompress(sxpar(blue_hdr[*,0], 'INSTRUME'), /rem)
if red_instrument eq 'MODS1R' then red_channel_code = 'm1r' $
else if red_instrument eq 'MODS2R' then red_channel_code = 'm2r' $
else begin
        print,red_instrument,' is unknown'
        stop
endelse
if blue_instrument eq 'MODS1B' then blue_channel_code = 'm1b' $
else if blue_instrument eq 'MODS2B' then blue_channel_code = 'm2b' $
else begin
        print,blue_instrument,' is unknown'
        stop
endelse
maskr = strcompress(sxpar(red_hdr[*,0], 'MASKNAME'), /rem)
maskb = strcompress(sxpar(blue_hdr[*,0], 'MASKNAME'), /rem)
if maskr ne maskb then begin
        print,'WARNING: These files are from different masks.'
        stop
endif
if maskr eq 'LS60x5' then apertures=[1]
if strmid(maskr,0,2) eq 'LS' then mask = 'LS' else mask = maskr

red_hdrcal = red_hdr
blue_hdrcal = blue_hdr

sxaddpar,red_hdrcal,'CUNIT1','Angstrom'
sxaddpar,red_hdrcal,'CTYPE1','Linear'
sxaddpar,red_hdrcal,'CRPIX1',1
sxaddpar,red_hdrcal,'CRVAL1',5500
sxaddpar,red_hdrcal,'CDELT1',0.5

sxaddpar,blue_hdrcal,'CUNIT1','Angstrom'
sxaddpar,blue_hdrcal,'CTYPE1','Linear'
sxaddpar,blue_hdrcal,'CRPIX1',1
sxaddpar,blue_hdrcal,'CRVAL1',3200
sxaddpar,blue_hdrcal,'CDELT1',0.5

; Calbration image
calR2 = mrdfits(red_cal,0,chd)
calB2 = mrdfits(blue_cal,0,chd)
;;;;;;;;;;;;;;;;;;;;;;;;;;
; MAKE THESE BE READ FROM THE CAL BLOCK HEADER
new_red_wave = fltarr(670)
;new_red_wave = fltarr(9001)
new_blue_wave = fltarr(151)
;new_blue_wave = fltarr(5201)
new_red_wave[0] = 5500
new_blue_wave[0] = 3200
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;for i=1,9001-1 do new_red_wave[i] = new_red_wave[i-1]+0.5
for i=1,670-1 do new_red_wave[i] = new_red_wave[i-1]+0.5
;for i=1,5201-1 do new_blue_wave[i] = new_blue_wave[i-1]+0.5
for i=1,151-1 do new_blue_wave[i] = new_blue_wave[i-1]+0.5
sset=bspline_iterfit(new_blue_wave,calB2,bkspace=6,nord=3)
calB = bspline_valu(new_blue_wave,sset)
sset=bspline_iterfit(new_red_wave,calR2,bkspace=6,nord=3)
calR = bspline_valu(new_red_wave,sset)
calR_coef = bspline_iterfit(new_red_wave,calR,bkspace=0.65,nord=3)
calB_coef = bspline_iterfit(new_blue_wave,calB,bkspace=0.65,nord=3)
extfile = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/mods_extiction.dat'
readcol,extfile,ext_wave,ext_val,/silent
ext_coef = bspline_iterfit(ext_wave,ext_val,bkspace=3,nord=3)
extR = bspline_valu(new_red_wave,ext_coef)
extB = bspline_valu(new_blue_wave,ext_coef)

; Wavelength image
if n_elements(wave_blue) EQ '' then begin
        wave_blue = 'wave-'+blue_channel_code+'_HgXeAr_'+mask+'.fits'
endif
if n_elements(wave_red) EQ '' then begin
        wave_red = 'wave-'+red_channel_code+'_NeXeAr_'+mask+'.fits'
endif
verify_red = file_test(wave_red)
verify_blue = file_test(wave_blue)
if verify_red then red_waveimg = mrdfits(wave_red,0,red_whdr) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid RED wave-image'
                print,wave_red,' does not exist.'
                print,"Use: wave_red='FILENAME'"
                stop
        endelse
if verify_blue then blue_waveimg = mrdfits(wave_blue,0,blue_whdr) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid BLUE wave-image'
                print,wave_blue,' does not exist.'
                print,"Use: wave_blue='FILENAME'"
                stop
        endelse

;open a window for plotting
window,1,retain=2;,XSize=1000,YSize=1000
if (not keyword_set(noPlotAp)) then window,2,retain=2

dims=size(red_sciimg,/dim)
nx=dims[0]
ny=dims[1]
mms = strcompress(sxpar(red_hdr[*,0], 'MASKINFO')+'.mms', /rem)
print,'mms file is ',mms

;create the slit structures
  if n_elements(red_slits) eq '' then red_slits = 'slits-'+red_channel_code+'_ill_'+mask+'.fits'
  if n_elements(blue_slits) eq '' then blue_slits = 'slits-'+blue_channel_code+'_ill_'+mask+'.fits'

  verify_red = file_test(red_slits)
  verify_blue = file_test(blue_slits)
  if verify_red then red_slits = mrdfits(red_slits,1,r_head) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid RED slit-image'
                print,red_slits,' does not exist.'
                print,"Use: red_slits='FILENAME'"
                stop
        endelse
   if verify_blue then blue_slits = mrdfits(blue_slits,1,b_head) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid BLUE slit-image'
                print,blue_slits,' does not exist.'
                print,"Use: blue_slits='FILENAME'"
                stop
        endelse

   dim_t = red_slits[0].dims
   print,'red_slits dimensions = ', dim_t
   nslit=size(red_slits.xx1,/dim)
   nslit=nslit[1]
   nx_t = dims[0]
   ny_t = dims[1]

   red_xx1 = red_slits[0].xx1
   red_xx2 = red_slits[1].xx2

   red_objtrace = red_xx1*0
   red_objspec = fltarr(ny_t)+1
   red_objwave = fltarr(ny_t)+1
   red_varspec = fltarr(ny_t)+1
   red_skyspec = fltarr(ny_t)+1

   blue_xx1 = blue_slits[0].xx1
   blue_xx2 = blue_slits[1].xx2

   blue_objtrace = blue_xx1*0
   blue_objspec = fltarr(ny_t)+1
   blue_objwave = fltarr(ny_t)+1
   blue_varspec = fltarr(ny_t)+1
   blue_skyspec = fltarr(ny_t)+1

   blue_objspec_cal = calB2
   red_objspec_cal = calR2
   blue_objspec_calb = calB2
   red_objspec_calb = calR2
   blue_varspec_cal = fltarr(5201)
   red_varspec_cal = fltarr(9001)
   blue_skyspec_cal = fltarr(5201)
   red_skyspec_cal = fltarr(9001)

;extraction struct
red_extract = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, $
	ra_deg: double(red_slits[0].tarra2), dec_deg: double(red_slits[0].tardec2), $
	ra_sex: red_slits[0].tarra, dec_sex: red_slits[0].tardec}
blue_extract = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, $
	ra_deg: double(blue_slits[0].tarra2), dec_deg: double(blue_slits[0].tardec2), $
	ra_sex: blue_slits[0].tarra, dec_sex: blue_slits[0].tardec}

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                      fit each slit on its own
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(apertures) EQ '' then apertures = indgen(nslit)+1
ndo = size(apertures,/dim)
ndo=ndo[0]
object=1

for ii=0,ndo-1 do begin
	LeObject: print,'Initializing'
	ap=apertures[ii]
	;currently, returning the center of the slit for ra and dec 
	red_extract_obj = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, $
		ra_deg: red_slits[0].apra2[ap-1], dec_deg: red_slits[0].apdec2[ap-1], $
		ra_sex: red_slits[0].apra[ap-1], dec_sex: red_slits[0].apdec[ap-1]}
	blue_extract_obj = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, $
		ra_deg: blue_slits[0].apra2[ap-1], dec_deg: blue_slits[0].apdec2[ap-1], $
		ra_sex: blue_slits[0].apra[ap-1], dec_sex: blue_slits[0].apdec[ap-1]}
	blank_lowb = [0]
	blank_lowr = [0]
	blank_hib = [0]
	blank_hir = [0]
	print,'APERTURE: ',ap
	print,''
	n_blank = 1
	center,red_sciimg,red_waveimg,ap,rlow_cent,rhi_cent,$
		red_result_x,red_result_y,red_good,red_wav1d_trim,$
		red_stack,red_stack2,red_xx1,red_xx2,trim_rb,trim_rt
        center,blue_sciimg,blue_waveimg,ap,blow_cent,bhi_cent,$
		blue_result_x,blue_result_y,blue_good,blue_wav1d_trim,$
		blue_stack,blue_stack2,blue_xx1,blue_xx2,trim_bb,trim_bt
        lower = 2.
        upper = 2.

	LeSigh: print,'Extracting Spectrum'
        red_spec1d = fltarr(ny)
        blue_spec1d = fltarr(ny)
        red_var1d = fltarr(ny)
        blue_var1d = fltarr(ny)
	red_sky1d = fltarr(ny)
	blue_sky1d = fltarr(ny)
        red_spec1db = fltarr(ny)
        blue_spec1db = fltarr(ny)
        red_wav1d = fltarr(ny)
        blue_wav1d = fltarr(ny)

	modsadcalc,blue_scifile,blue_wav1d_trim[blue_result_x],blue_wav1d_trim,bxdelt,bydelt
	modsadcalc,red_scifile,red_wav1d_trim[red_result_x],red_wav1d_trim,rxdelt,rydelt

	for iii=0,ny-1 do red_objtrace[iii,ap-1] = red_xx1[iii,ap-1] + trim_rb + red_result_y
        for iii=0,ny-1 do blue_objtrace[iii,ap-1] = blue_xx1[iii,ap-1] + trim_bb + blue_result_y
	
        ; added in for method demonstration
        blue_orig = blue_objtrace[*,ap-1]

	; trim for prism
        red_xx1_ap = red_xx1[*,ap-1]+trim_rb
        red_xx2_ap = red_xx2[*,ap-1]-trim_rt
        red_trim = where(red_xx1_ap gt 15)
	red_x_start = min(red_trim)
	red_x_end = max(red_trim)	
	;

        red_objtrace[red_trim,ap-1] = red_objtrace[red_trim,ap-1] - rydelt/0.123
        red_wav1d[red_trim] = red_wav1d_trim
       for iii=red_x_start,red_x_end do begin
;	for iii=0,ny-1 do begin
		bot = red_objtrace[iii,ap-1]-lower/0.123
                top = red_objtrace[iii,ap-1]+upper/0.123
                if bot lt 0 then bot = 0
		botint = total(bot,/integer)
		topint = total(top,/integer)
		red_spec1d[iii] = total(red_sciimg[bot:top,iii]) $
			+(1-(bot-botint))*total(red_sciimg[bot,iii]) $
			+(top-topint)*total(red_sciimg[top+1,iii])
		red_wav1d[iii] = median(red_waveimg[bot:top,iii])
		red_var1d[iii] = sqrt(total(red_varimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(red_varimg[bot,iii]) $
                        +(top-topint)*total(red_varimg[top+1,iii]))
; the sky image is non-existent...
		red_sky1d[iii] = sqrt(total(red_skyimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(red_skyimg[bot,iii]) $
                        +(top-topint)*total(red_skyimg[top+1,iii]))
; ......
;		red_spec1db[iii] = total(red_sciimg2[bot:top,iii]) $
;                        +(1-(bot-botint))*total(red_sciimg2[bot,iii]) $
;                        +(top-topint)*total(red_sciimg2[top+1,iii])
	endfor
        ; trim for prism
        blue_xx1_ap = blue_xx1[*,ap-1]+trim_bb
        blue_xx2_ap = blue_xx2[*,ap-1]-trim_bt
        blue_trim = where(blue_xx1_ap gt 15)
        blue_x_start = min(blue_trim)
        blue_x_end = max(blue_trim)
	;
        blue_objtrace[blue_trim,ap-1] = blue_objtrace[blue_trim,ap-1] - bydelt/0.12
        blue_wav1d[blue_trim] = blue_wav1d_trim

        for iii=blue_x_start,blue_x_end do begin
                bot = blue_objtrace[iii,ap-1]-lower/0.12
                top = blue_objtrace[iii,ap-1]+upper/0.12
		if bot lt 0 then bot = 0
                botint = total(bot,/integer)
                topint = total(top,/integer)
                blue_spec1d[iii] = total(blue_sciimg[bot:top,iii])$
                        +(1-(bot-botint))*total(blue_sciimg[bot,iii]) $
                        +(top-topint)*total(blue_sciimg[top+1,iii])
		blue_wav1d[iii] = median(blue_waveimg[bot:top,iii])
		blue_var1d[iii] = sqrt(total(blue_varimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(blue_varimg[bot,iii]) $
                        +(top-topint)*total(blue_varimg[top+1,iii]))
; sky img is broken....
		blue_sky1d[iii] = sqrt(total(blue_skyimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(blue_skyimg[bot,iii]) $
                        +(top-topint)*total(blue_skyimg[top+1,iii]))
; ..............
;		blue_spec1db[iii] = total(blue_sciimg2[bot:top,iii])$
;                        +(1-(bot-botint))*total(blue_sciimg2[bot,iii]) $
;                        +(top-topint)*total(blue_sciimg2[top+1,iii])
	endfor

	; plot slit and traces
	if (not keyword_set(noPlotAp)) then begin
		wset,2
	        ncolors = !D.Table_Size
	        device,DECOMPOSED=0
	        LoadCT,39
	        bmin = min(blue_xx1[blue_trim,ap-1])+trim_bb
	        bmax = max(blue_xx2[blue_trim,ap-1])-trim_bt
	        rmin = min(red_xx1[red_trim,ap-1])+trim_rb
	        rmax = max(red_xx2[red_trim,ap-1])-trim_rt


	bscale_min = -5.
	bscale_max = 5*mean(blue_sciimg[bmin:bmax,blue_trim])
	rscale_min = -5.
	rscale_max = 15*mean(red_sciimg[rmin:rmax,red_trim])

        if keyword_set(man_scale) then begin
	                Print,'Blue Min Scale for 2D plot (',bscale_min,')'
	                readnum,bscale_min
	                Print,'Blue Max Scale for 2D plot (',bscale_max,')'
	                readnum,bscale_max
	                Print,'Red Min Scale for 2D plot (',rscale_min,')'
	                readnum,rscale_min
	                Print,'Red Max Scale for 2D plot (',rscale_max,')'
	                readnum,rscale_max
        endif

	;extract for 2D plot
	        slitimb = transpose(BytScl(blue_sciimg[bmin:bmax,blue_trim],$
			min=bscale_min,max=bscale_max,/nan))
	        slitimr = transpose(BytScl(red_sciimg[rmin:rmax,red_trim],$
			min=rscale_min,max=rscale_max,/nan))
	        TVImage,slitimb,Position=[0.02,0.08,0.99,0.49],/erase
	        TVImage,slitimr,Position=[0.02,0.52,0.99,0.99]
		pathb = BLUE_OBJTRACE[blue_trim,ap-1] - bmin
	        pathr = RED_OBJTRACE[red_trim,ap-1] - rmin
	        blue_orig = blue_orig - bmin

		xpix_r = indgen(n_elements(pathr))
		xpix_b = indgen(n_elements(pathb))

		plot,xpix_b,pathb,/noerase,Position=[0.02,0.08,0.99,0.49],$
			xstyle=1,xrange=[0,n_elements(pathb)-1],ystyle=1,yrange=[0,bmax-bmin],$
			linestyle=1,color=cgColor('Dodger Blue')
	        plot,xpix_b,pathb-lower/0.12,/noerase,Position=[0.02,0.08,0.99,0.49],$
			xstyle=1,xrange=[0,n_elements(pathb)-1],ystyle=1,yrange=[0,bmax-bmin]
	        plot,xpix_b,pathb+upper/0.12,/noerase,Position=[0.02,0.08,0.99,0.49],$
			xstyle=1,xrange=[0,n_elements(pathb)-1],ystyle=1,yrange=[0,bmax-bmin]
	        plot,xpix_r,pathr,/noerase,Position=[0.02,0.52,0.99,0.99],$
	                xstyle=1,xrange=[0,n_elements(pathr)-1],ystyle=1,yrange=[0,rmax-rmin],$
			linestyle=1,color=cgColor('Dodger Blue')
	        plot,xpix_r,pathr-lower/0.123,/noerase,Position=[0.02,0.52,0.99,0.99],$
	                xstyle=1,xrange=[0,n_elements(pathr)-1],ystyle=1,yrange=[0,rmax-rmin]
	        plot,xpix_r,pathr+upper/0.123,/noerase,Position=[0.02,0.52,0.99,0.99],$
	                xstyle=1,xrange=[0,n_elements(pathr)-1],ystyle=1,yrange=[0,rmax-rmin]
	endif
	;Plot up the centering and extractions
	wset,1

	extractionplot,blue_wav1d[blue_trim],red_wav1d[red_trim],blue_good,red_good,$
	       blue_stack,red_stack,blue_result_x,red_result_x,$
	       blue_stack2,red_stack2,blue_result_y,red_result_y,$
	       blue_spec1d[blue_trim],red_spec1d[red_trim],lower,upper

	Ok=0

	while Ok ne 1 do begin
	        print,'Centered on lines:'
	        print,'Blue  - ',blue_wav1d[blue_result_x]
	        print,'Red   - ', red_wav1d[red_result_x]
	        print,'Peaks found at:' 
	        print,'Blue  - ',blue_result_y
	        print,'Red   - ',red_result_y
	        print,'Currently using: lower = ',lower,'  upper = ',upper
		print,''
		print,' NOTE: LOWER AND UPPER ARE IN ARCSECONDS'
		print,'      10 pixels is ~1.2"'
	        print,'Recenter extraction? (y/n/bye)'
	        redo = ''
	        read,': ',redo
	        if redo eq 'y' then begin
	                Print,'New center blue:'
	                read,': ',blue_result_y
	                Print,'New center red:'
	                read,': ',red_result_y
	                GOTO,LeSigh
	        endif else if redo eq 'n' then begin
	                print,'ReSize extraction? (y/n/bye)'
	                redo = ''
	                read,': ',redo
	                if redo eq 'y' then begin
	                        Print,'New Lower:'
	                        read,': ',lower
	                        Print,'New Upper:'
	                        read,': ',upper
	                        GOTO,LeSigh
	                endif else if redo eq 'n' then begin
	                        print,'Huzzah!'

                                blue_extract_obj.object = object
	                        blue_extract_obj.slit = ap
	                        blue_extract_obj.center = blue_result_y
	                        blue_extract_obj.lower = lower
	                        blue_extract_obj.upper = upper

                                pos_shift = blue_result_y * 0.12
                                if pos_shift LE blue_slits[0].length/2.D then $
                                        hyp = blue_slits[0].length/2.D - pos_shift else $
                                        hyp = pos_shift - blue_slits[0].length/2.D
                                delRA = hyp * sin(blue_slits[0].posang * !dpi/180) * $
                                        cos(blue_extract_obj.dec_deg * !dpi/180) / 3600.D
                                delDEC = hyp * cos(blue_slits[0].posang * !dpi/180) / 3600.D
                                if pos_shift LE blue_slits[0].length/2.D then begin
                                        objra_deg = blue_extract_obj.ra_deg - delRA
                                        objdec_deg = blue_extract_obj.dec_deg - delDEC
                                endif else begin
                                        objra_deg = blue_extract_obj.ra_deg + delRA
                                        objdec_deg = blue_extract_obj.dec_deg + delDEC
                                endelse
				blue_extract_obj.ra_deg = objra_deg
				blue_extract_obj.dec_deg = objdec_deg
				radec,objra_deg,objdec_deg,ihr,imin,xsec,ideg,imn,xsc
				blue_extract_obj.ra_sex = strtrim(string(ihr),2)+':' $
					+strtrim(string(imin),2)+':'+strtrim(string(xsec),2)
				blue_extract_obj.dec_sex = strtrim(string(ideg),2)+':' $
					+strtrim(string(imn),2)+':'+strtrim(string(xsc),2)
	
        	                red_extract_obj.object = object
				red_extract_obj.slit = ap
        	                red_extract_obj.center = red_result_y
        	                red_extract_obj.lower = lower
        	                red_extract_obj.upper = upper

                                pos_shift = red_result_y * 0.123
				if pos_shift LE red_slits[0].length/2.D then $
					hyp = red_slits[0].length/2.D - pos_shift else $
					hyp = pos_shift - red_slits[0].length/2.D
				delRA = hyp * sin(red_slits[0].posang * !dpi/180) * $
					cos(red_extract_obj.dec_deg * !dpi/180) / 3600.D
                                delDEC = hyp * cos(red_slits[0].posang * !dpi/180) / 3600.D
				;red_extract_obj.ra_sex = 'tmp' 
				;red_extract_obj.dec_sex = 'tmp'
				if pos_shift LE red_slits[0].length/2.D then begin
					objra_deg = red_extract_obj.ra_deg - delRA
					objdec_deg = red_extract_obj.dec_deg - delDEC
				endif else begin
					objra_deg = red_extract_obj.ra_deg + delRA
					objdec_deg = red_extract_obj.dec_deg + delDEC
				endelse
				red_extract_obj.ra_deg = objra_deg
				red_extract_obj.dec_deg = objdec_deg
                                radec,objra_deg,objdec_deg,ihr,imin,xsec,ideg,imn,xsc
                                red_extract_obj.ra_sex = strtrim(string(ihr),2)+':' $
					+strtrim(string(imin),2)+':'+strtrim(string(xsec),2)
                                red_extract_obj.dec_sex = strtrim(string(ideg),2)+':' $
					+strtrim(string(imn),2)+':'+strtrim(string(xsc),2)
	                        blue_extract = [[blue_extract],[blue_extract_obj]]
	                        red_extract = [[red_extract],[red_extract_obj]]
	                        object ++
	                        Ok = 1
        	        endif else if redo eq 'bye' then begin
	                        stop
	                        endif else print,'I do not understand ',redo
                endif else if redo eq 'bye' then begin
                        stop
	        endif else begin
	                print,'I do not understand ',redo
	        endelse
	endwhile

        ;read in the exptime and airmass from the headers
        time = strcompress(sxpar(red_hdr[*, 0], 'EXPTIME'), /rem)
        airmassR = strcompress(sxpar(red_hdr[*, 0], 'AIRMASS'), /rem)
        airmassB = strcompress(sxpar(blue_hdr[*, 0], 'AIRMASS'), /rem)

        ;Rebin the spectrum to the region covered by the standards in actuallity
;        mods_flexurecorr,blue_wav1d,blue_spec1d,blue_sky1d,blue_var1d,'blue',z=z,/plot_prog
;        mods_flexurecorr_p,red_wav1d,red_spec1d,red_sky1d,red_var1d,'red',z=z,/plot_prog

        ;Rebin the spectrum to the region covered by the standards in actuallity
        x_specrebin,blue_wav1d[blue_trim],blue_spec1d[blue_trim],new_blue_wave,corflux_blue
        x_specrebin,red_wav1d[red_trim],red_spec1d[red_trim],new_red_wave,corflux_red

;	x_specrebin,blue_wav1d,blue_spec1db,new_blue_wave,corflux_blueb
;	x_specrebin,red_wav1d,red_spec1db,new_red_wave,corflux_redb

	x_specrebin,blue_wav1d[blue_trim],blue_var1d[blue_trim],new_blue_wave,corvar_blue
	x_specrebin,red_wav1d[red_trim],red_var1d[red_trim],new_red_wave,corvar_red

	x_specrebin,blue_wav1d[blue_trim],blue_sky1d[blue_trim],new_blue_wave,corsky_blue
        x_specrebin,red_wav1d[red_trim],red_sky1d[red_trim],new_red_wave,corsky_red

        ;The extinction correction is given by the factor
        extBcor = 10. ^ (0.4 * airmassB * extB[*])
        extRcor = 10. ^ (0.4 * airmassR * extR[*])

        ;Flux the spectrum including the LBT extinction curve
        calfluxB = (corflux_blue*extBcor)/(time*10^(calB/2.5))
        calfluxR = (corflux_red*extRcor)/(time*10^(calR/2.5))

;	calfluxBb = (corflux_blueb*extBcor)/(time*10^(calB/2.5))
;	calfluxRb = (corflux_redb*extRcor)/(time*10^(calR/2.5))
	
	calvarB = (corvar_blue*extBcor)/(time*10^(calB/2.5))
	calvarR = (corvar_red*extRcor)/(time*10^(calR/2.5))

	calskyB = (corsky_blue*extBcor)/(time*10^(calB/2.5))
	calskyR = (corsky_red*extRcor)/(time*10^(calR/2.5))

	print,'Saving aperture: ',ap
        ;save this aperture
        blue_objspec = [[blue_objspec],[blue_spec1d]]
        blue_objwave = [[blue_objwave],[blue_wav1d]]
	blue_varspec = [[blue_varspec],[blue_var1d]]
	blue_skyspec = [[blue_skyspec],[blue_sky1d]]

        red_objspec = [[red_objspec],[red_spec1d]]
        red_objwave = [[red_objwave],[red_wav1d]]
	red_varspec = [[red_varspec],[red_var1d]]
	red_skyspec = [[red_skyspec],[red_sky1d]]

;        blue_objspec_cal = [[blue_objspec_cal],[calfluxB]]
;        red_objspec_cal = [[red_objspec_cal],[calfluxR]]
;	blue_objspec_calb = [[blue_objspec_calb],[calfluxBb]]
;	red_objspec_calb = [[red_objspec_calb],[calfluxRb]]
;	blue_varspec_cal = [[blue_varspec_cal],[calvarB]]
;	red_varspec_cal = [[red_varspec_cal],[calvarR]]
;	blue_skyspec_cal = [[blue_skyspec_cal],[calskyB]]
;	red_skyspec_cal = [[red_skyspec_cal],[calskyR]]

help,red_objspec
        ; Save the output
        mwrfits,red_objspec,redout,red_hdr,/create
	mwrfits,red_varspec,redout
	mwrfits,red_skyspec,redout
	mwrfits,red_objwave,redout
        mwrfits,red_extract,redout

;       mwrfits,red_objspec_cal,redoutcal,red_hdrcal,/create
;	mwrfits,red_varspec_cal,redoutcal
;	mwrfits,red_skyspec_cal,redoutcal
;	mwrfits,red_extract,redoutcal

        mwrfits,blue_objspec,blueout,blue_hdr,/create
	mwrfits,blue_varspec,blueout
	mwrfits,blue_skyspec,blueout
	mwrfits,blue_objwave,blueout
        mwrfits,blue_extract,blueout

;       mwrfits,blue_objspec_cal,blueoutcal,blue_hdrcal,/create
;	mwrfits,blue_varspec_cal,blueoutcal
;	mwrfits,blue_skyspec_cal,blueoutcal
;	mwrfits,blue_extract,blueoutcal

	; Add more than region in a slit
	Print,'Extract an additional Region? (y/n/bye)'
	read,': ',redo
	if redo eq 'y' then begin
		GOTO,LeObject
	endif else if redo eq 'bye' then begin
              stop
	endif else if redo eq 'n' then begin
		print,'Moving on'
	endif

endfor

print,'Extracted ',object-1,' objects from ',ap,' apertures'
end
