;+
; NAME: 
;    mods_extract1d
;
; PURPOSE:
;    Extract 1D row stacked spectra from MODS MOS observations
;
; CALLING SEPUENCE:
;    mods_skyfit2d, red_image,blue_image
;
; INPUTS:
;    red_image  - red channel science image
;    blue_image - blue channel science image
;
; OPTIONAL INUTS:
;    red_cal   - name of the RED standard file you wish to use 
;    blue_cal  - name of the BLUE standard file you wish to use
;    force     - forace a given flexure correction
;    skyslit   - Indicate which slit is a sky slit.  This will cause that slit 
;                to be fit first. The sky slit will be used to fill in the emission
;                under the brightest emission lines.  A full model of the sky-slit 
;                can also be applied to slits with no clean samples of sky.
;    z         - redshift of the object.  More important for emission line objects
;                when you wish to center on specific lines.
;    apertures - apertures to work on.  If not given
;                all apertures are processed. This is given as numbers 
;                in an array
;    wave_red  - red wave image.  If not supplied, it assumes
;                'wave-red_NeXeAr.fits' exists and is appropriate.
;    wave_blue - blue wave image.  if not suplied, it assumes
;                'wave-blue_NeXeAr.fits' exists and is appropriate.
;    outname   - string for use in naming the outfile
;    boxcar    - to use a boxcar smooth rather than an elliptical gaussian
;    convbeam  - Convolution beam used to smooth the sky. Given as an array
;                [dispersion direction,crossdispersion] in pixels. The default 
;                is [2,1].
;    centerLine- boundaries for where in the spectrum you want to search to center on 
;                the spectrum.  Default: [4850,5070, 6650,6800], i.e. [OIII] 5007 in 
;                blue and [SII] in red.
;    clobber   - Clobber existing images.  Otherwise they will not be overwritten.
;    noPlotAp  - turn off plotting of the aperture you are working on.
;    lw        - line-width in Angstroms to use in blanking strong lines when a skyslit 
;                is used. Default is 5A.
;    red_slits - red slit structure file
;    blue_slits - blue slit structure file
;    trim_??   - regions to trim at the top and bottom of every slit when fitting the sky.
;                rt = red top, rb = red bottom, bt = blue top, bb = blue bottom
;    ill_red   - red illumination file name
;    ill_blue  - blue illumination file name
;    illumination_corr - turn on illumination correction.  Not fully implemented.
;
; COMMENTS:
;
; EXAMPLE:
;   mods_extract1d,'Science/2Dsub_red_F2_Fin.fits','Science/2Dsub_blue_F2_Fin.fits',z=0.001544,outname='F2'
; BUGS:
;
; PROCEDURES CALLED:
;    bspline_iterfit
;    mrdfits
;    mwrfits
;
; Author:
;    K. Croxall, OSU Astronomy Dept.
;    croxall@astronomy.ohio-state.edu
;    2013 Jul 23
;------------------------------------------------------------------------------

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
        ap_im = x_ordrectify(file,xx1_ap[*],xx2_ap[*],/nocorrect)
        ap_wave = x_ordrectify(wave,xx1_ap[*],xx2_ap[*],/nocorrect)
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
	result_y = result_y[0]
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
plot,bluew[blueg],bluest>(-5),Position=[0.14,0.05,0.50,0.28],XRANGE=[xb0,xb1],xstyle=1,$
	ytitle='Av Spectrum (entire slit)',charsize=1.9
line_x = [bluew[bluex],bluew[bluex]]
line_x1 = [bluew[bluex-9],bluew[bluex-9]]
line_x2 = [bluew[bluex+10],bluew[bluex+10]]
minval = !Y.CRANGE[0]
maxval = !Y.CRANGE[1]
line_y = [minval,maxval]
oplot,line_x,line_y,Color=cgColor('orange')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')
plot,redw[redg],redst>(-5),Position=[0.59,0.05,0.95,0.28],$
        /noerase, XRANGE=[xr0,xr1],xstyle=1,charsize=1.9
line_x = [redw[redx],redw[redx]]
line_x1 = [redw[redx-9],redw[redx-9]]
line_x2 = [redw[redx+10],redw[redx+10]]
minval = !Y.CRANGE[0]
maxval = !Y.CRANGE[1]
line_y = [minval,maxval]
oplot,line_x,line_y,Color=cgColor('orange')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')

;middle pannel
stacksz = size(bluest2,/dim)
bluest2 = bluest2 +10.D
plot,indgen(stacksz),bluest2,XRANGE=[0,stacksz],$
	position=[0.14,0.35,0.50,0.62],ytitle='Slit Profile',$
        /noerase,yrange=[min(bluest2)>(0.1),max(bluest2)],xstyle=1, $
        charsize=1.9,ystyle=1,/ylog
line_x = [bluey,bluey]
minval = 10^(!Y.CRANGE[0])
maxval = 10^(!Y.CRANGE[1])
line_y = [minval,maxval]
line_x1 = [bluey+upper/0.12,bluey+upper/0.12]
line_x2 = [bluey-lower/0.12,bluey-lower/0.12]
niter = size(blank_lowb,/dim)
niter = niter[0]-1
for i=0,niter do polyfill,[blank_lowb[i]>0,blank_lowb[i]>0,blank_hib[i]<stacksz,blank_hib[i]<stacksz], $
	[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data
plot,indgen(stacksz),bluest2,XRANGE=[0,stacksz],/ylog,$
        position=[0.14,0.35,0.50,0.62],charsize=1.9,ystyle=1,$
        /noerase,yrange=[min(bluest2)>(0.1),max(bluest2)],xstyle=1,xtitle='Pixels'
oplot,line_x,line_y,Color=cgColor('orange')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')
axis,xaxis=1,XRANGE=[0,stacksz*0.12],xtitle='ARCSEC'

stacksz = size(redst2,/dim)
redst2 = redst2+10.D
plot,indgen(stacksz),redst2,XRANGE=[0,stacksz],$
	position=[0.59,0.35,0.95,0.62],$
        /noerase,yrange=[min(redst2)>(0.1),max(redst2)],xstyle=1, $
        charsize=1.9,ystyle=1,/ylog
line_x = [redy,redy]
minval = 10^(!Y.CRANGE[0])
maxval = 10^(!Y.CRANGE[1])
line_y = [minval,maxval]
line_x1 = [redy+upper/0.123,redy+upper/0.123]
line_x2 = [redy-lower/0.123,redy-lower/0.123]
for i=0,niter do polyfill,[blank_lowr[i],blank_lowr[i],blank_hir[i],blank_hir[i]],$
	[minval,maxval,maxval,minval],color=cgColor('firebrick')
oplot,line_x,line_y,Color=cgColor('orange')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')
plot,indgen(stacksz),redst2,XRANGE=[0,stacksz],/ylog,$
        position=[0.59,0.35,0.95,0.62],charsize=1.9,ystyle=1,$
        /noerase,yrange=[min(redst2)>(0.1),max(redst2)],xstyle=1,xtitle='Pixels'
axis,xaxis=1,XRANGE=[0,stacksz*0.123],xtitle='ARCSEC'

;top panel
plot,bluew,bluesp,position=[0.14,0.72,0.45,0.99],title='BLUE',$
        /noerase,XRANGE=[xb0-100,xb1+50],ytitle='Extracted Spectrum', $
        charsize=1.9,ystyle=1;,YRANGE=[200,900]
plot,redw,redsp,position=[0.59,0.72,0.95,0.99],title='RED',$
	/noerase,XRANGE=[xr0-100,xr1+100],$;ytitle='Av. Sky spectrum', $
        charsize=1.9,ystyle=1;,YRANGE=[200,900]

return
end

;=========================================================================
; mods_extract1d - extract and calibrate 1D dual channel MODS spectra
;			from 2D spectra
;=========================================================================
pro mods_extract1d,red_scifile,blue_scifile, $
	red_cal = red_cal, $
        blue_cal=blue_cal, $
        z = z, $
	apertures=apertures, $
        outname=outname, $
        wave_blue=wave_blue,wave_red=wave_red, $
	centerLine=centerLine, $
        trim_rt = trim_rt, trim_rb = trim_rb, trim_bb = trim_bb, trim_bt = trim_bt,$
	noPlotAp=noPlotAp, $
        red_slits = red_slits, blue_slits = blue_slits, $
        force_blue=force_blue, force_red=force_red, $
	ill_red=illred, ill_blue=ill_blue, $
        illumination_corr = illumination_corr, $
	clobber = clobber, $
	scale_minmax=scale_minmax

if not keyword_set(blue_cal) then blue_cal = 'bluesensfunc.fits'
if not keyword_set(red_cal) then red_cal = 'redsensfunc.fits'

if (n_elements(redill) EQ '') and (keyword_set(illumination_corr)) then begin
	redill = mrdfits('Proc/red_ill.fits',0,rilhd)
	redill = transpose(redill)
endif
if (n_elements(blueill) EQ '') and (keyword_set(illumination_corr)) then begin
	blueill = mrdfits('Proc/blue_ill.fits',0,bilhd)
	blueill = transpose(blueill)
endif

;deal with parameters that were not entered
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
	rlow_cent = 6650*z
	rhi_cent = 6800*z
	blow_cent = 4890*z
	bhi_cent =  5070*z
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

; Science image
verify_red = file_test(red_scifile)
verify_blue = file_test(blue_scifile)
if verify_red then begin
	red_sciimg = mrdfits(red_scifile,0,red_hdr)
	red_skyimg = mrdfits(red_scifile,3,red_hdrsky)
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
	blue_skyimg = mrdfits(blue_scifile,3,blue_hdrsky)
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

;outfile
if n_elements(outname) EQ '' then begin
        redstub = STRSPLIT(red_scifile,'_',/extract)
        redout = 'Science/1Dsub_'+red_channel_code+'_'  + redstub[1] + '.fits'
        redoutcal = 'Science/1Dsub_'+red_channel_code+'_'  + redstub[1] + '_cal.fits'
        bluestub = STRSPLIT(blue_scifile,'_',/extract)
        blueout = 'Science/1Dsub_'+blue_channel_code+'_' + bluestub[1] + '.fits'
        blueoutcal = 'Science/1Dsub_'+blue_channel_code+'_' + bluestub[1] + '_cal.fits'
endif else begin
        redout = 'Science/1Dsub_'+red_channel_code+'_'  + outname + '.fits'
        redoutcal = 'Science/1Dsub_'+red_channel_code+'_' + outname + '_cal.fits'
        blueout = 'Science/1Dsub_'+blue_channel_code+'_' + outname + '.fits'
        blueoutcal = 'Science/1Dsub_'+blue_channel_code+'_' + outname + '_cal.fits'
endelse
verify_red = file_test(redout)
verify_blue = file_test(blueout)
if (verify_red or verify_blue) and (not keyword_set(clobber)) then begin
        print,'=========================================================='
        if verify_red then print,redout,' already exists. If you want to overwite it please use /clobber'
        if verify_blue then print,blueout,' already exists. If you want to overwite it please use /clobber'
        stop
endif
verify_red = file_test(redoutcal)
verify_blue = file_test(blueoutcal)
if (verify_red or verify_blue) and (not keyword_set(clobber)) then begin
        print,'=========================================================='
        if verify_red then print,redoutcal,' already exists. If you want to overwite it please use /clobber'
        if verify_blue then print,blueoutcal,' already exists. If you want to overwite it please use /clobber'
        stop
endif

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

; Calibration images
verify_red = file_test(red_cal)
verify_blue = file_test(blue_cal)
if verify_red then calR = mrdfits(red_cal,1,chd) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid RED calibration-image'
                print,red_cal,' does not exist.'
                print,"Use: red_cal='FILENAME'"
                stop
        endelse
if verify_blue then calB = mrdfits(blue_cal,1,chd) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid BLUE calibration-image'
                print,blue_cal,' does not exist.'
                print,"Use: blue_cal='FILENAME'"
                stop
        endelse

; Extinction file
extfile = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/mods_extiction.dat'
readcol,extfile,ext_wave,ext_val,/silent
ext_coef = bspline_iterfit(ext_wave,ext_val,bkspace=3,nord=3)

;establish a uniform grid
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

new_red_wave = fltarr(9001) & new_blue_wave = fltarr(5201)
new_red_wave[0] = 5500 & new_blue_wave[0] = 3200
  for i=1,9001-1 do new_red_wave[i] = new_red_wave[i-1]+0.5
  for i=1,5201-1 do new_blue_wave[i] = new_blue_wave[i-1]+0.5

red_sensfunc = bspline_valu(new_red_wave,calR)
blue_sensfunc = bspline_valu(new_blue_wave,calB)

extR = bspline_valu(new_red_wave,ext_coef)
extB = bspline_valu(new_blue_wave,ext_coef)

;open a window for plotting
window,1,retain=2;,XSize=1000,YSize=1000
if (not keyword_set(noPlotAp)) then window,2,retain=2

dims=size(red_sciimg,/dim)
nx=dims[0]
ny=dims[1]
red_instrument = strcompress(sxpar(red_hdr[*, 0], 'INSTRUME'), /rem)
blue_instrument = strcompress(sxpar(blue_hdr[*, 0], 'INSTRUME'), /rem)
mms = strcompress(sxpar(red_hdr[*,0], 'MASKINFO')+'.mms', /rem)

;create the slit structures
  if n_elements(red_slits) eq '' then red_slits = 'slits-'+red_channel_code+'_ill_'+mask+'.fits'
  if n_elements(blue_slits) eq '' then blue_slits = 'slits-'+blue_channel_code+'_ill_'+mask+'.fits'
  verify_red = file_test(red_slits)
  verify_blue = file_test(blue_slits)
  if verify_red then begin
	red_apmask = mrdfits(red_slits,0,r_mask_head)
	red_slits = mrdfits(red_slits,1,r_head)
  endif else begin
        print,'------------------------------------------------------'
        print,'You need a valid RED slit-image'
        print,red_slits,' does not exist.'
        print,"Use: red_slits='FILENAME'"
        stop
   endelse
   if verify_blue then begin
	blue_apmask = mrdfits(blue_slits,0,b_mask_head)
	blue_slits = mrdfits(blue_slits,1,b_head)
   endif else begin
        print,'------------------------------------------------------'
        print,'You need a valid BLUE slit-image'
        print,blue_slits,' does not exist.'
        print,"Use: blue_slits='FILENAME'"
        stop
   endelse

   dim_t = red_slits[0].dims
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
   red_illspec = fltarr(ny_t)+1

   blue_xx1 = blue_slits[0].xx1
   blue_xx2 = blue_slits[1].xx2

   blue_objtrace = blue_xx1*0
   blue_objspec = fltarr(ny_t)+1
   blue_objwave = fltarr(ny_t)+1
   blue_varspec = fltarr(ny_t)+1
   blue_skyspec = fltarr(ny_t)+1
   blue_illspec = fltarr(ny_t)+1

   blue_objspec_cal = blue_sensfunc
   red_objspec_cal = red_sensfunc
   blue_varspec_cal = fltarr(5201)
   red_varspec_cal = fltarr(9001)
   blue_skyspec_cal = fltarr(5201)
   red_skyspec_cal = fltarr(9001)

   blue_illspec_cal = fltarr(5201)
   red_illspec_cal = fltarr(9001)


;extraction struct
red_extract = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D], $
	ra_deg: double(red_slits[0].tarra2), dec_deg: double(red_slits[0].tardec2), $
	ra_sex: red_slits[0].tarra, dec_sex: red_slits[0].tardec}
blue_extract = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D], $
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
	red_extract_obj = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D], $
		ra_deg: red_slits[0].apra2[ap-1], dec_deg: red_slits[0].apdec2[ap-1], $
		ra_sex: red_slits[0].apra[ap-1], dec_sex: red_slits[0].apdec[ap-1]}
	blue_extract_obj = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D], $
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
		red_result_x,red_result_y,red_good,red_wav1d,$
		red_stack,red_stack2,red_xx1,red_xx2,trim_rb,trim_rt
        center,blue_sciimg,blue_waveimg,ap,blow_cent,bhi_cent,$
		blue_result_x,blue_result_y,blue_good,blue_wav1d,$
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

	red_ill1d = fltarr(ny) + 1.D
        blue_ill1d = fltarr(ny) + 1.D

	modsadcalc,blue_scifile,blue_wav1d[blue_result_x],blue_wav1d,bxdelt,bydelt
	modsadcalc,red_scifile,red_wav1d[red_result_x],red_wav1d,rxdelt,rydelt

	for iii=0,ny-1 do red_objtrace[iii,ap-1] = red_xx1[iii,ap-1] + trim_rb + red_result_y
        for iii=0,ny-1 do blue_objtrace[iii,ap-1] = blue_xx1[iii,ap-1] + trim_bb + blue_result_y
	
        ; added in for method demonstration
        blue_orig = blue_objtrace[*,ap-1]
	
	blue_objtrace[*,ap-1] = blue_objtrace[*,ap-1] - bydelt/0.12
        red_objtrace[*,ap-1] = red_objtrace[*,ap-1] - rydelt/0.123

	for iii=0,ny-1 do begin
		bot = red_objtrace[iii,ap-1]-lower/0.123
                top = red_objtrace[iii,ap-1]+upper/0.123

		botint = total(bot,/integer)
		topint = total(top,/integer)
		red_spec1d[iii] = total(red_sciimg[bot:top,iii]) $
			+(1-(bot-botint))*total(red_sciimg[bot,iii]) $
			+(top-topint)*total(red_sciimg[top+1,iii])
		red_wav1d[iii] = median(red_waveimg[bot:top,iii])
		red_var1d[iii] = sqrt(total(red_varimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(red_varimg[bot,iii]) $
                        +(top-topint)*total(red_varimg[top+1,iii]))
		red_sky1d[iii] = sqrt(total(red_skyimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(red_skyimg[bot,iii]) $
                        +(top-topint)*total(red_skyimg[top+1,iii]))
		if keyword_set(illumination_corr) then $
                        red_ill1d[iii] = median(redill[bot:top,iii]) 


                bot = blue_objtrace[iii,ap-1]-lower/0.12
                top = blue_objtrace[iii,ap-1]+upper/0.12
                botint = total(bot,/integer)
                topint = total(top,/integer)
                blue_spec1d[iii] = total(blue_sciimg[bot:top,iii])$
                        +(1-(bot-botint))*total(blue_sciimg[bot,iii]) $
                        +(top-topint)*total(blue_sciimg[top+1,iii])
		blue_wav1d[iii] = median(blue_waveimg[bot:top,iii])
		blue_var1d[iii] = sqrt(total(blue_varimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(blue_varimg[bot,iii]) $
                        +(top-topint)*total(blue_varimg[top+1,iii]))
		blue_sky1d[iii] = sqrt(total(blue_skyimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(blue_skyimg[bot,iii]) $
                        +(top-topint)*total(blue_skyimg[top+1,iii]))
		if keyword_set(illumination_corr) then $
			blue_ill1d[iii] = median(blueill[bot:top,iii]) 
	endfor
	; plot slit and traces
	if (not keyword_set(noPlotAp)) then begin
		wset,2
	        ncolors = !D.Table_Size
	        device,DECOMPOSED=0
	        LoadCT,39
	        bmin = min(blue_xx1[*,ap-1])+trim_bb
	        bmax = max(blue_xx2[*,ap-1])-trim_bt
	        rmin = min(red_xx1[*,ap-1])+trim_rb
	        rmax = max(red_xx2[*,ap-1])-trim_rt

	        blue_slitpos = where(blue_apmask eq ap)
		red_slitpos = where(red_apmask eq ap)
		if keyword_set(scale_minmax) then begin
			bscale_min =  min(blue_sciimg[blue_slitpos])
			bscale_max =  max(blue_sciimg[blue_slitpos])
			rscale_min =  min(red_sciimg[red_slitpos])
			rscale_max =  max(red_sciimg[red_slitpos])
		endif else begin
                        bscale_min = -5
                        bscale_max = 60*median(blue_sciimg[blue_slitpos])>10.
                        rscale_min = -5
                        rscale_max = 60*median(red_sciimg[red_slitpos])>10.
		endelse

	        slitimb = transpose(BytScl(blue_sciimg[bmin:bmax,*], $
			min=bscale_min,max=bscale_max,/nan))
	        slitimr = transpose(BytScl(red_sciimg[rmin:rmax,*], $
			min=rscale_min,max=rscale_max,/nan))
	        TVImage,slitimb,Position=[0.02,0.08,0.99,0.49],/erase
		info = 'BLUE - aperture '+strcompress(ap)
                xyouts,300,bmax-bmin-20,info,CharThick=2,Size=2.3
	        TVImage,slitimr,Position=[0.02,0.52,0.99,0.99]
		info = 'RED - aperture '+strcompress(ap)
                xyouts,300,rmax-rmin-20,info,CharThick=2,Size=2.3
		pathb = BLUE_OBJTRACE[*,ap-1] - bmin
	        pathr = RED_OBJTRACE[*,ap-1] - rmin
	        blue_orig = blue_orig - bmin
		xpix = indgen(8192)

		plot,xpix,pathb,/noerase,Position=[0.02,0.08,0.99,0.49],$
			xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin],linestyle=2,color=cgColor('orange')
	        plot,xpix,pathb-lower/0.12,/noerase,Position=[0.02,0.08,0.99,0.49],$
			xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin]
	        plot,xpix,pathb+upper/0.12,/noerase,Position=[0.02,0.08,0.99,0.49],$
			xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin]
	        plot,xpix,pathr,/noerase,Position=[0.02,0.52,0.99,0.99],$
	                xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax-rmin],linestyle=2,color=cgColor('orange')
	        plot,xpix,pathr-lower/0.123,/noerase,Position=[0.02,0.52,0.99,0.99],$
	                xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax-rmin]
	        plot,xpix,pathr+upper/0.123,/noerase,Position=[0.02,0.52,0.99,0.99],$
	                xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax-rmin]
	endif

	;Plot up the centering and extractions
	wset,1
	extractionplot,blue_wav1d,red_wav1d,blue_good,red_good,$
	       blue_stack,red_stack,blue_result_x,red_result_x,$
	       blue_stack2,red_stack2,blue_result_y,red_result_y,$
	       blue_spec1d,red_spec1d,lower,upper

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
	                readnum,blue_result_y
	                Print,'New center red:'
	                readnum,red_result_y
	                GOTO,LeSigh
	        endif else if redo eq 'n' then begin
	                print,'ReSize extraction? (y/n/bye)'
	                redo = ''
	                read,': ',redo
	                if redo eq 'y' then begin
	                        Print,'New Lower:'
	                        readnum,lower
	                        Print,'New Upper:'
	                        readnum,upper
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

        ;Flexure correction
	wset,2
	if (keyword_set(force_blue) and keyword_set(force_red)) then begin
	   mods_flexurecorr,blue_wav1d,blue_spec1d,blue_sky1d,blue_var1d,'blue',$
		z=z,/plot_prog,force=force_blue,delta=bdelta,slope=bslope
	   mods_flexurecorr,red_wav1d,red_spec1d,red_sky1d,red_var1d,'red',$
		z=z,/plot_prog,force=force_red,delta=rdelta,slope=rslope
	endif else begin
	   mods_flexurecorr,blue_wav1d,blue_spec1d,blue_sky1d,blue_var1d,'blue',$
		z=z,/plot_prog,delta=bdelta,slope=bslope
	   mods_flexurecorr,red_wav1d,red_spec1d,red_sky1d,red_var1d,'red',$
		z=z,/plot_prog,delta=rdelta,slope=rslope
	endelse
	blue_extract[object-1].delta = bdelta
	red_extract[object-1].delta = rdelta
        blue_extract[object-1].slope = bslope
        red_extract[object-1].slope = rslope


;;;;;;;;; native resolution calibrations  No differences seen ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 	red_sensfunc_norebin = bspline_valu(red_wav1d,calR)                                 ;;
;;	blue_sensfunc_norebin = bspline_valu(blue_wav1d,calB)                               ;;
;;	extR_norebin = bspline_valu(red_wav1d,ext_coef)                                     ;;
;;	extB_norebin = bspline_valu(blue_wav1d,ext_coef)                                    ;;
;;                                                                                          ;;
;;       ;The extinction correction is given by the factor                                  ;;
;;       extBcor = 10. ^ (0.4 * airmassB * extB_norebin[*])                                 ;;
;;       extRcor = 10. ^ (0.4 * airmassR * extR_norebin[*])                                 ;;
;;	calfluxR_norebin = (red_spec1d*extRcor)/(time) * 10.^(red_sensfunc_norebin/2.5)     ;;
;;	calfluxB_norebin = (blue_spec1d*extBcor)/(time) * 10.^(blue_sensfunc_norebin/2.5)   ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        ;Rebin the spectrum to the region covered by the standards in actuallity
        x_specrebin,blue_wav1d,blue_spec1d,new_blue_wave,corflux_blue,/FLAMBDA
        x_specrebin,red_wav1d,red_spec1d,new_red_wave,corflux_red,/FLAMBDA

	x_specrebin,blue_wav1d,blue_var1d,new_blue_wave,corvar_blue,/FLAMBDA
	x_specrebin,red_wav1d,red_var1d,new_red_wave,corvar_red,/FLAMBDA

	x_specrebin,blue_wav1d,blue_sky1d,new_blue_wave,corsky_blue,/FLAMBDA
        x_specrebin,red_wav1d,red_sky1d,new_red_wave,corsky_red,/FLAMBDA

        x_specrebin,blue_wav1d,blue_ill1d,new_blue_wave,corill_blue,/FLAMBDA
        x_specrebin,red_wav1d,red_ill1d,new_red_wave,corill_red,/FLAMBDA

        sset=bspline_iterfit(new_blue_wave,corill_blue,bkspace=25,nord=3)
        corill_blue = bspline_valu(new_blue_wave,sset)
	sset=bspline_iterfit(new_red_wave,corill_red,bkspace=25,nord=3)
	corill_red = bspline_valu(new_red_wave,sset)

        ;The extinction correction is given by the factor
        extBcor = 10. ^ (0.4 * airmassB * extB[*])
        extRcor = 10. ^ (0.4 * airmassR * extR[*])

        ;Flux the spectrum including the LBT extinction curve
        calfluxB = (corflux_blue*extBcor)/(time) * 10.^(blue_sensfunc/2.5)
        calfluxR = (corflux_red*extRcor)/(time) * 10.^(red_sensfunc/2.5)
	
        calvarB = (corvar_blue*extBcor)/(time) * 10.^(blue_sensfunc/2.5)
	calvarR = (corvar_red*extRcor)/(time) * 10.^(red_sensfunc/2.5)

	calskyB = (corsky_blue*extBcor)/(time) * 10.^(blue_sensfunc/2.5)
	calskyR = (corsky_red*extRcor)/(time) * 10.^(red_sensfunc/2.5)

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

        blue_objspec_cal = [[blue_objspec_cal],[calfluxB]]
        blue_varspec_cal = [[blue_varspec_cal],[calvarB]]
        blue_skyspec_cal = [[blue_skyspec_cal],[calskyB]]
        blue_illspec_cal = [[blue_illspec_cal],[corill_blue]]

        red_objspec_cal = [[red_objspec_cal],[calfluxR]]
	red_varspec_cal = [[red_varspec_cal],[calvarR]]
	red_skyspec_cal = [[red_skyspec_cal],[calskyR]]
        red_illspec_cal = [[red_illspec_cal],[corill_red]]

        ; Save the output
        mwrfits,red_objspec,redout,red_hdr,/create
	mwrfits,red_varspec,redout
	mwrfits,red_skyspec,redout
	mwrfits,red_objwave,redout
        mwrfits,red_extract,redout

        mwrfits,red_objspec_cal,redoutcal,red_hdrcal,/create
	mwrfits,red_varspec_cal,redoutcal
	mwrfits,red_skyspec_cal,redoutcal
	mwrfits,red_extract,redoutcal
	mwrfits,red_illspec_cal,redoutcal

        mwrfits,blue_objspec,blueout,blue_hdr,/create
	mwrfits,blue_varspec,blueout
	mwrfits,blue_skyspec,blueout
	mwrfits,blue_objwave,blueout
        mwrfits,blue_extract,blueout

        mwrfits,blue_objspec_cal,blueoutcal,blue_hdrcal,/create
	mwrfits,blue_varspec_cal,blueoutcal
	mwrfits,blue_skyspec_cal,blueoutcal
	mwrfits,blue_extract,blueoutcal
	mwrfits,blue_illspec_cal,blueoutcal

	; Add more than region in a slit
	Print,'Extract an additional Region? (y/n/bye)'
	read,': ',redo
	if redo eq 'y' then begin
		GOTO,LeObject
	endif else if redo eq 'bye' then begin
                print,' '
                GOTO,LeEnd
	endif else if redo eq 'n' then begin
		print,'Moving on'
		print,' '
	endif

endfor
LeEnd: print,'Summary --------------------------------------------------'
print,'Extracted ',object-1,' objects from ',ap,' apertures'

print,'Average Flexure in blue: ', mean(blue_extract[1:object-1].delta),'+/-', $
	stddev(blue_extract[1:object-1].delta)/(n_elements(blue_extract[1:object-1].delta))
print,'Average Flexure in red: ',mean(red_extract[1:object-1].delta),'+/-', $
	stddev(red_extract[1:object-1].delta)/(n_elements(red_extract[1:object-1].delta))

print,'Average linear component in blue: ', mean(blue_extract[1:object-1].slope[0]),'+/-', $
        stddev(blue_extract[1:object-1].slope[0])/(n_elements(blue_extract[1:object-1].slope[0]))
print,'                          ', mean(blue_extract[1:object-1].slope[1]),'+/-', $
        stddev(blue_extract[1:object-1].slope[1])/(n_elements(blue_extract[1:object-1].slope[1]))

print,'Average linear component in red: ', mean(red_extract[1:object-1].slope[0]),'+/-', $
        stddev(red_extract[1:object-1].slope[0])/(n_elements(red_extract[1:object-1].slope[0]))
print,'                          ', mean(red_extract[1:object-1].slope[1]),'+/-', $
        stddev(red_extract[1:object-1].slope[1])/(n_elements(red_extract[1:object-1].slope[1]))


print,'If extracting again, please consider using:  force_blue=[',$
	strcompress(mean(blue_extract[1:object-1].delta)),',',$
	strcompress(mean(blue_extract[1:object-1].slope[0])),',',$
        strcompress(mean(blue_extract[1:object-1].slope[1])),']'
print,'If extracting again, please consider using:  force_red=[',$
        strcompress(mean(red_extract[1:object-1].delta)),',',$
        strcompress(mean(red_extract[1:object-1].slope[0])),',',$
        strcompress(mean(red_extract[1:object-1].slope[1])),']'

print,'------------------------------------------------------------------'

end
