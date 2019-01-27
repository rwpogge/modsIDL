;+
; NAME: 
;    mods_skyfit2d
;
; PURPOSE:
;    Perform 2D sky subtraction on DUAL mode MODS GRATING data
;    using a 2D B-spline.  Allows for interactive masking of
;    multiple emission sources withing the slit.
;
; CALLING SEPUENCE:
;    mods_skyfit2d, red_image,blue_image
;
; INPUTS:
;    red_image  - red channel science image
;    blue_image - blue channel science image
;
; OPTIONAL INUTS:
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
;
; COMMENTS:
;
; EXAMPLE:
;   mods_skyfit2d,'Science/std-mods1r.BD33.20120430.fits.gz',$
;                 'Science/std-mods1b.BD33.20120430.fits.gz',$
;                 outname='BD33',aperture=[1],$
;                 red_slits='slits-mods1r.Feige34.20120430.fits',$
;                 blue_slits='slits-mods1b.Feige34.20120430.fits',/clobber
;
; BUGS:
;
; PROCEDURES CALLED:
;    bspline_iterfit
;    mrdfits
;    filter_image
;    mwrfits
;
; Author:
;    K. Croxall, OSU Astronomy Dept.
;    croxall@astronomy.ohio-state.edu
;    2013 Jul 23
;------------------------------------------------------------------------------

;=========================================================================
;  center -- center on line emission in the slit
;=========================================================================
pro center,file,wave,slit_num,wav_min,wav_max,result_x,$
	result_y,goodwav,wav1d,stack,stack2,xx1,xx2,$
	trim1,trim2,skyslitfit = skyslitfit

if (NOT keyword_set(SKYSLITFIT)) then skyslitfit = 0 

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
	if skyslitfit then result_y=2 $;where(stack2 eq median(stack2)) $
       		else result_y=where(stack2 eq max(stack2))
return
end

;=========================================================================
;  extractionplot -- plot the extraction cuts
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
  plot,bluew[blueg],bluest,Position=[0.14,0.05,0.50,0.28],XRANGE=[xb0,xb1],xstyle=1,$
	ytitle='Av Spectrum (entire slit)',charsize=1.9
  line_x = [bluew[bluex],bluew[bluex]]
  line_x1 = [bluew[bluex-9],bluew[bluex-9]]
  line_x2 = [bluew[bluex+10],bluew[bluex+10]]
  minval = !Y.CRANGE[0]
  maxval = !Y.CRANGE[1]
  line_y = [minval,maxval]
  oplot,line_x,line_y,Color=cgColor('dodger blue')
  oplot,line_x1,line_y,Color=cgColor('beige')
  oplot,line_x2,line_y,Color=cgColor('beige')
  plot,redw[redg],redst,Position=[0.59,0.05,0.95,0.28],$
        /noerase, XRANGE=[xr0,xr1],xstyle=1,charsize=1.9
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
	position=[0.14,0.33,0.50,0.62],ytitle='Slit Profile',$
        /noerase,yrange=[min(bluest2),max(bluest2)],xstyle=1, $
	charsize=1.9,/ylog,ystyle=1
  line_x = [bluey,bluey]
  minval = 10^(!Y.CRANGE[0])
  maxval = 10^(!Y.CRANGE[1])
  line_y = [minval,maxval]
  line_x1 = [bluey+upper,bluey+upper]
  line_x2 = [bluey-lower,bluey-lower]
  niter = size(blank_lowb,/dim)
  niter = niter[0]-1
  for i=0,niter do polyfill,[blank_lowb[i]>0,blank_lowb[i]>0,blank_hib[i]<stacksz,blank_hib[i]<stacksz], $
	[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data
  polyfill,[0,0,1,1],[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data
  polyfill,[stacksz-1,stacksz-1,stacksz,stacksz],[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data

  plot,indgen(stacksz),bluest2,XRANGE=[0,stacksz],$
        position=[0.14,0.33,0.50,0.62],charsize=1.9,/ylog,ystyle=1,$
        /noerase,yrange=[min(bluest2),max(bluest2)],xstyle=1
  oplot,line_x,line_y,Color=cgColor('dodger blue')
  oplot,line_x1,line_y,Color=cgColor('beige')
  oplot,line_x2,line_y,Color=cgColor('beige')

  print,'---------------------------------------------------------------'
  info='Blue has '+strcompress(stacksz) + ' pixels across the slit'
  print,info
  stacksz = size(redst2,/dim)
  info = 'Red has '+strcompress(stacksz)+' pixels across the slit'
  print,info

  plot,indgen(stacksz),redst2,XRANGE=[0,stacksz],$
	position=[0.59,0.33,0.95,0.62],$;ytitle='Slit Profile',$
        /noerase,yrange=[min(redst2),max(redst2)],xstyle=1, $
	charsize=1.9,/ylog,ystyle=1
   line_x = [redy,redy]
   minval = 10^(!Y.CRANGE[0])
   maxval = 10^(!Y.CRANGE[1])
   line_y = [minval,maxval]
   line_x1 = [redy+upper,redy+upper]
   line_x2 = [redy-lower,redy-lower]
   for i=0,niter do polyfill,[blank_lowr[i]>0,blank_lowr[i]>0,blank_hir[i]<stacksz,blank_hir[i]<stacksz],$
	[minval,maxval,maxval,minval],color=cgColor('firebrick')
   polyfill,[0,0,1,1],[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data
   polyfill,[stacksz-1,stacksz-1,stacksz,stacksz],[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data
   oplot,line_x,line_y,Color=cgColor('dodger blue')
   oplot,line_x1,line_y,Color=cgColor('beige')
   oplot,line_x2,line_y,Color=cgColor('beige')

  plot,indgen(stacksz),redst2,XRANGE=[0,stacksz],$
        position=[0.59,0.33,0.95,0.62],charsize=1.9,/ylog,ystyle=1,$
        /noerase,yrange=[min(redst2),max(redst2)],xstyle=1

;top panel
  plot,bluew,bluesp,position=[0.14,0.66,0.50,0.97],title='BLUE',$
        /noerase,XRANGE=[xb0-100,xb1+50],ytitle='Av. Sky spectrum', $
	charsize=1.9,ystyle=1
  plot,redw,redsp,position=[0.59,0.66,0.95,0.97],title='RED',$
	/noerase,XRANGE=[xr0-100,xr1+100],$;ytitle='Av. Sky spectrum', $
        charsize=1.9,ystyle=1

return
end

;=========================================================================
;  mods_skyfit2d -- Subtract sky emission from DUAL mode MODS grating data 
;=========================================================================
pro mods_skyfit2d,red_scifile,blue_scifile, $
	skyslit=skyslit, $
	z=z, $
	convbeam=convbeam, $
	apertures=apertures, $
	outname=outname, $
	wave_blue=wave_blue,wave_red=wave_red, $
	centerLine=centerLine, $
	trim_rt = trim_rt, trim_rb = trim_rb, trim_bb = trim_bb, trim_bt = trim_bt, $
	boxcar = boxcar, $
	lw=lw, $
	red_slits = red_slits, blue_slits = blue_slits, $
	noPlotAp=noPlotAp, $
	clobber = clobber

;deal with parameters that were not entered
if n_elements(trim_rt) EQ 0 then trim_rt = 0
if n_elements(trim_bt) EQ 0 then trim_bt = 0

if n_elements(trim_rb) EQ 0 then trim_rb = 0
if n_elements(trim_bb) EQ 0 then trim_bb = 0

if n_elements(lw) EQ 0 then lw = 5. ; line width for skyslit blanking of strong lines

if n_elements(z) EQ '' then begin
   print,''
   print, 'You gave me no redshift.  I will assume (1+z)=1'
   print,''
   wait,0.5
   z=0
endif
z = z+1 ;later code uses z+1 so increment up

if n_elements(skyslit) EQ '' then begin
   print,''
   print, 'You did not list a skyslit.  I assume that means none was cut in the mask (or you have a long-slit observation).'
   print,''
   wait,0.5
   skyslit=0
endif

if n_elements(centerLine) EQ '' then begin
	rlow_cent = 6650*z
	rhi_cent = 6800*z
	blow_cent = 4850*z
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

; Load the images
verify_red = file_test(red_scifile)
verify_blue = file_test(blue_scifile)
if verify_red then red_sciimg = mrdfits(red_scifile,0,red_hdr) $
	else begin
                print,'------------------------------------------------------'
		print,'You need a valid RED science image'
		print,red_scifile,' does not exist. Check your Science folder'
		stop
	endelse
if verify_blue then blue_sciimg = mrdfits(blue_scifile,0,blue_hdr) $
        else begin
                print,'------------------------------------------------------'
		print,'You need a valid BLUE science image'
                print,blue_scifile,' does not exist. Check your Science folder'
                stop
        endelse

red_varimg = mrdfits(red_scifile,1,redv_hdr)
blue_varimg = mrdfits(blue_scifile,1,bluev_hdr)

; Get Mask and Instrument Info
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

if n_elements(outname) EQ '' then begin
        redstub = STRSPLIT(red_scifile,'-',/extract)
        redout = 'Science/2Dsub_'+red_channel_code+'_' + redstub[1]
        bluestub = STRSPLIT(blue_scifile,'-',/extract)
        blueout = 'Science/2Dsub_'+blue_channel_code+'_' + bluestub[1]
endif else begin
        redout = 'Science/2Dsub_' +red_channel_code+'_' + outname + '.fits'
        blueout = 'Science/2Dsub_' +blue_channel_code+'_' + outname + '.fits'
endelse
verify_red = file_test(redout)
verify_blue = file_test(blueout)
if (verify_red or verify_blue) and (not keyword_set(clobber)) then begin
        if verify_red then print,redout,' already exists. If you want to overwite it please use /clobber'
        if verify_blue then print,blueout,' already exists. If you want to overwite it please use /clobber'
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

; Gaussian convolution?
if n_elements(convbeam) EQ '' then convbeam=[1.0,2.0]
if n_elements(boxcar) NE '' then begin
	print,'USING A ',boxcar,' PIXEL BOXCAR TO SMOOTH'
	red_sciimg = smooth(red_sciimg,boxcar)
	blue_sciimg = smooth(blue_sciimg,boxcar)
	red_waveimg = smooth(red_waveimg,boxcar)
	blue_waveimg = smooth(blue_waveimg,boxcar)
endif else begin
	print,'USING A ',convbeam,' PIXEL GAUSSIAN TO SMOOTH'
	red_sciimg = filter_image(red_sciimg,FWHM=convbeam,/ALL)
	blue_sciimg = filter_image(blue_sciimg,FWHM=convbeam,/ALL)
	red_waveimg = filter_image(red_waveimg,FWHM=convbeam,/ALL)
	blue_waveimg = filter_image(blue_waveimg,FWHM=convbeam,/ALL)
endelse

;open a window for plotting
window,1,retain=2;,XSize=1000,YSize=1000
if (not keyword_set(noPlotAp)) then window,2,retain=2

;create a new file to fill with sky-subtracted data
red_fit = red_sciimg
blue_fit = blue_sciimg

dims=size(red_sciimg,/dim)
nx=dims[0]
ny=dims[1]
red_instrument = strcompress(sxpar(red_hdr[*, 0], 'INSTRUME'), /rem)
blue_instrument = strcompress(sxpar(blue_hdr[*, 0], 'INSTRUME'), /rem)
mms = strcompress(sxpar(red_hdr[*,0], 'MASKINFO')+'.mms', /rem)
grat_RED = strcompress(sxpar(red_hdr[*, 0], 'GRATNAME'), /rem)
grat_BLUE = strcompress(sxpar(blue_hdr[*, 0], 'GRATNAME'), /rem)


;update the headers
   header_line = 'mods_skyfit2d beta version 2014-04-15'
   sxaddpar,blue_hdr,'HISTORY',header_line
   sxaddpar,red_hdr,'HISTORY',header_line

   ;wavelength
   header_line = wave_red
   sxaddpar,red_hdr,'WAVEFILE',header_line
   header_line = wave_blue
   sxaddpar,blue_hdr,'WAVEFILE',header_line

   header_line = z
   sxaddpar,red_hdr,'REDSHFT',header_line
   sxaddpar,blue_hdr,'REDSHFT',header_line

   header_line = rlow_cent
   sxaddpar,red_hdr,'LOWCENT',header_line
   header_line = rhi_cent
   sxaddpar,red_hdr,'HICENT',header_line

   header_line = blow_cent
   sxaddpar,blue_hdr,'LOWCENT',header_line
   header_line = bhi_cent
   sxaddpar,blue_hdr,'HICENT',header_line

   ;convbeam
   if n_elements(boxcar) NE '' then begin
	header_line = 'mods_skyfit2d_ - using a boxcar'
	sxaddpar,red_hdr,'HISTORY',header_line
	sxaddpar,blue_hdr,'HISTORY',header_line
   endif
   header_line = 'mods_skyfit2d_ - smoothed for skysub using ['$
        +strtrim(string(convbeam[0]),2)+','+strtrim(string(convbeam[1]),2)+']'
   sxaddpar,red_hdr,'CONVBEAM',header_line
   sxaddpar,blue_hdr,'CONVBEAM',header_line

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
   nslit=size(red_slits.xx1,/dim)
   nslit=nslit[1]
   nx_t = dims[0]
   ny_t = dims[1]
   red_xx1 = red_slits[0].xx1
   red_xx2 = red_slits[1].xx2
   red_objtrace = red_xx1*0
   red_objspec = red_xx1*0
   red_objwave = red_xx1*0
   red_varmask = red_xx1*0

   blue_xx1 = blue_slits[0].xx1
   blue_xx2 = blue_slits[1].xx2
   blue_objtrace = blue_xx1*0
   blue_objspec = blue_xx1*0
   blue_objwave = blue_xx1*0
   blue_varmask = blue_xx1*0

   blue_objspec_cal = fltarr(5201,nslit)
   red_objspec_cal = fltarr(9001,nslit)

;extraction struct
   red_extract = {slit: 0, zones: 0, trim_bot: 0, trim_top: 0, $
	lower1: 0, upper1: 0, lower2: 999, upper2: 999, lower3: 999, upper3: 999, lower4: 999, upper4: 999}
   red_extract = replicate(red_extract,nslit)

   blue_extract = {slit: 0, zones: 0, trim_bot: 0, trim_top: 0, $
	lower1: 0, upper1: 0, lower2: 999, upper2: 999, lower3: 999, upper3: 999, lower4: 999, upper4: 999}
   blue_extract = replicate(blue_extract,nslit)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                   fit the best sky slit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if skyslit ne 0 then begin
	center,red_sciimg,red_waveimg,skyslit,rlow_cent,rhi_cent,$
	      red_result_x,red_result_y,red_good,red_wav1d,$
	      red_stack,red_stack2,red_xx1,red_xx2,trim_rb,trim_rt,skyslitfit=1
	center,blue_sciimg,blue_waveimg,skyslit,blow_cent,bhi_cent,$
	       blue_result_x,blue_result_y,blue_good,blue_wav1d,$
	       blue_stack,blue_stack2,blue_xx1,blue_xx2,trim_bb,trim_bt,skyslitfit=1
	
	lower = 5
	upper = 5
        n_blank = 1

	print,''
	LeSky: print,'Extracting Reference SKY'
	print,''
	red_spec1d = fltarr(ny)
	blue_spec1d = fltarr(ny)
	red_apmask = 0*red_sciimg
	blue_apmask = 0*blue_sciimg
	
	for iii=0,ny-1 do red_objtrace[iii,skyslit-1] = red_xx1[iii,skyslit-1] + red_result_y
	for iii=0,ny-1 do blue_objtrace[iii,skyslit-1] = blue_xx1[iii,skyslit-1] + blue_result_y
	
	for iii=0,ny-1 do begin
		for iiii=red_xx1[iii,skyslit-1]+trim_rb,red_xx2[iii,skyslit-1]-trim_rt do $
	              red_apmask[iiii,iii] = 1
	        for iiii=blue_xx1[iii,skyslit-1]+trim_bb,blue_xx2[iii,skyslit-1]-trim_bt do $
	              blue_apmask[iiii,iii] = 1
	endfor
	
	red_slitpos = where(red_apmask eq 1)
	blue_slitpos = where(blue_apmask eq 1)
	
	for iii=0,ny-1 do begin
		 bot = red_objtrace[iii,skyslit-1]-lower
	         top = red_objtrace[iii,skyslit-1]+upper
	         red_apmask[bot:top,iii] = 0
                 remaining_sky = where(red_apmask[*,iii] eq 1)
                 red_spec1d[iii] = median(red_sciimg[remaining_sky,iii])
                 red_wav1d[iii] = median(red_waveimg[remaining_sky,iii])

	         bot = blue_objtrace[iii,skyslit-1]-lower
	         top = blue_objtrace[iii,skyslit-1]+upper
	         blue_apmask[bot:top,iii] = 0
                 remaining_sky = where(blue_apmask[*,iii] eq 1)
                 blue_spec1d[iii] = median(blue_sciimg[remaining_sky,iii])
                 blue_wav1d[iii] = median(blue_waveimg[remaining_sky,iii])
	endfor

        ; plot slit and traces
        if (not keyword_set(noPlotAp)) then begin
                wset,2
                ncolors = !D.Table_Size
                device,DECOMPOSED=0
                LoadCT,39
                bmin = min(blue_xx1[*,skyslit-1])+trim_bb
                bmax = max(blue_xx2[*,skyslit-1])-trim_bt
                rmin = min(red_xx1[*,skyslit-1])+trim_rb
                rmax = max(red_xx2[*,skyslit-1])-trim_rt
                slitimb = transpose(BytScl(blue_sciimg[bmin:bmax,*],min=-1,max=100,/nan))
                slitimr = transpose(BytScl(red_sciimg[rmin:rmax,*],min=-1,max=500,/nan))
                TVImage,slitimb,Position=[0.02,0.08,0.99,0.49],/erase
                TVImage,slitimr,Position=[0.02,0.52,0.99,0.99]
                pathb = BLUE_OBJTRACE[*,skyslit-1] - bmin
                pathr = RED_OBJTRACE[*,skyslit-1] - rmin
                xpix = indgen(8192)
                plot,xpix,pathb,/noerase,Position=[0.02,0.08,0.99,0.49],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin],Color=cgColor('dodger blue')
                xyouts,300,bmax-bmin-20,'BLUE',CharThick=2,Size=2.3
                plot,xpix,pathb-lower,/noerase,Position=[0.02,0.08,0.99,0.49],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin]
                plot,xpix,pathb+upper,/noerase,Position=[0.02,0.08,0.99,0.49],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin]
                plot,xpix,pathr,/noerase,Position=[0.02,0.52,0.99,0.99],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax-rmin],Color=cgColor('dodger blue')
                xyouts,300,rmax-rmin-20,'RED',CharThick=2,Size=2.3
                plot,xpix,pathr-lower,/noerase,Position=[0.02,0.52,0.99,0.99],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax-rmin]
                plot,xpix,pathr+upper,/noerase,Position=[0.02,0.52,0.99,0.99],$
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
		print,''
		print,'Indicate the region to mask in the sky-slit.'
	        print,'Centered on lines:'
		print,'Blue  - ',blue_wav1d[blue_result_x]
		print,'Red   - ', red_wav1d[red_result_x]
		print,''
	        print,'Peaks found at:'
		print,'Blue  - ',blue_result_y
		print,'Red   - ',red_result_y
		print,''
	        print,'Currently using: lower = ',lower,'  upper = ',upper
	        print,'Re-Center mask? (y/n/bye)'
	        redo = ''
	        read,': ',redo
	        if redo eq 'y' then begin
	                Print,'New center blue:'
	                read,': ',blue_result_y
	                Print,'New center red:'
	                read,': ',red_result_y
	                GOTO,LeSky
	        endif else if redo eq 'n' then begin
		        print,'Re-Size Mask? (y/n/bye)'
	        	redo = ''
	        	read,': ',redo
	        	if redo eq 'y' then begin
	                	Print,'New Lower:'
	                	read,': ',lower
	                	Print,'New Upper:'
	                	read,': ',upper
	                	GOTO,LeSky
			endif else if redo eq 'n' then begin
				print,'Huzzah!'
				Ok = 1
	                endif else if redo eq 'bye' then begin
        	                stop
			endif else print,'I do not understand ',redo
	        endif else if redo eq 'n' then begin
	                print,'Huzzah!'
	                Ok = 1
		endif else if redo eq 'bye' then begin
			stop
	        endif else begin
	                print,'I do not understand ',redo
	        endelse
	endwhile
;
        Print,'Add an additional exclusion zone? (y/n/bye)'
            read,': ',redo
            if redo eq 'y' then begin
                 if n_blank gt 1 then begin
                        redaddmask = where(red_apmask_prev eq 0)
                        blueaddmask = where(blue_apmask_prev eq 0)
                        red_apmask[redaddmask] = 0
                        blue_apmask[blueaddmask] = 0
                 endif
                 n_blank = n_blank + 1
                 red_apmask_prev = red_apmask
                 blue_apmask_prev = blue_apmask
                 ;common exclude1, blank_lowb, blank_lowr, blank_hib, blank_hir
                 blank_lowb = [blank_lowb,blue_result_y-lower]
                 blank_lowr = [blank_lowr,red_result_y-lower]
                 blank_hib = [blank_hib,blue_result_y+upper]
                 blank_hir = [blank_hir,red_result_y+upper]
                 Print,'New center blue:'
                 read,': ',blue_result_y
                 Print,'New center red:'
                 read,': ',red_result_y
                 Print,'New Lower:'
                 read,': ',lower
                 Print,'New Upper:'
                 read,': ',upper

                 GOTO,LeSky
         endif else if redo eq 'bye' then begin
                 stop
         endif else if redo eq 'n' then begin
                 print,'You blanked ',n_blank,' region(s) in this slit.'
                 if n_blank gt 1 then begin
                         redaddmask = where(red_apmask_prev eq 0)
                         blueaddmask = where(blue_apmask_prev eq 0)
                         red_apmask[redaddmask] = 0
                         blue_apmask[blueaddmask] = 0
                 endif
         endif
;	
	;fit the pure sky
	red_sset_sky = bspline_iterfit(red_waveimg[red_slitpos],red_sciimg[red_slitpos],$
	                bkspace=.65,nord=4,invvar=red_apmask[red_slitpos])
	blue_sset_sky = bspline_iterfit(blue_waveimg[blue_slitpos],blue_sciimg[blue_slitpos],$
	                bkspace=.95,nord=2,invvar=blue_apmask[blue_slitpos])
;MARK
	;make a pure skyslit fit
	fullfit_blue = bspline_valu(blue_waveimg,blue_sset_sky)
	fullfit_red = bspline_valu(red_waveimg,red_sset_sky)
	fullfit_blue2 = (blue_sciimg - fullfit_blue)>(-100)
        fullfit_red2 = (red_sciimg - fullfit_red)>(-100)
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                      fit each slit on its own
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(apertures) EQ '' then apertures = indgen(nslit)+1
ndo = size(apertures,/dim)
ndo=ndo[0]

for ii=0,ndo-1 do begin
	blank_lowb = [0]
	blank_lowr = [0]
	blank_hib = [0]
	blank_hir = [0]
	ap=apertures[ii]
	print,'APERTURE: ',ap
	print,''
	n_blank = 1
	NOSKY = 0
        center,red_sciimg,red_waveimg,ap,rlow_cent,rhi_cent,$
		red_result_x,red_result_y,red_good,red_wav1d,$
		red_stack,red_stack2,red_xx1,red_xx2,trim_rb,trim_rt
        center,blue_sciimg,blue_waveimg,ap,blow_cent,bhi_cent,$
		blue_result_x,blue_result_y,blue_good,blue_wav1d,$
		blue_stack,blue_stack2,blue_xx1,blue_xx2,trim_bb,trim_bt

        lower = 10
        upper = 10

	LeSigh: print,'Extracting Spectrum'

        red_apmask = 0*red_sciimg
        blue_apmask = 0*blue_sciimg
        red_spec1d = fltarr(ny)
        blue_spec1d = fltarr(ny)

	for iii=0,ny-1 do red_objtrace[iii,ap-1] = red_xx1[iii,ap-1]+trim_rb + red_result_y
        for iii=0,ny-1 do blue_objtrace[iii,ap-1] = blue_xx1[iii,ap-1] +trim_bb + blue_result_y
	
	for iii=0,ny-1 do begin
		for iiii=red_xx1[iii,ap-1]+trim_rb,red_xx2[iii,ap-1]-trim_rt do $
			red_apmask[iiii,iii] = 1
                for iiii=blue_xx1[iii,ap-1]+trim_bb,blue_xx2[iii,ap-1]-trim_bt do $
			blue_apmask[iiii,iii] = 1
	endfor
	red_slitpos = where(red_apmask eq 1)
	blue_slitpos = where(blue_apmask eq 1)

	for iii=0,ny-1 do begin
		;blank thin strips along top and bottom
		bot = red_xx1[iii,ap-1]-4
		top = red_xx1[iii,ap-1]+2
		red_apmask[bot:top,iii] = 0
		bot = red_xx2[iii,ap-1]-2
		top = red_xx2[iii,ap-1]+4
		red_apmask[bot:top,iii] = 0
		;blank the slected region in the aperture masks
		bot = red_objtrace[iii,ap-1]-lower>0
                top = red_objtrace[iii,ap-1]+upper<3088
		red_apmask[bot:top,iii] = 0
		;extract the remaining sky spectrum?
		remaining_sky = where(red_apmask[*,iii] eq 1)
		red_spec1d[iii] = median(red_sciimg[remaining_sky,iii])
		red_wav1d[iii] = median(red_waveimg[remaining_sky,iii])

		;and for blue
                bot = blue_xx1[iii,ap-1]-4
                top = blue_xx1[iii,ap-1]+2
                blue_apmask[bot:top,iii] = 0
                bot = blue_xx2[iii,ap-1]-2
                top = blue_xx2[iii,ap-1]+4
                blue_apmask[bot:top,iii] = 0
                bot = blue_objtrace[iii,ap-1]-lower>0
                top = blue_objtrace[iii,ap-1]+upper<3088
                blue_apmask[bot:top,iii] = 0
                remaining_sky = where(blue_apmask[*,iii] eq 1)
                blue_spec1d[iii] = median(blue_sciimg[remaining_sky,iii])
		blue_wav1d[iii] = median(blue_waveimg[remaining_sky,iii])
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
                slitimb = transpose(BytScl(blue_sciimg[bmin:bmax,*],min=1,max=10*median(blue_sciimg[bmin:bmax,*]),/nan))
                slitimr = transpose(BytScl(red_sciimg[rmin:rmax,*],min=-1,max=10*median(red_sciimg[rmin:rmax,*]),/nan))
                TVImage,slitimb,Position=[0.02,0.08,0.99,0.49],/erase
                TVImage,slitimr,Position=[0.02,0.52,0.99,0.99]
                pathb = BLUE_OBJTRACE[*,ap-1] - bmin
                pathr = RED_OBJTRACE[*,ap-1] - rmin
                xpix = indgen(8192)
                plot,xpix,pathb,/noerase,Position=[0.02,0.08,0.99,0.49],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin-1],Color=cgColor('dodger blue')
                xyouts,300,bmax-bmin-20,'BLUE',CharThick=2,Size=2.3
                plot,xpix,pathb-lower,/noerase,Position=[0.02,0.08,0.99,0.49],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin-1]
                plot,xpix,pathb+upper,/noerase,Position=[0.02,0.08,0.99,0.49],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin-1]
                plot,xpix,pathr,/noerase,Position=[0.02,0.52,0.99,0.99],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax-rmin-1],Color=cgColor('dodger blue')
                xyouts,300,rmax-rmin-20,'RED',CharThick=2,Size=2.3
                plot,xpix,pathr-lower,/noerase,Position=[0.02,0.52,0.99,0.99],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax-rmin-1]
                plot,xpix,pathr+upper,/noerase,Position=[0.02,0.52,0.99,0.99],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax-rmin-1]
		plot,xpix,red_xx1[*,ap-1]-rmin,/noerase,Position=[0.02,0.52,0.99,0.99],$
			xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax-rmin-1],linestyle=2
                plot,xpix,red_xx2[*,ap-1]-rmin,/noerase,Position=[0.02,0.52,0.99,0.99],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax-rmin-1],linestyle=2
                plot,xpix,blue_xx1[*,ap-1]-bmin,/noerase,Position=[0.02,0.08,0.99,0.49],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin-1],linestyle=2
                plot,xpix,blue_xx2[*,ap-1]-bmin,/noerase,Position=[0.02,0.08,0.99,0.49],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin-1],linestyle=2
        endif

        ;Plot up the centering and extractions
        wset,1
	extractionplot,blue_wav1d,red_wav1d,blue_good,red_good,$
	       blue_stack,red_stack,blue_result_x,red_result_x,$
	       blue_stack2,red_stack2,blue_result_y,red_result_y,$
	       blue_spec1d,red_spec1d,lower,upper

	Ok=0

	while Ok ne 1 do begin
		print,''
	        print,'Centered on lines:'
	        print,'Blue  - ',blue_wav1d[blue_result_x]
	        print,'Red   - ', red_wav1d[red_result_x]
		print,''
	        print,'Peaks found at:' 
	        print,'Blue  - ',blue_result_y
	        print,'Red   - ',red_result_y
		print,''
	        print,'Currently using: lower = ',lower,'  upper = ',upper
	        print,'Re-Center Mask? (y/n/NOSKY/bye)'
	        redo = ''
	        read,': ',redo
	        if redo eq 'y' then begin
	                Print,'New center blue:'
	                read,': ',blue_result_y
	                Print,'New center red:'
	                read,': ',red_result_y
	                GOTO,LeSigh
	        endif else if redo eq 'n' then begin
	                print,'Re-Size Mask? (y/n/bye)'
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
	                        Ok = 1
        	        endif else if redo eq 'bye' then begin
	                        stop
	                        endif else print,'I do not understand ',redo
	        endif else if redo eq 'n' then begin
	                print,'Huzzah!'
	                Ok = 1
                endif else if redo eq 'bye' then begin
                        stop
                endif else if redo eq 'NOSKY' then begin
			NOSKY=1
			Ok = 1
	        endif else begin
	                print,'I do not understand ',redo
	        endelse
	endwhile

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;  Add the ability to include more than one non-sky region in a slit
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        if NOSKY eq 0 then begin
		Print,'Add an additional exclusion zone? (y/n/bye)'
	        read,': ',redo
		if redo eq 'y' then begin
			if n_blank gt 1 then begin
				redaddmask = where(red_apmask_prev eq 0)
				blueaddmask = where(blue_apmask_prev eq 0)
				red_apmask[redaddmask] = 0
				blue_apmask[blueaddmask] = 0
			endif
			n_blank = n_blank + 1
			red_apmask_prev = red_apmask
			blue_apmask_prev = blue_apmask
			;common exclude1, blank_lowb, blank_lowr, blank_hib, blank_hir
			blank_lowb = [blank_lowb,blue_result_y-lower]
			blank_lowr = [blank_lowr,red_result_y-lower]
       	        	blank_hib = [blank_hib,blue_result_y+upper]
        	        blank_hir = [blank_hir,red_result_y+upper]
			Print,'New center blue:'
         	      	read,': ',blue_result_y
         	        Print,'New center red:'
	                read,': ',red_result_y
			Print,'New Lower:'
	                read,': ',lower
	                Print,'New Upper:'
	                read,': ',upper
	 	
			GOTO,LeSigh
	        endif else if redo eq 'bye' then begin
	                stop
		endif else if redo eq 'n' then begin
                        blank_lowb = [blank_lowb,blue_result_y-lower]
                        blank_lowr = [blank_lowr,red_result_y-lower]
                        blank_hib = [blank_hib,blue_result_y+upper]
                        blank_hir = [blank_hir,red_result_y+upper]

			print,'You blanked ',n_blank,' region(s) in this slit.'
			blue_extract[ii].slit = ap
			blue_extract[ii].zones = n_blank
                        red_extract[ii].slit = ap
			red_extract[ii].zones = n_blank

			blue_extract[ii].trim_bot = trim_bb
                        blue_extract[ii].trim_top = trim_bt
                        red_extract[ii].trim_bot = trim_rb
                        red_extract[ii].trim_top = trim_rt

			blue_extract[ii].lower1 = blank_lowb[1]
			blue_extract[ii].upper1 = blank_hib[1]
			red_extract[ii].lower1 = blank_lowr[1]
			red_extract[ii].upper1 = blank_hir[1]
			if n_blank GE 2 then begin 
				blue_extract[ii].lower2 = blank_lowb[2]
				blue_extract[ii].upper2 = blank_hib[2]
                        	red_extract[ii].lower2 = blank_lowr[2]
                        	red_extract[ii].upper2 = blank_hir[2]
			endif
                        if n_blank GE 3 then begin 
                                blue_extract[ii].lower3 = blank_lowb[3]
                                blue_extract[ii].upper3 = blank_hib[3]
                                red_extract[ii].lower3 = blank_lowr[3]
                                red_extract[ii].upper3 = blank_hir[3]
                        endif
                        if n_blank GE 4 then begin 
                                blue_extract[ii].lower4 = blank_lowb[4]
                                blue_extract[ii].upper4 = blank_hib[4]
                                red_extract[ii].lower4 = blank_lowr[4]
                                red_extract[ii].upper4 = blank_hir[4]
                        endif

			if n_blank gt 1 then begin
	                        redaddmask = where(red_apmask_prev eq 0)
	                        blueaddmask = where(blue_apmask_prev eq 0)
	                        red_apmask[redaddmask] = 0
	                        blue_apmask[blueaddmask] = 0
			endif
		endif
	endif
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	; FIT THE SKY IN THE SELETCED APERTURE
        if NOSKY eq 0 then begin
	        sset=bspline_iterfit(red_waveimg[red_slitpos],red_sciimg[red_slitpos],$
			bkspace=.99,nord=3,invvar=red_apmask[red_slitpos])
		extract=red_sciimg[red_slitpos]
		extract_wave=red_waveimg[red_slitpos]
		lw = lw+1
	        red_line = where($;(extract_wave lt (6563+lw)*z and extract_wave gt (6563-lw)*z) or $ ; Halpha
	                (extract_wave lt (6584+lw)*z and extract_wave gt (6584-lw)*z) or $ ; NII 6584
	                (extract_wave lt (6716+lw)*z and extract_wave gt (6716-lw)*z) or $ ; SII 6716
	                (extract_wave lt (6731+lw)*z and extract_wave gt (6731-lw)*z) $ ; SII 6731
	               ; (extract_wave lt (9071+lw)*z and extract_wave gt (9067-lw)*z) or $ ; SIII 9069
	               ; (extract_wave lt (9532+lw)*z and extract_wave gt (9532-lw)*z) $ ; SIII 9532
	                )
		red_line2 = where((extract_wave lt (6548+lw)*z and extract_wave gt (6548-lw)*z) or $  ;NII 6548
			(extract_wave lt (6563+lw+2)*z and extract_wave gt (6563-lw+2)*z)) ;Halpha
                red_scalereg = where((extract_wave lt 6670*z and extract_wave gt 6620*z))
		red_scalereg2 = where((extract_wave lt 6530*z and extract_wave gt 6450*z))

	        yfit = bspline_valu(red_waveimg[red_slitpos],sset)

		if skyslit ne 0 then begin
			yfit_sky = bspline_valu(red_waveimg[red_slitpos],red_sset_sky)
                        scale = median(yfit[red_scalereg]/yfit_sky[red_scalereg])
			scale2 = median(yfit[red_scalereg2]/yfit_sky[red_scalereg2])
			yfit[red_line] = scale*yfit_sky[red_line]
			yfit[red_line2] = scale2*yfit_sky[red_line2]
		endif
	
		red_subimg = red_sciimg
		red_subimg[red_slitpos] = extract-yfit
		red_fit[red_slitpos] = yfit	

;window,3,retain=2
;plot,red_waveimg[red_slitpos],red_sciimg[red_slitpos],xrange=[8800,9000],psym=4,yrange=[-10,200],$
;plot,red_waveimg[red_slitpos],red_sciimg[red_slitpos],xrange=[6500,7000],psym=4,yrange=[-10,200],$
;plot,red_waveimg[1850:2280,*],red_sciimg[1850:2280,*]+50,xrange=[9000,9600],psym=4,yrange=[-10,200],$
;	title='Red Sky Subtraction',xtitle='Wavelength [A]',ytitle='Counts [ADU]',charsize=2,charthick=2,$
;	xthick=3,ythick=3
;oplot,red_waveimg[red_slitpos],yfit+50,color=cgColor('Dodger blue'),psym=6,symsize=0.5
;oplot,red_waveimg[1929,*],red_subimg[1929,*],color=cgColor('Orange'),psym=10,thick=3
;
;oplot,red_waveimg[1903,*],red_sciimg[1903,*],color=cgColor('red'),psym=2
;oplot,red_waveimg[2204,*],red_sciimg[2204,*],color=cgColor('dark green'),psym=1
;oplot,red_waveimg[2204,*],red_subimg[2204,*],color=cgColor('green');,psym=1
;oplot,red_waveimg[1903,*],red_subimg[1903,*],color=cgColor('lime green'),psym=10
;oplot,red_waveimg[2204,*],red_fit[2204,*],color=cgColor('violet')
;stop

		lw = lw-1
	        sset=bspline_iterfit(blue_waveimg[blue_slitpos],blue_sciimg[blue_slitpos],$
	                bkspace=.95,nord=2,invvar=blue_apmask[blue_slitpos])
	        extract=blue_sciimg[blue_slitpos]
	        extract_wave=blue_waveimg[blue_slitpos]
	        blue_line = where((extract_wave lt (5007+lw)*z and extract_wave gt (5007-lw)*z) or $ ;OIII 5007
	                (extract_wave lt (4959+lw)*z and extract_wave gt (4959-lw)*z) or $ ; OIII 4959
	                (extract_wave lt (4861+lw)*z and extract_wave gt (4861-lw)*z)) ; Hbeta
	        blue_line5 = where((extract_wave lt (4102+lw)*z and extract_wave gt (4102-lw)*z)); Hdelta
	        blue_line2 = where((extract_wave lt (3727+lw)*z and extract_wave gt (3727-lw)*z)); OII 3727
                blue_line4 = where((extract_wave lt (4340+lw-1.5)*z and extract_wave gt (4340-lw)*z)) ; Hgamma

		blue_scalereg = where((extract_wave lt 4800*z and extract_wave gt 4700*z))
		blue_scalereg5 = where((extract_wave lt 4200*z and extract_wave gt 4106*z))
                blue_scalereg2 = where((extract_wave lt 3720*z and extract_wave gt 3700*z))
                blue_scalereg4 = where((extract_wave lt 4339*z and extract_wave gt 4300*z))

	        yfit = bspline_valu(blue_waveimg[blue_slitpos],sset)
		if skyslit ne 0 then begin
			yfit_sky = bspline_valu(blue_waveimg[blue_slitpos],blue_sset_sky)
			scale = median(yfit[blue_scalereg]/yfit_sky[blue_scalereg])
                        scale2 = median(yfit[blue_scalereg2]/yfit_sky[blue_scalereg2])
                        scale5 = median(yfit[blue_scalereg5]/yfit_sky[blue_scalereg5])
			scale4 = median(yfit[blue_scalereg4]/yfit_sky[blue_scalereg4])
	        	yfit[blue_line] = scale*yfit_sky[blue_line]
                        yfit[blue_line2] = scale2*yfit_sky[blue_line2]
			yfit[blue_line4] = scale4*yfit_sky[blue_line4]
                        yfit[blue_line5] = scale5*yfit_sky[blue_line5]
		endif
	        blue_subimg = blue_sciimg
;add in some extra extensions for alterate sky sub products?
	        blue_subimg[blue_slitpos] = extract-yfit
	        blue_fit[blue_slitpos] = yfit
	endif else begin
		if skyslit EQ 0 then begin
			print, 'You did not list a skyslit.'
			print, ' Please re-run MODS_SKYFIT2D with a skyslit indicated.'
   			stop
		endif
		print,'Applying the sky solution from Aperture ',skyslit,' to aperture ',ap
                extract=blue_sciimg[blue_slitpos]
		yfit = bspline_valu(blue_waveimg[blue_slitpos],blue_sset_sky)
		blue_subimg = blue_sciimg
                blue_subimg[blue_slitpos] = extract-yfit
                blue_fit[blue_slitpos] = yfit

                extract=red_sciimg[red_slitpos]
		yfit = bspline_valu(red_waveimg[red_slitpos],red_sset_sky)
                red_subimg = red_sciimg
                red_subimg[red_slitpos] = extract-yfit
                red_fit[red_slitpos] = yfit
	endelse

	red_subimg = red_subimg>(-100)
	blue_subimg = blue_subimg>(-100)
	red_sciimg = red_subimg
	blue_sciimg = blue_subimg

	; Save the output
		;0 - sky-subtracted image
		;1 - variance image
		;2 - extraction table
		;3 - sky fit
		;4 - sky fit using only skyslit
		;5 - sky-subtracted image using only skyslit

	mwrfits,red_subimg,redout,red_hdr,/create
	mwrfits,red_varimg,redout
	mwrfits,red_extract,redout
	mwrfits,red_fit,redout
	
	mwrfits,blue_subimg,blueout,blue_hdr,/create
	mwrfits,blue_varimg,blueout
	mwrfits,blue_extract,blueout
	mwrfits,blue_fit,blueout
endfor

	IF skyslit NE 0 then mwrfits,fullfit_red,redout
	IF skyslit NE 0 then mwrfits,fullfit_blue,blueout
        IF skyslit NE 0 then mwrfits,fullfit_red2,redout
        IF skyslit NE 0 then mwrfits,fullfit_blue2,blueout
end
