;+
; NAME: 
;    mods_extract_dualchan
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
;  plot_2d_spec -- plot the 2D emission for a given aperture/slit
;=========================================================================
pro plot_2d_spec,red_sciimg,blue_sciimg,$
                red_xx1,red_xx2,blue_xx1,blue_xx2,$
                red_center_trace=red_center_trace,$
		blue_center_trace=blue_center_trace,$
                lower=lower,upper=upper,$
		red_platescale=red_platescale,$
		blue_platescale=blue_platescale,$
		man_scale = man_scale, $
		title = title

	wset,2
        ncolors = !D.Table_Size
        device,DECOMPOSED=0
        LoadCT,39 ;0 - B&W ;39 - Color
        rmin_y = min(red_xx1[*]) & bmin_y = min(blue_xx1[*])
        rmax_y = max(red_xx2[*]) & bmax_y = max(blue_xx2[*])

;	if keyword_set(scale_minmax) then begin
;               bscale_min =  min(blue_sciimg[blue_slitpos])
;               bscale_max =  max(blue_sciimg[blue_slitpos])
;               rscale_min =  min(red_sciimg[red_slitpos])
;		rscale_max =  max(red_sciimg[red_slitpos])
;        endif else begin
                bscale_min = -1
                bscale_max = 50*abs(median(blue_sciimg[bmin_y:bmax_y,*]))
                rscale_min = -1
                rscale_max = 50*abs(median(red_sciimg[rmin_y:rmax_y,*]))
;       endelse
;
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

	red_plot_position = [0.05,0.58,0.99,0.95]
	blue_plot_position = [0.05,0.08,0.99,0.45]

        blue_slitim = transpose(BytScl(blue_sciimg[bmin_y:bmax_y,*],$
		min=bscale_min,max=bscale_max,/nan))
	red_slitim = transpose(BytScl(red_sciimg[rmin_y:rmax_y,*],$
		min=rscale_min,max=rscale_max,/nan))
        TVImage,blue_slitim,Position=blue_plot_position,/erase
        info = 'BLUE'; - aperture '+strcompress(ap)
	xyouts,300,bmax_y-bmin_y-20,info,CharThick=2,Size=2.3
	TVImage,red_slitim,Position=red_plot_position
	info = 'RED'; - aperture '+strcompress(ap)
	xyouts,300,rmax_y-rmin_y-20,info,CharThick=2,Size=2.3
        xpix = indgen(8192)

        plot,xpix,blue_xx1[*]-bmin_y,/noerase,Position=blue_plot_position,$
                xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax_y-bmin_y-1], $
                linestyle=2,thick=2,xtitle='Pixels',charsize=1.2,title=title[0]
        plot,xpix,blue_xx2[*]-bmin_y,/noerase,Position=blue_plot_position,$
                xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax_y-bmin_y-1], $
                linestyle=2,thick=2,xtitle='Pixels',charsize=1.2
        plot,xpix,red_xx1[*]-rmin_y,/noerase,Position=red_plot_position,$
                xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax_y-rmin_y-1], $
                linestyle=2,thick=2,xtitle='Pixels',charsize=1.2,title=title[1]
        plot,xpix,red_xx2[*]-rmin_y,/noerase,Position=red_plot_position,$
                xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,rmax_y-rmin_y-1], $
                linestyle=2,thick=2,xtitle='Pixels',charsize=1.2


        if keyword_set(blue_center_trace) then $
                plot,xpix,blue_center_trace[*]-bmin_y,/noerase, $
                        Position=blue_plot_position,xstyle=5,xrange=[0,8191],$
                        ystyle=5,yrange=[0,bmax_y-bmin_y-1],Color=cgColor('red'),$
                        linestyle=0,thick=2,xtitle='Pixels',charsize=1.2
	if keyword_set(red_center_trace) then $
		plot,xpix,red_center_trace[*]-rmin_y,/noerase, $
			Position=red_plot_position,xstyle=5,xrange=[0,8191],$
			ystyle=5,yrange=[0,rmax_y-rmin_y-1],Color=cgColor('red'),$
			linestyle=0,thick=2,xtitle='Pixels',charsize=1.2

        if keyword_set(lower) then begin
                plot,xpix,blue_center_trace[*]-bmin_y-lower/blue_platescale,/noerase, $
                        Position=blue_plot_position,xstyle=5,xrange=[0,8191],$
                        ystyle=5,yrange=[0,bmax_y-bmin_y-1],Color=cgColor('red'),$
                        linestyle=2,thick=2,xtitle='Pixels',charsize=1.2
		plot,xpix,red_center_trace[*]-rmin_y-lower/red_platescale,/noerase, $
			Position=red_plot_position,xstyle=5,xrange=[0,8191],$
			ystyle=5,yrange=[0,rmax_y-rmin_y-1],Color=cgColor('red'),$
			linestyle=2,thick=2,xtitle='Pixels',charsize=1.2
	endif
        if keyword_set(upper) then begin
                plot,xpix,blue_center_trace[*]-bmin_y+upper/blue_platescale,/noerase, $
                        Position=blue_plot_position,xstyle=5,xrange=[0,8191],$
                        ystyle=5,yrange=[0,bmax_y-bmin_y-1],Color=cgColor('red'),$
                        linestyle=2,thick=2,xtitle='Pixels',charsize=1.2
		plot,xpix,red_center_trace[*]-rmin_y+upper/red_platescale,/noerase, $
			Position=red_plot_position,xstyle=5,xrange=[0,8191],$
			ystyle=5,yrange=[0,rmax_y-rmin_y-1],Color=cgColor('red'),$
			linestyle=2,thick=2,xtitle='Pixels',charsize=1.2
	endif
return
end


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
	trim1,trim2,lower=lower,upper=upper,$
        em_line=em_line,centersum=centersum,subcont=subcont

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
	if keyword_set(em_line) then result=where(stack eq max(stack)) else result=Value_Locate(stack,mean(stack))
        result_x = goodwav[result]

        ;find which part of the slit is at the median... mininimum pulled out the slit edge
	tmp = ap_im[result_x-(centersum/2-1):result_x+(centersum/2),*]
        stack2 = smooth(total(tmp>0,1),2)
	if keyword_set(subcont) then begin
		tmp_cont = ap_im[result_x-4+subcont:result_x+5+subcont,*]
		;plot,stack2,/ylog,yrange=[1,1000]
		stack2_cont = smooth(total(tmp_cont>0,1),2)
		stack2 = stack2 - stack2_cont
		;oplot,stack2_cont,color=cgcolor('red')
		;oplot,stack2,color=cgcolor('green')
		;stop
	endif
	result_y=where(stack2[10:n_elements(stack2)-10] eq max(stack2[10:n_elements(stack2)-10]))
	result_y = result_y[0]+10
return
end

;=========================================================================
; extractionplot - plot the variosu cuts and extractions
;			for the selected aperture
;=========================================================================
pro extractionplot,bluew,redw,$
	blueg,redg, $
	bluest,redst,$
	bluex,redx,$
	bluest2,redst2,$
	bluey,redy, $
	bluesp,redsp,$
	lower,upper,$
	subcont=subcont,$
	edge_mask = edge_mask, $
	centersum=centersum, $
	logprofile = logprofile, $
	title = title

common exclude1, blank_lowb, blank_lowr, blank_hib, blank_hir, blue_platescale, red_platescale

;Plot up the centering and extractions
xb0 = bluew[blueg[0]]-20    &  xr0 = redw[redg[0]]-20
xb1 = max(bluew[blueg])+20  &  xr1 = max(redw[redg])+20

;bottom pannel
stacksz = size(bluest2,/dim)
bluest2 = bluest2 +10.D

if keyword_set(logprofile) then minval = min(bluest2/max(bluest2))/5>0.0001 else minval = 0.
plot,indgen(stacksz),bluest2/max(bluest2),XRANGE=[0,stacksz],$
	position=[0.14,0.09,0.50,0.48],ytitle='Slit Profile',$
        yrange=[minval,1],xstyle=9, $
        charsize=1.5,ystyle=1,ylog=logprofile
line_x = [bluey,bluey]
maxval = 1;10^(!Y.CRANGE[1])
line_y = [minval,maxval]
line_x1 = [bluey+upper/0.12,bluey+upper/0.12]
line_x2 = [bluey-lower/0.12,bluey-lower/0.12]
niter = size(blank_lowb,/dim)
niter = niter[0]-1
for i=0,niter do polyfill,[blank_lowb[i]>0,blank_lowb[i]>0,blank_hib[i]<stacksz,blank_hib[i]<stacksz], $
	[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data
plot,indgen(stacksz),bluest2/max(bluest2),XRANGE=[0,stacksz],$
        position=[0.14,0.09,0.50,0.48],charsize=1.5,ystyle=1,$
        /noerase,yrange=[minval,1],xstyle=9,xtitle='Pixels',ylog=logprofile
oplot,line_x,line_y,Color=cgColor('dodger blue')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')
axis,xaxis=1,XRANGE=[-(line_x[0]*blue_platescale),(stacksz-line_x[0])*blue_platescale],xtitle='ARCSEC',charsize=1.3,xstyle=1

stacksz = size(redst2,/dim)
redst2 = redst2+10.D
if keyword_set(logprofile) then minval = min(redst2/max(redst2))/5>0.0001 else minval = 0.
plot,indgen(stacksz),redst2/max(redst2),XRANGE=[0,stacksz],$
	position=[0.59,0.09,0.95,0.48],$
        /noerase,yrange=[minval,1],xstyle=9, $
        charsize=1.5,ystyle=1,ylog=logprofile
line_x = [redy,redy]
maxval = 1;10^(!Y.CRANGE[1])
line_y = [minval,maxval]
line_x1 = [redy+upper/0.123,redy+upper/0.123]
line_x2 = [redy-lower/0.123,redy-lower/0.123]
for i=0,niter do polyfill,[blank_lowr[i],blank_lowr[i],blank_hir[i],blank_hir[i]],$
	[minval,maxval,maxval,minval],color=cgColor('firebrick')
oplot,line_x,line_y,Color=cgColor('dodger blue')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')
plot,indgen(stacksz),redst2/max(redst2),XRANGE=[0,stacksz],$
        position=[0.59,0.09,0.95,0.48],charsize=1.5,ystyle=1,$
        /noerase,yrange=[minval,1],xstyle=9,xtitle='Pixels',ylog=logprofile
axis,xaxis=1,XRANGE=[-(line_x[0]*red_platescale),(stacksz-line_x[0])*red_platescale],xtitle='ARCSEC',charsize=1.3,xstyle=1

;top panel
plot,bluew,bluesp,position=[0.14,0.62,0.50,0.97],title=title[0],$
        /noerase,XRANGE=[xb0-100,xb1+50],ytitle='Extracted Spectrum', $
        charsize=1.2,ystyle=1,xtitle='Wavelength',/ylog
oplot,[bluew[bluex],bluew[bluex]],[10^!Y.CRANGE[0],10^!Y.CRANGE[1]],color=cgColor('firebrick'),$
        thick=2
oplot,bluew,bluesp
plot,redw,redsp,position=[0.59,0.62,0.95,0.97],title=title[1],$
	/noerase,XRANGE=[xr0-100,xr1+100],$;ytitle='Av. Sky spectrum', $
        charsize=1.2,ystyle=1,xtitle='Wavelength',/ylog
oplot,[redw[redx],redw[redx]],[10^!Y.CRANGE[0],10^!Y.CRANGE[1]],color=cgColor('firebrick'),$
        thick=2
oplot,redw,redsp

return
end

;=========================================================================
; mods_extract_dualchan - extract and calibrate 1D dual channel MODS spectra
;			from 2D spectra
;=========================================================================
pro mods_extract_dualchan,blue_scifile,red_scifile, 	$
	red_cal = red_cal, 				$
        blue_cal=blue_cal, 				$
        z = z, 						$	
	apertures=apertures, 				$
        outname=outname, 				$
        wave_blue=wave_blue,wave_red=wave_red, 		$
	centerLine_blue=centerLine_blue, 		$
	centerLine_red=centerLine_red, 			$
        trim_rt = trim_rt, trim_rb = trim_rb, 		$
	trim_bb = trim_bb, trim_bt = trim_bt,		$
	centersum_red = centersum_red,                  $
	centersum_blue = centersum_blue,                $
	noPlotAp=noPlotAp, 				$
        red_slits = red_slits, blue_slits = blue_slits, $
        force_blue=force_blue, force_red=force_red, 	$
	ill_red=illred, ill_blue=ill_blue, 		$
        illumination_corr = illumination_corr, 		$
	clobber = clobber, 				$
	scale_minmax=scale_minmax, 			$
	acquisition_centering_correction = acquisition_centering_correction, $
	man_scale = man_scale, 				$
	subcont = subcont, 				$
	skyslit = skyslit,				$
	verbose=verbose, 				$
        em_line=em_line,              			$
	logprofile = logprofile,		        $
        extraction_reference = extraction_reference

if not keyword_set(verbose) then !Quiet=1

if not keyword_set(red_cal) and keyword_set(blue_cal) then begin
        print,'To produce a calibrated sectrum two response functions are needed.'
        print,'no calibration file(s) provided.'
        stop
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

common exclude1, blank_lowb, blank_lowr, blank_hib, blank_hir, blue_platescale, red_platescale
blank_lowb = [0]
blank_lowr = [0]
blank_hib = [0]
blank_hir = [0]

; Science image
verify_red = file_test(red_scifile)
verify_blue = file_test(blue_scifile)
if verify_red then begin
        red_sciimg = mrdfits(red_scifile,0,red_hdr,/silent)
        red_skyimg = mrdfits(red_scifile,3,red_hdrsky,/silent)
        red_ivarimg = mrdfits(red_scifile,1,rheadvar,/silent)
        red_varimg = sqrt(2/red_ivarimg)
        if keyword_set(skyslit) then red_sciimg = mrdfits(red_scifile,5,tmpred_hdr,/silent)
 endif else begin
                print,'------------------------------------------------------'
                print,'You need a valid RED science image'
                print,red_scifile,' does not exist. Check your Science folder'
                stop
 endelse
if verify_blue then begin
        blue_sciimg = mrdfits(blue_scifile,0,blue_hdr,/silent)
        blue_skyimg = mrdfits(blue_scifile,3,blue_hdrsky,/silent)
        blue_ivarimg = mrdfits(blue_scifile,1,bheadvar,/silent)
        blue_varimg = sqrt(2/blue_ivarimg)
        if keyword_set(skyslit) then blue_sciimg = mrdfits(blue_scifile,5,tmpblue_hdr,/silent)
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
if strmid(maskr,0,2) eq 'LS' then mask = 'ls' else mask = maskr

if n_elements(centerLine_red) EQ '' then begin
        rlow_cent = 6630*z
        rhi_cent = 6800*z
endif else begin
        rlow_cent = centerLine_red[0]*z
        rhi_cent  = centerLine_red[1]*z
endelse

if n_elements(centerLine_blue) EQ '' then begin
	blow_cent = 4890*z
	bhi_cent =  5070*z
endif else begin
        blow_cent = centerLine_blue[0]*z
        bhi_cent =  centerLine_blue[1]*z
endelse

;subtract continuum
if keyword_set(subcont) then begin
	Print,'Continuum offset (pixels):'
        read,': ',subcont
endif else subcont=0

;check that these should be extrated together
obj = strcompress(sxpar(red_hdr[*,0], 'object'), /rem)
obj2 = strcompress(sxpar(blue_hdr[*,0], 'object'), /rem)
if obj ne obj2 then begin
	print,'These appear to be RED and BLUE spectra of different objects'
	stop
endif

;outfile
if n_elements(outname) EQ '' then begin
        name_len=STRLEN(red_scifile)
	redout = STRMID(red_scifile, 0, name_len-8) + 'r1d.fits'
	redoutcal = STRMID(red_scifile, 0, name_len-8) + 'x1d.fits'
        name_len=STRLEN(blue_scifile)
        blueout = STRMID(blue_scifile, 0, name_len-8) + 'r1d.fits'
	blueoutcal = STRMID(blue_scifile, 0, name_len-8) + 'x1d.fits' 
endif else begin
	stub = STRSPLIT(red_scifile,'-',/extract)
        redout = stub[0] + '-' + red_channel_code+'_' + outname + '_r1d.fits'
        redoutcal = stub[0] + '-' + red_channel_code+'_' + outname + '_x1d.fits'
        blueout = stub[0] +'-' +blue_channel_code+'_' + outname + '_r1d.fits' 
        blueoutcal = stub[0] +'-' +blue_channel_code+'_' + outname + '_x1d.fits'
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
if n_elements(wave_blue) EQ '' then wave_blue = 'wave-comp_'+blue_channel_code+'_'+mask+'.fits'
if n_elements(wave_red) EQ '' then wave_red = 'wave-comp_'+red_channel_code+'_'+mask+'.fits'
verify_red = file_test(wave_red)
verify_blue = file_test(wave_blue)
if verify_red then red_waveimg = mrdfits(wave_red,0,red_whdr,/silent) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid RED wave-image'
                print,wave_red,' does not exist.'
                print,"Use: wave_red='FILENAME'"
                stop
		red_waveimg = mrdfits(wave_red,0,red_whdr,/silent)
        endelse
if verify_blue then blue_waveimg = mrdfits(wave_blue,0,blue_whdr,/silent) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid BLUE wave-image'
                print,wave_blue,' does not exist.'
                print,"Use: wave_blue='FILENAME'"
                stop
		blue_waveimg = mrdfits(wave_blue,0,blue_whdr,/silent)
        endelse

; Calibration images
verify_red = file_test(red_cal)
verify_blue = file_test(blue_cal)
if verify_red then calR = mrdfits(red_cal,1,chd,/silent) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid RED calibration-image'
                print,red_cal,' does not exist.'
                print,"Use: red_cal='FILENAME'"
                stop
		calR = mrdfits(red_cal,1,chd,/silent)
        endelse
if verify_blue then calB = mrdfits(blue_cal,1,chd,/silent) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid BLUE calibration-image'
                print,blue_cal,' does not exist.'
                print,"Use: blue_cal='FILENAME'"
                stop
		calB = mrdfits(blue_cal,1,chd,/silent)
        endelse

name_len=STRLEN(red_cal)
red_cal1d = STRMID(red_cal, 0, name_len-8) + 'x1d.fits'
name_len=STRLEN(blue_cal)
blue_cal1d = STRMID(blue_cal, 0, name_len-8) + 'x1d.fits'
verify_red = file_test(red_cal1d)
verify_blue = file_test(blue_cal1d)
if verify_red then calR_1d = mrdfits(red_cal1d,4,chd1d,/silent) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid RED standard'
                print,red_cal1d,' does not exist.'
                print,"Use: red_cal='FILENAME'"
                stop
                calR_1d = mrdfits(red_cal,4,chd,/silent)
        endelse
if verify_blue then calB_1d = mrdfits(blue_cal1d,4,chd1d,/silent) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid BLUE standard'
                print,blue_cal1d,' does not exist.'
                print,"Use: blue_cal='FILENAME'"
                stop
                calB_1d = mrdfits(blue_cal,4,chd,/silent)
        endelse
;altlen = STRLEN(red_cal)
;alt_red_cal = STRMID(red_cal, 0, altlen-5) + '2.fits'
;alt_calR = mrdfits(alt_red_cal,1,chd,/silent)
;altlen = STRLEN(blue_cal)
;alt_blue_cal = STRMID(blue_cal, 0, altlen-5) + '2.fits'
;alt_calB = mrdfits(alt_blue_cal,1,chd,/silent)


; Flat image
line = 'ls Proc/flat*' + blue_channel_code + '*fits*'
spawn,line,blue_flat_filenames
line = 'ls Proc/flat*' + red_channel_code + '*fits*'
spawn,line,red_flat_filenames

nfile_red = n_elements(red_flat_filenames)   &   nfile_blue = n_elements(blue_flat_filenames)

if nfile_red gt 0 then begin
        if nfile_red eq 1 then red_flatfile = red_flat_filenames
        if nfile_red gt 1 then begin
        	print,'Multpile flat files found.  Please advise.'
		print,"red_flatfile = '<filename>'"
		print,red_flat_filenames
		stop
        endif
endif else begin
        print,'No flat files found.  Please be careful with fluxing.'
        stop
endelse
flat_red_im = mrdfits(red_flatfile,0,fr_hd,/silent)
flat_red_im = transpose(flat_red_im)
if nfile_blue gt 0 then begin
        if nfile_blue eq 1 then blue_flatfile = blue_flat_filenames
        if nfile_blue gt 1 then begin
                print,'Multpile flat files found.  Please advise.'
		print,"blue_flatfile = '<filename>'"
		print,blue_flat_filenames
                stop
        endif
endif else begin
        print,'No flat files found.  Please be careful with fluxing.'
        stop
endelse
flat_blue_im = mrdfits(blue_flatfile,0,fb_hd,/silent)
flat_blue_im = transpose(flat_blue_im)

;if keyword_set(flat_blue) then begin
;	flat_blue_im = mrdfits(flat_blue,0,fb_hd,/silent)
;	flat_blue_im = transpose(flat_blue_im)

;establish a uniform grid
red_hdrcal = red_hdr
blue_hdrcal = blue_hdr

sxaddpar,red_hdrcal,'HISTORY','mods_extract_dualchan: beta version'
sxaddpar,red_hdrcal,'CUNIT1','Angstrom'
sxaddpar,red_hdrcal,'CTYPE1','Linear'
sxaddpar,red_hdrcal,'CRPIX1',1
sxaddpar,red_hdrcal,'CRVAL1',5500
sxaddpar,red_hdrcal,'CDELT1',0.5
sxaddpar,red_hdrcal,'BUNIT','erg/cm^2/s/A'

sxaddpar,blue_hdrcal,'HISTORY','mods_extract_dualchan: beta version'
sxaddpar,blue_hdrcal,'CUNIT1','Angstrom'
sxaddpar,blue_hdrcal,'CTYPE1','Linear'
sxaddpar,blue_hdrcal,'CRPIX1',1
sxaddpar,blue_hdrcal,'CRVAL1',3000
sxaddpar,blue_hdrcal,'CDELT1',0.5
sxaddpar,blue_hdrcal,'BUNIT','erg/cm^2/s/A'

new_red_wave = indgen(9001)*0.5 + 5500.D
new_blue_wave= indgen(5601)*0.5 + 3000.D

red_sensfunc = bspline_valu(new_red_wave,calR)
blue_sensfunc = bspline_valu(new_blue_wave,calB)
;red_altSensfunc = bspline_valu(new_red_wave,alt_calR)
;blue_altSensfunc = bspline_valu(new_blue_wave,alt_calB)

; Extinction file
extfile = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/mods_extiction.dat'
readcol,extfile,ext_wave,ext_val,/silent
ext_coef = bspline_iterfit(ext_wave,ext_val,bkspace=3,nord=3)
extR = bspline_valu(new_red_wave,ext_coef)
extB = bspline_valu(new_blue_wave,ext_coef)

;open a window for plotting
title = 'Extraction and Slit Profile'; for ' + blue_scifile + ' and ' + red_scifile
if not WindowAvailable(1) then window,1,retain=2,title=title
title = '2D Spectrum'; (' + blue_scifile + ' and ' + red_scifile + ')'
if (not keyword_set(noPlotAp)) and (not WindowAvailable(2)) then window,2,retain=2,title=title

dims=size(red_sciimg,/dim)
nx=dims[0]
ny=dims[1]
red_instrument = strcompress(sxpar(red_hdr[*, 0], 'INSTRUME'), /rem)
blue_instrument = strcompress(sxpar(blue_hdr[*, 0], 'INSTRUME'), /rem)
mms = strcompress(sxpar(red_hdr[*,0], 'MASKINFO')+'.mms', /rem)

;create the slit structures
  if n_elements(red_slits) eq '' then red_slits = 'slits-flat_'+red_channel_code+'_'+mask+'.fits'
  if n_elements(blue_slits) eq '' then blue_slits = 'slits-flat_'+blue_channel_code+'_'+mask+'.fits'
  verify_red = file_test(red_slits)
  verify_blue = file_test(blue_slits)
  if verify_red then begin
	red_apmask = mrdfits(red_slits,0,r_mask_head,/silent)
	red_slits = mrdfits(red_slits,1,r_head,/silent)
  endif else begin
        print,'------------------------------------------------------'
        print,'You need a valid RED slit-image'
        print,red_slits,' does not exist.'
        print,"Use: red_slits='FILENAME'"
        stop
	red_apmask = mrdfits(red_slits,0,r_mask_head,/silent)
        red_slits = mrdfits(red_slits,1,r_head,/silent)
   endelse
   if verify_blue then begin
	blue_apmask = mrdfits(blue_slits,0,b_mask_head,/silent)
	blue_slits = mrdfits(blue_slits,1,b_head,/silent)
   endif else begin
        print,'------------------------------------------------------'
        print,'You need a valid BLUE slit-image'
        print,blue_slits,' does not exist.'
        print,"Use: blue_slits='FILENAME'"
        stop
        blue_apmask = mrdfits(blue_slits,0,b_mask_head,/silent)
        blue_slits = mrdfits(blue_slits,1,b_head,/silent)
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
;   red_illspec = fltarr(ny_t)+1
   red_flatspec = fltarr(ny_t)+1

   blue_xx1 = blue_slits[0].xx1
   blue_xx2 = blue_slits[1].xx2

   blue_objtrace = blue_xx1*0
   blue_objspec = fltarr(ny_t)+1
   blue_objwave = fltarr(ny_t)+1
   blue_varspec = fltarr(ny_t)+1
   blue_skyspec = fltarr(ny_t)+1
;   blue_illspec = fltarr(ny_t)+1
   blue_flatspec = fltarr(ny_t)+1

   blue_objspec_cal = blue_sensfunc
   red_objspec_cal = red_sensfunc
   blue_varspec_cal = fltarr(size(new_blue_wave,/dim))
   red_varspec_cal = fltarr(size(new_red_wave,/dim))
   blue_skyspec_cal = fltarr(size(new_blue_wave,/dim))
   red_skyspec_cal = fltarr(size(new_red_wave,/dim))
;   blue_illspec_cal = fltarr(5601)
;   red_illspec_cal = fltarr(9001)
   blue_flatspec_cal = fltarr(size(new_blue_wave,/dim))
   red_flatspec_cal = fltarr(size(new_red_wave,/dim))


;extraction struct
red_extract = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D], $
	ra_deg: double(red_slits[0].tarra2), dec_deg: double(red_slits[0].tardec2), $
	ra_sex: red_slits[0].tarra, dec_sex: red_slits[0].tardec}
blue_extract = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D], $
	ra_deg: double(blue_slits[0].tarra2), dec_deg: double(blue_slits[0].tardec2), $
	ra_sex: blue_slits[0].tarra, dec_sex: blue_slits[0].tardec}

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                      EXTRACTION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(apertures) EQ '' then apertures = indgen(nslit)+1
ndo = size(apertures,/dim)
ndo=ndo[0]
object=1

for ii=0,ndo-1 do begin
	LeObject: 
	ap=apertures[ii]

        ;get the traces
	red_xx1_ap = red_xx1[*,ap-1]+trim_rb
        red_xx2_ap = red_xx2[*,ap-1]-trim_rt
        blue_xx1_ap = blue_xx1[*,ap-1]+trim_bb
        blue_xx2_ap = blue_xx2[*,ap-1]-trim_bt

        lower = 2.
        upper = 2.
        blank_lowb = [0]
        blank_lowr = [0]
        blank_hib = [0]
        blank_hir = [0]
        n_blank = 1

	;currently, returning the center of the slit for ra and dec 
	red_extract_obj = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D], $
		ra_deg: red_slits[0].apra2[ap-1], dec_deg: red_slits[0].apdec2[ap-1], $
		ra_sex: red_slits[0].apra[ap-1], dec_sex: red_slits[0].apdec[ap-1]}
	blue_extract_obj = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D], $
		ra_deg: blue_slits[0].apra2[ap-1], dec_deg: blue_slits[0].apdec2[ap-1], $
		ra_sex: blue_slits[0].apra[ap-1], dec_sex: blue_slits[0].apdec[ap-1]}
	

	center,red_sciimg,red_waveimg,ap,rlow_cent,rhi_cent,$
                red_result_x,red_result_y,red_good,red_wav1d,$
                red_stack,red_stack2,red_xx1,red_xx2,trim_rb,trim_rt,$
		subcont=subcont,lower=lower,upper=upper,$
                centersum=centersum_red,em_line=em_line
        center,blue_sciimg,blue_waveimg,ap,blow_cent,bhi_cent,$
		blue_result_x,blue_result_y,blue_good,blue_wav1d,$
		blue_stack,blue_stack2,blue_xx1,blue_xx2,trim_bb,trim_bt,$
		subcont=subcont,lower=lower,upper=upper,$
		centersum=centersum_blue,em_line=em_line

;       modsadcalc,blue_scifile,blue_wav1d[blue_result_x],blue_wav1d,bxdelt,bydelt
;       modsadcalc,red_scifile,red_wav1d[red_result_x],red_wav1d,rxdelt,rydelt
	modsadcalc,blue_scifile,6300.,blue_wav1d,bxdelt,bydelt
	modsadcalc,red_scifile,6300.,red_wav1d,rxdelt,rydelt

	LeSigh: print,'Extracting Spectrum from aperture #',ap 
        red_spec1d = fltarr(ny)  &  blue_spec1d = fltarr(ny)
        red_var1d = fltarr(ny)   &  blue_var1d = fltarr(ny)
	red_sky1d = fltarr(ny)   &  blue_sky1d = fltarr(ny)
        red_spec1db = fltarr(ny) &  blue_spec1db = fltarr(ny)
        red_flat1d = fltarr(ny)  &  blue_flat1d = fltarr(ny)

        for iii=0,ny-1 do red_objtrace[iii,ap-1] = red_xx1[iii,ap-1] + trim_rb + red_result_y
        for iii=0,ny-1 do blue_objtrace[iii,ap-1] = blue_xx1[iii,ap-1] + trim_bb + blue_result_y
	
        ; added in for method demonstration
        blue_orig = blue_objtrace[*,ap-1]

	if (red_channel_code eq 'm1r') or (red_channel_code eq 'm2r') then red_platescale = 0.123
	if (blue_channel_code eq 'm1b') or (blue_channel_code eq 'm2b') then blue_platescale = 0.120	
	blue_objtrace[*,ap-1] = blue_objtrace[*,ap-1] - bydelt/blue_platescale
        red_objtrace[*,ap-1] = red_objtrace[*,ap-1] - rydelt/red_platescale

	for iii=0,ny-1 do begin
		bot = red_objtrace[iii,ap-1]-lower/red_platescale
                top = red_objtrace[iii,ap-1]+upper/red_platescale
		botint = total(bot,/integer)
		topint = total(top,/integer)
		red_spec1d[iii] = total(red_sciimg[bot:top,iii]) $
			+(1-(bot-botint))*total(red_sciimg[bot,iii]) $
			+(top-topint)*total(red_sciimg[top+1,iii])
		red_wav1d[iii] = median(red_waveimg[bot:top,iii]); + rxdelt[iii]*0.123
		red_var1d[iii] = sqrt(total(red_varimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(red_varimg[bot,iii]) $
                        +(top-topint)*total(red_varimg[top+1,iii]))
		red_sky1d[iii] = sqrt(total(red_skyimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(red_skyimg[bot,iii]) $
                        +(top-topint)*total(red_skyimg[top+1,iii]))
		red_flat1d[iii] = mean(flat_red_im[bot:top,iii])

                bot = blue_objtrace[iii,ap-1]-lower/blue_platescale
                top = blue_objtrace[iii,ap-1]+upper/blue_platescale
                botint = total(bot,/integer)
                topint = total(top,/integer)
                blue_spec1d[iii] = total(blue_sciimg[bot:top,iii])$
                        +(1-(bot-botint))*total(blue_sciimg[bot,iii]) $
                        +(top-topint)*total(blue_sciimg[top+1,iii])
		blue_wav1d[iii] = median(blue_waveimg[bot:top,iii]); + bxdelt[iii]*0.12
		blue_var1d[iii] = sqrt(total(blue_varimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(blue_varimg[bot,iii]) $
                        +(top-topint)*total(blue_varimg[top+1,iii]))
		blue_sky1d[iii] = sqrt(total(blue_skyimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(blue_skyimg[bot,iii]) $
                        +(top-topint)*total(blue_skyimg[top+1,iii]))
		blue_flat1d[iii] = mean(flat_blue_im[bot:top,iii])
	endfor

	; plot slit and traces
	if (not keyword_set(noPlotAp)) then plot_2d_spec,red_sciimg,blue_sciimg,$
		red_xx1_ap,red_xx2_ap,blue_xx1_ap,blue_xx2_ap,$
		red_center_trace=red_objtrace[*,ap-1],blue_center_trace=blue_objtrace[*,ap-1],$
		lower=lower,upper=upper,red_platescale=red_platescale,blue_platescale=blue_platescale,$
		man_scale = man_scale,title=[blue_scifile,red_scifile]

	;Plot up the centering and extractions
	wset,1
	extractionplot,blue_wav1d,red_wav1d,blue_good,red_good,$
	       blue_stack,red_stack,blue_result_x,red_result_x,$
	       blue_stack2,red_stack2,blue_result_y,red_result_y,$
	       blue_spec1d,red_spec1d,lower,upper,subcont=subcont,$
	       centersum=centersum,logprofile = logprofile,title=[blue_scifile,red_scifile]

	Ok=0
	while Ok ne 1 do begin
		print,'RED Slit profile taken at: ', blue_wav1d[blue_result_x[0]],' Angstroms'
		print,'BLUE Slit profile taken at: ', red_wav1d[red_result_x[0]],' Angstroms'
	        print,'Blue Peak found at pixel: ',blue_result_y
	        print,'Red Peak found at pixel: ',red_result_y
	        print,'Current extraction bounds: lower = ',lower,'  upper = ',upper
		print,'      10 pixels is ~1.2"'
	        print,'Recenter extraction? (y/n/bye)'
	        redo = ''
		redo = strlowcase(redo)
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
			redo = strlowcase(redo)
	                if redo eq 'y' then begin
	                        Print,'New Lower (in arcseconds):'
	                        readnum,lower
	                        Print,'New Upper (in arcseconds):'
	                        readnum,upper
				lower = abs(lower)
				upper = abs(upper)
	                        GOTO,LeSigh
	                endif else if redo eq 'n' then begin
                                blue_extract_obj.object = object
	                        blue_extract_obj.slit = ap
	                        blue_extract_obj.center = blue_result_y
	                        blue_extract_obj.lower = lower
                 		blue_extract_obj.upper = upper

                                pos_shift = blue_result_y * blue_platescale
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

                                pos_shift = red_result_y * red_platescale
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
	if keyword_set(force_blue) then $
           mods_flexurecorr,blue_wav1d,blue_spec1d,blue_sky1d,blue_var1d,'blue',$
                z=z,/plot_prog,force=force_blue,delta=bdelta,slope=bslope $
	   else mods_flexurecorr,blue_wav1d,blue_spec1d,blue_sky1d,blue_var1d,'blue',$
                z=z,/plot_prog,delta=bdelta,slope=bslope,/noslope

	if keyword_set(force_red) then $
	   mods_flexurecorr,red_wav1d,red_spec1d,red_sky1d,red_var1d,'red',$
		z=z,/plot_prog,force=force_red,delta=rdelta,slope=rslope $
	   else mods_flexurecorr,red_wav1d,red_spec1d,red_sky1d,red_var1d,'red',$
		z=z,/plot_prog,delta=rdelta,slope=rslope

	blue_extract[object-1].delta = bdelta
	red_extract[object-1].delta = rdelta
        blue_extract[object-1].slope = bslope
        red_extract[object-1].slope = rslope

        if not keyword_set(acquisition_centering_correction) then acquisition_centering_correction = 0.D
	blue_wav1d = blue_wav1d+acquisition_centering_correction
	red_wav1d = red_wav1d-acquisition_centering_correction	

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

	x_specrebin,blue_wav1d,blue_flat1d,new_blue_wave,blue_corFlat,/FLAMBDA
	x_specrebin,red_wav1d,red_flat1d,new_red_wave,red_corFlat,/FLAMBDA

        ;The extinction correction is given by the factor
        extBcor = 10. ^ (0.4 * airmassB * extB[*])
        extRcor = 10. ^ (0.4 * airmassR * extR[*])

        ;Flux the spectrum including the LBT extinction curve
        calfluxB = 1.0d-17*(corflux_blue*extBcor)/(time) * 10.^(blue_sensfunc/2.5)
        calfluxR = 1.0d-17*(corflux_red*extRcor)/(time) * 10.^(red_sensfunc/2.5)
	
        calvarB = 1.0d-17*(corvar_blue*extBcor)/(time) * 10.^(blue_sensfunc/2.5)
	calvarR = 1.0d-17*(corvar_red*extRcor)/(time) * 10.^(red_sensfunc/2.5)

	calskyB = 1.0d-17*(corsky_blue*extBcor)/(time) * 10.^(blue_sensfunc/2.5)
	calskyR = 1.0d-17*(corsky_red*extRcor)/(time) * 10.^(red_sensfunc/2.5)

;	corflux2_blue = corflux_blue/(blue_corFlat/mean(blue_corFlat))
;	corflux2_red  = corflux_red/(red_corFlat/mean(red_corFlat))
;        calfluxB2 = 1.0d-17*(corflux2_blue*extBcor)/(time) * 10.^(blue_altSensfunc/2.5)
;        calfluxR2 = 1.0d-17*(corflux2_red*extRcor)/(time) * 10.^(red_altSensfunc/2.5)

fit_coef = bspline_iterfit(new_blue_wave,((blue_corFlat/calB_1d)/mean(blue_corFlat/calB_1d)),bkspace=250,nord=3)
smoothfitB = bspline_valu(new_blue_wave,fit_coef)
fit_coef = bspline_iterfit(new_red_wave,((red_corFlat/calR_1d)/mean(red_corFlat/calR_1d)),bkspace=250,nord=3)
smoothfitR = bspline_valu(new_red_wave,fit_coef)

;plot,new_blue_wave,calfluxB,/ylog,/ystyle,/xstyle,xrange=[3100,5700],yrange=[3e-16,3e-15]
;oplot,new_blue_wave,calfluxB2,color=cgcolor('dodger blue')
;oplot,new_blue_wave,((blue_corFlat/calB_1d)/mean(blue_corFlat/calB_1d))*1e-15,color=cgcolor('orange')
;oplot,new_blue_wave,smoothfitB*1e-15,color=cgcolor('green')
;oplot,new_blue_wave,((blue_corFlat/calB_1d)/mean(blue_corFlat/calB_1d))/smoothfitB*1e-15,color=cgcolor('dark green')
;oplot,new_blue_wave,calfluxB/(((blue_corFlat/calB_1d)/mean(blue_corFlat/calB_1d))/smoothfitB),color=cgcolor('red')

	calfluxB = calfluxB/(((blue_corFlat/calB_1d)/mean(blue_corFlat/calB_1d))/smoothfitB)
	calfluxR = calfluxR/(((red_corFlat/calR_1d)/mean(red_corFlat/calR_1d))/smoothfitR)

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

        red_objspec_cal = [[red_objspec_cal],[calfluxR]]
	red_varspec_cal = [[red_varspec_cal],[calvarR]]
	red_skyspec_cal = [[red_skyspec_cal],[calskyR]]

        ; Save the output
        mwrfits,red_objspec,redout,red_hdr,/create
	mwrfits,red_varspec,redout,/silent
	mwrfits,red_skyspec,redout,/silent
	mwrfits,red_objwave,redout,/silent
        mwrfits,red_extract,redout,/silent

        mwrfits,red_objspec_cal,redoutcal,red_hdrcal,/create
	mwrfits,red_varspec_cal,redoutcal,/silent
	mwrfits,red_skyspec_cal,redoutcal,/silent
	mwrfits,red_extract,redoutcal,/silent
;	mwrfits,red_illspec_cal,redoutcal,/silent

        mwrfits,blue_objspec,blueout,blue_hdr,/create
	mwrfits,blue_varspec,blueout,/silent
	mwrfits,blue_skyspec,blueout,/silent
	mwrfits,blue_objwave,blueout,/silent
        mwrfits,blue_extract,blueout,/silent

        mwrfits,blue_objspec_cal,blueoutcal,blue_hdrcal,/create
	mwrfits,blue_varspec_cal,blueoutcal,/silent
	mwrfits,blue_skyspec_cal,blueoutcal,/silent
	mwrfits,blue_extract,blueoutcal,/silent
;	mwrfits,blue_illspec_cal,blueoutcal,/silent

	; Add more than region in a slit
	Print,'Extract an additional Region? (y/n/bye)'
	read,': ',redo
        Ok=0
        while Ok ne 1 do begin
		if redo eq 'y' then begin
			GOTO,LeObject
		endif else if redo eq 'bye' then begin
                	print,' '
                	GOTO,LeEnd
		endif else if redo eq 'n' then begin
			print,'Moving on'
			print,' '
			Ok = 1
		endif else begin
			print,'I do not understand ',redo,', please try again.'
			read,': ',redo
		endelse
	endwhile

endfor
LeEnd: print,'Summary --------------------------------------------------'
print,'Extracted ',object-1,' objects from ',ap,' apertures'
if (object-1) gt 3 then begin
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
endif
print,'------------------------------------------------------------------'
;NEED TO ADD THE FORCE COMMANDS TO THE HEADER...
;sxaddpar,red_hdrcal,'BUNIT','erg/cm^2/s/A'


end
