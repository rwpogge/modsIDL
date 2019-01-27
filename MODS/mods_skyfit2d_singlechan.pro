;Nov 6 2014 KVC OSU
;2019 Jan 22 - patch to mods_skyfit2d_singlechan [rwp/osu]
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
;  plot_2d_spec_skyfit -- plot the 2D emission for a given aperture/slit
;=========================================================================
pro plot_2d_spec_skyfit,sciimg,xx1,xx2, $
	lower1=lower1, upper1=upper1, $
	lower_arr=lower_arr, upper_arr=upper_arr,$
        man_scale = man_scale, $
	title = title

	wset,2
	ncolors = !D.Table_Size
	device,DECOMPOSED=0
	LoadCT,39;0
	min_y = min(xx1[*])
	max_y = max(xx2[*])

        scale_min = -1
        scale_max = 10*median(sciimg[min_y:max_y,*])
        if keyword_set(man_scale) then begin
                Print,'Min Scale for 2D plot (',scale_min,')'
                readnum,scale_min
                Print,'Max Scale for 2D plot (',scale_max,')'
                readnum,scale_max
        endif

	plot_position = [0.05,0.12,0.99,0.95]
	slitim = transpose(BytScl(sciimg[min_y:max_y,*],min=scale_min,max=scale_max,/nan))
	TVImage,slitim,Position=plot_position,/erase
	xpix = indgen(8192)

	plot,xpix,xx1[*]-min_y,/noerase,Position=plot_position,$
		xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,max_y-min_y-1], $
		linestyle=2,thick=2,xtitle='Pixels',charsize=1.2,title=title
	plot,xpix,xx2[*]-min_y,/noerase,Position=plot_position,$
		xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,max_y-min_y-1], $
		linestyle=2,thick=2,xtitle='Pixels',charsize=1.2

	if keyword_set(lower1) then $
		plot,xpix,xx1[*]-min_y+lower1,/noerase, $
			Position=plot_position,xstyle=5,xrange=[0,8191],$
			ystyle=5,yrange=[0,max_y-min_y-1],Color=cgColor('red'),$
			linestyle=2,thick=2,xtitle='Pixels',charsize=1.2
	if keyword_set(upper1) then $
		plot,xpix,xx1[*]-min_y+upper1,/noerase, $
			Position=plot_position,xstyle=5,xrange=[0,8191],$
			ystyle=5,yrange=[0,max_y-min_y-1],Color=cgColor('red'),$
                        linestyle=2,thick=2,xtitle='Pixels',charsize=1.2

	if keyword_set(lower_arr) and n_elements(lower_arr) gt 1 then for j=1,n_elements(lower_arr)-1 do $
		plot,xpix,xx1[*]-min_y+lower_arr[j],/noerase, $
                        Position=plot_position,xstyle=5,xrange=[0,8191],$
                        ystyle=5,yrange=[0,max_y-min_y-1],Color=cgColor('firebrick'),$
                        linestyle=2,thick=2,xtitle='Pixels',charsize=1.2
	if keyword_set(upper_arr) and n_elements(upper_arr) gt 1 then for j=1,n_elements(upper_arr)-1 do $
		plot,xpix,xx1[*]-min_y+upper_arr[j],/noerase, $
                        Position=plot_position,xstyle=5,xrange=[0,8191],$
                        ystyle=5,yrange=[0,max_y-min_y-1],Color=cgColor('firebrick'),$
                        linestyle=2,thick=2,xtitle='Pixels',charsize=1.2


return
end

;=========================================================================
;  extractionplot_skyfit -- plot the extraction cuts 
;=========================================================================
pro extractionplot_skyfit, wave, wave_good, $
	spectrum, midwave, slit_profile, $
	lower, upper, lower_arr, upper_arr, $
	centersum=centersum, $
        edge_mask = edge_mask, $
	linemask = linemask, lw=lw,z=z, $
	logprofile = logprofile,       $
	slit_profile_mean = slit_profile_mean,result_xmean=result_xmean, $
	title = title

common skyfit_exclude1, blank_low, blank_hi

;Plot up the centering and extractions
x0 = wave[wave_good[0]]-20
x1 = max(wave[wave_good])+20

;bottom panel
print,'---------------------------------------------------------------'
profile_size = size(slit_profile,/dim)
info = strcompress(profile_size)+' pixels across the slit'
;print,info

bottom_plot_pos = [0.13,0.09,0.95,0.45]
top_plot_pos    = [0.13,0.56,0.95,0.95]
if keyword_set(logprofile) then minval = min(slit_profile/max(slit_profile))/5 else minval = 0.
plot,indgen(profile_size),slit_profile/max(slit_profile),XRANGE=[0,profile_size],$
	position=bottom_plot_pos,ytitle='Slit Profile',xtitle='Pixels',$
        yrange=[minval,1],xstyle=1, $
	charsize=1.5,ystyle=1,ylog=logprofile, $
	title = 'Solid line = maxium cut, Dashed line = mean cut'
   minval = 0;10^(!Y.CRANGE[0])
   maxval = 1;10^(!Y.CRANGE[1])
   line_y = [minval,maxval]
   line_x1 = [upper,upper]
   line_x2 = [lower,lower]
   niter = size(blank_low,/dim)
   niter = niter[0]-1
;   for i=0,niter do polyfill,[blank_low[i]>0,blank_low[i]>0,blank_hi[i]<profile_size,blank_hi[i]<profile_size],$
;	[minval,maxval,maxval,minval],color=cgColor('firebrick')
;   polyfill,[0,0,edge_mask,edge_mask],[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data
;   polyfill,[profile_size-edge_mask,profile_size-edge_mask,profile_size,profile_size],[minval,maxval,maxval,minval],$
;	color=cgColor('firebrick'),/data
   if n_elements(lower_arr) gt 1 then for j=1,n_elements(lower_arr)-1 do $
	polyfill,[lower_arr[j],lower_arr[j],upper_arr[j],upper_arr[j]], $
		[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data
   polyfill,[lower,lower,upper,upper],[minval,maxval,maxval,minval],color=cgColor('red'),/data
   oplot,line_x1,line_y,Color=cgColor('beige')
   oplot,line_x2,line_y,Color=cgColor('beige')

if keyword_set(logprofile) then minval = min(slit_profile/max(slit_profile))/5 else minval = 0.
plot,indgen(profile_size),slit_profile/max(slit_profile),XRANGE=[0,profile_size],$
        position=bottom_plot_pos,charsize=1.5,ystyle=1,$
        /noerase,yrange=[minval,1],xstyle=1,ylog=logprofile
plot,indgen(profile_size),slit_profile_mean/max(slit_profile_mean),XRANGE=[0,profile_size],$
        position=bottom_plot_pos,charsize=1.5,ystyle=1,$
        /noerase,yrange=[minval,1],xstyle=1,linestyle=2,color=cgcolor('gold'),thick=2,ylog=logprofile

;top panel
plot,wave,spectrum,position=top_plot_pos,thick=2,$
	/noerase,XRANGE=[x0,x1],ytitle='Summed Spectrum', xstyle=1,$
        charsize=1.5,ystyle=16;,title='Vertical lines highlight region being summed to create the slit profile.'
line_x = [wave[midwave],wave[midwave]]
line_xmean = [wave[result_xmean],wave[result_xmean]]
line_x1 = [wave[midwave-(centersum/2-1)],wave[midwave-(centersum/2-1)]]
line_x2 = [wave[midwave+(centersum/2)],wave[midwave+(centersum/2)]]
minval = !Y.CRANGE[0]
maxval = !Y.CRANGE[1]
line_y = [minval,maxval]

if keyword_set(linemask) then begin
	for i=0,n_elements(linemask)-1 do begin
		if (linemask[i] gt !X.CRANGE[0]) and (linemask[i] lt !X.CRANGE[1]) then $
			polyfill,[(linemask[i]-lw[i])*z,(linemask[i]-lw[i])*z,$
				  (linemask[i]+lw[i])*z,(linemask[i]+lw[i])*z],$
				  [!Y.CRANGE[0],!Y.CRANGE[1],!Y.CRANGE[1],!Y.CRANGE[0]],$
				  color=cgColor('dodger blue'),/data
	endfor
endif
oplot,line_x,line_y,Color=cgColor('firebrick'),thick=2
oplot,line_xmean,line_y,Color=cgColor('firebrick'),linestyle=2,thick=2
;oplot,line_x1,line_y,Color=cgColor('beige')
;oplot,line_x2,line_y,Color=cgColor('beige')
plot,wave,spectrum,/noerase,position=top_plot_pos,XRANGE=[x0,x1],ystyle=16,charsize=1.5,xstyle=1,$
	thick=2,title=title

;re-establish the cut as the data reference
if keyword_set(logprofile) then minval = min(slit_profile/max(slit_profile))/5 else minval = 0.
plot,indgen(profile_size),slit_profile/max(slit_profile),XRANGE=[0,profile_size],$
        position=bottom_plot_pos,charsize=1.5,ystyle=1,$
        /noerase,yrange=[minval,1],xstyle=1,thick=2,ylog=logprofile

return
end

;=========================================================================
;  mods_skyfit2d_singlechan -- Subtract sky emission from single channel
;		MODS grating spectra 
;=========================================================================
pro mods_skyfit2d_singlechan,scifile, $
	skyslit=skyslit, $
	z=z, $
	convbeam=convbeam, $
	apertures=apertures, $
	outname=outname, $
	wave=wave, $
	centerLine=centerLine, $
	trim_t = trim_t, trim_b = trim_b, $
	boxcar = boxcar, $
	slits = slits, $
	noPlotAp=noPlotAp, $
	man_scale = man_scale, $
        clobber = clobber,verbose=verbose, $
	linemaskfile = linemaskfile, $
	centersum = centersum, $
	logprofile = logprofile, $
	edge_mask = edge_mask, $
	broad_lines = broad_lines

if not keyword_set(verbose) then !Quiet=1

;deal with parameters that were not entered
if n_elements(trim_t) EQ 0 then trim_t = 0
if n_elements(trim_b) EQ 0 then trim_b = 0
if n_elements(edge_mask) EQ 0 then edge_mask = 4

if keyword_set(linemaskfile) then begin	;masking strong emission lines
	if strtrim(string(linemaskfile),2) EQ '1' then linemask_list = GETENV('XIDL_DIR') + $
		'/Spec/Longslit/pro/LBT/MODS/Calib_Lib/linemask.lst' $
		else linemask_list = linemask
	readcol,linemask_list,linemask,lw,scale_reg1,scale_reg2,/silent
	if keyword_set(broad_lines) then lw = lw +2.
endif

if n_elements(z) EQ '' then begin
   if keyword_set(verbose) then print, 'You gave me no redshift.  I will assume (1+z)=1'
   if keyword_set(verbose) then wait,0.5
   z=0
endif
z = z+1 ;later code uses z+1 so increment up

if n_elements(skyslit) EQ '' then begin
   if keyword_set(verbose) then print, 'You did not list a skyslit.  I assume that means none was cut in the mask (or you have a long-slit observation).'
   if keyword_set(verbose) then wait,0.5
   skyslit=0
endif

common skyfit_exclude1, blank_low, blank_hi
blank_low = [0]
blank_hi = [0]

; Load the images
verify = file_test(scifile)
if verify then sciimg = mrdfits(scifile,0,hdr,/silent) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid science image'
                print,scifile,' does not exist. Check your Science folder'
                stop
        endelse

varimg = mrdfits(scifile,1,v_hdr,/silent)

; Get Mask and Instrument Info
instrument = strcompress(sxpar(hdr[*,0], 'INSTRUME'), /rem)
if instrument eq 'MODS1R' then channel_code = 'm1r' $
else if instrument eq 'MODS2R' then channel_code = 'm2r' $
else if instrument eq 'MODS1B' then channel_code = 'm1b' $
else if instrument eq 'MODS2B' then channel_code = 'm2b' $
else begin
        print,instrument,' is unknown'
        stop
endelse
mask = strcompress(sxpar(hdr[*,0], 'MASKNAME'), /rem)
if mask eq 'LS60x5' then apertures=[1]
if strmid(mask,0,2) eq 'LS' then mask = 'ls' 

if n_elements(centerLine) EQ '' then begin
        if (channel_code eq 'm1b') or (channel_code eq 'm2b') then low_cent = floor(4850*z)
        if (channel_code eq 'm1b') or (channel_code eq 'm2b') then hi_cent = ceil(5070*z)
        if (channel_code eq 'm1r') or (channel_code eq 'm2r') then low_cent = floor(6650*z)
        if (channel_code eq 'm1r') or (channel_code eq 'm2r') then hi_cent = ceil(6800*z)
endif else begin
        low_cent = floor(centerLine[0])
        hi_cent = ceil(centerLine[1])
endelse

; Out file
if n_elements(outname) EQ '' then begin
	name_len=STRLEN(scifile)
        stub = STRSPLIT(scifile,'-',/extract)
	out = STRMID(scifile, 0, name_len-8) + '_s2d.fits'
endif else begin
        stub = STRSPLIT(scifile,'-',/extract)
        out = stub[0] +'-' + channel_code + '_' + outname + '_s2d.fits'
endelse
verify = file_test(out)
if (verify) and (not keyword_set(clobber)) then begin
        if verify then print,out,' already exists. If you want to overwite it please use /clobber'
	return
        stop
endif

; Wavelength image
if (n_elements(wave) EQ '') then wave = 'wave-comp_'+ channel_code+'_'+mask+'.fits'
verify = file_test(wave)
if verify then waveimg = mrdfits(wave,0,whdr,/silent) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid wave-image'
                print,wave,' does not exist.'
                stop
		waveimg = mrdfits(blue,0,whdr,/silent)
        endelse

; Gaussian convolution?
if n_elements(convbeam) EQ '' then convbeam=[1.0,2.0]
if n_elements(boxcar) NE '' then begin
	if keyword_set(verbose) then print,'USING A ',boxcar,' PIXEL BOXCAR TO SMOOTH'
	sciimg = smooth(sciimg,boxcar)
	waveimg = smooth(waveimg,boxcar)
endif else begin
	if keyword_set(verbose) then print,'USING A ',convbeam,' PIXEL GAUSSIAN TO SMOOTH'
	sciimg = filter_image(sciimg,FWHM=convbeam,/ALL)
	waveimg = filter_image(waveimg,FWHM=convbeam,/ALL)
endelse

;open a window for plotting
title = 'Extraction and Slit Profile'
if not WindowAvailable(1) then window,1,retain=2,title=title
title = '2D Spectrum'
if (not keyword_set(noPlotAp)) and (not WindowAvailable(2)) then window,2,retain=2,title=title

;create a new file to fill with sky-subtracted data
fit = sciimg

dims=size(sciimg,/dim)
nx=dims[0]
ny=dims[1]
instrument = strcompress(sxpar(hdr[*, 0], 'INSTRUME'), /rem)
mms = strcompress(sxpar(hdr[*,0], 'MASKINFO')+'.mms', /rem)
grat = strcompress(sxpar(hdr[*, 0], 'GRATNAME'), /rem)

;update the headers
   header_line = 'mods_skyfit2d: beta version 2014-09-02'
   sxaddpar,hdr,'HISTORY',header_line

   ;wavelength
   header_line = wave
   sxaddpar,hdr,'WAVEFILE',header_line
   header_line = z
   sxaddpar,hdr,'REDSHFT',header_line
   header_line = skyslit
   sxaddpar,hdr,'SKYSLIT',header_line
   header_line = low_cent
   sxaddpar,hdr,'LOWCENT',header_line
   header_line = hi_cent
   sxaddpar,hdr,'HICENT',header_line

   ;convbeam
   if n_elements(boxcar) NE '' then begin
        header_line = 'mods_skyfit2d_ - using a boxcar'
        sxaddpar,hdr,'HISTORY',header_line
   endif
   header_line = 'mods_skyfit2d: smoothed for skysub using ['$
        +strtrim(string(convbeam[0]),2)+','+strtrim(string(convbeam[1]),2)+']'
   sxaddpar,hdr,'CONVBEAM',header_line

;create the slit structures 
  if n_elements(slits) eq '' then slits = 'slits-flat_'+channel_code+'_'+mask+'.fits'
  verify = file_test(slits)
  if verify then slits = mrdfits(slits,1,tset_head,/silent) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid slit-image'
                print,slits,' does not exist.'
                stop
                slits = mrdfits(slits,1,tset_head,/silent)
        endelse

   dim_t = slits[0].dims
   nslit=size(slits.xx1,/dim)
   nslit=nslit[1]
   nx_t = dims[0]
   ny_t = dims[1]

   xx1 = slits[0].xx1
   xx2 = slits[1].xx2

   objtrace = xx1*0
   objspec = xx1*0
   objwave = xx1*0
   varmask = xx1*0

   objspec_cal = fltarr(9001,nslit)

;extraction struct
extract_tab = {slit: 0, zones: 0, trim_bot: 0, trim_top: 0, lower: 'test', upper: 'test', lower2: 999, upper2: 999, lower3: 999, upper3: 999, lower4: 999, upper4: 999}
extract_tab = replicate(extract_tab,nslit)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                   fit the best sky slit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if skyslit ne 0 then begin
        print,'Performing SKY SUBTRACTION on APERTURE ',skyslit,' in ',scifile
	print,'This is the currently designated PRIMARY SKYSLIT'
	n_selection = 1
	spec1d = fltarr(ny)

        ;get the traces
        xx1_ap = xx1[*,skyslit-1]+trim_b
        xx2_ap = xx2[*,skyslit-1]-trim_t

        lower = 0.
        upper = 0.

        ;define the aperture/slit
        apmask = 0*sciimg
        for iii=0,ny-1 do begin
                for iiii=xx1_ap[iii],xx2_ap[iii] do apmask[iiii,iii] = 1
        endfor
        slitpos = where(apmask eq 1)

        ; get the profile
        ap_im = x_ordrectify(sciimg,xx1_ap[*],xx2_ap[*],/nocorrect)
        ap_wave = x_ordrectify(waveimg,xx1_ap[*],xx2_ap[*],/nocorrect)
        ap_im = transpose(ap_im)
        ap_wave = transpose(ap_wave)

        slit_dims=size(ap_im,/dim) 
        slit_nx=slit_dims[0]
        slit_ny=slit_dims[1]
    
        wav1d = ap_wave[*,slit_ny/2]
        spec1d = wav1d*0
        spec1d_tmp = wav1d*0 
        wav1d_tmp = wav1d
    
        goodwav = where(wav1d lt hi_cent and wav1d gt low_cent)
        tmp = ap_im[goodwav,*]
        total_spec = total(tmp,2)
        result=where(total_spec[*] eq max(total_spec))
        result_x = goodwav[result]

	sorted = total_spec[sort(total_spec)]
	loc_mean = value_locate(sorted,median(total_spec))
        result_mean=where(total_spec[*] eq sorted[loc_mean])
	result_xmean = goodwav[result_mean[0]]
        tmp = ap_im[result_xmean-(centersum/2-1):result_xmean+(centersum/2),*]
	slit_profile_mean = total(tmp,1)

        tmp = ap_im[result_x-(centersum/2-1):result_x+(centersum/2),*]
        slit_profile = total(tmp,1)

        ;get the initial 1D specturm
        for iii=0,ny-1 do begin
                ;blank thin strips along top and bottom
                bot = xx1_ap[iii]-edge_mask
                top = xx1_ap[iii]+edge_mask
                apmask[bot:top,iii] = 0
                bot = xx2_ap[iii]-edge_mask
                top = xx2_ap[iii]+edge_mask
                apmask[bot:top,iii] = 0
                ;extract the remaining mean spectrum
                remaining_sky = where(apmask[*,iii] eq 1)
                spec1d[iii] = total(sciimg[remaining_sky,iii])
                wav1d[iii] = median(waveimg[remaining_sky,iii])
        endfor
        apmask = 0*sciimg

        lower_arr = [0]
        upper_arr = [0]
        ;define sky regions
        defineskyslit:
        ; plot slit and traces
        lower = 0. & upper = 0.
        if (not keyword_set(noPlotAp)) then plot_2d_spec_skyfit,sciimg,xx1_ap,xx2_ap, $
                lower1=lower,upper1=upper,lower_arr=lower_arr,upper_arr=upper_arr,man_scale=man_scale,title=scifile
        wset,1
        extractionplot_skyfit,wav1d,goodwav,spec1d,result_x,$
               slit_profile,lower,upper,lower_arr,upper_arr,centersum=centersum,edge_mask = edge_mask,$
                linemask = linemask, lw=lw, z=z,slit_profile_mean=slit_profile_mean,result_xmean=result_xmean,logprofile = logprofile,title=scifile
        print,'Please click where you would like to place the lower bound for a sky block.'
        while !MOUSE.button ne 1 do cursor,xloc,yloc,/data,wait=3
        lower = xloc
        wait,0.3
        !MOUSE.button = 0
        print,'Left click your mouse on the upper edge of the sky block.'
        while !MOUSE.button ne 1 do cursor,xloc,yloc,/data,wait=3
        upper = xloc
        wait,0.3
        !MOUSE.button = 0
        if lower gt upper then begin
                tmp = lower
                lower = upper
                upper = tmp
                print,'Your lower bound was greater than your upper bound.  They have been reversed'
        endif
        for iii=0,ny-1 do begin
                 bot = xx1_ap[iii] + lower
                 top = xx1_ap[iii] + upper
                 apmask[bot:top,iii] = 1
                 remaining_sky_tmp = where(apmask[*,iii] eq 1)
                 spec1d_tmp[iii] = total(sciimg[remaining_sky_tmp,iii])
                 wav1d_tmp[iii] = median(waveimg[remaining_sky_tmp,iii])
         endfor
        extractionplot_skyfit,wav1d_tmp,goodwav,spec1d_tmp,result_x,$
               slit_profile,lower,upper,lower_arr,upper_arr,centersum=centersum,edge_mask = edge_mask, $
                linemask = linemask, lw=lw, z=z,slit_profile_mean=slit_profile_mean,result_xmean=result_xmean,logprofile = logprofile
        if (not keyword_set(noPlotAp)) then plot_2d_spec_skyfit,sciimg,xx1_ap,xx2_ap, $
                lower1=lower,upper1=upper,lower_arr=lower_arr,upper_arr=upper_arr,man_scale=man_scale,title=scifile

        Ok=0
        while Ok ne 1 do begin
                print,'Is this the region you wish to designate as sky for modeling purposes? (y/n/abort)'
                redo = ''
                read,': ',redo
                if redo eq 'y' then begin
                        for iii=0,ny-1 do begin
                                bot = xx1_ap[iii] + lower
                                top = xx1_ap[iii] + upper
                                apmask[bot:top,iii] = 1
                                remaining_sky = where(apmask[*,iii] eq 1)
                                spec1d[iii] = total(sciimg[remaining_sky,iii])
                                wav1d[iii] = median(waveimg[remaining_sky,iii])
                        endfor
                        lower_arr = [lower_arr,lower]
                        upper_arr = [upper_arr,upper]
                        Ok = 1
                endif else if redo eq 'n' then begin
                        GOTO,defineskyslit
                endif else if redo eq 'abort' then begin
                        stop
                endif else begin
                        print,'I do not understand ',redo
                endelse
        endwhile
        Print,'Add an additional sky region for modeling? (y/n/abort)'
        Ok=0
        while Ok ne 1 do begin
                 read,': ',redo
                 if redo eq 'y' then begin
                         n_selection += 1
                         GOTO,defineskyslit
                 endif else if redo eq 'abort' then begin
                         stop
                 endif else if redo eq 'n' then begin
			ok=1
 		 endif else begin
                        print,'I do not understand ',redo
                 endelse
        endwhile

        ;fit the pure sky
        if (channel_code eq 'm1b') or (channel_code eq 'm2b') then $
                sset_sky = bspline_iterfit(waveimg[slitpos],sciimg[slitpos],$
                        bkspace=.95,nord=2,invvar=apmask[slitpos])

        if (channel_code eq 'm1r') or (channel_code eq 'm2r') then $
                sset_sky = bspline_iterfit(waveimg[slitpos],sciimg[slitpos],$
                        bkspace=.65,nord=4,invvar=apmask[slitpos])

        ;make a pure skyslit fit
        fullfit = bspline_valu(waveimg,sset_sky)
        fullfit2 = sciimg - fullfit
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                      fit each slit on its own
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(apertures) EQ '' then apertures = indgen(nslit)+1
ndo = size(apertures,/dim)
ndo=ndo[0]

for ii=0,ndo-1 do begin
	ap=apertures[ii]
	print,'Performing SKY SUBTRACTION on APERTURE ',ap,' in ',scifile
	n_selection = 1
	NOSKY = 0

	;get the traces
	xx1_ap = xx1[*,ap-1]+trim_b
	xx2_ap = xx2[*,ap-1]-trim_t

        lower = 0.
        upper = 0.

	;define the aperture/slit
	apmask = 0*sciimg
	for iii=0,ny-1 do begin
                for iiii=xx1_ap[iii],xx2_ap[iii] do apmask[iiii,iii] = 1
        endfor
	slitpos = where(apmask eq 1)

	; get the profile
	ap_im = x_ordrectify(sciimg,xx1_ap[*],xx2_ap[*],/nocorrect)
	ap_wave = x_ordrectify(waveimg,xx1_ap[*],xx2_ap[*],/nocorrect)
	ap_im = transpose(ap_im)
	ap_wave = transpose(ap_wave)

	slit_dims=size(ap_im,/dim)
	slit_nx=slit_dims[0]
	slit_ny=slit_dims[1]

	wav1d = ap_wave[*,slit_ny/2]
	spec1d = wav1d*0
	spec1d_tmp = wav1d*0
	wav1d_tmp = wav1d

	goodwav = where(wav1d lt hi_cent and wav1d gt low_cent)
	tmp = ap_im[goodwav,*]
	total_spec = total(tmp,2)
        result=where(total_spec[*] eq max(total_spec))
	result_x = goodwav[result]

        sorted = total_spec[sort(total_spec)]
        loc_mean = value_locate(sorted,median(total_spec))
        result_mean=where(total_spec[*] eq sorted[loc_mean])
        result_xmean = goodwav[result_mean[0]]
        tmp = ap_im[result_xmean-(centersum/2-1):result_xmean+(centersum/2),*]
        slit_profile_mean = total(tmp,1)

	tmp = ap_im[result_x-(centersum/2-1):result_x+(centersum/2),*]
	slit_profile = total(tmp,1)

	;get the initial 1D specturm
	for iii=0,ny-1 do begin
		;blank thin strips along top and bottom
		bot = xx1_ap[iii]-edge_mask
		top = xx1_ap[iii]+edge_mask
                apmask[bot:top,iii] = 0
                bot = xx2_ap[iii]-edge_mask
                top = xx2_ap[iii]+edge_mask
                apmask[bot:top,iii] = 0
                ;extract the remaining mean spectrum
                remaining_sky = where(apmask[*,iii] eq 1)
                spec1d[iii] = total(sciimg[remaining_sky,iii])
                wav1d[iii] = median(waveimg[remaining_sky,iii])
	endfor
	apmask = 0*sciimg

	lower_arr = [0]
	upper_arr = [0]
	;define sky regions
	defineregion:
        ; plot slit and traces
        lower = 0. & upper = 0.
        if (not keyword_set(noPlotAp)) then plot_2d_spec_skyfit,sciimg,xx1_ap,xx2_ap, $
		lower1=lower,upper1=upper,lower_arr=lower_arr,upper_arr=upper_arr,man_scale=man_scale,title=scifile
	wset,1
	extractionplot_skyfit,wav1d,goodwav,spec1d,result_x,$
               slit_profile,lower,upper,lower_arr,upper_arr,centersum=centersum,edge_mask = edge_mask,$
		linemask = linemask, lw=lw, z=z,slit_profile_mean=slit_profile_mean,result_xmean=result_xmean,$
		logprofile = logprofile,title=scifile
	print,'Please click where you would like to place one edge of a sky block.'
	while !MOUSE.button ne 1 do cursor,xloc,yloc,/data,wait=3
	lower = xloc
	wait,0.3
	!MOUSE.button = 0
	print,'Left click your mouse on the other edge of the sky block.'	
	while !MOUSE.button ne 1 do cursor,xloc,yloc,/data,wait=3
	upper = xloc
	wait,0.3
	!MOUSE.button = 0
	if lower gt upper then begin
		tmp = lower
		lower = upper
		upper = tmp
	endif
	for iii=0,ny-1 do begin
                 bot = xx1_ap[iii] + lower
                 top = xx1_ap[iii] + upper
                 apmask[bot:top,iii] = 1
                 remaining_sky_tmp = where(apmask[*,iii] eq 1)
                 spec1d_tmp[iii] = total(sciimg[remaining_sky_tmp,iii])
                 wav1d_tmp[iii] = median(waveimg[remaining_sky_tmp,iii])
         endfor
	extractionplot_skyfit,wav1d_tmp,goodwav,spec1d_tmp,result_x,$
               slit_profile,lower,upper,lower_arr,upper_arr,centersum=centersum,edge_mask = edge_mask, $
		linemask = linemask, lw=lw, z=z,slit_profile_mean=slit_profile_mean,result_xmean=result_xmean,$
		logprofile = logprofile,title=scifile
        if (not keyword_set(noPlotAp)) then plot_2d_spec_skyfit,sciimg,xx1_ap,xx2_ap, $
                lower1=lower,upper1=upper,lower_arr=lower_arr,upper_arr=upper_arr,man_scale=man_scale,title=scifile
	
	Ok=0
	while Ok ne 1 do begin
		print,'Is this the region you wish to designate as sky for modeling purposes? (y/n/nosky/abort)' 
		redo = ''
		read,': ',redo
		redo = strlowcase(redo)
		if redo eq 'y' then begin
			for iii=0,ny-1 do begin
				bot = xx1_ap[iii] + lower
				top = xx1_ap[iii] + upper
				apmask[bot:top,iii] = 1
				remaining_sky = where(apmask[*,iii] eq 1)
				spec1d[iii] = total(sciimg[remaining_sky,iii])
				wav1d[iii] = median(waveimg[remaining_sky,iii])
			endfor
			lower_arr = [lower_arr,lower]
			upper_arr = [upper_arr,upper]
			Ok = 1
		endif else if redo eq 'n' then begin
			GOTO,defineregion
		endif else if redo eq 'nosky' then begin
			NOSKY=1
			Ok = 1
		endif else if redo eq 'abort' then begin
			stop
		endif else begin
			print,'I do not understand ',redo
		endelse
	endwhile
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;  Add the ability to include more than one sky region in a slit
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if NOSKY eq 0 then begin
		Print,'Add an additional sky region? (y/n/abort)'
	        Ok=0
		while Ok ne 1 do begin
			read,': ',redo
			if redo eq 'y' then begin
				n_selection += 1
				GOTO,defineregion
			endif else if redo eq 'abort' then begin
				stop
			endif else if redo eq 'n' then begin
				print,'You marked ',n_selection,' region(s) in this slit.'
				extract_tab[ii].slit = ap
				extract_tab[ii].zones = n_selection
				extract_tab[ii].trim_bot = trim_b
				extract_tab[ii].trim_top = trim_t
				
				lower_str = '' & upper_str = ''
				for j=1,n_selection do lower_str = lower_str + string(lower_arr[j])
				for j=1,n_selection do upper_str = upper_str + string(upper_arr[j])
				extract_tab[ii].lower = lower_str
				extract_tab[ii].upper = upper_str
				Ok = 1	
			endif else begin
				print,'I do not understand ',redo
			endelse
		endwhile	
	endif

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; FIT THE SKY IN THE SELETCED APERTURE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        if NOSKY eq 0 then begin
		if (channel_code eq 'm1b') or (channel_code eq 'm2b') then begin
		        sset=bspline_iterfit(waveimg[slitpos],sciimg[slitpos],$
				bkspace=.99,nord=3,invvar=apmask[slitpos])
			extract=sciimg[slitpos]
			extract_wave=waveimg[slitpos]
			yfit = bspline_valu(waveimg[slitpos],sset)
			if keyword_set(linemask) then begin
				valid_line = where((linemask lt max(waveimg[slitpos])) and $
						   (linemask gt min((waveimg[slitpos]))) and $
						   (linemask lt 6000.))
				for i=0,n_elements(valid_line)-1 do begin
					;linemask[valid_line[i]]						
					line = where((extract_wave lt (linemask[valid_line[i]]+lw[valid_line[i]])*z) and $
						     (extract_wave gt (linemask[valid_line[i]]-lw[valid_line[i]])*z))
					scalereg = where(extract_wave lt scale_reg2[valid_line[i]]*z and $
						         extract_wave gt scale_reg1[valid_line[i]]*z)
					if skyslit ne 0 then begin
						yfit_sky = bspline_valu(waveimg[slitpos],sset_sky)
;						scale = median(yfit[scalereg]/yfit_sky[scalereg])
;						yfit[line] = scale*yfit_sky[line]
						scale = median(yfit[scalereg]-yfit_sky[scalereg])
						yfit[line] = scale + yfit_sky[line]
					endif
				endfor
			endif
			subimg = sciimg
			subimg[slitpos] = extract-yfit
			fit[slitpos] = yfit	
		endif
                if (channel_code eq 'm1r') or (channel_code eq 'm2r') then begin
	                sset=bspline_iterfit(waveimg[slitpos],sciimg[slitpos],$
	                        bkspace=.99,nord=3,invvar=apmask[slitpos])
	                extract=sciimg[slitpos]
	                extract_wave=waveimg[slitpos]
			yfit = bspline_valu(waveimg[slitpos],sset)
			if keyword_set(linemask) then begin
                                valid_line = where((linemask lt max(waveimg[slitpos])) and $
                                                   (linemask gt min((waveimg[slitpos]))) and $
						   (linemask gt 5500))
                                for i=0,n_elements(valid_line)-1 do begin
                                        ;linemask[valid_line[i]]                                                
                                        line = where((extract_wave lt (linemask[valid_line[i]]+lw[valid_line[i]])*z) and $
                                                     (extract_wave gt (linemask[valid_line[i]]-lw[valid_line[i]])*z))
                                        scalereg = where(extract_wave lt scale_reg2[valid_line[i]]*z and $
                                                         extract_wave gt scale_reg1[valid_line[i]]*z)
                                        if skyslit ne 0 then begin
                                                yfit_sky = bspline_valu(waveimg[slitpos],sset_sky)
;                                                scale = median(yfit[scalereg]/yfit_sky[scalereg])
						scale = median(yfit[scalereg]-yfit_sky[scalereg])
;                                                yfit[line] = scale*yfit_sky[line]
						yfit[line] = scale + yfit_sky[line]
                                        endif
                                endfor
                        endif
	                subimg = sciimg
	                subimg[slitpos] = extract-yfit
	                fit[slitpos] = yfit
                endif

	endif else begin
		if skyslit EQ 0 then begin
			print, 'You did not list a skyslit.'
			print, ' Please re-run MODS_SKYFIT2D_SINGLECHAN with a skyslit indicated in the parameter file.'
   			stop
		endif
		print,'Applying the sky solution from Aperture ',skyslit,' to aperture ',ap
                extract=sciimg[slitpos]
		yfit = bspline_valu(waveimg[slitpos],sset_sky)
                subimg = sciimg
                subimg[slitpos] = extract-yfit
                fit[slitpos] = yfit
	endelse

        subimg = subimg>(-50)
	sciimg = subimg

	; Save the output
		;0 - sky-subtracted image
		;1 - variance image
		;2 - extraction table
		;3 - sky fit
		;4 - sky fit using only skyslit
		;5 - sky-subtracted image using only skyslit

	mwrfits,subimg,out,hdr,/create
	mwrfits,varimg,out
	mwrfits,extract_tab,out
	mwrfits,fit,out
endfor

	IF skyslit NE 0 then mwrfits,fullfit,out
        IF skyslit NE 0 then mwrfits,fullfit2,out
;print,"MODS_EXTRACT COMMAND: mods_extract_singlechan,'",out,"'"
end
