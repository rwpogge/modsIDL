;Jan 2 2015 KVC OSU
;2019 Jan 22 - mods_fluxstand_singlechan - [rwp/osu]
;=========================================================================
;  plot_2d_spec_fluxstand -- plot the 2D emission for a given aperture/slit
;=========================================================================
pro plot_2d_spec_fluxstand,sciimg,xx1,xx2, $
        center_trace = center_trace, $
	upper = upper, $
	lower = lower, $
        title = title

	plot_position = [0.05,0.12,0.99,0.95]

        wset,2
        ncolors = !D.Table_Size
        device,DECOMPOSED=0
        LoadCT,0 			; 0 - B&W   39 - Color
        min_y = min(xx1[*])
        max_y = max(xx2[*])
	scale_min = -1
        scale_max = 50*abs(median(sciimg[min_y:max_y,*]))

        slitim = transpose(BytScl(sciimg[min_y:max_y,*],min=scale_min,max=scale_max,/nan))

        TVImage,slitim,Position=plot_position,/erase
        xpix = indgen(8192)

        plot,xpix,xx1[*]-min_y,/noerase,Position=plot_position,$
                xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,max_y-min_y-1], $
                linestyle=2,thick=2,xtitle='Pixels',charsize=1.2,title=title
        plot,xpix,xx2[*]-min_y,/noerase,Position=plot_position,$
                xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,max_y-min_y-1], $
                linestyle=2,thick=2,xtitle='Pixels',charsize=1.2

	if keyword_set(center_trace) then $
		plot,xpix,center_trace[*]-min_y,/noerase, $
                        Position=plot_position,xstyle=5,xrange=[0,8191],$
                        ystyle=5,yrange=[0,max_y-min_y-1],Color=cgColor('red'),$
                        linestyle=0,thick=2,xtitle='Pixels',charsize=1.2
        if keyword_set(lower) then $
                plot,xpix,center_trace[*]-min_y-lower,/noerase, $
                        Position=plot_position,xstyle=5,xrange=[0,8191],$
                        ystyle=5,yrange=[0,max_y-min_y-1],Color=cgColor('red'),$
                        linestyle=2,thick=2,xtitle='Pixels',charsize=1.2
        if keyword_set(upper) then $
                plot,xpix,center_trace[*]-min_y+upper,/noerase, $
                        Position=plot_position,xstyle=5,xrange=[0,8191],$
                        ystyle=5,yrange=[0,max_y-min_y-1],Color=cgColor('red'),$
                        linestyle=2,thick=2,xtitle='Pixels',charsize=1.2
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
;  center -- center on line emission in the slit
;=========================================================================
pro center,file,wave,slit_num,wav_min,wav_max,result_x,$
	result_y,goodwav,wav1d,stack,stack2,xx1,xx2,$
	trim1,trim2,lower=lower,upper=upper,$
        emm_line=emm_line,centersum=centersum,subcont=subcont

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
        if keyword_set(emm_line) then result=where(stack eq max(stack)) else result=Value_Locate(stack,mean(stack))
        result_x = goodwav[result]

        ;find which part of the slit is at the median... mininimum pulled out the slit edge
        tmp = ap_im[result_x-(centersum/2-1):result_x+(centersum/2),*]
        stack2 = smooth(total(tmp>0,1),2)
	if keyword_set(subcont) then begin
		tmp_cont = ap_im[result_x-4+subcont:result_x+5+subcont,*]
		stack2_cont = smooth(total(tmp_cont>0,1),2)
		stack2 = stack2 - stack2_cont
	endif
	result_y=where(stack2[10:n_elements(stack2)-10] eq max(stack2[10:n_elements(stack2)-10]))
	result_y = result_y[0]+10

	;set width to include 95% of Flux
        limit = 0.10*max(stack2[10:n_elements(stack2)-10])
        below_lim = where(stack2 lt limit)
        lower_loc = max(where(below_lim lt result_y[0]))
        upper_loc = min(where(below_lim gt result_y[0]))
        lower = result_y[0] - below_lim[lower_loc]<result_y
        upper = below_lim[upper_loc] - result_y[0]<(n_elements(stack2)-result_y)
        lower = lower[0]
        upper = upper[0]
return
end

;=========================================================================
;  extractionplot_fluxstand -- plot the extraction cuts 
;=========================================================================
pro extractionplot_fluxstand,w,g,$
	st,x,st2,y, $
	sp,lower,upper,$
        centersum=centersum, $
        edge_mask = edge_mask,subcont=subcont, $
	title = title

common flux_exclude1, blank_low, blank_hi
	bottom_plot_pos = [0.13,0.09,0.95,0.45]
	top_plot_pos    = [0.13,0.59,0.95,0.95]

;Plot up the centering and extractions
;bottom panel
x0 = w[g[0]]-20
x1 = max(w[g])+20
line_x = [w[x],w[x]]
line_x1 = [w[x-(centersum/2-1)],w[x-(centersum/2-1)]]
line_x2 = [w[x+(centersum/2)],w[x+(centersum/2)]]

stacksz = size(st2,/dim)
st2 = st2 +10.D
plot,indgen(stacksz),st2,XRANGE=[0,stacksz],$
	position=bottom_plot_pos,ytitle='Slit Profile',$
        yrange=[min(st2)>(0.1),max(st2)],xstyle=9, $
        charsize=1.5,ystyle=1,/ylog
line_x = [y,y]
minval = 10^(!Y.CRANGE[0])
maxval = 10^(!Y.CRANGE[1])
line_y = [minval,maxval]
line_x1 = [y+upper/0.12,y+upper/0.12]
line_x2 = [y-lower/0.12,y-lower/0.12]
niter = size(blank_low,/dim)
niter = niter[0]-1
for i=0,niter do polyfill,[blank_low[i]>0,blank_low[i]>0,blank_hi[i]<stacksz,blank_hi[i]<stacksz], $
	[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data
plot,indgen(stacksz),st2,XRANGE=[0,stacksz],/ylog,$
        position=bottom_plot_pos,charsize=1.5,ystyle=1,$
        /noerase,yrange=[min(st2)>(0.1),max(st2)],xstyle=9,xtitle='Pixels'
oplot,line_x,line_y,Color=cgColor('dodger blue')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')
axis,xaxis=1,XRANGE=[-(line_x[0]*0.12),(stacksz-line_x[0])*0.12],xtitle='ARCSEC',charsize=1.3,xstyle=1

;top panel
plot,w,sp,position=top_plot_pos,$
        /noerase,XRANGE=[x0-100,x1+50],ytitle='Extracted Spectrum', $
        charsize=1.5,ystyle=1,xtitle='Wavelength',title=title

;re-establish the cut as the data reference
plot,indgen(stacksz),st2,XRANGE=[0,stacksz],$
        position=bottom_plot_pos,ytitle='Slit Profile',$
        /noerase,yrange=[min(st2)>(0.1),max(st2)],xstyle=9, $
        charsize=1.5,ystyle=1,/ylog
axis,xaxis=1,XRANGE=[-(line_x[0]*0.12),(stacksz-line_x[0])*0.12],xtitle='ARCSEC',charsize=1.3,xstyle=1

return
end

;=========================================================================
; mods_fluxstand - extract and create sensitivity functions for
;			 1D dual channel MODS spectra from 2D spectra
;=========================================================================
pro mods_fluxstand_singlechan,scifile, std_name=std_name, $
	outname=outname,wave=wave, $
	centerLine=centerLine, $
	verbose=verbose, $
	noPlotAp=noPlotAp, slits = slits,force=force, $
        centersum = centersum, $
	clobber=clobber, $
	no_tell_mask=no_tell_mask

if not keyword_set(verbose) then !Quiet=1	

if not keyword_set(std_name) then begin
	print,'Please give me the name of a standard in my library'
	print,"AKA: std_name = 'NAME'"
	print,''
	print,'Recommended standard files:'
        print,'         g191b2b_mod_005   g191b2b_stisnic_002  g191b2b_10a'
        print,'         gd71_mod_006      gd71_stisnic_002 '
        print,'         feige34_stis_001  feige34_10a'
        print,'         feige66_002'
        print,'         feige67_002'
        print,'         gd153_mod_005     gd153_stisnic_002  '
        print,'         hz43_mod_005      hz43_stis_001 '
        print,'         hz44_stis_001     '
        print,'         bd_33d2642_fos_003  '
        print,'         bd_28d4211_stis_001 '
        print,'         feige110_stisnic_002'
	stop
endif

;deal with parameters that were not entered
  z=0
  z = z+1 ;later code uses z+1 so increment up
  apertures=[1]

  if n_elements(trim_t) EQ 0 then trim_t = 0
  if n_elements(trim_b) EQ 0 then trim_b = 0

common flux_exclude1, blank_low, blank_hi
blank_low = [0]
blank_hi = [0]

; Load the images
; Science image
verify = file_test(scifile)
if verify then sciimg = mrdfits(scifile,0,hdr,/silent) $
        else begin
                print,'You need a valid science image'
                print,scifile,'does not exist. Check your Science folder'
                stop
        endelse
skyimg = mrdfits(scifile,3,hdrsky,/silent)
varimg = mrdfits(scifile,1,headvar,/silent)
varimg = skyimg * varimg

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
if strmid(mask,0,2) eq 'LS' then mask = 'LS'

if n_elements(centerLine) EQ '' then begin
        if (channel_code eq 'm1b') or (channel_code eq 'm2b') then low_cent = floor(4850*z)
        if (channel_code eq 'm1b') or (channel_code eq 'm2b') then hi_cent = ceil(5070*z)
        if (channel_code eq 'm1r') or (channel_code eq 'm2r') then low_cent = floor(6650*z)
        if (channel_code eq 'm1r') or (channel_code eq 'm2r') then hi_cent = ceil(6800*z)
endif else begin
        low_cent = floor(centerLine[0])
        hi_cent = ceil(centerLine[1])
endelse

hdrcal = hdr
sxaddpar,hdrcal,'CUNIT1','Angstrom'
sxaddpar,hdrcal,'CTYPE1','Linear'
sxaddpar,hdrcal,'CRPIX1',1
if (channel_code eq 'm1r') or (channel_code eq 'm2r') then sxaddpar,hdrcal,'CRVAL1',5500
if (channel_code eq 'm1b') or (channel_code eq 'm2b') then sxaddpar,hdrcal,'CRVAL1',3000
sxaddpar,hdrcal,'CDELT1',0.5
sxaddpar,hdrcal,'FLUX_UNITS','erg/cm^2/s/A'

if (channel_code eq 'm1r') or (channel_code eq 'm2r') then new_wave = 0.5*indgen(9001) + 5500
if (channel_code eq 'm1b') or (channel_code eq 'm2b') then new_wave = 0.5*indgen(5601) + 3000

extfile = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/mods_extiction.dat'
readcol,extfile,ext_wave,ext_val,/silent
ext_coef = bspline_iterfit(ext_wave,ext_val,bkspace=3,nord=3)
ext = bspline_valu(new_wave,ext_coef)

; Out file
if n_elements(outname) EQ '' then begin
        name_len=STRLEN(scifile)
        out = STRMID(scifile, 0, name_len-8) + 'r1d.fits'
	outcal = STRMID(scifile, 0, name_len-8) + 'x1d.fits'
endif else begin
        stub = STRSPLIT(scifile,'-',/extract)
        out = stub[0] +'-' + channel_code + '_' + outname + '_r1d.fits'
	outcal = stub[0] +'-' + channel_code + '_' + outname + '_x1d.fits'
endelse
verify = file_test(out)
if (verify) and (not keyword_set(clobber)) then begin
        print,'=========================================================='
        if verify then print,out,' already exists. If you want to overwite it please use /clobber'
        stop
endif

;;; Slit Flats
line = 'ls Proc/flat*' + channel_code + '*fits*'
spawn,line,flat_filenames
nfile = n_elements(flat_filenames)
if nfile gt 0 then begin
	if nfile eq 1 then flatfile = flat_filenames
	if nfile gt 1 then begin
		select_LS = where( strmatch(flat_filenames,'*ls*' ,/fold_case) eq 1 ,flatcount)	
		if flatcount gt 1 then begin
			print,'Multpile flat files found.'
			print,'SCIENCE FILE: ',scifile
			print,'INDEX NUMBERS: ',select_LS
			print,'FLAT FILES: ',flat_filenames[select_LS]
			print,'Please enter the desired index numer.'
			readnum,select_LS
		endif
		flatfile = flat_filenames[select_LS]
	endif
	flatimg = mrdfits(flatfile)
	flatimg = transpose(flatimg)
endif else begin
	print,'No flat files found.  Please be careful with fluxing.'
	stop
endelse

; Wavelength image
if (n_elements(wave) EQ '') then wave = 'wave-comp_'+ channel_code+'_'+mask+'.fits'
verify = file_test(wave)
if verify then waveimg = mrdfits(wave,0,whdr,/silent) $
        else begin
	        print,'=========================================================='
                print,'You need a valid wave-image'
                print,wave,' does not exist.'
                print,"Use: wave='FILENAME'"
                stop
		waveimg = mrdfits(wave,0,whdr,/silent)
        endelse

;open a window for plotting
if not WindowAvailable(1) then window,1,retain=2,title='Extraction and Slit Profile'
if (not keyword_set(noPlotAp)) and (not WindowAvailable(2)) then window,2,retain=2,title='2D Spectrum'

dims=size(sciimg,/dim)
nx=dims[0]
ny=dims[1]
instrument = strcompress(sxpar(hdr[*, 0], 'INSTRUME'), /rem)
mms = strcompress(sxpar(hdr[*,0], 'MASKINFO')+'.mms', /rem)

;create the slit structures
  if n_elements(slits) eq '' then slits = 'slits-flat_'+channel_code+'_'+mask+'.fits'
  verify = file_test(slits)
  if verify then slits = mrdfits(slits,1,tset_head,/silent) $
        else begin
                print,'=========================================================='
                print,'You need a valid slit-image'
                print,slits,' does not exist.'
                print,"Use: slits='FILENAME'"
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
   objspec = fltarr(ny_t)+1
   objwave = fltarr(ny_t)+1
   varspec = fltarr(ny_t)+1
   skyspec = fltarr(ny_t)+1
   flatspec = fltarr(ny_t)+1

   objspec_cal = fltarr(size(new_wave,/dim))
   varspec_cal = fltarr(size(new_wave,/dim))
   skyspec_cal = fltarr(size(new_wave,/dim))
   illspec_cal = fltarr(size(new_wave,/dim))

;extraction struct
extract = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D], $
        ra_deg: double(slits[0].tarra2), dec_deg: double(slits[0].tardec2), $
        ra_sex: slits[0].tarra, dec_sex: slits[0].tardec}

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
        xx1_ap = xx1[*,ap-1]+trim_b
        xx2_ap = xx2[*,ap-1]-trim_t

        lower = 6.
        upper = 6.
        blank_low = [0]
        blank_hi = [0]
        n_blank = 1

	;currently, returning the center of the slit for ra and dec 
        extract_obj = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D],$
                ra_deg: slits[0].apra2[ap-1], dec_deg: slits[0].apdec2[ap-1], $
                ra_sex: slits[0].apra[ap-1], dec_sex: slits[0].apdec[ap-1]}

        center,sciimg,waveimg,ap,low_cent,hi_cent,$
                result_x,result_y,good,wav1d,$
                stack,stack2,xx1,xx2,trim_b,trim_t,lower=lower,upper=upper,$
                centersum=centersum

	modsadcalc,scifile,wav1d[result_x],wav1d,xdelt,ydelt

	LeSigh: print,'Extracting APERTURE: ',ap
        spec1d = fltarr(ny)
        var1d = fltarr(ny)
	sky1d = fltarr(ny)
	flat1d = fltarr(ny)
        spec1db = fltarr(ny)

        for iii=0,ny-1 do objtrace[iii,ap-1] = xx1[iii,ap-1] + trim_b + result_y
	
        ; added in for method demonstration
        orig = objtrace[*,ap-1]
	
	if (channel_code eq 'm1r') or (channel_code eq 'm2r') then platescale = 0.123
	if (channel_code eq 'm1b') or (channel_code eq 'm2b') then platescale = 0.120
        objtrace[*,ap-1] = objtrace[*,ap-1] - ydelt/platescale

	for iii=0,ny-1 do begin
                bot = objtrace[iii,ap-1]-lower/platescale
                top = objtrace[iii,ap-1]+upper/platescale
                botint = total(bot,/integer)
                topint = total(top,/integer)
                spec1d[iii] = total(sciimg[bot:top,iii])$
                        +(1-(bot-botint))*total(sciimg[bot,iii]) $
                        +(top-topint)*total(sciimg[top+1,iii])
                wav1d[iii] = median(waveimg[bot:top,iii])
                var1d[iii] = sqrt(total(varimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(varimg[bot,iii]) $
                        +(top-topint)*total(varimg[top+1,iii]))
                sky1d[iii] = sqrt(total(skyimg[bot:top,iii]) $
                        +(1-(bot-botint))*total(skyimg[bot,iii]) $
                        +(top-topint)*total(skyimg[top+1,iii]))
		flat1d[iii] = mean(flatimg[bot:top,iii])
	endfor

	; plot slit and traces
        if (not keyword_set(noPlotAp)) then plot_2d_spec_fluxstand,sciimg,xx1_ap,xx2_ap,center_trace=objtrace[*,ap-1],$
			lower=lower/platescale,upper=upper/platescale,title=scifile 
                
	;Plot up the centering and extractions
	wset,1
        extractionplot_fluxstand,wav1d,good,stack,result_x,$
               stack2,result_y, spec1d,lower,upper, $
	       centersum=centersum,title=scifile

	Ok=0
	while Ok ne 1 do begin
                ;print,'Centered on lines: ', wav1d[result_x]
                ;print,'Peak found at (pixels): ',result_y
                print,'Currently using: lower = ',lower,'  upper = ',upper
		print,' NOTE: LOWER AND UPPER ARE IN ARCSECONDS'
		print,'      10 pixels is ~1.2"'
	        print,'Recenter extraction? (y/n/abort)'
	        redo = ''
	        read,': ',redo
	        if redo eq 'y' then begin
			print,'Please click where you would like to center the extraction.'
		        while !MOUSE.button ne 1 do cursor,xloc,yloc,/data,wait=3
		        result_y = xloc
			result_y = where(max(stack2(result_y-5:result_y+5)) eq stack2[*])
		        wait,0.3
		        !MOUSE.button = 0
	                GOTO,LeSigh
	        endif else if redo eq 'n' then begin
	                print,'ReSize extraction? (y/n/bye)'
	                redo = ''
	                read,': ',redo
	                if redo eq 'y' then begin
	                        Print,'New Lower:'
	                        readnum,lower
				lower = abs(lower)
	                        Print,'New Upper:'
	                        readnum,upper
				upper = abs(upper)
	                        GOTO,LeSigh
	                endif else if redo eq 'n' then begin
                                extract_obj.object = object
	                        extract_obj.slit = ap
	                        extract_obj.center = result_y
	                        extract_obj.lower = lower
	                        extract_obj.upper = upper

                                pos_shift = result_y * platescale
                                if pos_shift LE slits[0].length/2.D then $
                                        hyp = slits[0].length/2.D - pos_shift else $
                                        hyp = pos_shift - slits[0].length/2.D
                                delRA = hyp * sin(slits[0].posang * !dpi/180) * $
                                        cos(extract_obj.dec_deg * !dpi/180) / 3600.D
                                delDEC = hyp * cos(slits[0].posang * !dpi/180) / 3600.D
                                if pos_shift LE slits[0].length/2.D then begin
                                        objra_deg = extract_obj.ra_deg - delRA
                                        objdec_deg = extract_obj.dec_deg - delDEC
                                endif else begin
                                        objra_deg = extract_obj.ra_deg + delRA
                                        objdec_deg = extract_obj.dec_deg + delDEC
                                endelse
				extract_obj.ra_deg = objra_deg
				extract_obj.dec_deg = objdec_deg
				radec,objra_deg,objdec_deg,ihr,imin,xsec,ideg,imn,xsc
				extract_obj.ra_sex = strtrim(string(ihr),2)+':' $
					+strtrim(string(imin),2)+':'+strtrim(string(xsec),2)
				extract_obj.dec_sex = strtrim(string(ideg),2)+':' $
					+strtrim(string(imn),2)+':'+strtrim(string(xsc),2)
	                        extract = [[extract],[extract_obj]]
	                        object ++
	                        Ok = 1
        	        endif else if redo eq 'bye' then begin
	                        stop
	                        endif else print,'I do not understand ',redo
                endif else if redo eq 'abort' then begin
                        stop
	        endif else begin
	                print,'I do not understand ',redo
	        endelse
	endwhile

        ;read in the exptime and airmass from the headers
        time = strcompress(sxpar(hdr[*, 0], 'EXPTIME'), /rem)
        airmass = strcompress(sxpar(hdr[*, 0], 'AIRMASS'), /rem)

        ;Rebin the spectrum to the region covered by the standards in actuallity
        x_specrebin,wav1d,spec1d,new_wave,corflux,/FLAMBDA
	x_specrebin,wav1d,var1d,new_wave,corvar,/FLAMBDA
	x_specrebin,wav1d,sky1d,new_wave,corsky,/FLAMBDA
	x_specrebin,wav1d,flat1d,new_wave,corFlat,/FLAMBDA

        ;The extinction correction is given by the factor
        extcor = 10. ^ (0.4 * airmass * ext[*])

	corflux2 = corflux/(corFlat/mean(corFlat))
	std_spectrum = {wave: new_wave, flux: corflux, ivar: corvar, spec_flat: corFlat, flux2: corflux2}
	if (channel_code eq 'm1r') or (channel_code eq 'm2r') then BEGIN
		BALM_MASK_WID=10.
		NRESLN=7.
	endif
	if (channel_code eq 'm1b') or (channel_code eq 'm2b') then begin
		BALM_MASK_WID=8
		NRESLN=7.5
	endif
        name_len=STRLEN(scifile)
        name = STRMID(scifile, 0, name_len-8) + 'ssf.fits'

;plot,std_spectrum.wave,std_spectrum.flux,/ylog,yrange=[1e4,2e5],/ystyle
;oplot,std_spectrum.wave,corFlat*(median(std_spectrum.flux)/median(corFlat)),color=cgcolor('red')
;oplot,std_spectrum.wave,std_spectrum.flux/(corFlat/mean(corFlat))/3.,color=cgcolor('dodger blue')
;print,'Pause for plots...'
;stop

	sensfunc = mods_sensfunc(std_spectrum,name $
		, SCHHDR = hdr, std_name = std_name $
		, /MSK_BALM, BALM_MASK_WID=BALM_MASK_WID $
		, NRESLN=NRESLN,airmass = airmass,exptime = time $
		, no_tell_mask=no_tell_mask,alt=alt_sensfunc)

	sensfunc = bspline_valu(new_wave,sensfunc)
	alt_sensfunc = bspline_valu(new_wave,alt_sensfunc)

	calflux = 1.0d-17*(corflux*extcor)/(time) * 10.^(sensfunc/2.5)
	calflux2 = 1.0d-17*(corflux2*extcor)/(time) * 10.^(alt_sensfunc/2.5)

	calvar  = 1.0d-17*(corvar)/time *10.^(sensfunc/2.5)
	calsky  = 1.0d-17*(corsky)/time *10.^(sensfunc/2.5) 

	print,'Saving aperture: ',ap
        ;save this aperture
        objspec = [[objspec],[spec1d]]
        objwave = [[objwave],[wav1d]]
	varspec = [[varspec],[var1d]]
	skyspec = [[skyspec],[sky1d]]

        objspec_cal = [[sensfunc],[calflux]]
	varspec_cal = [[varspec_cal],[calvar]]
	skyspec_cal = [[skyspec_cal],[calsky]]
	flatspec_cal = [corFlat]

        ; Save the output
        mwrfits,objspec,out,hdr,/create
	mwrfits,varspec,out
	mwrfits,skyspec,out
	mwrfits,objwave,out
        mwrfits,extract,out

        mwrfits,objspec_cal,outcal,hdrcal,/create
	mwrfits,varspec_cal,outcal
	mwrfits,skyspec_cal,outcal
	mwrfits,extract,outcal
	mwrfits,flatspec_cal,outcal

endfor

end
