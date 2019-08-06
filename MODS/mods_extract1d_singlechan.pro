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
pro extractionplot,w,g, $
	st,x,st2,y, $
	sp,lower,upper

common exclude1, blank_low, blank_hi,  platescale

;Plot up the centering and extractions
;bottom panel
x0 = w[g[0]]-20
x1 = max(w[g])+20
plot,w[g],st>(-5),Position=[0.14,0.05,0.99,0.28],XRANGE=[x0,x1],xstyle=1,$
	ytitle='Av Spectrum (entire slit)',charsize=1.9
line_x = [w[x],w[x]]
line_x1 = [w[x-9],w[x-9]]
line_x2 = [w[x+10],w[x+10]]
minval = !Y.CRANGE[0]
maxval = !Y.CRANGE[1]
line_y = [minval,maxval]
oplot,line_x,line_y,Color=cgColor('dodger blue')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')

;middle pannel
stacksz = size(st2,/dim)
st2 = st2 +10.D
plot,indgen(stacksz),st2,XRANGE=[0,stacksz],$
	position=[0.14,0.35,0.99,0.62],ytitle='Slit Profile',$
        /noerase,yrange=[min(st2)>(0.1),max(st2)],xstyle=1, $
        charsize=1.9,ystyle=1,/ylog
line_x = [y,y]
minval = 10^(!Y.CRANGE[0])
maxval = 10^(!Y.CRANGE[1])
line_y = [minval,maxval]
line_x1 = [y+upper/platescale,y+upper/platescale]
line_x2 = [y-lower/platescale,y-lower/platescale]
niter = size(blank_lowb,/dim)
niter = niter[0]-1
for i=0,niter do polyfill,[blank_low[i]>0,blank_low[i]>0,blank_hi[i]<stacksz,blank_hi[i]<stacksz], $
	[minval,maxval,maxval,minval],color=cgColor('firebrick'),/data
plot,indgen(stacksz),st2,XRANGE=[0,stacksz],/ylog,$
        position=[0.14,0.35,0.99,0.62],charsize=1.9,ystyle=1,$
        /noerase,yrange=[min(st2)>(0.1),max(st2)],xstyle=1,xtitle='Pixels'
oplot,line_x,line_y,Color=cgColor('dodger blue')
oplot,line_x1,line_y,Color=cgColor('beige')
oplot,line_x2,line_y,Color=cgColor('beige')
axis,xaxis=1,XRANGE=[0,stacksz*0.12],xtitle='ARCSEC'

;top panel
plot,w,sp,position=[0.14,0.72,0.99,0.99],/noerase,XRANGE=[x0-100,x1+50], $
	ytitle='Extracted Spectrum', charsize=1.9,ystyle=1;,YRANGE=[200,900]

return
end

;=========================================================================
; mods_extract1d_singlechan - extract and calibrate 1D single channel MODS spectra
;			from 2D spectra
;=========================================================================
pro mods_extract1d_singlechan,scifile, $
	cal=cal, z = z, clobber=clobber, $
	apertures=apertures,outname=outname,wave=wave, $
	centerLine=centerLine, trim_b = trim_b, trim_t = trim_t,$
	noPlotAp=noPlotAp, slits = slits,force=force, $
	ill=ill, illumination_corr = illumination_corr

if not keyword_set(cal) then cal = 'sensfunc.fits'

if (n_elements(ill) EQ '') and (keyword_set(illumination_corr)) then begin
	ill = mrdfits('Proc/ill.fits',0,ilhd)
	ill = transpose(ill)
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

if n_elements(trim_t) EQ 0 then trim_t = 1
if n_elements(trim_b) EQ 0 then trim_b = 1

if n_elements(centerLine) EQ '' then begin
	if keyword_set(red) then begin
		low_cent = 6650*z
		hi_cent = 6800*z
	endif
	if keyword_set(blue) then begin
		low_cent = 4890*z
		hi_cent =  5070*z
	endif
endif else begin
        low_cent = centerLine[0]*z
        hi_cent =  centerLine[1]*z
endelse

common exclude1, blank_low, blank_hi, platescale
blank_low = [0]
blank_hi = [0]

; Science image
verify = file_test(scifile)
if verify then sciimg = mrdfits(scifile,0,hdr) $
        else begin
                print,'You need a valid science image'
                print,scifile,'does not exist. Check your Science folder'
                stop
        endelse
skyimg = mrdfits(scifile,3,hdrsky)
varimg = mrdfits(scifile,1,headvar)
varimg = skyimg * varimg

; Get Mask and Instrument Info
instrument = strcompress(sxpar(hdr[*,0], 'INSTRUME'), /rem)
if instrument eq 'MODS1R' then channel_code = 'm1r' $
else if instrument eq 'MODS2R' then red_channel_code = 'm2r' $
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
        if (channel_code eq 'm1b') or (channel_code eq 'm2b') then low_cent = 4850*z
        if (channel_code eq 'm1b') or (channel_code eq 'm2b') then hi_cent = 5070*z
        if (channel_code eq 'm1r') or (channel_code eq 'm2r') then low_cent = 6650*z
        if (channel_code eq 'm1r') or (channel_code eq 'm2r') then hi_cent = 6800*z
endif else begin
        low_cent = centerLine[0]*z
        hi_cent = centerLine[1]*z
endelse

; Out file
if n_elements(outname) EQ '' then begin
        stub = STRSPLIT(scifile,'_',/extract)
        out = 'Science/1D_' + channel_code + '_' + stub[2] + '.fits'
        outcal = 'Science/1D_' + channel_code + '_' + stub[2] + '_cal.fits'
endif else begin
        out = 'Science/1D_' + channel_code + '_' + outname + '_cal.fits'
        outcal = 'Science/1D_' + channel_code + '_' + outname + '_cal.fits'
endelse
verify = file_test(out)
if (verify) and (not keyword_set(clobber)) then begin
        if verify then print,out,' already exists. If you want to overwite it please use /clobber'
        stop
endif

; Wavelength image
if ((n_elements(wave) EQ '') and ((channel_code eq 'm1b') or (channel_code eq 'm2b'))) then $
        wave = 'wave-'+ channel_code+'_HgXeAr_'+mask+'.fits'
if ((n_elements(wave) EQ '') and ((channel_code eq 'm1r') or (channel_code eq 'm2r'))) then $
        wave = 'wave-'+ channel_code+'_NeXeAr_'+mask+'.fits'
verify = file_test(wave)
if verify then waveimg = mrdfits(wave,0,whdr) $
        else begin
                print,'You need a valid wave-image'
                print,wave,' does not exist.'
                stop
        endelse

; Calibration images
cal = mrdfits(cal,1,chd)

; Extinction file
extfile = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/mods_extiction.dat'
readcol,extfile,ext_wave,ext_val,/silent
ext_coef = bspline_iterfit(ext_wave,ext_val,bkspace=3,nord=3)

;establish a uniform grid
hdrcal = hdr

sxaddpar,hdrcal,'CUNIT1','Angstrom'
sxaddpar,hdrcal,'CTYPE1','Linear'
sxaddpar,hdrcal,'CRPIX1',1
if (channel_code eq 'm1r') or (channel_code eq 'm2r') then sxaddpar,hdrcal,'CRVAL1',5500
if (channel_code eq 'm1b') or (channel_code eq 'm2b') then sxaddpar,hdrcal,'CRVAL1',3200
sxaddpar,hdrcal,'CDELT1',0.5

if (channel_code eq 'm1r') or (channel_code eq 'm2r') then new_wave = 0.5*indgen(9001) + 5500
if (channel_code eq 'm1b') or (channel_code eq 'm2b') then new_wave = 0.5*indgen(5201) + 3200

sensfunc = bspline_valu(new_wave,cal)
ext = bspline_valu(new_wave,ext_coef)

;open a window for plotting
window,1,retain=2;,XSize=1000,YSize=1000
if (not keyword_set(noPlotAp)) then window,2,retain=2

dims=size(sciimg,/dim)
nx=dims[0]
ny=dims[1]
instrument = strcompress(sxpar(hdr[*, 0], 'INSTRUME'), /rem)
mms = strcompress(sxpar(hdr[*,0], 'MASKINFO')+'.mms', /rem)

;create the slit structures
  if n_elements(slits) eq '' then slits = 'slits-'+channel_code+'_ill_'+mask+'.fits'
  verify = file_test(slits)
  if verify then slits = mrdfits(slits,1,tset_head) $
        else begin
                print,'You need a valid slit-image'
                print,slits,' does not exist.'
                stop
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
   illspec = fltarr(ny_t)+1


   objspec_cal = sensfunc
   varspec_cal = fltarr(size(new_wave,/dim))
   skyspec_cal = fltarr(size(new_wave,/dim))
   illspec_cal = fltarr(size(new_wave,/dim))


;extraction struct
extract = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D], $
	ra_deg: double(slits[0].tarra2), dec_deg: double(slits[0].tardec2), $
	ra_sex: slits[0].tarra, dec_sex: slits[0].tardec}

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
	extract_obj = {object: 0, slit: 0, center: 0, lower: 0.D, upper: 0.D, delta: 0.D, slope: [0.D,0.D],$
		ra_deg: slits[0].apra2[ap-1], dec_deg: slits[0].apdec2[ap-1], $
		ra_sex: slits[0].apra[ap-1], dec_sex: slits[0].apdec[ap-1]}
	blank_low = [0]
	blank_hi = [0]
	print,'APERTURE: ',ap
	print,''
	n_blank = 1
        center,sciimg,waveimg,ap,low_cent,hi_cent,$
		result_x,result_y,good,wav1d,$
		stack,stack2,xx1,xx2,trim_b,trim_t

        lower = 2.
        upper = 2.

	LeSigh: print,'Extracting Spectrum'
        spec1d = fltarr(ny)
        var1d = fltarr(ny)
	sky1d = fltarr(ny)
        spec1db = fltarr(ny)
        ill1d = fltarr(ny) + 1.D

	modsadcalc,scifile,wav1d[result_x],wav1d,xdelt,ydelt

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
		if keyword_set(illumination_corr) then $
			ill1d[iii] = median(ill[bot:top,iii]) 
	endfor
	; plot slit and traces
	if (not keyword_set(noPlotAp)) then begin
		wset,2
	        ncolors = !D.Table_Size
	        device,DECOMPOSED=0
	        LoadCT,39
	        min_ap = min(xx1[*,ap-1])+trim_b
	        max_ap = max(xx2[*,ap-1])-trim_t
	        slitim = transpose(BytScl(sciimg[min_ap:max_ap,*],min=-5,max=60*median(sciimg[min_ap:max_ap,*])>10,/nan))
	        TVImage,slitim,Position=[0.02,0.08,0.99,0.99],/erase
		path = OBJTRACE[*,ap-1] - min_ap
	        orig = orig - min_ap
		xpix = indgen(8192)

		plot,xpix,path,/noerase,Position=[0.02,0.08,0.99,0.99],xstyle=1,xrange=[0,8191],$
			ystyle=1,yrange=[0,max_ap-min_ap],linestyle=1,color=cgColor('Dodger Blue')
	        plot,xpix,path-lower/0.12,/noerase,Position=[0.02,0.08,0.99,0.99],$
			xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,max_ap-min_ap]
	        plot,xpix,path+upper/0.12,/noerase,Position=[0.02,0.08,0.99,0.99],$
			xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,max_ap-min_ap]
	endif

	;Plot up the centering and extractions
	wset,1
	extractionplot,wav1d,good,$
	       stack,result_x,$
	       stack2,result_y,$
	       spec1d,lower,upper

	Ok=0

	while Ok ne 1 do begin
	        print,'Centered on lines: ', wav1d[result_x]
	        print,'Peak found at: ',result_y
	        print,'Currently using: lower = ',lower,'  upper = ',upper
		print,''
		print,' NOTE: LOWER AND UPPER ARE IN ARCSECONDS'
		print,'      10 pixels is ~1.2"'
	        print,'Recenter extraction? (y/n/bye)'
	        redo = ''
	        read,': ',redo
	        if redo eq 'y' then begin
	                Print,'New center:'
	                read,': ',result_y
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
                endif else if redo eq 'bye' then begin
                        stop
	        endif else begin
	                print,'I do not understand ',redo
	        endelse
	endwhile

        ;read in the exptime and airmass from the headers
        time = strcompress(sxpar(hdr[*, 0], 'EXPTIME'), /rem)
        airmass = strcompress(sxpar(hdr[*, 0], 'AIRMASS'), /rem)

        ;Flexure correction
	wset,2
	if (channel_code eq 'm1b') or (channel_code eq 'm2b') then color = 'blue'
	if (channel_code eq 'm1r') or (channel_code eq 'm2r') then color = 'red'
	if keyword_set(force) then begin
	   mods_flexurecorr,wav1d,spec1d,sky1d,var1d,color,$
		z=z,/plot_prog,force=force,delta=delta,slope=slope
	endif else begin
	   mods_flexurecorr,wav1d,spec1d,sky1d,var1d,color,$
		z=z,/plot_prog,delta=delta,slope=slope
	endelse
	extract[object-1].delta = delta
	extract[object-1].slope = slope

        ;Rebin the spectrum to the region covered by the standards in actuallity
        x_specrebin,wav1d,spec1d,new_wave,corflux,/FLAMBDA
	x_specrebin,wav1d,var1d,new_wave,corvar,/FLAMBDA
	x_specrebin,wav1d,sky1d,new_wave,corsky,/FLAMBDA
        x_specrebin,wav1d,ill1d,new_wave,corill,/FLAMBDA

        sset=bspline_iterfit(new_wave,corill,bkspace=25,nord=3)
        corill = bspline_valu(new_wave,sset)

        ;The extinction correction is given by the factor
        extcor = 10. ^ (0.4 * airmass * ext[*])

        ;Flux the spectrum including the LBT extinction curve
        calflux = (corflux*extcor)/(time) * 10.^(sensfunc/2.5)
        calvar = (corvar*extcor)/(time) * 10.^(sensfunc/2.5)
	calsky = (corsky*extcor)/(time) * 10.^(sensfunc/2.5)

	print,'Saving aperture: ',ap
        ;save this aperture
        objspec = [[objspec],[spec1d]]
        objwave = [[objwave],[wav1d]]
	varspec = [[varspec],[var1d]]
	skyspec = [[skyspec],[sky1d]]

        objspec_cal = [[objspec_cal],[calflux]]
        varspec_cal = [[varspec_cal],[calvar]]
        skyspec_cal = [[skyspec_cal],[calsky]]
        illspec_cal = [[illspec_cal],[corill]]

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
	mwrfits,illspec_cal,outcal

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
print,'Average Flexure Shift: ', mean(extract[1:object-1].delta),'+/-', $
	stddev(extract[1:object-1].delta)/(n_elements(extract[1:object-1].delta))
print,'Average linear component: ', mean(extract[1:object-1].slope[0]),'+/-', $
        stddev(extract[1:object-1].slope[0])/(n_elements(extract[1:object-1].slope[0]))
print,'                          ', mean(extract[1:object-1].slope[1]),'+/-', $
        stddev(extract[1:object-1].slope[1])/(n_elements(extract[1:object-1].slope[1]))
print,'If extracting again, please consider using:  force=[',$
	strcompress(mean(extract[1:object-1].delta)),',',$
        strcompress(mean(extract[1:object-1].slope[0])),',', $
        strcompress(mean(extract[1:object-1].slope[1])),']'
print,'------------------------------------------------------------------'

end
