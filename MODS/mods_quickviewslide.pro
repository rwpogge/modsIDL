;=========================================================================
; mods_quickview - View MODS spectra 2D spectra interactively
;=========================================================================

pro mods_quickviewslide,scifile,slitfile,outname=outname, $
	apertures=apertures, flexure=flexure,color=color,outcal=outcal

; Science image
  sciimg = mrdfits(scifile,1,hdr)
  varimg = mrdfits(scifile,2,headvar)
  waveimg =mrdfits(scifile,4,whdr)
  skyimg = mrdfits(scifile,3,hdrsky)
  tset_slits = mrdfits(scifile,6,tset_head)
  slitim = mrdfits(scifile,5,shdr)

  sciimg = transpose(sciimg)
  varimg = transpose(varimg)
  waveimg = transpose(waveimg)
  skyimg = transpose(skyimg)
  slitim = transpose(slitim)

;image size
  dims=size(sciimg,/dim)
  nx=dims[0] & nx_t = dims[0]
  ny=dims[1] & ny_t = dims[1]

;slit info
  dim_t = tset_slits[0].dims
  nslit=size(tset_slits.xx1,/dim)
  nslit=nslit[1]
  xx1 = tset_slits[0].xx1
  xx2 = tset_slits[1].xx2

  objspec = fltarr(ny_t)+1
  objwave = fltarr(ny_t)+1
  skyspec = fltarr(ny_t)+1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;      plot each slit on its own
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(apertures) EQ '' then apertures = indgen(nslit)+1
ndo = size(apertures,/dim)
ndo=ndo[0]
object=1

for ii=0,ndo-1 do begin
	ap=apertures[ii]
	;for now extract the whole aperture

        spec1d = fltarr(ny)
        sky1d = fltarr(ny)
        wav1d = fltarr(ny)

	 ; plot slit and traces
	screen = GET_SCREEN_SIZE()
	screen = 0.75*screen
;	window,1,retain=2,xsize=screen[0],ysize=screen[1]
       	ncolors = !D.Table_Size
        device,DECOMPOSED=0
        LoadCT,39
	objtrace = xx1[*,ap-1]
        im_min = min(xx1[*,ap-1])
        im_max = max(xx2[*,ap-1])
        slitimb = transpose(BytScl(sciimg[im_min:im_max,*],$
	min=-5,max=100*median(sciimg[im_min:im_max,*])>10,/nan))
	slit_dims = size(slitimb,/dim)
	slit_dim_x = slit_dims[0]
	slit_dim_y = slit_dims[1]

        for iii=0,ny-1 do begin
                mask = where(slitim[*,iii] eq ap)
                spec1d[iii] = total(sciimg[mask,iii] >(-1))
                wav1d[iii] = median(waveimg[mask,iii])
                sky1d[iii] = total(skyimg[mask,iii])
        endfor


        !mouse.button = 1
	xloc = 0 & yloc = 0
	wave_val = 0.0 & value = 0.0
        while (!MOUSE.button ne 4) do begin
		;2D image
		title='Aperture '+string(ap)

		main = widget_base(title=title,/column,xsize=0.7*screen[0],ysize=1.2*(size(slitimb))[2])
		wavebase = widget_base(main,/row,/frame)

		draw = widget_draw(main,x_scroll_size=0.7*screen[0],y_scroll_size=(size(slitimb))[2],$
			xsize=(size(slitimb))[1],ysize=0.9*(size(slitimb))[2],/scroll,retain=2)
		widget_control,main,/realize
		Widget_control,draw,get_value=drawid
		wset,drawID
		tvImage,slitimb
		xmanager,'draw_scroll',main,/no_block
stop

                slide_Image,slitimb,title=title,retain=2,yvisible=slit_dim_y<screen[1],xvisible=0.6*screen[0]
stop
		; plot the spectrum
		xyouts,0.04*screen[0],0.45*screen[1],'Wavelength: ',/device,charsize=2
		xyouts,0.5*screen[0],0.45*screen[1],'ADU: ',/device,charsize=2
		xyouts,0.1*screen[0],0.45*screen[1],Wave_val,/device,charsize=2
		xyouts,0.6*screen[0],0.45*screen[1],value,/device,charsize=2	
                xyouts,xloc,yloc,'X',charthick=3,color=cgcolor('black'),charsize=1.5

;                cursor,xloc,yloc,/data,/nowait
                cursor,xloc,yloc,/data,wait=3
;		if !MOUSE.button eq 1 then begin
			Plot,IndGen(slit_dim_x),IndGen(slit_dim_y),Position=[0.04,0.52,0.97,0.99],$
                        /nodata,Color=FCS2,xstyle=1,ystyle=1,charsize=1.5
			xyouts,xloc,yloc,'X'
			value = sciimg[yloc+im_min,xloc]
			Wave_val = waveimg[yloc+im_min,xloc]
		
			delta = objtrace[xloc] - min(objtrace)	
			new_trace = objtrace + yloc -delta	
			for iii=0,ny-1 do begin
			     wid = 15.
        		     spec1d[iii] = total(sciimg[(new_trace[iii]-wid):(new_trace[iii]+wid),iii] >(-3))
                	     wav1d[iii] = median(waveimg[(new_trace[iii]-wid):(new_trace[iii]+wid),iii])
                	     sky1d[iii] = total(skyimg[(new_trace[iii]-wid):(new_trace[iii]+wid),iii])
        		endfor

;		endif
	endwhile

	!mouse.button = 0
stop
        path = objtrace[*,ap-1] - im_min
        blue_orig = blue_orig - im_min
                xpix = indgen(8192)

                plot,xpix,pathb,/noerase,Position=[0.02,0.08,0.99,0.49],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin],linestyle=1,color=cgColor('Dodger Blue')
                plot,xpix,pathb-lower/0.12,/noerase,Position=[0.02,0.08,0.99,0.49],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin]
                plot,xpix,pathb+upper/0.12,/noerase,Position=[0.02,0.08,0.99,0.49],$
                        xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,bmax-bmin]



        spec1d = fltarr(ny)
        sky1d = fltarr(ny)
	wav1d = fltarr(ny)
	for iii=0,ny-1 do begin
		mask = where(slitim[*,iii] eq ap)
		spec1d[iii] = total(sciimg[mask,iii] >(-1))
		wav1d[iii] = median(waveimg[mask,iii])
		sky1d[iii] = total(skyimg[mask,iii])
	endfor

        ;read in the exptime and airmass from the headers
        time = strcompress(sxpar(hdr[*, 0], 'EXPTIME'), /rem)
        airmass = strcompress(sxpar(hdr[*, 0], 'AIRMASS'), /rem)

        ;Rebin the spectrum to the region covered by the standards in actuallity
        x_specrebin,wav1d,spec1d,new_wave,corflux,/FLAMBDA
        x_specrebin,wav1d,sky1d,new_wave,corsky,/FLAMBDA

        ;The extinction correction is given by the factor
        extcor = 10. ^ (0.4 * airmass * ext[*])

        ;Flux the spectrum including the LBT extinction curve
        calflux = (corflux*extcor)/(time*10^(cal/2.5))
	calsky = (corsky*extcor)/(time*10^(cal/2.5))

	print,'Saving aperture: ',ap
        ;save this aperture
        objspec = [[objspec],[spec1d]]
        objwave = [[objwave],[wav1d]]
	skyspec = [[skyspec],[sky1d]]

        objspec_cal = [[objspec_cal],[calflux]]
	skyspec_cal = [[skyspec_cal],[calsky]]

if (color eq 'RED') then xrange=[5700,9999]
if (color eq 'BLUE') then xrange=[3600,5500]
plot,objwave[*,ap],objspec[*,ap],xstyle=1,xrange=xrange
oplot,objwave[*,ap],skyspec[*,ap],color=cgColor('Dodger Blue')
plot,new_wave,objspec_cal[*,ap],xstyle=1,xrange=xrange,Charsize=2,$
	xtitle='Wavelength [A]',ytitle='Flux',title='Automated Quicklook Spectrum'
oplot,new_wave,skyspec_cal[*,ap],color=cgColor('Dodger Blue')

        ; Save the output
        mwrfits,objspec,outname,hdr,/create
	mwrfits,skyspec,outname
	mwrfits,objwave,outname

        mwrfits,objspec_cal,outcal,hdrcal,/create
	mwrfits,skyspec_cal,outcal
endfor
end
