;+
; NAME: 
;    mods_skyfitprism
;
; PURPOSE:
;    Perform 2D sky subtraction on DUAL mode MODS PRISM data
;    using a 2D B-spline.  Allows for interactive masking of
;    slit edges.
;
; CALLING SEPUENCE:
;    mods_skyfitprism, red_image,blue_image,red_slits,blue_slits, $
;          apertures = apertures, wave_red = wave_red, wave_blue = wave_blue
;
; INPUTS:
;    red_image 	- red channel science image
;    blue_image	- blue channel science image
;    red_slits	- red slit structure file
;    blue_slits	- blue slit structure file
;
; OPTIONAL INUTS:
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
;    clobber   - Clobber existing images.  Otherwise they will not be overwritten.
;
; COMMENTS:
;
; EXAMPLE:
;    mods_skyfitprism,'Science/sci-mods1r.20130504.0034_otfc.fits.gz', $
;                     'Science/sci-mods1b.20130504.0025_otfc.fits.gz', $
;                     'slits-red.fits','slits-blue.fits', $
;                     apertures=[12], $
;                     wave_red='wave-mods1r.20130504.0034_otfc.fits',wave_blue='wave-mods1b.20130504.0025_otfc.fits'
;
; BUGS:
;
; PROCEDURES CALLED:
;    prism_fits (included)
;    bspline_iterfit
;    mrdfits
;    filter_image
;    mwrfits
;
; Author:
;    K. Croxall, OSU Astronomy Dept.
;    croxall@astronomy.ohio-state.edu
;    2014 Mar 7
;------------------------------------------------------------------------------

;=========================================================================
; prismfits - Run the fitting iteratively allowing for mask alteration
;=========================================================================
pro prismfits,sciimg,slitpos,waveimg,fit,bkspace,top,bot,ap,xx1,xx2,subimg,color=color

  fit_full = fit                        ; initalize a temporary fit frame

  ncolors = !D.Table_Size               ; load the color table
  device,DECOMPOSED=0
  LoadCT,39

  xx1_ap = xx1[*,ap-1]                  ; initialize the lower trace
  xx2_ap = xx2[*,ap-1]                  ; initialize the upper trace

  xloc = where(xx1_ap ne 0.0)           ; find the extent of the slit in the dispersion direction
  min = min(xx1_ap[xloc])+1 > 5.
  max = max(xx2_ap[xloc])-1 < 3087
  xpix = (indgen(size(xloc,/dim))) + min(xloc)


;full slit sub
  sset_full = bspline_iterfit(waveimg[slitpos],sciimg[slitpos],bkspace=bkspace,nord=2)
  extract=sciimg[slitpos]
  extract_wave=waveimg[slitpos]
  yfit_full = bspline_valu(waveimg[slitpos],sset_full)

  subimg_full = sciimg
  subimg_full[slitpos] = extract-yfit_full
  fit_full[slitpos] = yfit_full

  REPEAT BEGIN
        accept_subtraction = 0                  ;initialize flag
  	apmask = 0*sciimg                     ; create a mask for the skyfit
  	apmask[slitpos] = 1                   ; validate the region for the givn slit

   	;draw the full slit subtraction
	maxscale = max(subimg_full[min:max,min(xloc):max(xloc)])+5
	presubim = transpose(BytScl(subimg_full[min:max,min(xloc):max(xloc)],min=-5,max=maxscale,/nan))
	tvImage,presubim,Position=[0.02,0.53,0.99,0.73],/erase
	plot,xpix,xx1_ap[xloc]-min,/noerase,Position=[0.02,0.53,0.99,0.73],$
        	xstyle=1,xrange=[min(xloc),max(xloc)],ystyle=1,yrange=[0,max-min-1]
        oplot,xpix,xx1_ap[xloc]-min,color=cgcolor('red'),thick=2
        oplot,xpix,xx2_ap[xloc]-min,color=cgcolor('red'),thick=2
        xyouts,min(xloc)+10,10,'Full Slit Model Subtraction'

        ;draw the original data  
  	maxscale = max(sciimg[min:max,min(xloc):max(xloc)])
        slitim = transpose(BytScl(sciimg[min:max,min(xloc):max(xloc)],min=0,max=maxscale,/nan))
        tvImage,slitim,Position=[0.02,0.76,0.99,0.96],/noerase
        plot,xpix,xx1_ap[xloc]-min,/noerase,Position=[0.02,0.76,0.99,0.96],$
              xstyle=1,xrange=[min(xloc),max(xloc)],ystyle=1,yrange=[0,max-min-1]
        oplot,xpix,xx1[xloc]-min,color=cgcolor('red'),thick=2
        oplot,xpix,xx2[xloc]-min,color=cgcolor('red'),thick=2
        xyouts,min(xloc)+10,10,'Un-Subtracted Data'
	title = 'Aperture ' + strtrim(ap,2) + ' - ' + color
        xyouts,min(xloc)+10,20,title
        
        ;initial masking
        for j=min(xloc),max(xloc) do apmask[xx1_ap[j]-4:xx1_ap[j]+bot,j] = 0
        for j=min(xloc),max(xloc) do apmask[xx2_ap[j]-top:xx2_ap[j]+4,j] = 0
        oplot,xpix,xx1_ap[xloc]-min+bot,color=cgcolor('magenta'),linestyle=1,thick=2
        oplot,xpix,xx2_ap[xloc]-min-top,color=cgcolor('magenta'),linestyle=1,thick=2

        ;fit the masked slit  -- bkspace 1 - 5 seems best depending on quadrant issues
        sset=bspline_iterfit(waveimg[slitpos],sciimg[slitpos],bkspace=bkspace,nord=2,invvar=apmask[slitpos])
        extract=sciimg[slitpos]
        extract_wave=waveimg[slitpos]

        yfit = bspline_valu(waveimg[slitpos],sset)
        subimg = sciimg
        subimg[slitpos] = extract-yfit
        fit[slitpos] = yfit

        ;draw the sky fit
  	maxscale = max(fit[min:max,min(xloc):max(xloc)])
        fitim = transpose(BytScl(fit[min:max,min(xloc):max(xloc)],min=0,max=maxscale,/nan))
        tvImage,fitim,Position=[0.02,0.29,0.99,0.49],/noerase
        plot,xpix,xx1_ap[xloc]-min,/noerase,Position=[0.02,0.29,0.99,0.49],$
                xstyle=1,xrange=[min(xloc),max(xloc)],ystyle=1,yrange=[0,max-min-1]
        oplot,xpix,xx1_ap[xloc]-min,color=cgcolor('red')
        oplot,xpix,xx2_ap[xloc]-min,color=cgcolor('red')
        xyouts,min(xloc)+10,10,'Fit to masked sky'

        subim = transpose(BytScl(subimg[min:max,min(xloc):max(xloc)],min=-5,max=50,/nan))
        tvImage,subim,Position=[0.02,0.06,0.99,0.26],/noerase
        plot,xpix,xx1[xloc]-min,/noerase,Position=[0.02,0.06,0.99,0.26],$
                xstyle=1,xrange=[min(xloc),max(xloc)],ystyle=1,yrange=[0,max-min-1]
        oplot,xpix,xx1_ap[xloc]-min,color=cgcolor('red')
        oplot,xpix,xx2_ap[xloc]-min,color=cgcolor('red')
        xyouts,min(xloc)+10,10,'Masked Subtraction'

	print,'Currently using: top = ',strtrim(top,2),', bot = ',strtrim(bot,2),', bkspace = ',strtrim(bkspace,2)
	Print,'Try new Parameters? (y/n/bye)'
        redo = ''
        read,': ',redo
        if redo eq 'y' then begin
                Print,'New TOP mask size (Current value = ',top,'):'
                read,': ',top
                Print,'New BOTTOM maks size:(Current value = ',bot,')'
                read,': ',bot
                Print,'New BKSPACE (Current value = ',bkspace,', 1--5 currently recommended):'
                read,': ',bkspace
        endif else if redo eq 'bye' then begin
                stop
        endif else if redo eq 'n' then accept_subtraction = 1

  ENDREP UNTIL accept_subtraction

END

;=========================================================================
; mods_skyfitprism - Subtract sky emission from Dual mode MODS Prism spectra
;=========================================================================
pro mods_skyfitprism,red_scifile,blue_scifile, $
	red_slits=red_slits,blue_slits=blue_slits,$
	convbeam=convbeam, $
	apertures=apertures,outname=outname, $
	wave_blue=wave_blue,wave_red=wave_red, $
	boxcar = boxcar, $
	clobber = clobber

;deal with parameters that were not entered
if n_elements(outname) EQ '' then begin
	redstub = STRSPLIT(red_scifile,'-',/extract)
	redout = 'Science/2Dsub_' + redstub[1]
	bluestub = STRSPLIT(blue_scifile,'-',/extract)
	blueout = 'Science/2Dsub_' + bluestub[1]
endif else begin
	redout = 'Science/2Dsub_red_' + outname + '.fits'
	blueout = 'Science/2Dsub_blue_' + outname + '.fits'
endelse
verify_red = file_test(redout)
verify_blue = file_test(blueout)
if (verify_red or verify_blue) and (not keyword_set(clobber)) then begin
        if verify_red then print,redout,' already exists. If you want to overwite it please use /clobber'
        if verify_blue then print,blueout,' already exists. If you want to overwite it please use /clobber'
        stop
endif

; Load the images
red_sciimg = mrdfits(red_scifile,0,red_hdr)
blue_sciimg = mrdfits(blue_scifile,0,blue_hdr)
verify_red = file_test(red_scifile)
verify_blue = file_test(blue_scifile)
if verify_red then red_sciimg = mrdfits(red_scifile,0,red_hdr) $
        else begin
                print,'You need a valid RED science image'
                print,red_scifile,'does not exist. Check your Science folder'
                stop
        endelse
if verify_blue then blue_sciimg = mrdfits(blue_scifile,0,blue_hdr) $
        else begin
                print,'You need a valid BLUE science image'
                print,blue_scifile,'does not exist. Check your Science folder'
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
                print,'You need a valid RED wave-image'
                print,wave_red,'does not exist.'
                stop
        endelse
if verify_blue then blue_waveimg = mrdfits(wave_blue,0,blue_whdr) $
        else begin
                print,'You need a valid BLUE wave-image'
                print,wave_blue,'does not exist.'
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
;window,1,retain=2;,XSize=1000,YSize=1000
;if (not keyword_set(noPlotAp)) then window,2,retain=2

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

;create the slit structures 
  if n_elements(red_slits) eq '' then red_slits = 'slits-'+red_channel_code+'_ill_'+mask+'.fits'
  if n_elements(blue_slits) eq '' then blue_slits = 'slits-'+blue_channel_code+'_ill_'+mask+'.fits'
  verify_red = file_test(red_slits)
  verify_blue = file_test(blue_slits)
  if verify_red then red_slits_im = mrdfits(red_slits,0,r_slit_head) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid RED slit-image'
                print,red_slits,' does not exist.'
                print,"Use: red_slits='FILENAME'"
                stop
        endelse
   if verify_blue then blue_slits_im = mrdfits(blue_slits,0,b_slit_head) $
        else begin
                print,'------------------------------------------------------'
                print,'You need a valid BLUE slit-image'
                print,blue_slits,' does not exist.'
                print,"Use: blue_slits='FILENAME'"
                stop
        endelse
   red_tset_slits = mrdfits(red_slits,1,r_tset_head)
   blue_tset_slits = mrdfits(blue_slits,1,b_tset_head)


   dim_t = red_tset_slits[0].dims
   nslit=size(red_tset_slits.xx1,/dim)
   nslit=nslit[1]
   nx_t = dims[0]
   ny_t = dims[1]

   red_xx1 = red_tset_slits[0].xx1
   red_xx2 = red_tset_slits[1].xx2

   red_objtrace = red_xx1*0
   red_objspec = red_xx1*0
   red_objwave = red_xx1*0
   red_varmask = red_xx1*0

   blue_xx1 = blue_tset_slits[0].xx1
   blue_xx2 = blue_tset_slits[1].xx2

   blue_objtrace = blue_xx1*0
   blue_objspec = blue_xx1*0
   blue_objwave = blue_xx1*0
   blue_varmask = blue_xx1*0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                      fit each slit on its own
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(apertures) EQ '' then apertures = indgen(nslit)+1
ndo = size(apertures,/dim)
ndo=ndo[0]

for ii=0,ndo-1 do begin
	ap=apertures[ii]
	print,'APERTURE: ',ap
	print,''

	red_slitpos = where(red_slits_im eq ap)
	blue_slitpos = where(blue_slits_im eq ap)

	; FIT THE SKY IN THE SELETCED APERTURE
	bkspace = 3 & top = 5 & bot = 5
	prismfits,red_sciimg,red_slitpos,red_waveimg,red_fit,bkspace,top,bot,ap,$
		red_xx1,red_xx2,red_subimg,color='Red'

       prismfits,blue_sciimg,blue_slitpos,blue_waveimg,blue_fit,bkspace,top,bot,ap,$
                blue_xx1,blue_xx2,blue_subimg,color='Blue'

	red_sciimg = red_subimg
	blue_sciimg = blue_subimg

	; Save the output
		;0 - sky-subtract0d image
		;1 - variance image
		;3 - sky fit

	mwrfits,red_subimg,redout,red_hdr,/create
	mwrfits,red_varimg,redout
	mwrfits,red_fit,redout
	
	mwrfits,blue_subimg,blueout,blue_hdr,/create
	mwrfits,blue_varimg,blueout
	mwrfits,blue_fit,blueout
endfor

end
