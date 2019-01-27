
;=========================================================================
;  mods_skyfit2d -- 
;=========================================================================
pro mods_quickskyfit2d, scifile, wavefile, slitfile, color=color,$
	convbeam=convbeam, apertures=apertures, outname=outname

;deal with parameters that were not entered
if n_elements(outname) EQ '' then begin
	stub = STRSPLIT(scifile,'-',/extract)
	outname = '2Dautosub_' + stub[1] 
endif 

; Load the images
sciimg = mrdfits(scifile,0,scihdr)
varimg = mrdfits(scifile,1,varhdr)
orig = sciimg

; Wavelength image
waveimg = mrdfits(wavefile,0,whdr)

; Gaussian convolution?
if n_elements(convbeam) EQ '' then convbeam=[1.0,2.0]
print,'USING A ',convbeam,' PIXEL GAUSSIAN TO SMOOTH for sky subtraction'
sciimg = filter_image(sciimg,FWHM=convbeam,/ALL)
waveimg = filter_image(waveimg,FWHM=convbeam,/ALL)

header_line = 'mods_quickreduce - smoothed for skysub using ['$
	+strtrim(string(convbeam[0]),2)+','+strtrim(string(convbeam[1]),2)+']'
sxaddpar,scihdr,'HISTORY',header_line

;create a new file to fill with sky-subtracted data
fit = sciimg

dims=size(sciimg,/dim)
nx=dims[0]
ny=dims[1]
instrument = strcompress(sxpar(scihdr[*, 0], 'INSTRUME'), /rem)
mms = strcompress(sxpar(scihdr[*,0], 'MASKINFO')+'.mms', /rem)
grat = strcompress(sxpar(scihdr[*, 0], 'GRATNAME'), /rem)

;create the slit structures
slitimg = mrdfits(slitfile,0,slithd) 
tset_slits = mrdfits(slitfile,1,tset_head)
   dim_t = tset_slits[0].dims
   nslit=size(tset_slits.xx1,/dim)
   nslit=nslit[1]
   nx_t = dims[0]
   ny_t = dims[1]

   xx1 = tset_slits[0].xx1
   xx2 = tset_slits[1].xx2

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

        apmask = 0*sciimg
        spec1d = fltarr(ny)
	for iii=0,ny-1 do begin
		for iiii=xx1[iii,ap-1],xx2[iii,ap-1] do apmask[iiii,iii] = 1
	endfor
	slitpos = where(apmask eq 1)

	for iii=0,ny-1 do begin
		;blank thin strips along top and bottom
		delta = (xx2[iii,ap-1]-xx1[iii,ap-1])*0.1 <10.
		bot = xx1[iii,ap-1]-delta
		top = xx1[iii,ap-1]+delta
		apmask[bot:top,iii] = 0
		bot = xx2[iii,ap-1]-delta
		top = xx2[iii,ap-1]+delta
		apmask[bot:top,iii] = 0
		; blank the center
		delta = (xx2[iii,ap-1]-xx1[iii,ap-1])*0.25
		bot = xx1[iii,ap-1]+delta
		top = xx2[iii,ap-1]-delta
		apmask[bot:top,iii] = 0
	endfor

	; FIT THE SKY IN THE SELETCED APERTURE
	if (color eq 'RED') then begin
		print,'RED'
		sset=bspline_iterfit(waveimg[slitpos],sciimg[slitpos],$
		bkspace=1.0,nord=3,invvar=apmask[slitpos])
	endif else if (color eq 'BLUE') then begin
		print,'BLUE'
		sset=bspline_iterfit(waveimg[slitpos],sciimg[slitpos],$
                bkspace=0.90,nord=3,invvar=apmask[slitpos])
	endif else stop
 
	extract = sciimg[slitpos]
	extract_wave = waveimg[slitpos]
        yfit = bspline_valu(waveimg[slitpos],sset)

        subimg = sciimg
        subimg[slitpos] = extract-yfit
        fit[slitpos] = yfit
        sciimg = subimg

	; Save the output
	;0 - original OTF image
	;1 - sky-subtracted image
	;2 - variance image
	;3 - sky fit
	;4 - wavelength image
	;5 - slit image
	;6 - slit structure
	header_line = 'mods_quickreduce Ext 0 - OTF'
	sxaddpar,scihdr,'HISTORY',header_line
        header_line = 'mods_quickreduce Ext 1 - sky subtracted'
        sxaddpar,scihdr,'HISTORY',header_line
        header_line = 'mods_quickreduce Ext 2 - variance'
        sxaddpar,scihdr,'HISTORY',header_line
        header_line = 'mods_quickreduce Ext 3 - sky'
        sxaddpar,scihdr,'HISTORY',header_line
        header_line = 'mods_quickreduce Ext 4 - wavelength image'
        sxaddpar,scihdr,'HISTORY',header_line
        header_line = 'mods_quickreduce Ext 5 - slit mask'
        sxaddpar,scihdr,'HISTORY',header_line
        header_line = 'mods_quickreduce Ext 6 - Slit BinTable'
        sxaddpar,scihdr,'HISTORY',header_line

	;normailzed masking
	orig = orig>(-100)
	normmask = slitimg/slitimg
	subimg = (subimg * normmask)>(-100)
	varimg = varimg * normmask
	fit = fit * normmask < max(orig)
	
	mwrfits,transpose(orig),outname,scihdr,/create
	mwrfits,transpose(subimg),outname
	mwrfits,transpose(varimg),outname
        mwrfits,transpose(fit),outname
	mwrfits,transpose(waveimg),outname
	mwrfits,transpose(slitimg),outname
	mwrfits,tset_slits,outname
endfor

end
