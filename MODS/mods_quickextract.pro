;=========================================================================
; mods_extract1d - extract and calibrate 1D dual channel MODS spectra
;			from 2D spectra
;=========================================================================
pro mods_quickextract,scifile,slitfile,wavefile,calfile,outname=outname, $
	apertures=apertures, flexure=flexure,color=color,outcal=outcal

;variance spectrum hard coded for now.... add into the skysub as another layer
skyimg = mrdfits(scifile,2,hdrsky)
varimg = mrdfits(scifile,1,headvar)
varimg = skyimg * varimg
; Science image
sciimg = mrdfits(scifile,0,hdr)
hdrcal = hdr
waveimg = mrdfits(wavefile,0,whdr)
tset_slits = mrdfits(slitfile,1,tset_head)
slitim = mrdfits(slitfile,0,shdr)
; Calbration image
cal = mrdfits(calfile,0,chd)

; set the outname
if n_elements(outname) EQ '' then begin
	stub = STRSPLIT(scifile,'_',/extract)
	out = '1Dauto_' + stub[1] + '.fits'
	outcal = '1Dauto_' + stub[1] + '_cal.fits'
endif

; set re-binned cal-wavelength
sxaddpar,hdrcal,'CUNIT1','Angstrom'
sxaddpar,hdrcal,'CTYPE1','Linear'
sxaddpar,hdrcal,'CRPIX1',1
sxaddpar,hdrcal,'CDELT1',0.5
new_wave = fltarr(n_elements(cal))
if (color eq 'RED') then begin
	sxaddpar,hdrcal,'CRVAL1',5500
	new_wave[0] = 5500
endif
if (color eq 'BLUE') then begin
	sxaddpar,hdrcal,'CRVAL1',3200
	new_wave[0] = 3200
endif
for i=1,n_elements(cal)-1 do new_wave[i] = new_wave[i-1]+0.5

;set the response function
sset=bspline_iterfit(new_wave,cal,bkspace=6,nord=3)
cal = bspline_valu(new_wave,sset)
cal_coef = bspline_iterfit(new_wave,cal,bkspace=0.65,nord=3)
;extinction correction
extfile = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/mods_extiction.dat'
readcol,extfile,ext_wave,ext_val,/silent
ext_coef = bspline_iterfit(ext_wave,ext_val,bkspace=3,nord=3)
ext = bspline_valu(new_wave,ext_coef)

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

objtrace = xx1*0
objspec = fltarr(ny_t)+1
objwave = fltarr(ny_t)+1
skyspec = fltarr(ny_t)+1
objspec_cal = cal
skyspec_cal = cal
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;       fit each slit on its own
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(apertures) EQ '' then apertures = indgen(nslit)+1
ndo = size(apertures,/dim)
ndo=ndo[0]
object=1

print,apertures
for ii=0,ndo-1 do begin
	ap=apertures[ii]
	;for now extract the whole aperture
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
