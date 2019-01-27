pro mods_apply_new_sensfunc,scifile,cal,outname=outname,clobber=clobber

; Science image
verify = file_test(scifile)
if verify then sciimg = mrdfits(scifile,0,hdr) $
        else begin
                print,'You need a valid science image'
                print,scifile,'does not exist. Check your Science folder'
                stop
        endelse
skyimg = mrdfits(scifile,2,hdrsky)
varimg = mrdfits(scifile,1,headvar)
wavimg = mrdfits(scifile,3,head3)
extractstruct = mrdfits(scifile,4,head4)

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
if (channel_code eq 'm1b') or (channel_code eq 'm2b') then sxaddpar,hdrcal,'CRVAL1',3000
sxaddpar,hdrcal,'CDELT1',0.5

if (channel_code eq 'm1r') or (channel_code eq 'm2r') then new_wave = 0.5*indgen(9001) + 5500
if (channel_code eq 'm1b') or (channel_code eq 'm2b') then new_wave = 0.5*indgen(5601) + 3000

sensfunc = bspline_valu(new_wave,cal)
ext = bspline_valu(new_wave,ext_coef)
time = strcompress(sxpar(hdr[*, 0], 'EXPTIME'), /rem)
airmass = strcompress(sxpar(hdr[*, 0], 'AIRMASS'), /rem)


ndims = size(sciimg,/dim)
ndo = ndims[1]

objspec_cal = sensfunc
varspec_cal = fltarr(size(new_wave,/dim))
skyspec_cal = fltarr(size(new_wave,/dim))

for ii=1,ndo-1 do begin
	wav1d = wavimg[*,ii]
	spec1d = sciimg[*,ii]
	var1d = varimg[*,ii]
	sky1d = skyimg[*,ii]

        ;Rebin the spectrum to the region covered by the standards in actuallity
        x_specrebin,wav1d,spec1d,new_wave,corflux,/FLAMBDA
        x_specrebin,wav1d,var1d,new_wave,corvar,/FLAMBDA
        x_specrebin,wav1d,sky1d,new_wave,corsky,/FLAMBDA

        ;The extinction correction is given by the factor
        extcor = 10. ^ (0.4 * airmass * ext[*])

        ;Flux the spectrum including the LBT extinction curve
        calflux = (corflux*extcor)/(time) * 10.^(sensfunc/2.5)
        calvar = (corvar*extcor)/(time) * 10.^(sensfunc/2.5)
        calsky = (corsky*extcor)/(time) * 10.^(sensfunc/2.5)

        objspec_cal = [[objspec_cal],[calflux]]
        varspec_cal = [[varspec_cal],[calvar]]
        skyspec_cal = [[skyspec_cal],[calsky]]
endfor

; Save the output
mwrfits,objspec_cal,outcal,hdrcal,/create
mwrfits,varspec_cal,outcal
mwrfits,skyspec_cal,outcal
mwrfits,extractstruct,outcal

end
