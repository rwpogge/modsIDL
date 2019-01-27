; NAME:
;   mods_wavesolve
;
; PURPOSE:
;  Full processing of an arc image.  This includes fitting the lines
;  and generating a wavelength image.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;  wavefile -- Name of output file for wavelength image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;  long_proc
;  long_waveimg
;   
; REVISION HISTORY:
;   10-Mar-2005  Written by J. Hennawi (UCB) and Scott Burles (MIT)
;   2013         Re-Written for MODS data K. Croxall
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
pro mods_wavesolve, filename, wavefile, slitfile = slitfile $
                    , biasfile = biasfile, pixflatfile = pixflatfile $
                    , verbose = verbose, LINELIST=linelist, CHK = CHK $
                    , REID_FILE= reid_file, BIN_RATIO=bin_ratio $
                    , ARC_INTER=arc_inter, PROCARC=procarc, TWEAK_ARC=tweak_arc

t0 = systime(1)

splog, 'Computing wavelength solution from file ', filename
splog, 'Using linelist ', linelist
;;----------
;; Read the arc image
if not keyword_set(PROCARC) then begin
;    long_proc, filename, arcimg, arcivar, hdr = hdr $
;               , biasfile = biasfile, verbose = verbose, pixflatfile = pixflatfile
	hdr = xheadfits(filename)
	if (size(hdr, /tname) NE 'STRING') then begin
	    splog, 'Invalid FITS header for file ', filename
	    flux = 0
	    invvar = 0
	    return
	endif
	arcimg = xmrdfits(filename, 0, hdr, /fscale,/silent)
	rnoise = 2.5
	arcivar = 0.0 * arcimg
	arcivar = 1.0/(abs(arcimg - sqrt(2.0)*rnoise) +rnoise^2)
	arcimg = transpose(arcimg)
	arcivar = transpose(arcivar)
endif else begin
    arcimg = xmrdfits(filename, 0 ,/silent)
    arcivar = xmrdfits(filename, 1, hdr,/silent)
;    restore, filename
endelse

dims = size(arcimg, /dimen) 
nx = dims[0]
ny = dims[1]

;; Parse the header information to set paramters for the wavelength structure
wstruct = mods_wstruct(hdr, LINELIST = linelist, REID_FILE = reid_file $
                       , BIN_RATIO = bin_ratio)

qafile = repstr(wavefile, '.fits', '.ps')
savefile = repstr(wavefile, '.fits', '.sav')

;----------
; Read in slit structure 
tset_slits = xmrdfits(slitfile, 1, silent = (keyword_set(verbose) EQ 0))

;----------
; Compute wavelength solution
xfit = mods_waveimg(arcimg, arcivar, tset_slits, wstruct, savefile $
                    , fwhmset = fwhmset, qafile = qafile, arc_inter=arc_inter,$
                    tweak_arc=tweak_arc)

; Now construct pixel wavelengths and apply wavelength solutions
pixset = mods_wavepix(arcimg, tset_slits, fwhm = fwhmset.MEDIAN $
                      , box_radius = wstruct.radius $
                      , sig_thresh = wstruct.sig_wpix $
                      , pkwdth = wstruct.pkwdth $
                      , TOLER = wstruct.TOLER, CHK = CHK)
piximg = mods_wpix2image(pixset, tset_slits, XFIT = xfit $
                         , waveimg = waveimg)
;--------------                                
;;  Write output to wavefile
splog, 'Writing output file'
sxdelpar, hdr, 'XTENSION'
sxdelpar, hdr, 'NAXIS2'
sxdelpar, hdr, 'NAXIS1'
sxdelpar, hdr, 'NAXIS'
sxdelpar, hdr, 'O_BZERO'
sxaddpar, hdr, 'BITPIX', -64

mwrfits, waveimg, wavefile, hdr, /create 
mwrfits, pixset, wavefile 
mwrfits, fwhmset, wavefile 

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

return
end
;------------------------------------------------------------------------------
