;+
; NAME:
;   mods_quickreduce
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
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
;   
; REVISION HISTORY:
;   03-Dec-2013  KVC
;-  


;-----------------------------------------------------------------------------
; The code that follows is the science frame reduction code

PRO mods_quickreduce, filename $
          , apertures=apertures $
	  , calwave = calwave $
	  , noclean=noclean
;options to add later
; slitfile = slitfile, wavefile = wavefile $
;          , biasfile = biasfile, pixflatfile = pixflatfile $
;          , illumflatfile = illumflatfile $

if n_elements(filename) eq 0 then begin
	doc_library,'mods_quickreduce'
	stop
endif
  
t0 = systime(1)
;;----------
;; Read the raw science image
scihdr = xheadfits(filename)
if (size(scihdr, /tname) NE 'STRING') then begin
     splog, 'Invalid FITS header for file ', filename
     sciimg = 0
     sciivar = 0
     return
endif

;;  What telescope are we using?
telescope = strcompress(sxpar(scihdr[*, 0], 'TELESCOP'), /rem)
gratname  = strcompress(sxpar(scihdr[*, 0], 'GRATNAME'), /rem)
channel   = strcompress(sxpar(scihdr[*, 0], 'CHANNEL'), /rem)
mms       = strcompress(sxpar(scihdr[*, 0], 'MASKINFO')+'.mms', /rem)
detector  = strcompress(sxpar(scihdr[*, 0], 'DETECTOR'), /rem)
object    = strcompress(sxpar(scihdr[*, 0], 'object'), /rem)

if STRMID(mms,11,12,/REVERSE_OFFSET) eq 'LongSlit.mms' then mms = '1.0arcsecSegmentedLongSlit.mms'

;;----------
;; Set up the calibration library files
; pixel flat
redpix = 'rpixFlat_current.fits'
bluepix = 'bpixFlat_current.fits'

if (channel eq 'RED') then begin
   spawnline = 'modsProc -b ' + $
	filename + ' ' + GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/Calib_Lib/' + redpix
   SPAWN, spawnline
endif else if (channel eq 'BLUE') then begin
   spawnline = 'modsProc -b ' + $
        filename + ' ' + GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/Calib_Lib/' + bluepix
   SPAWN, spawnline	
endif else begin
   print,channel,' is an Unknown Channel, pixflats do not exist.'
   stop
endelse

; update the filename to the modsProc output
name_len=STRLEN(filename)
filename = STRMID(filename, 0, name_len-5)+ '_otf.fits'
sciimg = xmrdfits(filename, 0, hdr, /fscale)

; create a slit structure
;if STRMID(mms,0,30) eq '1.0arcsecSegmentedLongSlit.mms' then begin
if mms eq '1.0arcsecSegmentedLongSlit.mms' then begin
   print,'-----------------------------------------------------------------'
   print,'You are using a MODS longslit, I will create an MMS file for you.'
   print,'-----------------------------------------------------------------'
   mkmms_ls,filename
endif

;if STRMID(mms,0,36) eq '5arcsecSpectrophotometryWideSlit.mms' then begin
if mms eq '5arcsecSpectrophotometryWideSlit.mms' then begin
   print,'-----------------------------------------------------------------'
   print,'You are using the MODS wide-slit, I will create an MMS file for you.'
   print,'-----------------------------------------------------------------'
   mkmms_ws,filename
   apertures = [1]
endif

verify = file_test(mms)
if not verify then begin
   print,'I need the MMS file to run properly.'
   print,'You are using a MOS field.'
   print,'Please place ',mms,' in this directory and re-run mods_quickreduce.pro'
   stop
endif

name_len=STRLEN(filename)
slitfile = STRMID(filename, 0, name_len - 9)+ '_slits.fits' ;assume a standard naming convention
mods_mask,filename,slitfile,INTERACTIVE_SLITS=0,AUTOTUNE_SLITS=0

; make a 2D wavelength image
wavefile = STRMID(filename, 0, name_len-9)+ '_wave.fits' ;assume a standard naming convention

;for speed hardcode in the LS wave solutions
if keyword_set(calwave) then begin
	if ((mms eq '1.0arcsecSegmentedLongSlit.mms') OR $
	   (mms eq '5arcsecSpectrophotometryWideSlit.mms')) then begin
		if (channel eq 'RED') then begin
			wavefile = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/Calib_Lib/' + 'wave_LS_red.fits'
			wave = mrdfits(wavefile,0,whd)
		endif else if (channel eq 'BLUE') then begin
			wavefile = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/Calib_Lib/' + 'wave_LS_blue.fits'
			wave = mrdfits(wavefile,0,whd)
		endif else begin
		   print,channel,' is an Unknown Channel, no wavelength info known.'
		   stop
		endelse
	endif else begin
		print,'You are processing data from MOS observations.'
		print,'Standard wave calibration files only exist for LS mode'
		stop
	endelse
endif else begin
	wave = mods_empwave(slitfile,filename=wavefile, apertures=apertures)
endelse

wave=float(wave)

;create a variance image
rnoise = 2.5
sciivar = 0.0 * sciimg
sciivar = 1.0/(abs(sciimg - sqrt(2.0)*rnoise) +rnoise^2)

; XIDL works with transposed images
sciimg = transpose(sciimg)
sciivar = transpose(sciivar)

; Save the pre-Skysub images
name_len=STRLEN(filename)
presub_outname = STRMID(filename, 0, name_len-9)+ '_otfw.fits' ;assume a standard naming convention
header_line = 'mods_quickreduce beta test version 2014-02-18'
sxaddpar,hdr,'HISTORY',header_line
if (channel eq 'RED') then header_line = 'mods_quickreduce - pixel flat used - '+redpix $
	else header_line = 'mods_quickreduce - pixel flat used - '+bluepix
sxaddpar,hdr,'HISTORY',header_line
if keyword_set(calwave) then header_line = 'mods_quickreduce - wavelength solution from Calib_lib' $
	else header_line = 'mods_quickreduce - wavelength solution from sudoku solution'
sxaddpar,hdr,'HISTORY',header_line
if keyword_set(calwave) then header_line = wavefile
sxaddpar,hdr,'HISTORY',header_line

mwrfits,sciimg,presub_outname,hdr,/create
mwrfits,sciivar,presub_outname
mwrfits,wave,presub_outname

;skysub
name_len=STRLEN(filename)
sub_outname = STRMID(filename, 0, name_len-9)+ '_x2d.fits' ;assume a standard naming convention
mods_quickskyfit2d,presub_outname,wavefile,slitfile, outname=sub_outname, color=channel, apertures=apertures

;clean up
if not keyword_set(noclean) then begin
	spawnline='rm -f '+presub_outname
	spawn,'rm -f coords.dat'
	spawn,'rm -f fit.dat'
	if not keyword_set(calwave) then begin
		spawn,spawnline
        	spawnline='rm -f '+wavefile
	endif
        spawn,spawnline
        spawnline='rm -f '+slitfile
        spawn,spawnline
        spawnline='rm -f '+filename
        spawn,spawnline
endif

;extract
;OLD EXTRACTION... USES OLD STYLE SENSFUNC FROM IRAF
;if (channel eq 'RED') then begin
;	respfile = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/Calib_Lib/' + 'resp_red_new.fits'
;endif else if (channel eq 'BLUE') then begin
;        respfile = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/Calib_Lib/' + 'resp_blue_new.fits'
;endif
;outname = '1dauto_' + object + '.fits'
;outcal = '1dauto_' + object + '_cal.fits'
;mods_quickextract,'skysub.fits',slitfile,wavefile,respfile,outname=outname,$
;	color=channel,apertures=apertures,outcal=outcal

print,'Elapsed time = ', systime(1)-t0, ' sec'
;print,'Launching mods_quickview:'
;print,"    mods_quickviewslide,'",sub_outname,"', color='",channel,"',apertures=",apertures
;mods_quickviewslide,sub_outname, color=channel,apertures=apertures

END
