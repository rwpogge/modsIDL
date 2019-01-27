;+
; NAME:
;   long_wstruct
;
; PURPOSE:
;   Create a structure which guides the data reduction process.  Parts
;   of the structure are initialized according to the instrument and
;   its configuration.
;
; CALLING SEQUENCE:
;
; INPUTS:
;   hdr  -- Image header
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;  Returns a structure describing the instrument configuration which
;  is used to guide the reduction steps.
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
;   10-Mar-2005  Written by S. Burles (MIT), David Schlegel (LBL), and
;                Joe Hennawi (UC Berkeley)
;   22-Jan-2019  MODS2 hooks [rwp/osu]
;-
;------------------------------------------------------------------------------
FUNCTION MODS_WSTRUCT, hdr, LINELIST=linelist, REID_FILE=REID_FILE $
                       , BIN_RATIO = bin_ratio

;;-----------
;;  Line list directory
calib_path = GETENV('LONGSLIT_DIR') + '/calib/linelists/'
line_path = GETENV('XIDL_DIR') + '/Spec/Arcs/Lists/'
;;----------
;; Initial guess for central wavelength and dispersion
telescope = strtrim(sxpar(hdr, 'TELESCOP'))
instrument = strtrim(sxpar(hdr, 'INSTRUME'))
detector = strtrim(sxpar(hdr, 'DETECTOR'))
binning = long(strsplit(sxpar(hdr, 'BINNING'), ',', /extract))
mask    = strtrim(sxpar(hdr, 'SLITNAME'))
mswave = double(strtrim(sxpar(hdr, 'MSWAVE')))
telid     = strtrim(sxpar(hdr, 'TELID'), 2)


wstruct = create_struct('INSTRUMENT', instrument $
                        , 'REID', 0L             $
                        , 'AUTOID', 0L           $
                        , 'CRUDEID', 0L            $
                        , 'PKWDTH', 0.0          $
                        , 'TOLER', 0.0           $
                        , 'NFIND', 0L            $
                        , 'DISP_GUESS', 0.0D     $
                        , 'FUNC', ' '            $
                        , 'NORD_CRUDE', 0L       $
                        , 'SIGREJ_CRUDE', 0.0    $
                        , 'PSIG_CRUDE', 0.0      $
                        , 'PSIG', fltarr(6)      $
                        , 'MXOFF', fltarr(6)     $
                        , 'SIGREJ', fltarr(6)    $
                        , 'NORD', fltarr(6)      $
                        , 'FLG_QUAL', lonarr(6)      $
                        , 'PSIG_REID', fltarr(6)     $
                        , 'MXOFF_REID', fltarr(6)    $
                        , 'SIGREJ_REID', fltarr(6)   $
                        , 'NORD_REID', lonarr(6)     $
                        , 'FLG_QUAL_REID', lonarr(6) $
                        , 'LINELIST_AP', ' '     $
                        , 'LINELIST',    ' '     $
                        , 'NPANIC', 0L           $
                        , 'NORD_PANIC', 0L       $
                        , 'DR_WAVE', [0.0, 0.0]  $
                        , 'MXSHIFT', 0.0  $
                        , 'REID_FILE', ' ' $
                        , 'BIN_RATIO',1. $
                        , 'RADIUS', 0.0 $
                        , 'SIG_WPIX', 0.0 $
                        , 'FWEIGHT',1L)

wstruct.radius = 3.0
wstruct.SIG_WPIX = 3.0

;; change from LBT-SX to LBT-* for both MODS instances [rwp/osu 2019-01-22]

IF strmatch(telescope,'LBT-*') THEN BEGIN
    binning = [1, 1]
    grating = strtrim(sxpar(hdr, 'GRATINFO'),2)
    grating = strmid(grating,0,7)
    mask = strtrim(sxpar(hdr, 'MASKNAME'),2)
    ;; Blue side
    CASE grating OF
       '400l/mm': BEGIN
        ;'G450L': BEGIN
           ;; parameters for arc peak finding
            wstruct.pkwdth        = 12.0D/binning[1]
            wstruct.TOLER         = 4.0D/binning[1]
            wstruct.psig[0:4]     = [12.0, 8.0, 7.0, 4.0, 4.0]
            wstruct.mxoff[0:4]    = [7.0, 7.0, 5.0, 5.0, 5.0]/binning[1]
            ;; parameters for wavelength soln fits
            wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 2.5, 2.5]
            wstruct.nord[0:4]     = [4L, 5L, 5L, 5L, 5L]
            wstruct.FLG_QUAL[0:4] = [1L, 1L, 1L, 2L, 2L]
            wstruct.LINELIST      = line_path+'/mods_blue_400kev.lst'
            wstruct.npanic        = 5L
            wstruct.nord_panic    = 2L
            wstruct.AUTOID        = 0L
            wstruct.REID          = 1 ; reidentify using archived solutions
            wstruct.REID_FILE     = calib_path+'/mods_blue_400ms.sav'
            wstruct.BIN_RATIO     = binning[1] ;arxiv soln unbinned spectrally
         END
       '250l/mm': BEGIN
        ;'G670L': BEGIN
           wstruct.pkwdth       = 6.0/binning[1]
           wstruct.TOLER        = 1.5/binning[1]
           wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
           wstruct.MXOFF[0:5]    = [5.0, 5.0, 5.0, 4.0, 3.0, 3.0]/binning[1]
           wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
           wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
           wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
           wstruct.LINELIST      = line_path+'/mods_red_670kev.lst'
           wstruct.npanic      = 20L
           wstruct.nord_panic  = 2L
           ;; AUTOID parameters
           wstruct.nfind        = 30L
           wstruct.disp_guess   = 0.847d0*binning[1]
           wstruct.dr_wave     = [0.9, 1.1]*(0.847d0)*binning[1]
           wstruct.AUTOID      = 0L
           wstruct.REID        = 1
           if STRMATCH(mask,'*LS*') then begin
               wstruct.REID_FILE   = calib_path+'/mods_red_670.sav'
           endif else begin
               wstruct.REID_FILE   = calib_path+'/mods_red_670ms.sav'
           endelse
           wstruct.BIN_RATIO   = binning[1]     ;;arxiv soln unbinned
           wstruct.FUNC  = 'CHEBY'
        END
       'FuSi+Al': BEGIN  ;BLUE Prism
           wstruct.pkwdth        = 6.0/binning[1]
           wstruct.TOLER         = 1.5/binning[1]
           wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
           wstruct.MXOFF[0:5]    = [5.0, 5.0, 5.0, 4.0, 3.0, 3.0]/binning[1]
           wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
           wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
           wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
           wstruct.LINELIST      = line_path+'/mods_blue_prism.lst'
           wstruct.npanic        = 20L
           wstruct.nord_panic    = 2L
           ;; AUTOID parameters
           wstruct.nfind       = 30L
           wstruct.disp_guess  = 0.847d0*binning[1]
           wstruct.dr_wave     = [0.9, 1.1]*(0.847d0)*binning[1]
           wstruct.AUTOID      = 0L
           wstruct.REID        = 1
           wstruct.REID_FILE   = calib_path+'/mods_blue_prism.sav'
           wstruct.BIN_RATIO   = binning[1] ;arxiv soln unbinned spectrally
           wstruct.FUNC  = 'CHEBY'

	   print,'Sooo..... you want to use a PRISM!'
	   stop
        END
        ELSE: message, 'ERROR: Unknown grating'
    ENDCASE
endif else BEGIN
    message, 'ERROR: Unknown instrument!'
ENDELSE

if keyword_set(LINELIST) then begin
   wstruct.LINELIST = line_path+LINELIST
   ;; Kludge for no Ne lamps
   ipos = strpos(LINELIST,'noNe')
   if ipos GT 0 then wstruct.npanic = 5L
endif

if keyword_set(REID_FILE) then wstruct.REID_FILE = calib_path+REID_FILE
if keyword_set(BIN_RATIO) then wstruct.BIN_RATIO = BIN_RATIO

splog,'LINELIST ',wstruct.LINELIST
RETURN, wstruct
END

