;+
; NAME:
;   mods_science
; PURPOSE:
;   Create continue processing sci files for MODS(s)
;
; CALLING SEQUENCE:
;   mods_science, [ fileexpr, indir, planfile= ]
;
; INPUTS:
;
; OUTPUT:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   15-Oct-2014  Written by KV Croxall, OSU.
;-
;------------------------------------------------------------------------------

;=========================================================================
; read_index - make sure the index is a valid number
;=========================================================================
pro read_index,value,maxnum=maxnum
        satisfactory = 0
	valid = 0
        input = ''
	repeat begin
           repeat begin
                read,': ',input
                if not valid_num(input) then print,'I need a number' else begin
                        value = float(input)
			if keyword_set(maxnum) and value le maxnum then valid = 1 $
				else print,fix(value),' is larger than allowed by the input array.'
                        satisfactory = 1
                endelse
           endrep until satisfactory
	endrep until valid
return
end


function plan_struct, nfile
        planstr = create_struct(name='SCI', $
            'FILENAME'    , '',   $ ;file to process
	    'INSTRUMENT'  , '',   $ ; which camera
            'WAVEFILE'    , '',   $ ; which wave file
            'SLITFILE'    , '',   $ ; which slitfile
;	    'OUTNAME'     , '',   $ ; name for out
	    'SSFFILE'     , '',   $ ; which calibration file
            'APERTURES'   , '"ALL"',   $ ; which apertures to process
            'SKYSLIT'     ,  0,   $ ; which slit to use for the sky-slit
            'Z'           , 0.,   $ ; redshift
	    'DUAL'        ,  0,   $ ; for use if you want to extract red and blue together
            'MASKNAME'    , ''  )
           return, replicate(planstr, nfile)
end

function param_struct,nfile
	paramstr = create_struct(name='PARAMS', $
;            'CONVBEAM'    , "[1.0,2.0]", $ ; convolution beam
            'FILENAME'    , '',   $ ;file to process
            'BOXCAR'      , 0,    $ ; use of boxcar y/n
            'CENTERLINE'  , 0,   $
            'CENTERWIDTH' , 500., $ ; size of the extraction window
            'CENTERSUM'   , 10. , $
            'TRIM_TOP'    , 2.  , $
            'TRIM_BOT'    , 2.  , $
            'MASK_LINES'  , 0  )
           return, replicate(paramstr, nfile)
end

function info_struct, nfile
        infostr = create_struct(name='INFO', $
            'FILENAME'    , '',   $ ;file to process
	    'DATE'        , '',   $ ; Date of obs
            'OBJECT'      , '',   $ ; which object
            'EXPTIME'     , '',   $ ; exposure time
            'OBSMODE'     , '',   $ ; which mode
            'MASKNAME'    , ''  )
           return, replicate(infostr, nfile)
end

;------------------------------------------------------------------------------

pro mods_science,                   $
	planfile = planfile,        $
	paramfile = paramfile,      $ 
	infofile = infofile,        $
	dual = dual,                $
	redonly = redonly,          $
	blueonly = blueonly,        $
	outname = outname,          $
	clobber = clobber,          $
	verbose = verbose,          $
	linemask = linemask,        $  
	edge_mask = edge_mask,      $
	no_sky_sub = no_sky_sub,    $
	no_tell_mask = no_tell_mask,$ 
	em_line = em_line,          $
        logprofile = logprofile,    $
	man_scale = man_scale,	    $
	planonly = planonly,        $
        force_blue = force_blue,    $
	force_red = force_red,      $
	subcont=subcont,	    $
	broad_lines = broad_lines

   if not (keyword_set(dual) or keyword_set(redonly) or keyword_set(blueonly)) then begin
	print,'No instrument mode set.  Please use one of the mode flags:'
	print,'       /dual, /blueonly, /redonly'
	stop
   endif
   if not keyword_set(verbose) then !Quiet=1
   if not keyword_set(paramfile) then paramfile = 'sciparam.par'
   if NOT keyword_set(planfile)  then planfile  = 'sci.par'
   if not keyword_set(infofile)  then infofile  = 'sciinfo.par'

;make plan file if none is given
   valid_plan = file_test(planfile)
   if valid_plan then print,planfile,' already exists - proceeding with the given file.'
   if not valid_plan then begin
	planonly = 1
	; find files	
	spawn,"ls Science/sci*fits* | grep -v 's2d' |grep -v 'x1d' | grep -v 'ssf' | grep -v 'r1d'",sci_filenames
	spawn,"ls Science/std*ssf.fits*",ssf_filenames
        nfile = n_elements(sci_filenames)

	planstr =  plan_struct(nfile)
	paramstr = param_struct(nfile)
	infostr = info_struct(nfile)

	for i = 0L, nfile-1L do begin
		planstr[i].FILENAME = sci_filenames[i]
                paramstr[i].FILENAME = sci_filenames[i]
		infostr[i].FILENAME = sci_filenames[i]
		hdr = xheadfits(sci_filenames[i])
		planstr[i].INSTRUMENT = strcompress(sxpar(hdr[*,0], 'INSTRUME'), /rem)
		planstr[i].SLITFILE = strcompress(sxpar(hdr[*,0], 'SLITFILE'), /rem)
		planstr[i].WAVEFILE = strcompress(sxpar(hdr[*,0], 'WAVEFILE'), /rem)
		planstr[i].MASKNAME = strcompress(sxpar(hdr[*,0], 'MASKNAME'), /rem)

		infostr[i].OBJECT   = strcompress(sxpar(hdr[*,0], 'OBJECT'), /rem)
		infostr[i].DATE   = strcompress(sxpar(hdr[*,0], 'DATE-OBS'), /rem)
		infostr[i].EXPTIME  = strcompress(sxpar(hdr[*,0], 'EXPTIME'), /rem)
		infostr[i].OBSMODE  = strcompress(sxpar(hdr[*,0], 'DICHNAME'), /rem)
		infostr[i].MASKNAME = strcompress(sxpar(hdr[*,0], 'MASKNAME'), /rem)

		if (planstr[i].INSTRUMENT eq 'MODS1B' or planstr[i].INSTRUMENT eq 'MODS2B') then begin
			loc_ssf = where(strmatch(ssf_filenames,'*_m1b_*',/fold_case) eq 1)
			if n_elements(loc_ssf) gt 1 then begin
				print,'There is more than 1 valid SSF file for ',sci_filenames[i],'.  Please choose one before .continue'
				print,"ssf_file = '<FILENAME>'"
				print,ssf_filenames[loc_ssf]
				stop
			endif else if loc_ssf eq (-1) then begin
				print,'No valid SSF files found, have you processed all flux standards?'
				print,ssf_filenames
				stop
			endif else ssf_file = ssf_filenames[loc_ssf]
			planstr[i].SSFFILE = ssf_file
		endif
		if (planstr[i].INSTRUMENT eq 'MODS1R' or planstr[i].INSTRUMENT eq 'MODS2R') then begin
			loc_ssf = where(strmatch(ssf_filenames,'*_m1r_*',/fold_case) eq 1)
                        if n_elements(loc_ssf) gt 1 then begin
                                print,'There is more than 1 valid SSF file for ',sci_filenames[i],'.  Please choose one before .continue'
                                print,"ssf_file = '<FILENAME>'"
                                print,ssf_filenames[loc_ssf]
                                stop
                        endif else if loc_ssf eq (-1) then begin
                                print,"No valid SSF files recognized, If you have processed all flux standards this could be due to adopting a different naming convention.  Please select a file from the following list using ssf_file = '<ssf file>'"
                                print,ssf_filenames
                                stop
                        endif else ssf_file = ssf_filenames[loc_ssf]
                        planstr[i].SSFFILE = ssf_file

		endif
		if (planstr[i].INSTRUMENT eq 'MODS1B' or planstr[i].INSTRUMENT eq 'MODS2B') then $ 
			paramstr[i].CENTERLINE = 4500.
		if (planstr[i].INSTRUMENT eq 'MODS1R' or planstr[i].INSTRUMENT eq 'MODS2R') then $
			paramstr[i].CENTERLINE = 7750.
		objname = strcompress(sxpar(hdr[*,0], 'OBJNAME'), /rem)
	endfor

	hdr = ["############################################################################################################################",$
                "##  indexType      filename       instrument    wavefile   slitfile    ssffile   apertures   skyslit   z    dual    maskname",$
                "############################################################################################################################"]
        yanny_write, planfile, ptr_new(planstr),/align,hdr=hdr
	hdr[1] = "##  indexType      filename       boxcar     centerline     centerwidth      centersum     trim_top     trim_bot     mask_lines"
        if not file_test(paramfile) then yanny_write, paramfile, ptr_new(paramstr),/align,hdr=hdr
	hdr[1] = "##  indexType      filename       obsDate     object     expTime     obsMode     maskName"
	if not file_test(infofile) then yanny_write,infofile, ptr_new(infostr),/align,hdr=hdr
   endif
 
  IF keyword_set(planonly) then BEGIN
	print,'---------------------------------------------------------------'
	print,'Three parameter files have been created'
	print,'Please check them to ensure the desired extraction parameters are set before continuing.'
        print,'---------------------------------------------------------------'
	stop
   ENDIF

   planstr = yanny_readone(planfile, /anonymous)
   paramstr = yanny_readone(paramfile, /anonymous) 
   infostr = yanny_readone(infofile, /anonymous)

   if not (ARRAY_EQUAL(planstr.FILENAME,infostr.FILENAME) and $
       ARRAY_EQUAL(planstr.FILENAME,paramstr.FILENAME) ) then begin
	print,'----------------------------------------------------------------------------------'
	print,'The parameter files have mismatched input.  Make sure they reference the same files'
	print,'----------------------------------------------------------------------------------'
	stop
   endif

; run skysubtraction
   if not keyword_set(no_sky_sub) then begin
	   for i = 0L, n_elements(planstr)-1L do begin
		if strlen(planstr[i].apertures) gt 0 then begin
			if planstr[i].apertures eq 'ALL'then undefine,apertures else begin
			     apertures = strmid(planstr[i].apertures,1,strlen(planstr[i].apertures)-2)
			     apertures = fix(strsplit(apertures,',',/extract))
			     endelse
		endif else begin
			print,'No apertures selected, thus no actions performed'
			stop
		endelse
  		   mods_skyfit2d_singlechan,planstr[i].FILENAME,			$
			wave= planstr[i].WAVEFILE, 					$
			slits=planstr[i].SLITFILE,					$
			trim_t = paramstr[i].TRIM_TOP, 					$
			trim_b = paramstr[i].TRIM_BOT, 					$
			centerline=[paramstr[i].CENTERLINE-paramstr[i].CENTERWIDTH/2.,	$
				paramstr[i].CENTERLINE+paramstr[i].CENTERWIDTH/2.],   	$
			apertures=apertures, 						$
			clobber=clobber, 						$
			skyslit = planstr[i].SKYSLIT,					$
			verbose=verbose,						$
                        logprofile = logprofile,                        $
			linemaskfile=linemask,						$
			centersum=paramstr[i].CENTERSUM,				$
			edge_mask = edge_mask, 						$
			z = planstr[i].Z, 						$
			outname = outname,						$
			broad_lines = broad_lines
	   endfor
   endif

; run extraction
   if keyword_set(dual) then begin
	for i=1,n_elements(infostr)/2. do begin
		ok=0
		blue_index = 0
		red_index  = 1
		print_struct,infostr
		while ok ne 1 do begin
	                print,'Please Enter the INDEX of the BLUE file (To the left of the filenames)'
			print,'Note: Indexes are hidden if there are only two files in that case they are 0 & 1'
	                read_index,blue_index,maxnum=fix(n_elements(infostr)-1)
	                print,'Please Enter the INDEX of the RED file (To the left of the filenames)'
	                read_index,red_index,maxnum=fix(n_elements(infostr)-1)
		        name_len=STRLEN(planstr[blue_index].FILENAME)
			skysub_name_blue = STRMID(planstr[blue_index].FILENAME, 0, name_len-8) + '_s2d.fits'
			name_len=STRLEN(planstr[red_index].FILENAME)
			skysub_name_red = STRMID(planstr[red_index].FILENAME, 0, name_len-8) + '_s2d.fits'

			print_struct,infostr,which=[blue_index,red_index]
			print,'Do you want to extract these images as a dual pair?'	

	                redo = ''
	                read,': ',redo
	                redo = strlowcase(redo)
	                if redo eq 'y' then begin
				ok = 1
                        endif else if redo eq 'n' then begin
				print_struct,infostr
	                endif else if redo eq 'abort' then begin
	                       stop
	                endif else begin
	                        print,'I do not understand ',redo
				print_struct,infostr
	                endelse
		endwhile

                if strlen(planstr[blue_index].apertures) gt 0 then begin
                        if planstr[blue_index].apertures eq 'ALL'then undefine,apertures else begin
                             apertures = strmid(planstr[blue_index].apertures,1,strlen(planstr[blue_index].apertures)-2)
                             apertures = fix(strsplit(apertures,',',/extract))
                             endelse
                endif else begin
                        print,'No apertures selected, thus no actions performed'
                        stop
                endelse

		mods_extract_dualchan,skysub_name_blue,skysub_name_red, $
			wave_blue = planstr[blue_index].WAVEFILE,       $
			wave_red  = planstr[red_index].WAVEFILE,        $
			red_cal   = planstr[red_index].ssffile,         $
			blue_cal  = planstr[blue_index].ssffile,        $
			red_slits = planstr[red_index].SLITFILE,        $
			blue_slits= planstr[blue_index].SLITFILE,       $
			trim_rt   = paramstr[red_index].TRIM_TOP,       $
			trim_rb   = paramstr[red_index].TRIM_BOT,       $
			trim_bt   = paramstr[blue_index].TRIM_TOP,      $
			trim_bb   = paramstr[blue_index].TRIM_BOT,      $
			centerline_red=[paramstr[red_index].CENTERLINE-paramstr[red_index].CENTERWIDTH/2.,$
				paramstr[red_index].CENTERLINE+paramstr[red_index].CENTERWIDTH/2.],   $
			centerline_blue=[paramstr[blue_index].CENTERLINE-paramstr[blue_index].CENTERWIDTH/2.,$
				paramstr[blue_index].CENTERLINE+paramstr[blue_index].CENTERWIDTH/2.],   $
			verbose=verbose,                                $
			clobber=clobber,                                $
			centersum_red = paramstr[red_index].CENTERSUM,  $
			centersum_blue = paramstr[blue_index].CENTERSUM,$ 
			apertures=apertures,                            $
			skyslit = skyslit,                              $
			z = planstr[blue_index].Z,                      $
	                em_line=em_line,                                $
			logprofile = logprofile,                        $
			extraction_reference = extraction_reference,    $
			force_blue = force_blue,                        $
        		force_red = force_red,				$
 			man_scale = man_scale, 				$
			subcont=subcont,				$
                        outname = outname

	endfor
   endif else begin
      for i = 0L, n_elements(planstr)-1L do begin
        if strlen(planstr[i].apertures) gt 0 then begin
                if planstr[i].apertures eq 'ALL'then undefine,apertures else begin
                        apertures = strmid(planstr[i].apertures,1,strlen(planstr[i].apertures)-2)
                	apertures = fix(strsplit(apertures,',',/extract))
                endelse
        endif else begin
        	print,'No apertures selected, thus no actions performed'
                stop
        endelse
   	name_len=STRLEN(planstr[i].FILENAME)
	skysub_name = STRMID(planstr[i].FILENAME, 0, name_len-8) + '_s2d.fits' 
        if (keyword_set(redonly) and $
               (planstr[i].INSTRUMENT eq 'MODS1R' or planstr[i].INSTRUMENT eq 'MODS2R')) or $
           (keyword_set(blueonly) and $
               (planstr[i].INSTRUMENT eq 'MODS1B' or planstr[i].INSTRUMENT eq 'MODS2B')) then begin
		if (keyword_set(redonly) and keyword_set(force_red)) then force = force_red
		if (keyword_set(blueonly) and keyword_set(force_blue)) then force = force_blue
	  	mods_extract_singlechan, skysub_name, $
			cal=planstr[i].ssffile,       $
	                wave= planstr[i].WAVEFILE,    $
	                slits=planstr[i].SLITFILE,    $
	                trim_t = paramstr[i].TRIM_TOP,$
			trim_b = paramstr[i].TRIM_BOT,$
	                centerline=[paramstr[i].CENTERLINE-paramstr[i].CENTERWIDTH/2.,$
	                       paramstr[i].CENTERLINE+paramstr[i].CENTERWIDTH/2.],   $
	                verbose=verbose,              $
			clobber=clobber,              $
	                centersum=paramstr[i].CENTERSUM, $
	                apertures=apertures,          $
	                z = planstr[i].Z,	      $
			em_line=em_line,              $
	                logprofile = logprofile,      $
			extraction_reference = extraction_reference,$
			force = force, $
			subcont=subcont,$
	                outname = outname
	endif
      endfor
   endelse

;force_blue=force_blue, force_red=force_red, $
;scale_minmax=scale_minmax, $
;acquisition_centering_correction = acquisition_centering_correction, $
;man_scale = man_scale, $
;subcont = subcont, $


end
