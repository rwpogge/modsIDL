;+
; NAME:
;   mods_standard
; PURPOSE:
;   Create continue processing std files for MODS(s)
;
; CALLING SEQUENCE:
;   mods_standard, [ fileexpr, indir, planfile= ]
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
;   27-Aug-2014  Written by KV Croxall, OSU.
;-
;------------------------------------------------------------------------------
function plan_struct, nfile
        planstr = create_struct(name='STD', $
            'FILENAME'    , '',   $ ;file to process
	    'INSTRUMENT'  , '',   $ ; which camera
            'WAVEFILE'    , '',   $ ; which wave file
            'SLITFILE'    , '',   $ ; which slitfile
;	    'OUTNAME'     , '',   $ ; name for out
	    'STD_NAME'    , '',   $ ; name of standard star
            'APERTURES'   , '[1]',$ ; which apertures to process
;            'Z'           , 0.,   $ ; redshift
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
            'TRIM_TOP'    , 1.  , $
            'TRIM_BOT'    , 1.  , $
            'MASK_LINES'  , 0  )
           return, replicate(paramstr, nfile)
end

;------------------------------------------------------------------------------


pro mods_standard,                $
	planfile=planfile,        $
	paramfile=paramfile,      $ ; dual=dual,single=single, $
	outname = outname,        $
        logprofile = logprofile,                        $
	clobber=clobber,          $
	verbose=verbose,          $
	linemask = linemask,      $
	edge_mask = edge_mask,    $
	no_sky_sub = no_sky_sub,  $
	no_tell_mask=no_tell_mask,$
	onlyplan = onlyplan

   if not keyword_set(verbose) then !Quiet=1
   if not keyword_set(paramfile) then paramfile = 'stdparam.par'
   if (NOT keyword_set(planfile)) then planfile = 'std.par'

;make plan file if none is given
   valid_plan = file_test(planfile)
   if not valid_plan then begin
	; find files	
	spawn,"ls Science/std*fits* | grep -v 's2d' |grep -v 'x1d' | grep -v 'ssf' | grep -v 'r1d'",std_filenames
        nfile = n_elements(std_filenames)

	planstr =  plan_struct(nfile)
	paramstr = param_struct(nfile)
	for i = 0L, nfile-1L do begin
		planstr[i].FILENAME = std_filenames[i]
		hdr = xheadfits(std_filenames[i])
		planstr[i].INSTRUMENT = strcompress(sxpar(hdr[*,0], 'INSTRUME'), /rem)
		planstr[i].SLITFILE = strcompress(sxpar(hdr[*,0], 'SLITFILE'), /rem)
		planstr[i].WAVEFILE = strcompress(sxpar(hdr[*,0], 'WAVEFILE'), /rem)
		planstr[i].MASKNAME = strcompress(sxpar(hdr[*,0], 'MASKNAME'), /rem)

		if (planstr[i].INSTRUMENT eq 'MODS1B' or planstr[i].INSTRUMENT eq 'MODS2B') then $ 
			paramstr[i].CENTERLINE = 4500.
		if (planstr[i].INSTRUMENT eq 'MODS1R' or planstr[i].INSTRUMENT eq 'MODS2R') then $
			paramstr[i].CENTERLINE = 7750.
		name_len=STRLEN(std_filenames[i])
		objname = strcompress(sxpar(hdr[*,0], 'OBJNAME'), /rem)
		if objname eq 'Feige34' then planstr[i].STD_NAME = 'feige34_stis_001' $
		else if objname eq 'Feige66'    then planstr[i].STD_NAME = 'feige66_002' $
                else if objname eq 'Feige67'    then planstr[i].STD_NAME = 'feige67_10a' $
                else if objname eq 'G191B2B'    then planstr[i].STD_NAME = 'g191b2b_mod_005' $
                else if objname eq 'G191-B2B'   then planstr[i].STD_NAME = 'g191b2b_mod_005' $
                else if objname eq 'Feige110'   then planstr[i].STD_NAME = 'feige110_stisnic_002' $
		else if objname eq 'Hz43'       then planstr[i].STD_NAME = 'hz43_mod_005' $
                else if objname eq 'Hz44'       then planstr[i].STD_NAME = 'hz44_stis_001' $
		else if objname eq 'GD71'       then planstr[i].STD_NAME = 'gd71_mod_006' $
		else if objname eq 'BD+33d2642' then planstr[i].STD_NAME = 'bd_33d2642_fos_003' $
		else if objname eq 'BD+284211' then planstr[i].STD_NAME = 'bd_28d4211_stis_001' $
		;NEED TO ADD IN THE OTHER STANDARDS NAMES AS SUPPLIED BY CAL SCRIPTS
		else begin
			print,'Not sure which standard to choose for this object'
			print,''
			print,'FILE: ',std_filenames[i]
			print,'OBJNAME: ',objname
			print,''
		        print,'Please give me the name of a standard in my library'
		        print,''
		        print,'Recommended standard files:'
		        print,'         g191b2b_mod_005   g191b2b_stisnic_002  g191b2b_10a'
		        print,'         gd71_mod_006      gd71_stisnic_002 '
		        print,'         feige34_stis_001  feige34_10a'
		        print,'         feige66_002'
		        print,'         feige67_10a       feige67_002'
		        print,'         gd153_mod_005     gd153_stisnic_002  '
		        print,'         hz43_mod_005      hz43_stis_001 '
		        print,'         hz44_stis_001     '
		        print,'         bd_33d2642_fos_003  '
		        print,'         bd_28d4211_stis_001 '
		        print,'         feige110_stisnic_002'
			stdname = ''
	                read,': ',stdname
			planstr[i].STD_NAME = stdname
		endelse
	endfor

        yanny_write, planfile, ptr_new(planstr),/align
        if not file_test(paramfile) then yanny_write, paramfile, ptr_new(paramstr),/align
   endif
   if keyword_set(onlyplan) then stop

   planstr = yanny_readone(planfile, /anonymous)
   paramstr = yanny_readone(paramfile, /anonymous) 

; run skysubtraction
   if not keyword_set(no_sky_sub) then begin
	   for i = 0L, n_elements(planstr)-1L do begin
		mods_skyfit2d_singlechan,planstr[i].FILENAME,$
			wave= planstr[i].WAVEFILE,     $
			slits=planstr[i].SLITFILE,     $
			trim_t = paramstr[i].TRIM_TOP, $
			trim_b = paramstr[i].TRIM_BOT, $
			centerline=[paramstr[i].CENTERLINE-paramstr[i].CENTERWIDTH/2.,$
				paramstr[i].CENTERLINE+paramstr[i].CENTERWIDTH/2.],   $
			boxcar=paramstr[i].BOXCAR,      $
			apertures=[1],clobber=clobber,  $
			verbose=verbose,                $
			linemaskfile=linemask,          $
                        logprofile = logprofile,        $
			centersum=paramstr[i].CENTERSUM,$
			edge_mask = edge_mask
	   endfor
   endif

; run fluxing
   for i = 0L, n_elements(planstr)-1L do begin
	name_len=STRLEN(planstr[i].FILENAME)
	skysub_name = STRMID(planstr[i].FILENAME, 0, name_len-8) + '_s2d.fits' 
	mods_fluxstand_singlechan, skysub_name, std_name=planstr[i].STD_NAME,$
                wave= planstr[i].WAVEFILE, $
                slits=planstr[i].SLITFILE,$
                centerline=[paramstr[i].CENTERLINE-paramstr[i].CENTERWIDTH/2.,$
                       paramstr[i].CENTERLINE+paramstr[i].CENTERWIDTH/2.],   $
                verbose=verbose,clobber=clobber, $
                centersum=paramstr[i].CENTERSUM,no_tell_mask=no_tell_mask
   endfor

end
