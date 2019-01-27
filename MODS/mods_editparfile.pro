;=========================================================================
;=========================================================================
pro readnum,value
        satisfactory = 0
        input = ''
        repeat begin
                read,': ',input
                if not valid_num(input) then print,'I need a number' else begin
                        value = float(input)
                        satisfactory = 1
                endelse
        endrep until satisfactory
return
end

;=========================================================================
;=========================================================================
pro mods_editparfile,parfile,$
        verbose=verbose

	old_parfile = parfile + '_old'
   	if not keyword_set(verbose) then !Quiet=1
	parameters = yanny_readone(parfile, /anonymous)
	yanny_write, old_parfile, ptr_new(parameters),/align ;save old version
	print_struct,parameters
	tnames = TAG_NAMES(parameters)

	rinseAndRepeat:
	print,''
	print,'Type then name of the parameter you wish to edit'
	field = ''
	chooseParameter:
	read,': ',field
	tindex = where(STRCMP(tnames,field,/fold_case) eq 1)
	if tindex eq (-1) then begin
		print,'The selected parameter name (',field,')is not recognized.'
		print,'Please choose from the following list: '
		print,tnames
		goto,chooseParameter
	endif	
	data = parameters.(tindex)

	n_cols = n_elements(data)
	chooseIndex:
	print,'Type the index of the line you wish to edit (0 - ',strtrim(n_cols-1,2),')'
	readnum,nindex
	if nindex gt n_cols-1 then begin
		print, 'There are not that many lines in the parameter file.'
		goto,chooseIndex
	endif
	print,'You have selected:'
	help,data[nindex]

	print,'Please enter new value for the selected parameter.'
	new = ''
	read,': ',new
	parameters[nindex].(tindex) = new
	
	print_struct,parameters
	print,'Would you like to edit another parameter? (Y/N)'
	Ok=0
	while Ok ne 1 do begin
		again = ''
		read,': ',again
		if STRCMP(again,'y',/fold_case) eq 1 then GOTO,rinseAndRepeat $
			else if STRCMP(again,'n',/fold_case) then Ok = 1 $
			else print,'I do not understand ',again
	endwhile

	yanny_write, parfile, ptr_new(parameters),/align
end
	
