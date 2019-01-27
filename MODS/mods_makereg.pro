pro mods_makereg,mms,out

openw,42,out,width=500
printf,42,'# Region file format: DS9 version 4.1'
printf,42,'global color=magenta dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
printf,42,'fk5'

readcol,mms,input,val,FORMAT='A,D',/silent

col1 = input
col2 = input
col3 = input
dim=size(input, /dimensions)
for n=0,dim[0]-1 do begin
        parts=str_sep(input[n],'.')
        col1[n]=parts[0]
        col2[n]=parts[1]
        col3[n]=parts[2]
endfor
for n=0,dim[0]-1  do begin
	if col2[n] eq 'ROT' then pangle = strtrim(string(val[n],Format='(F20.1)'),2)
	if col2[n] eq 'TARG' then begin
		if col3[n] eq 'ALPHA' then begin
			taralpha = strtrim(string(val[n],Format='(F20.3)'),2)
			tardelta = strtrim(string(val[n+1],Format='(F20.3)'),2)
		endif
	endif
endfor
radot = strpos(taralpha,'.')
ra1=strmid(taralpha,radot-6,2)
ra2=strmid(taralpha,radot-4,2)
ra3=strmid(taralpha,radot-2)
decdot = strpos(tardelta,'.')
dec1=strmid(tardelta,decdot-6,2)
dec2=strmid(tardelta,decdot-4,2)
dec3=strmid(tardelta,decdot-2)
tarRA_deg = 15 * (double(ra1) + double(ra2)/60 + double(ra3)/3600)
tarDEC_deg = double(dec1) + double(dec2)/60 + double(dec3)/3600
tarRA_sex = ra1 + ':' + ra2 + ':' + ra3
tarDEC_sex = dec1 + ':' + dec2 + ':' + dec3
printf,42,'circle(',tarRA_sex,',+',tarDEC_sex,',5") # text={CENTER}'

slitno=0
for n=0,dim[0]-1  do begin
        if col3[n] eq 'WIDMM' then begin
                WIDMM=val[n]
                LENMM=val[n+1]
                if WIDMM gt 0.6 or LENMM gt 0.6 and WIDMM ne LENMM then slitno+=1
        endif
endfor
print,'Found ',slitno,' slits in the mms file.'
slit=1
for n=0,dim[0]-1  do begin
    if col3[n] eq 'WIDMM' then begin
	WID = strtrim(string(val[n-5],Format='(F20.1)'),2)
        LEN = strtrim(string(val[n-4],Format='(F20.1)'),2)
        WIDMM=val[n]
        LENMM=val[n+1]
        XMM=val[n+2]
        YMM=val[n+3]
	ALPHA=string(val[n-2],Format='(F20.3)')
	DELTA=string(val[n-1],Format='(F20.3)')
        if WIDMM gt 0.6 or LENMM gt 0.6 and WIDMM ne LENMM then begin
	   radot = strpos(ALPHA,'.')
           ra1=strmid(ALPHA,radot-6,2)
           ra2=strmid(ALPHA,radot-4,2)
	   ra3=strmid(ALPHA,radot-2)
	   decdot = strpos(DELTA,'.')
	   dec1=strmid(DELTA,decdot-6,2)
	   dec2=strmid(DELTA,decdot-4,2)
	   dec3=strmid(DELTA,decdot-2)
	   apRA_deg = 15.D * (double(ra1) + double(ra2)/60.D + double(ra3)/3600.D)
	   apDEC_deg = double(dec1) + double(dec2)/60.D + double(dec3)/3600.D
	   apRA_sex = ra1 + ':' + ra2 + ':' + ra3
	   apDEC_sex = dec1 + ':' + dec2 + ':' + dec3
	   offRA = strtrim(string(((tarRA_deg - apRA_deg) * cos(apDEC_deg * !dpi/180))*3600.D,Format='(F20.2)'),2)
	   offDEC = strtrim(string((tarDEC_deg - apDEC_deg)*3600.D,Format='(F20.2)'),2)
	   ang_off = djs_diff_angle(tarRA_deg,tarDEC_deg,apRA_deg,apDEC_deg)
	   printf,42,'box(',apRA_sex,',+',apDEC_sex,',',WID,'",',LEN,'",',pangle,$
		') # text={#',strtrim(slit,2),'}'
	   printf,42,'# slit#',strtrim(slit,2),' has an offset of ',offRA,',',offDEC
	   slit=slit+1
	endif
    endif
endfor

close,42
end
