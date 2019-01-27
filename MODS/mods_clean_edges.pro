pro mods_clean_edges,list

  readcol,list,format='A',filelist
  num_ims = n_elements(filelist)
  for i=0,num_ims-1 do begin
	name_len=STRLEN(filelist[i])
	outname = STRMID(filelist[i], 0, name_len-5) + 'e.fits'
	image = mrdfits(filelist[i],0,hdr)
	image[*,0:112] = 0	;bottom
	image[0,*] = 0		;left 
	image[*,3062:3087] = 0  ;top
	mwrfits,image,outname,hdr,/create
  endfor

end
   
