pro mods_blank_nonslit,file,slitfile,out

  image = mrdfits(file,0,hdr)
  slits = mrdfits(slitfile,0,slitheader)

  sl_trans = transpose(slits)
  blank = where(sl_trans eq 0)
  image[blank] = 0

  mwrfits,image,out,hdr,/create

end
   
