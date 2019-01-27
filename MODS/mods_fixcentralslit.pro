pro mods_fixcentralslit,file,file2,slitfile,slit=slit,top=top,outname=outname,bot=bot
; file = file where the other slits are good
; file2 file where the other half of the central slit is good

  if n_elements(slit) eq 0 then slit=1
  if n_elements(top) eq  0 then top=0
  if n_elements(bot) eq  0 then top=1
  if n_elements(bot) eq  1 then top=0
  if n_elements(outname) eq  0 then outname = 'the_great_nameless_file.fits'

  im_0 = mrdfits(file,0,head)
  im_1 = mrdfits(file,1,hd)
  im_2 = mrdfits(file,2,hd)
  im_3 = mrdfits(file,3,hd)
;  im_4 = mrdfits(file,4,hd)
;  im_5 = mrdfits(file,5,hd)

  im_0b = mrdfits(file2,0,head2)

  slim = mrdfits(slitfile,0,shead)

  sl_loc = where(slim eq slit)
  tmpim = im_0*0
  tmpim[sl_loc] = im_0b[sl_loc]
  if (top eq 1) then tmpim[1545:3087,*] =0
  if (top eq 0) then tmpim[0:1544,*] =0
  new=where(tmpim ne 0)
  im_0[new]=im_0b[new]
  
  mwrfits,im_0,outname,head,/create
  mwrfits,im_1,outname
  mwrfits,im_2,outname
  mwrfits,im_3,outname
;  mwrfits,im_4,outname
;  mwrfits,im_5,outname
end

