;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  designed for 2D skysub images to replace apertures that needed tweaking
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mods_apreplace,infile_good,outfile_tofix,apertures,slit_mask

; Load the images
  good_img = mrdfits(infile_good,0,hdr_2)
  good_img1 = mrdfits(infile_good,1,hdr_2b)
  good_img2 = mrdfits(infile_good,2,hdr_2b)
  good_img3 = mrdfits(infile_good,3,hdr_2b)

  tofix_img = mrdfits(outfile_tofix,0,hdr_1)
  tofix_img1 = mrdfits(outfile_tofix,1,hdr_1b)
  tofix_img2 = mrdfits(outfile_tofix,2,hdr_1b)
  tofix_img3 = mrdfits(outfile_tofix,3,hdr_1b)

;create the slit structures 
  tset_slits = mrdfits(slit_mask,1,tsethead)
  slit_im = mrdfits(slit_mask,0,tset_head)

  ndo = size(apertures,/dim)
  ndo=ndo[0]
  for ii=0,ndo-1 DO BEGIN
	ap=apertures[ii]
	location = where(slit_im EQ ap)
	tofix_img[location] = good_img[location]
        tofix_img1[location] = good_img1[location]
        tofix_img3[location] = good_img3[location]
  ENDFOR

  mwrfits,tofix_img,outfile_tofix,hdr_1,/create
  mwrfits,tofix_img1,outfile_tofix
  mwrfits,tofix_img2,outfile_tofix
  mwrfits,tofix_img3,outfile_tofix
END

