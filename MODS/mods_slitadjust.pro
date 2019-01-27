;+
; NAME:
;   MODS_SLITADJUST
;
; VERSION:
;   1.1 (Jul, 2013)
;
; PURPOSE:
;   Adjust the upper and lower bounds of MODS slits. The slits may need
;   ajdustment due to a slight shift resulting in how the mask was placed
;   in the instrument; misadjustment due to neighboring slits.
;
; REFERENCE:
;   Croxall, Pogge, and Stoll in prep
;
; CALLING SEQUENCE:
;   mods_slitadjust, slit_file, image, [apertures=]
;
; INPUTS:
;   slit_file: XIDL slit file to be altered/checked
;
;   image: XIDL transposed image to compare against the slit traces
;
; OPTIONAL INPUTS:
;   apertures: Array of apertures to be adjusted.  If no apertures are passed
;      then all apertures are checked.
;
; OUTPUT:
;   Saves the new slit structure in the place of the old slit structure.
;
; EXAMPLE: 
;   mods_slitadjust,'slits-blue_ill.fits','Science/sci-ngc628_f1b_sum.fits.gz'
;
; History: K.V. Croxall Jul-2 2013
;    2013-07-10 (K. Croxall): Fixed slits that go off the detector.
;    2013-07-02 (K. Croxall): Written
;;;

pro mods_slitadjust,slit_structure_file,image_file,apertures=apertures,raw=raw

;open files & window
  window,1,retain=2
  image = mrdfits(image_file,0,hdr)
  slit_structure = mrdfits(slit_structure_file,1,tset_head)
  slit_image = mrdfits(slit_structure_file,0,tset_head)
  if keyword_set(raw) then image = transpose(image)

;get image size
  dims=size(image,/dim)
  nx=dims[0]
  ny=dims[1]

;get structure size
  dim_t = slit_structure[0].dims
  nslit=size(slit_structure.xx1,/dim)
  nslit=nslit[1]
  nx_t = dims[0]
  ny_t = dims[1]

  xx1 = slit_structure[0].xx1
  xx2 = slit_structure[1].xx2

;tweak apertures
  if n_elements(apertures) EQ '' THEN apertures = indgen(nslit)+1
  ndo = size(apertures,/dim)
  ndo=ndo[0]

;adjust the selected slits
  FOR ii=0,ndo-1 DO BEGIN
	ap=apertures[ii]
	print,'APERTURE: ',ap

	LeSigh: print,'Plotting Aperture ',ap

	ncolors = !D.Table_Size
        device,DECOMPOSED=0
        LoadCT,39
        ap_min = min(xx1[*,ap-1])-10 >0
;nonzero = where(xx1[*,ap-1] gt 0)
;ap_min = min(xx1[nonzero])-10 >0
;print,ap_min
;stop
        ap_max = max(xx2[*,ap-1])+10 <3078
	slit_image = image[ap_min:ap_max,*]
        slit_image = transpose(BytScl(slit_image,min=1,max=10*median(slit_image),/nan))
        TVImage,slit_image,Position=[0.02,0.08,0.99,0.98],/erase
        path_xx1 = xx1[*,ap-1] - ap_min
        path_xx2 = xx2[*,ap-1] - ap_min
        xpix = indgen(8192)
        plot,xpix,path_xx1,/noerase,Position=[0.02,0.08,0.99,0.98],$
             xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,ap_max-ap_min]
        plot,xpix,path_xx2,/noerase,Position=[0.02,0.08,0.99,0.98],$
             xstyle=1,xrange=[0,8191],ystyle=1,yrange=[0,ap_max-ap_min]

	;alter boundaries if needed
        Ok=0
        WHILE Ok NE 1 DO BEGIN
                print,''
                print,'Adjust place ment of Slit traces:'
                print,'1 - Tweak lower trace'
                print,'2 - Tweak upper trace'
		print,'3 - Accept Slit'
                print,'4 - STOP'
                redo = ''
                read,': ',redo
                IF redo EQ '1' THEN BEGIN
                        Print,'Pixel Shift:'
                        read,': ',delta_y
;			xx1[nonzero] = xx1[nonzero]+delta_y
			xx1[*,ap-1] = xx1[*,ap-1]+delta_y
                        GOTO,LeSigh
                ENDIF ELSE IF redo EQ '2' THEN BEGIN
                        Print,'Pixel Shift:'
                        read,': ',delta_y
;			xx2[nonzero] = xx2[nonzero]+delta_y
                        xx2[*,ap-1] = xx2[*,ap-1]+delta_y
                        GOTO,LeSigh
                ENDIF ELSE IF redo EQ '3' THEN BEGIN
			Ok = 1
                ENDIF ELSE IF redo EQ '4' THEN BEGIN
                        stop
                ENDIF ELSE BEGIN
                        print,'Not A valid option.',redo
                ENDELSE
        ENDWHILE
  	slit_structure[0].xx1 = xx1 > 0.
  	slit_structure[1].xx2 = xx2 < 3088.
	
	print,'Updating mask in ',slit_structure_file
	slit_image = mods_slits2mask(slit_structure, nslit = nslit)
  	mwrfits,slit_image,slit_structure_file,hdr,/create
	mwrfits,slit_structure,slit_structure_file
  ENDFOR
END
