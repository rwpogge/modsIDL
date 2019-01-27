;+
; NAME:
;   mods_empwave
;
; PURPOSE:
;   Generate 
;
; CALLING SEQUENCE:
;   mods_
;
; INPUTS:
;
; OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:

; REVISION HISTORY:
;   Jul-2013  Written by K. Croxall
;-
;------------------------------------------------------------------------------

;=========================================================================
;  mods_empwave -- Create the MODS-XIDL 
;=========================================================================
PRO MODS_EMPWAVE, slitfile,FILENAME=FILENAME, apertures=apertures

   if not keyword_set(verbose) then !Quiet=1

slitimage = mrdfits(slitfile,0,im_hd,/silent)
slit_struct = mrdfits(slitfile,1,struct_hd,/silent)

IF n_elements(FILENAME) EQ 0 then FILENAME = 'wave-empirical.fits'
n_slits=size(slit_struct.XX1,/dim)
n_slits=n_slits[1]
if n_slits eq 2 then begin
	print,'You have n_slits = 2'
	print,'If this is a standard star you only have 1 slit'
	print,'to continue with 2 slits use .continue'
	print,'to continue with 1 slit: n_slits = 1 then .continue'
	stop
endif

grat = strcompress(sxpar(im_hd[*,0], 'GRATNAME'), /rem)
instr= strcompress(sxpar(im_hd[*,0], 'INSTRUME'), /rem)

empfit=dblarr(3088,8192)

if n_slits eq 0 then n_slits = 1
if n_elements(apertures) EQ '' then apertures = indgen(n_slits)+1
ndo = size(apertures,/dim)
ndo=ndo[0]

;for n=0,1  do begin
for n=0,ndo-1 do begin
   print,'Generating 2D wavelength map for slit', n+1
   WIDMM =slit_struct[0].WIDMM_ARR[n]
   LENMM =slit_struct[0].LENMM_ARR[n]
   XMM = slit_struct[0].XMM_ARR[n]
   YMM = slit_struct[0].YMM_ARR[n]
   red_lam=fltarr(101)

   ;BLUE GRATING
   if GRAT eq 'G450L' then begin
      splog,'No G450L grating exists.'
      splog,'there may have been an issue in the header (this is known to have been mis-entered once in the past).'
      splog,'please check your data to be sure you are using the MODS Blue grating (G400L)'
      print,'To adopt G400L, you may .continue'
      stop
      GRAT = 'G400L'
   endif

   if GRAT eq 'G400L' then begin
      lam=fltarr(4601) & coord_pix0=fltarr(4601,2) & coord_pix1=fltarr(4601,2) & coord_pix2=fltarr(4601,2)
      coord_pix3=fltarr(4601,2) & coord_pix4=fltarr(4601,2) & coord_pix5=fltarr(4601,2) & coord_pix6=fltarr(4601,2)
      coord_pix7=fltarr(4601,2) & coord_pix8=fltarr(4601,2) & coord_pix9=fltarr(4601,2) & coord_pix10=fltarr(4601,2)
      ;calculate the slit trace
      for i=0,4600 do begin
          lam[i]=3000+0.65*i
          coord_pix0[i,*]=blue_map(instr,XMM,-(YMM-LENMM),lam[i])
          coord_pix1[i,*]=blue_map(instr,XMM,-(YMM-LENMM*4.0/5.0),lam[i])
          coord_pix2[i,*]=blue_map(instr,XMM,-(YMM-LENMM*3.0/5.0),lam[i])
          coord_pix3[i,*]=blue_map(instr,XMM,-(YMM-LENMM*2.0/5.0),lam[i])
          coord_pix4[i,*]=blue_map(instr,XMM,-(YMM-LENMM*1.0/5.0),lam[i])
          coord_pix5[i,*]=blue_map(instr,XMM,-(YMM+LENMM*0),lam[i])
          coord_pix6[i,*]=blue_map(instr,XMM,-(YMM+LENMM*1.0/5.0),lam[i])
          coord_pix7[i,*]=blue_map(instr,XMM,-(YMM+LENMM*2.0/5.0),lam[i])
          coord_pix8[i,*]=blue_map(instr,XMM,-(YMM+LENMM*3.0/5.0),lam[i])
          coord_pix9[i,*]=blue_map(instr,XMM,-(YMM+LENMM*4.0/5.0),lam[i])
          coord_pix10[i,*]=blue_map(instr,XMM,-(YMM+LENMM),lam[i])
      endfor
      xcords = [coord_pix0[*,1],coord_pix1[*,1],coord_pix2[*,1],coord_pix3[*,1],coord_pix4[*,1],$
                coord_pix5[*,1],coord_pix6[*,1],coord_pix7[*,1],coord_pix8[*,1],coord_pix9[*,1],$
                coord_pix10[*,1]]
      ycords = [coord_pix0[*,0],coord_pix1[*,0],coord_pix2[*,0],coord_pix3[*,0],coord_pix4[*,0],$
                coord_pix5[*,0],coord_pix6[*,0],coord_pix7[*,0],coord_pix8[*,0],coord_pix9[*,0],$
                coord_pix10[*,0]]
      wavestack = [lam,lam,lam,lam,lam,lam,lam,lam,lam,lam,lam]

   ;BLUE PRISM
   endif else if GRAT eq 'P450L' then begin
      lam=fltarr(2001) & coord_pix0=fltarr(2001,2) & coord_pix1=fltarr(2001,2) & coord_pix2=fltarr(2001,2)
      coord_pix3=fltarr(2001,2) & coord_pix4=fltarr(2001,2) & coord_pix5=fltarr(2001,2) & coord_pix6=fltarr(2001,2)
      coord_pix7=fltarr(2001,2) & coord_pix8=fltarr(2001,2) & coord_pix9=fltarr(2001,2) & coord_pix10=fltarr(2001,2)
      for i=0,2000 do begin
          lam[i]=3200+1.8*i
          coord_pix0[i,*]=blue_map_prism(XMM,-(YMM-LENMM),lam[i])
          coord_pix1[i,*]=blue_map_prism(XMM,-(YMM-LENMM*4.0/5.0),lam[i])
          coord_pix2[i,*]=blue_map_prism(XMM,-(YMM-LENMM*3.0/5.0),lam[i])
          coord_pix3[i,*]=blue_map_prism(XMM,-(YMM-LENMM*2.0/5.0),lam[i])
          coord_pix4[i,*]=blue_map_prism(XMM,-(YMM-LENMM*1.0/5.0),lam[i])
          coord_pix5[i,*]=blue_map_prism(XMM,-(YMM+LENMM*0),lam[i])
          coord_pix6[i,*]=blue_map_prism(XMM,-(YMM+LENMM*1.0/5.0),lam[i])
          coord_pix7[i,*]=blue_map_prism(XMM,-(YMM+LENMM*2.0/5.0),lam[i])
          coord_pix8[i,*]=blue_map_prism(XMM,-(YMM+LENMM*3.0/5.0),lam[i])
          coord_pix9[i,*]=blue_map_prism(XMM,-(YMM+LENMM*4.0/5.0),lam[i])
          coord_pix10[i,*]=blue_map_prism(XMM,-(YMM+LENMM),lam[i])
      endfor
      coord_pix0[*,0] += 2048 & coord_pix2[*,0] += 2048 & coord_pix3[*,0] += 2048 & coord_pix4[*,0] += 2048
       coord_pix5[*,0] += 2048 & coord_pix6[*,0] += 2048 & coord_pix7[*,0] += 2048 & coord_pix8[*,0] += 2048
       coord_pix9[*,0] += 2048 & coord_pix10[*,0] += 2048 & coord_pix1[*,0] += 2048
      xcords = [coord_pix0[*,1],coord_pix1[*,1],coord_pix2[*,1],coord_pix3[*,1],coord_pix4[*,1],$
                coord_pix5[*,1],coord_pix6[*,1],coord_pix7[*,1],coord_pix8[*,1],coord_pix9[*,1],$
                coord_pix10[*,1]]
      ycords = [coord_pix0[*,0],coord_pix1[*,0],coord_pix2[*,0],coord_pix3[*,0],coord_pix4[*,0],$
                coord_pix5[*,0],coord_pix6[*,0],coord_pix7[*,0],coord_pix8[*,0],coord_pix9[*,0],$
                coord_pix10[*,0]]
      wavestack = [lam,lam,lam,lam,lam,lam,lam,lam,lam,lam,lam]

   ;RED GRATING
   endif else if GRAT eq 'G670L' then begin
      lam=fltarr(4601) & coord_pix0=fltarr(4601,2) & coord_pix1=fltarr(4601,2) & coord_pix2=fltarr(4601,2)
      coord_pix3=fltarr(4601,2) & coord_pix4=fltarr(4601,2) & coord_pix5=fltarr(4601,2) & coord_pix6=fltarr(4601,2)
      coord_pix7=fltarr(4601,2) & coord_pix8=fltarr(4601,2) & coord_pix9=fltarr(4601,2) & coord_pix10=fltarr(4601,2)
      for i=0,4600 do begin
          lam[i]=5500+1.*i
          coord_pix0[i,*]=red_map(XMM,-(YMM-LENMM),lam[i])
          coord_pix1[i,*]=red_map(XMM,-(YMM-LENMM*4.0/5.0),lam[i])
          coord_pix2[i,*]=red_map(XMM,-(YMM-LENMM*3.0/5.0),lam[i])
          coord_pix3[i,*]=red_map(XMM,-(YMM-LENMM*2.0/5.0),lam[i])
          coord_pix4[i,*]=red_map(XMM,-(YMM-LENMM*1.0/5.0),lam[i])
          coord_pix5[i,*]=red_map(XMM,-(YMM+LENMM*0),lam[i])
          coord_pix6[i,*]=red_map(XMM,-(YMM+LENMM*1.0/5.0),lam[i])
          coord_pix7[i,*]=red_map(XMM,-(YMM+LENMM*2.0/5.0),lam[i])
          coord_pix8[i,*]=red_map(XMM,-(YMM+LENMM*3.0/5.0),lam[i])
          coord_pix9[i,*]=red_map(XMM,-(YMM+LENMM*4.0/5.0),lam[i])
          coord_pix10[i,*]=red_map(XMM,-(YMM+LENMM),lam[i])
      endfor
      xcords = [coord_pix0[*,1],coord_pix1[*,1],coord_pix2[*,1],coord_pix3[*,1],coord_pix4[*,1],$
                coord_pix5[*,1],coord_pix6[*,1],coord_pix7[*,1],coord_pix8[*,1],coord_pix9[*,1],$
                coord_pix10[*,1]]
      ycords = [coord_pix0[*,0],coord_pix1[*,0],coord_pix2[*,0],coord_pix3[*,0],coord_pix4[*,0],$
                coord_pix5[*,0],coord_pix6[*,0],coord_pix7[*,0],coord_pix8[*,0],coord_pix9[*,0],$
                coord_pix10[*,0]]
      wavestack = [lam,lam,lam,lam,lam,lam,lam,lam,lam,lam,lam]

   ;RED PRISM
   endif else if GRAT eq 'P700L' then begin
      lam=fltarr(2001) & coord_pix0=fltarr(2001,2) & coord_pix1=fltarr(2001,2) & coord_pix2=fltarr(2001,2)
      coord_pix3=fltarr(2001,2) & coord_pix4=fltarr(2001,2) & coord_pix5=fltarr(2001,2) & coord_pix6=fltarr(2001,2)
      coord_pix7=fltarr(2001,2) & coord_pix8=fltarr(2001,2) & coord_pix9=fltarr(2001,2) & coord_pix10=fltarr(2001,2)
      for i=0,2000 do begin
          lam[i]=5500+2.1*i	; changed the spacing to make sure the tested wavelengths fully
          coord_pix0[i,*]=red_map_prism(XMM,-(YMM-LENMM),lam[i])
          coord_pix1[i,*]=red_map_prism(XMM,-(YMM-LENMM*4.0/5.0),lam[i])
          coord_pix2[i,*]=red_map_prism(XMM,-(YMM-LENMM*3.0/5.0),lam[i])
          coord_pix3[i,*]=red_map_prism(XMM,-(YMM-LENMM*2.0/5.0),lam[i])
          coord_pix4[i,*]=red_map_prism(XMM,-(YMM-LENMM*1.0/5.0),lam[i])
          coord_pix5[i,*]=red_map_prism(XMM,-(YMM+LENMM*0),lam[i])
          coord_pix6[i,*]=red_map_prism(XMM,-(YMM+LENMM*1.0/5.0),lam[i])
          coord_pix7[i,*]=red_map_prism(XMM,-(YMM+LENMM*2.0/5.0),lam[i])
          coord_pix8[i,*]=red_map_prism(XMM,-(YMM+LENMM*3.0/5.0),lam[i])
          coord_pix9[i,*]=red_map_prism(XMM,-(YMM+LENMM*4.0/5.0),lam[i])
          coord_pix10[i,*]=red_map_prism(XMM,-(YMM+LENMM),lam[i])
      endfor
      coord_pix0[*,0] += 2048 & coord_pix2[*,0] += 2048 & coord_pix3[*,0] += 2048 & coord_pix4[*,0] += 2048
       coord_pix5[*,0] += 2048 & coord_pix6[*,0] += 2048 & coord_pix7[*,0] += 2048 & coord_pix8[*,0] += 2048
       coord_pix9[*,0] += 2048 & coord_pix10[*,0] += 2048 & coord_pix1[*,0] += 2048
      xcords = [coord_pix0[*,1],coord_pix1[*,1],coord_pix2[*,1],coord_pix3[*,1],coord_pix4[*,1],$
                coord_pix5[*,1],coord_pix6[*,1],coord_pix7[*,1],coord_pix8[*,1],coord_pix9[*,1],$
                coord_pix10[*,1]]
      ycords = [coord_pix0[*,0],coord_pix1[*,0],coord_pix2[*,0],coord_pix3[*,0],coord_pix4[*,0],$
                coord_pix5[*,0],coord_pix6[*,0],coord_pix7[*,0],coord_pix8[*,0],coord_pix9[*,0],$
                coord_pix10[*,0]]
      wavestack = [lam,lam,lam,lam,lam,lam,lam,lam,lam,lam,lam]
   endif

   ; write out test point wavelengths to file
   WRITE_CSV, 'coords.dat', xcords,ycords,wavestack,HEADER=['XCOORDS','YCOORDS','WAVELENGTH']
   ; use python code to fit 2D wavelength solution and print wavelengths for each pixel
   ;	within the vicinity of the slit

;     sky_line_file = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/uves_red.dat'
   spawnline  = 'python ' + GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/surfacefit.py'
   SPAWN, spawnline
   fitdata = READ_ASCII('fit.dat')
   fit = fitdata.(0)

   newt = fltarr(3088,8192)
   minx = round(min(xcords)) >0
   maxx = round(max(xcords)) <3088
   miny = round(min(ycords)) >0
   maxy = round(max(ycords)) <8192

   ; fit.dat coords are relative to the slit region -> transform to 8x3k coordinates
   for i=minx,maxx-1 do for j=miny,maxy-1 do newt[i,j] = fit[j-miny,i-minx] ;newt[i,j] = fit[i-minx,j-miny]

   aploc = where(slitimage eq n+1)
   empfit[aploc] = newt[aploc]
endfor

mwrfits,empfit,FILENAME,im_hd,/create
END
