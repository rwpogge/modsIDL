function gaussmods, x, p1,p2,p3, skew=skew, peak=peak, _EXTRA=extra

  sz = size(x)
  if sz(sz(0)+1) EQ 5 then smax = 26D else smax = 13.

  p = [p1,p2,p3]

  if n_elements(p) GE 3 then norm = p(2) else norm = x(0)*0 + 1

  u = ((x-p(0))/(abs(p(1)) > 1e-20))^2
  mask = u LT (smax^2)
  if NOT keyword_set(peak) then norm = norm / (sqrt(2.D * !dpi)*p(1))
  f = norm * mask * exp(-0.5*temporary(u) * mask)
  mask = 0

  if n_elements(skew) GT 0 then $
    f = (1.D + skew * (x-p(0))/p(1))*f

  return, f
end
;===============================================================

pro mods_flexurecorr,wave,spec,sky,var,colstr,delta=delta,slope = slope, $
	z=z,PLOT_PROGRESS=plot_prog,force=force,noslope=noslope

  plot_progress=keyword_set(plot_prog)

  if colstr eq 'blue' then col=1
  if colstr eq 'red' then col=0
  if ((col ne 1) and (col ne 0)) then begin
        print,'Not a valid color option.  Please use blue or red'
        stop
  endif

  pix_arr = indgen(n_elements(wave))

  ;trim to regions where we can find lines quickly
  if col then begin
          trim = where((wave le 5610.0) and (wave ge 4000))
  endif else begin
          trim = where((wave ge 5550.0) and (wave le 10000))
  endelse

  ;wave model
  wtrim = where((wave gt 2500.) and (wave lt 10000))
  expr = 'P[0] + P[1]*X'
  P1f = mpfitexpr(expr,pix_arr[wtrim],wave[wtrim],yfit=lin_fit,/quiet)
  sub = wave[wtrim] - lin_fit
;     plot,pix_arr[wtrim],sub,xstyle=1,position=[0.12,0.52,0.98,0.98]
  P3f = poly_fit(pix_arr[wtrim],sub,3,yfit=cub_fit)
     subsub = sub - cub_fit
;     oplot,pix_arr[wtrim],subsub,color=cgcolor('red')
;     test_pixel = 4000.D
;     test_wave = (P1f[0] + P1f[1]*test_pixel) + $
;	(P3f[0] + P3f[1]*test_pixel + P3f[2]*test_pixel^2 + P3f[3]*test_pixel^3)
;     print,test_pixel,test_wave

  expr = 'P[0] + P[1]*X'
  P1b = mpfitexpr(expr,wave[wtrim],pix_arr[wtrim],yfit=lin_fit,/quiet)
  sub = pix_arr[wtrim] - lin_fit
;     plot,wave[wtrim],sub,xstyle=1,position=[0.12,0.02,0.98,0.48],/noerase
  P3b = poly_fit(wave[wtrim],sub,3,yfit=cub_fit)
     subsub = sub - cub_fit
;     oplot,wave[wtrim],subsub,color=cgcolor('red')
;     test_pixel = (P1b[0] + P1b[1]*test_wave) + $
;        (P3b[0] + P3b[1]*test_wave + P3b[2]*test_wave^2 + P3b[3]*test_wave^3)
;     print,test_pixel,test_wave     
;     stop

     ;read in sky lines from UVES
  if col then begin
     sky_line_file = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/uves_blue.dat'
  endif else begin
     sky_line_file = GETENV('XIDL_DIR') + '/Spec/Longslit/pro/LBT/MODS/uves_red.dat'
  endelse

  readcol,sky_line_file,sky_lines

  ;if force, skip the fitting
  if keyword_set(force) then begin
        delta = force(0)
	slope = [force(1),force(2)]
        print,'Adopting a pixel shift of ',delta,' pixels as instructed'
        goto,Forced
  endif

  REPEAT BEGIN  ;iterate to clean out outlying lines
     ;determine the pixel positions of the sky lines
     orig_pixels = (P1b[0] + P1b[1]*sky_lines) + $
         (P3b[0] + P3b[1]*sky_lines + P3b[2]*sky_lines^2 + P3b[3]*sky_lines^3)

     ;make a model
     expr = 'P[0]' ; flat continuum
     off = 3
     ; Add lines
     for i=0,n_elements(sky_lines)-1 do begin
        segment = '+GAUSSMODS(X,P[' $
      	   + strcompress(string(off)) + '],P[1],P[' + strcompress(string(off+1)) + '])'
        expr = expr + segment
        off+=2
     endfor

     ;initial parameters
     start = fltarr(2*n_elements(sky_lines)+3)
     start[0] = min(sky[trim]) &  start[1] = 1.5D & start[2] = 0.D              ; continuum,FWHM,Vsys
     for i=0,n_elements(sky_lines)-1 do start[2*i+3] = orig_pixels[i] ; initialize wavelength
     for i=0,n_elements(sky_lines)-1 do start[2*i+4] = 500.D         ; initial flux
     pi = replicate({fixed: 0, limited: [0,0],limits: [0.D,0.D]},2*n_elements(sky_lines)+3)
;     pi(0).fixed = 1
     pi(1).limited = 1 & pi(1).limits = [1.,3]                       ; limits on the FWHM
     for i=0,n_elements(sky_lines)-1 do pi(2*i+4).limited[0,*] = 1
     for i=0,n_elements(sky_lines)-1 do pi(2*i+4).limits[0,*]=0.0D   ;non-negative flux
   
     ;fit strong emission lines
     result = MPFITEXPR(expr,pix_arr[trim],sky[trim],var[trim],start,yfit=yfit,parinfo=pi,/quiet)

     ;find the shift
     delta_arr = fltarr(n_elements(sky_lines))
     for i=0,n_elements(sky_lines)-1 do delta_arr[i] = orig_pixels[i]-result[2*i+3]
     delta = mean(delta_arr)
     residuals = delta_arr - delta

     ;linear component
     expr = 'P[0] + P[1]*X'
     start = [0,0]
     if keyword_set(noslope) then begin
	expr = 'P[0]'
	start = [0]
     endif
;     slope = mpfitexpr(expr,orig_pixels,delta_arr,/quiet,perror=slope_err)
     slope = mpfitexpr(expr,orig_pixels,residuals,/quiet,perror=slope_err,start=start)
     if keyword_set(noslope) then begin
	new_slope = [slope,0.D]
	slope = new_slope
     endif
     delta_arr_lin = delta_arr - (slope(0)+orig_pixels*slope(1))
     residuals_lin = delta_arr_lin - delta

     ;find outliers and trim the line list
     outliers = where(abs(residuals_lin) gt 2.6*stddev(residuals_lin))
     if outliers[0] ne -1 then begin
	   print,'Shift of ',delta,' pixels.  Linear components slope is ',slope[1]
	   print,'Outliers detected, iterating.'
	   print,'Rejecting ',n_elements(outliers),' lines.'
	   good = where(abs(residuals_lin) lt 2.5*stddev(residuals_lin))
	   sky_lines = sky_lines[good]
     endif

     if plot_progress then begin
         if col then begin
            w_max = 5600 & w_min = 3100
	    p_max = 7200 & p_min = 2000
         endif else begin
            w_max = 10000 & w_min = 5500
            p_max = 7200 & p_min = 1000
         endelse
;         plot,wave,spec,yrange=[-100,0.25*max(spec)], $
;            ystyle=1,Position=[0.06,0.06,0.70,0.33], xstyle=1,xrange=[w_min,w_max]
;         oplot,wave,sky/max(sky)*0.25*max(spec),color=cgColor('Dodger Blue')
;         for i =0,n_elements(sky_lines)-1 do oplot,[sky_lines[i],sky_lines[i]],[-100,1e6],linestyle=1

         plot,pix_arr[trim],sky[trim],psym=0,$;/noerase,$
             Position=[0.06,0.65,0.98,0.99],xstyle=1,xtitle='Pixels', $
	     title='FIT',XTickFormat = "(A1)",xrange=[p_min,p_max],/ylog,ystyle=1
         oplot,pix_arr[trim],yfit,color=cgColor('Dodger Blue')
         for i =0,n_elements(sky_lines)-1 do $
             oplot,[orig_pixels[i],orig_pixels[i]],[-100,1e6],linestyle=1
	 pinfo = 'Median ' + cgSymbol('Delta') + 'x = ' + string(delta) + $
		string(177b) + string(stddev(residuals)/sqrt(n_elements(sky_lines))) + ' pixels'
	 xyouts,pix_arr[trim[100]],2.5*median(sky[trim]),pinfo,charsize=1.5

         plot,orig_pixels,residuals_lin,psym=6,/noerase,xtitle='Pixels', ytitle='Residuals',$
               xstyle=1,Position=[0.06,0.40,0.98,0.65],xrange=[p_min,p_max]
         oplot,[0,1e5],[0,0],linestyle=1
	 pinfo = cgSymbol('sigma') + ' = ' + string(stddev(residuals_lin)) + ' pixels' 
	 xyouts,(p_min+500),0.5*stddev(residuals_lin),pinfo,charsize=2
;	 oplot,pix_arr[trim],subsub,color=cgcolor('red')
         pinfo = 'Rejecting ' + string(n_elements(outliers)) + ' line(s) and iterating'
         if outliers[0] ne -1 then xyouts,(p_min+500),-0.5*stddev(residuals),pinfo,charsize=2

     endif

  ENDREP UNTIL outliers[0] eq -1

; set the flexure and recast the wavelength array
  FORCED:print,'Flexure has caused a wavelength shift of ',delta,' pixel(s)'

;  new_pixels = indgen(n_elements(wave)) + delta
;  new_wave = (P1f[0] + P1f[1]*new_pixels) + $
;      (P3f[0] + P3f[1]*new_pixels + P3f[2]*new_pixels^2 + P3f[3]*new_pixels^3)

  new_pixels = indgen(n_elements(wave)) + delta + (slope(0)+indgen(n_elements(wave))*slope(1))
  new_wave = (P1f[0] + P1f[1]*new_pixels) + $
      (P3f[0] + P3f[1]*new_pixels + P3f[2]*new_pixels^2 + P3f[3]*new_pixels^3)

  if ((plot_progress) and (not keyword_set(force))) then begin
	if col then begin
            w_max=5600 & w_min = 5540
         endif else begin
            w_max=9910 & w_min = 9700
         endelse
         plot,wave,sky,position=[0.72,0.06,0.98,0.33], title='Blue is the CORRECTED sky', $
            ystyle=1,xstyle=1,xrange=[w_min,w_max],xtitle='Wavelength',/noerase
         for i =0,n_elements(sky_lines)-1 do oplot,[sky_lines[i],sky_lines[i]],[-100,1e6],linestyle=1
         oplot,new_wave,sky,color=cgColor('Dodger Blue')
;     oplot,new_wave,spec,color=cgColor('Dodger Blue')

        if col then begin
            w_max=4070 & w_min = 4020
         endif else begin
            w_max=5600 & w_min = 5550
         endelse
         plot,wave,sky,position=[0.05,0.06,0.18,0.33], title='Blue is the CORRECTED sky', $
            ystyle=1,xstyle=1,xrange=[w_min,w_max],xtitle='Wavelength',/noerase
         for i =0,n_elements(sky_lines)-1 do oplot,[sky_lines[i],sky_lines[i]],[-100,1e6],linestyle=1
         oplot,new_wave,sky,color=cgColor('Dodger Blue')

        if col then begin
            w_max=4380 & w_min = 4330
         endif else begin
            w_max=6350 & w_min = 6250
         endelse
         plot,wave,sky,position=[0.20,0.06,0.38,0.33], title='Blue is the CORRECTED sky', $
            ystyle=1,xstyle=1,xrange=[w_min,w_max],xtitle='Wavelength',/noerase
         for i =0,n_elements(sky_lines)-1 do oplot,[sky_lines[i],sky_lines[i]],[-100,1e6],linestyle=1
         oplot,new_wave,sky,color=cgColor('Dodger Blue')

        if col then begin
            w_max=5480 & w_min = 5440
         endif else begin
            w_max=7880 & w_min = 7700
         endelse
         plot,wave,sky,position=[0.41,0.06,0.70,0.33], title='Blue is the CORRECTED sky', $
            ystyle=1,xstyle=1,xrange=[w_min,w_max],xtitle='Wavelength',/noerase
         for i =0,n_elements(sky_lines)-1 do oplot,[sky_lines[i],sky_lines[i]],[-100,1e6],linestyle=1
         oplot,new_wave,sky,color=cgColor('Dodger Blue')
  endif

  Ok=0
  while Ok ne 1 do begin
  print,'Adopt flexure correction? (y/n/bye)'
  permit = ''
  permit = strlowcase(permit)
  read,': ',permit
    if permit eq 'y' then begin
	wave = new_wave
	Ok = 1
    endif else if permit eq 'n' then begin
	wave = wave
	delta = 0.
	slope = 0.
	Ok = 1
    endif else if permit eq 'bye' then begin
	stop
    endif else begin
	print,'I do not understand ',permit
    endelse
  endwhile
end


