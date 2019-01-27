pro mods2xy, tset, xpos, ypos, ignore_jump=ignore_jump

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - traceset2xy, tset, xpos, ypos'
      return
   endif

   if (tset.func EQ 'legendre') then begin
      ndim = size(tset.coeff, /n_dim)
      dims = size(tset.coeff, /dim)
      ncoeff = dims[0]
      nTrace = dims[1]
      nx = long(tset.xmax - tset.xmin + 1)
      xrange = tset.xmax - tset.xmin
      if (NOT keyword_set(xpos)) then $
        xpos = djs_laxisgen([nx, nTrace], iaxis=0) + tset.xmin
      if (median(tset.pix_end - tset.pix_beg) lt 1000) then begin
        xpos = dblarr(490,nTrace)
        for iTrace=0, nTrace-1 do begin
           steps = findgen(490)
	   xpos[*,iTrace] = tset.pix_beg[iTrace]
           xpos[*,iTrace] = xpos[*,iTrace] + steps[*]
        endfor
      endif
      ypos=0*xpos
      for iTrace=0, nTrace-1 do begin
	;now take slitend_x and figure out where to place this trace
	 ypos[*,iTrace]=polyleg(xpos[*,iTrace],tset.coeff[*,iTrace])
      endfor

;stop
   endif else begin
      error, 'Unknown function' + func
   endelse

   return
end

;------------------------------------------------------------------------------



