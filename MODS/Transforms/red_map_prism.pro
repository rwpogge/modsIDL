;+
; NAME: 
;   RED_MAP_PRISM
;
; VERSION:
;   1.0 (Jul, 2013)
;
; PURPOSE:
;   Convert MODS1 mask X,Y,lambda coordinates in mm and angstroms
;    to X,Y in pixels for the prism in the Red Channel.
;
; Reference:
;   Croxall, Pogge, and Stoll in prep.
;
; CALLING SEQUENCE:
;   coord=red_map_prism(X_MMS,Y_MMS,LAMBDA)
;
; INPUTS:
;   X_MMS: X position in milimeters from the MODS mask
;   Y_MMS: Y position in milimeters from the MODS mask
;   LAMBDA: The wavelength in angstroms
;
; OUTPUT:
;   The pixel coordinates are returned as a 2x1 array
;
; EXAMPLE:
;   coord_pix1[i,*]=red_map_prism(XMM,YMM,red_lam[i])
;
; HISTORY:
;   2013-07 (K Croxall): Written
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;=========================================================================
;   coef_red -- return the coefficient from the cubic fit
;               Each coefficient (CL##, CXnm, and CYnm) are 
;               functions of wavelength (l): 
;               C...(l) = c0 + c1*l + c2*l^2 + .. + c_deg*l^deg,  
;               where deg = degree of polynomial
;=========================================================================
function coef_red,coef_name,lambda
; Each coefficient (CL##, CXnm, and CYnm) are functions of wavelength (l): 
;   C...(l) = c0 + c1*l + c2*l^2 + .. + c_deg*l^deg,  where deg = degree of polynomial

coef=coef_name[0] + coef_name[1]*lambda + $
     coef_name[2]*lambda^2 + coef_name[3]*lambda^3 + $
     coef_name[4]*lambda^4
return,coef
end

;=========================================================================
;   red_map -- map maks positions from milimeters on mask 
;              to pixels on chip 
;              coeff for MODS1 RED PRISM
;=========================================================================
function red_map_prism,instrume,x_mms,y_mms,lambda

if instrume eq 'MODS1R' then begin
	; Cubic Fit Coefficients
	; Red Prism - 2013 Apr 7 [rwp]
	;id      c0             c1             c2             c3             c4
	CL11 = [ -4502.23,       2.63853,      -0.405848E-03,  0.290674E-07, -0.799361E-12]
	CL21 = [  1513.41,       0.196031E-01, -0.197078E-05,  0.698102E-10,  0.0]
	CL12 = [ -15.4620,       0.101764E-03, -0.608708E-08,  0.0,           0.0]
	CL22 = [ -0.131589,      0.484660E-05, -0.318392E-09,  0.0,           0.0]
	CL13 = [ -0.681868E-01, -0.396967E-05,  0.0,           0.0,           0.0]
	CL23 = [ 13.5016,        0.187593E-04, -0.833724E-09,  0.0,           0.0]
	CX00 = [  8.64803,      -0.250249E-03,  0.0,           0.0,           0.0]
	CY00 = [  0.491235,     -0.149508E-03,  0.857387E-08,  0.0,           0.0]
	CX10 = [  0.164200,     -0.423940E-04,  0.306330E-08,  0.0,           0.0]
	CY10 = [  0.154902E-01, -0.266857E-05,  0.742761E-10,  0.0,           0.0]
	CX20 = [ -0.163173E-02,  0.496125E-07,  0.0,           0.0,           0.0]
	CY20 = [ -0.505811E-05,  0.228396E-08,  0.0,           0.0,           0.0]
	CX30 = [ -0.121477E-04,  0.276148E-08, -0.174295E-12,  0.0,           0.0]
	CY30 = [ -0.997458E-06,  0.132022E-09,  0.0,           0.0,           0.0]
	CX01 = [ -0.286303E-01,  0.330329E-05,  0.0,           0.0,           0.0]
	CY01 = [ -0.121261E-01,  0.184746E-05, -0.957170E-10,  0.0,           0.0]
	CX11 = [ -0.616070E-04,  0.544271E-08,  0.0,           0.0,           0.0]
	CY11 = [  0.525682E-03, -0.122516E-06,  0.604853E-11,  0.0,           0.0]
	CX21 = [  0.629422E-06, -0.740222E-10,  0.0,           0.0,           0.0]
	CY21 = [  9.091423E-07,  0.0,           0.0,           0.0,           0.0]
	CX02 = [ -0.131047E-02,  0.537496E-07, -0.266836E-11,  0.0,           0.0]
	CY02 = [  3.115546E-05,  0.0,           0.0,           0.0,           0.0]
	CX12 = [ -1.058861E-06,  0.0,           0.0,           0.0,           0.0]
	CY12 = [ -0.140033E-05,  0.200057E-09,  0.0,           0.0,           0.0]
	CX03 = [ -0.223800E-06,  0.413940E-10,  0.0,           0.0,           0.0]
	CY03 = [  3.915830E-07,  0.0,           0.0,           0.0,           0.0]
endif else if instrume eq 'MODS2R' then begin
        print,'Coefficients not yet derived for MODS2 Prism Mode.'
        stop
endif else begin
        print,instrume,' is not a recognized instrument for the current optical  mapping.'
endelse


; Linear Part:
   cx1 = coef_red(CL11,lambda)
   cx2 = coef_red(CL12,lambda)
   cx3 = coef_red(CL13,lambda)

   cy1 = coef_red(CL21,lambda)
   cy2 = coef_red(CL22,lambda)
   cy3 = coef_red(CL23,lambda)

   x_lin  = cx1 + cx2*x_mms + cx3*y_mms
   y_lin  = cy1 + cy2*x_mms + cy3*y_mms

; Polynomial Part (the l dependence omitted for clarity):
   cpx0 = coef_red(CX00,lambda)
   cpx1 = coef_red(CX10,lambda)
   cpx2 = coef_red(CX20,lambda)
   cpx3 = coef_red(CX30,lambda)
   cpx4 = coef_red(CX01,lambda)
   cpx5 = coef_red(CX11,lambda)
   cpx6 = coef_red(CX21,lambda)
   cpx7 = coef_red(CX02,lambda)
   cpx8 = coef_red(CX12,lambda)
   cpx9 = coef_red(CX03,lambda)

   cpy0 = coef_red(CY00,lambda)
   cpy1 = coef_red(CY10,lambda)
   cpy2 = coef_red(CY20,lambda)
   cpy3 = coef_red(CY30,lambda)
   cpy4 = coef_red(CY01,lambda)
   cpy5 = coef_red(CY11,lambda)
   cpy6 = coef_red(CY21,lambda)
   cpy7 = coef_red(CY02,lambda)
   cpy8 = coef_red(CY12,lambda)
   cpy9 = coef_red(CY03,lambda)

   x_poly = cpx0 + cpx1*x_mms   + cpx2*x_mms^2       + cpx3*x_mms^3 $
                 + cpx4*y_mms   + cpx5*x_mms*y_mms   + cpx6*x_mms^2*y_mms $
                 + cpx7*y_mms^2 + cpx8*x_mms*y_mms^2 + cpx9*y_mms^3 

   y_poly = cpy0 + cpy1*x_mms   + cpy2*x_mms^2       + cpy3*x_mms^3 $
                 + cpy4*y_mms   + cpy5*x_mms*y_mms   + cpy6*x_mms^2*y_mms $
                 + cpy7*y_mms^2 + cpy8*x_mms*y_mms^2 + cpy9*y_mms^3 

coord=[x_lin+x_poly,y_lin+y_poly]

return,coord
end
