;+
; NAME: 
;   BLUE_MAP_PRISM
;
; VERSION:
;   1.0 (Jul, 2013)
;
; PURPOSE:
;   Convert MODS1 mask X,Y,lambda coordinates in mm and angstroms
;    to X,Y in pixels for the prism in the Blue Channel.
;
; Reference:
;   Croxall, Pogge, and Stoll in prep.
;
; CALLING SEQUENCE:
;   coord=blue_map_prism(X_MMS,Y_MMS,LAMBDA)
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
;   coord_pix1[i,*]=blue_map_prism(XMM,YMM,blue_lam[i])
;
; HISTORY:
;   2013-06 (K Croxall): Written
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;=========================================================================
;   coef_blue -- return the coefficient from the cubic fit
;               Each coefficient (CL##, CXnm, and CYnm) are 
;               functions of wavelength (l): 
;               C...(l) = c0 + c1*l + c2*l^2 + .. + c_deg*l^deg,  
;               where deg = degree of polynomial
;=========================================================================
function coef_blue,coef_name,lambda
coef=coef_name[0] + coef_name[1]*lambda + $
     coef_name[2]*lambda^2 + coef_name[3]*lambda^3 + $
     coef_name[4]*lambda^4
return,coef
end

;=========================================================================
;   blue_map_prism -- map maks positions from milimeters on mask 
;              to pixels on chip 
;              coeff for MODS1 BLUE PRISM
;=========================================================================
function blue_map_prism,instrume,x_mms,y_mms,lambda

if instrume eq 'MODS1B' then begin
	; Cubic Fit Coefficients
	; Blue Prism - 2013 Feb 20 [rwp]
	;id      c0             c1             c2             c3             c4
	CL11 = [-1751.90,       2.15974,      -0.472961E-03,  0.478632E-07, -0.183714E-11]
	CL21 = [ 1496.23,       0.585859E-01, -0.951949E-05,  0.534407E-09,  0.0]
	CL12 = [ 15.1790,      -0.155977E-04,  0.0,           0.0,           0.0]
	CL22 = [ 0.719607E-01, -0.255662E-05,  0.0,           0.0,           0.0]
	CL13 = [-4.177400E-02,  0.0,           0.0,           0.0,           0.0]
	CL23 = [ 13.7901,       0.459839E-04, -0.386222E-08,  0.0,           0.0]
	CX00 = [ 5.32895,       0.0,           0.0,           0.0,           0.0]
	CY00 = [ 0.338659,     -0.100309E-03,  0.0,           0.0,           0.0]
	CX10 = [-1.653185E-02,  0.0,           0.0,           0.0,           0.0]
	CY10 = [-0.835183E-01,  0.316848E-04, -0.294624E-08,  0.0,           0.0]
	CX20 = [-0.130417E-02,  0.516041E-07,  0.0,           0.0,           0.0]
	CY20 = [-0.815897E-04,  0.157349E-07,  0.0,           0.0,           0.0]
	CX30 = [ 0.280322E-05, -0.347145E-09,  0.0,           0.0,           0.0]
	CY30 = [ 0.582325E-06, -0.131084E-09,  0.0,           0.0,           0.0]
	CX01 = [-0.138940E-01,  0.327814E-05,  0.0,           0.0,           0.0]
	CY01 = [-0.273931E-01,  0.420186E-05,  0.0,           0.0,           0.0]
	CX11 = 	[-0.200845E-03,  0.445991E-07,  0.0,           0.0,           0.0]
	CY11 = [-0.232001E-03,  0.409042E-07,  0.0,           0.0,           0.0]
	CX21 = [ 0.119599E-05, -0.202171E-09,  0.0,           0.0,           0.0]
	CY21 = [ 0.683855E-05, -0.111251E-08,  0.0,           0.0,           0.0]
	CX02 = [-0.725448E-03, -0.245101E-07,  0.0,           0.0,           0.0]
	CY02 = [-0.877680E-04,  0.274468E-07,  0.0,           0.0,           0.0]
	CX12 = [ 2.574787E-07,  0.0,           0.0,           0.0,           0.0]
	CY12 = [ 0.140328E-05, -0.289865E-09,  0.0,           0.0,           0.0]
	CX03 = [ 0.130341E-05, -0.267100E-09,  0.0,           0.0,           0.0]
	CY03 = [ 0.214046E-05, -0.311006E-09,  0.0,           0.0,           0.0]
endif else if instrume eq 'MODS2B' then begin
        print,'Coefficients not yet derived for MODS2 Prism Mode.'
        stop
endif else begin
        print,instrume,' is not a recognized instrument for the current optical  mapping.'
endelse

; Linear Part:
   cx1 = coef_blue(CL11,lambda)
   cx2 = coef_blue(CL12,lambda)
   cx3 = coef_blue(CL13,lambda)

   cy1 = coef_blue(CL21,lambda)
   cy2 = coef_blue(CL22,lambda)
   cy3 = coef_blue(CL23,lambda)

   x_lin  = cx1 + cx2*x_mms + cx3*y_mms
   y_lin  = cy1 + cy2*x_mms + cy3*y_mms

; Polynomial Part (the l dependence omitted for clarity):
   cpx0 = coef_blue(CX00,lambda)
   cpx1 = coef_blue(CX10,lambda)
   cpx2 = coef_blue(CX20,lambda)
   cpx3 = coef_blue(CX30,lambda)
   cpx4 = coef_blue(CX01,lambda)
   cpx5 = coef_blue(CX11,lambda)
   cpx6 = coef_blue(CX21,lambda)
   cpx7 = coef_blue(CX02,lambda)
   cpx8 = coef_blue(CX12,lambda)
   cpx9 = coef_blue(CX03,lambda)

   cpy0 = coef_blue(CY00,lambda)
   cpy1 = coef_blue(CY10,lambda)
   cpy2 = coef_blue(CY20,lambda)
   cpy3 = coef_blue(CY30,lambda)
   cpy4 = coef_blue(CY01,lambda)
   cpy5 = coef_blue(CY11,lambda)
   cpy6 = coef_blue(CY21,lambda)
   cpy7 = coef_blue(CY02,lambda)
   cpy8 = coef_blue(CY12,lambda)
   cpy9 = coef_blue(CY03,lambda)

   x_poly = cpx0 + cpx1*x_mms   + cpx2*x_mms^2       + cpx3*x_mms^3 $
                 + cpx4*y_mms   + cpx5*x_mms*y_mms   + cpx6*x_mms^2*y_mms $
                 + cpx7*y_mms^2 + cpx8*x_mms*y_mms^2 + cpx9*y_mms^3 

   y_poly = cpy0 + cpy1*x_mms   + cpy2*x_mms^2       + cpy3*x_mms^3 $
                 + cpy4*y_mms   + cpy5*x_mms*y_mms   + cpy6*x_mms^2*y_mms $
                 + cpy7*y_mms^2 + cpy8*x_mms*y_mms^2 + cpy9*y_mms^3 

coord=[x_lin+x_poly,y_lin+y_poly]

return,coord
end
