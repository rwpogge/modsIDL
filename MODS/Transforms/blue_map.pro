;+
; NAME: 
;   BLUE_MAP
;
; VERSION:
;   1.0 (Nov, 2012)
;
; PURPOSE:
;   Convert MODS1 mask X,Y,lambda coordinates in mm and angstroms
;    to X,Y in pixels for the grating in the Blue Channel.
;
; Reference:
;   Croxall, Pogge, and Stoll in prep.
;
; CALLING SEQUENCE:
;   coord=blue_map(X_MMS,Y_MMS,LAMBDA)
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
;   coord_pix1[i,*]=blue_map(XMM,YMM,blue_lam[i])
;
; HISTORY:
;   2012-11 (K Croxall): Written
;   2019-01-15 (R Pogge): Added MODS2 Blue Grating Transforms
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
;   blue_map -- map maks positions from milimeters on mask 
;              to pixels on chip 
;              coeff for MODS1 BLUE GRATING
;=========================================================================
function blue_map,instrume,x_mms,y_mms,lambda
   if not keyword_set(verbose) then !Quiet=1

if instrume eq 'MODS1B' then begin
	; Cubic Fit Coefficients
	; Blue Grating - 2012 Nov 23 [rwp]
	;id      c0             c1             c2             c3            c4
	CL11 = [-3911.33,	1.32579	,      0.240512E-03, -0.400448E-07, 0.242569E-11]
	CL21 = [1510.13,	0.484731E-01, -0.696644E-05,  0.456536E-09, 0]
	CL12 = [28.0541,       -0.120617E-01,  0.376213E-05, -0.532120E-09, 0.287103E-13]
	CL22 = [0.564021E-01,  -0.146439E-04,  0.171790E-08,  0,            0]
	CL13 = [-0.696389E-01,	0.301566E-05,  0.101929E-08,  0,            0]
	CL23 = [14.8234,       -0.505605E-03,  0.865401E-07, -0.418654E-11, 0]
	CX00 = [21.5439,       -0.860033E-02,  0.971777E-06, -0.442110E-10, 0]
	CY00 = [ 9.27603,      -0.375665E-02,  0.371047E-06,  0,            0]
	CX10 = [ 0.332801,     -0.157520E-03,  0.181682E-07,  0,            0]
	CY10 = [-0.413914E-02,	0.102603E-05,  0,             0,            0]
	CX20 = [-0.296153E-02,	0.724482E-06,  0,             0,            0]
	CY20 = [-0.112912E-02,	0.437223E-06, -0.419629E-10,  0,            0]
	CX30 = [ 0.230899E-05, -0.351801E-09,  0,             0,            0]
	CY30 = [ 0.500026E-05, -0.213028E-08,  0.222690E-12,  0,            0]
	CX01 = [ 6.899079E-04,	0           ,  0,             0,            0]
	CY01 = [-0.985898E-02,	0.458420E-06,  0,             0,            0]
	CX11 = [ 0.235346E-04,	0.196159E-09,  0,             0,            0]
	CY11 = [-0.288128E-02,	0.776892E-06, -0.282843E-10,  0,            0]
	CX21 = [ 0.811930E-06, -0.183961E-09,  0,             0,            0]
	CY21 = [ 0.416787E-05, -0.125967E-08,  0.118710E-12,  0,            0]
	CX02 = [-0.916052E-03,	0.253423E-06,  0,             0,            0]
	CY02 = [-0.153048E-02,	0.627725E-06, -0.613007E-10,  0,            0]
	CX12 = [ 0.189239E-05, -0.192270E-09,  0,             0,            0]
	CY12 = [-1.332816E-07,	0,             0,             0,            0]
	CX03 = [-8.338054E-09,	0,             0,             0,            0]
	CY03 = [ 1.044642E-06,	0,             0,             0,            0]
endif else if instrume eq 'MODS2B' then begin
	CL11 = [    -4989.9187,       2.2336563,  -6.2773959e-05,   4.3958280e-09,   0]
	CL12 = [     16.272753,   -0.0014755819,   2.1445927e-07,  -8.4209753e-12,   0]
	CL13 = [   -0.11235339,   2.8399523e-05,  -2.0112681e-09,   0,               0]
	CL21 = [     1609.0538,    -0.085802248,   2.2864547e-05,  -1.5894026e-09,   0]
	CL22 = [   0.058939785,  -1.7968558e-05,   1.4973963e-09,   0,               0]
	CL23 = [     15.123112,  -0.00071810602,   1.3042698e-07,  -7.1169874e-12,   0]
	CX00 = [     6.0871915,    0.0018401508,  -1.2558305e-06,   1.0302182e-10,   0]
	CX10 = [ -0.0031576554,  -4.4028047e-07,   0,               0,               0]
	CX20 = [ -0.0033896522,    9.126697e-07,  -2.0914868e-11,   0,               0]
	CX30 = [ 5.6047985e-06,    -2.10439e-09,   2.0399361e-13,   0,               0]
	CX01 = [  -0.017202856,    9.791781e-06,  -1.2368746e-09,   0,               0]
	CX11 = [-0.00040657296,   1.8799025e-07,  -1.9700917e-11,   0,               0]
	CX21 = [-3.1963620e-07,   0,               0,               0,               0]
	CX02 = [-0.00094075154,   2.5813338e-07,   0,               0,               0]
	CX12 = [ 1.8253681e-06,  -2.9209507e-10,   0,               0,               0]
	CX03 = [ 4.2997443e-06,  -2.0790681e-09,   2.3318142e-13,   0,               0]
	CY00 = [     2.8955724,   -0.0011987761,   1.0534182e-07,   0,               0]
	CY10 = [   0.046310285,  -1.8928747e-05,   1.9443019e-09,   0,               0]
	CY20 = [-0.00045651284,   1.8237831e-07,  -1.7139235e-11,   0,               0]
	CY30 = [-3.9561653e-06,   1.5822761e-09,  -1.5817791e-13,   0,               0]
	CY01 = [  0.0039123931,   -2.651954e-06,   0,               0,               0]
	CY11 = [  -0.002324546,   5.2278450e-07,   0,               0,               0]
	CY21 = [ 2.2060578e-07,   9.2581236e-11,   0,               0,               0]
	CY02 = [-0.00049938045,   2.2147321e-07,  -1.9869726e-11,   0,               0]
	CY12 = [-2.6707218e-06,     1.27038e-09,  -1.5934599e-13,   0,               0]
	CY03 = [ 9.2817317e-07,   0,               0,               0,               0]

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
