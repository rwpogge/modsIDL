;+
; NAME: 
;   RED_MAP
;
; VERSION:
;   1.0 (Nov, 2012)
;
; PURPOSE:
;   Convert MODS1 mask X,Y,lambda coordinates in mm and angstroms
;    to X,Y in pixels for the grating in the Red Channel.
;
; Reference:
;   Croxall, Pogge, and Stoll in prep.
;
; CALLING SEQUENCE:
;   coord=red_map(X_MMS,Y_MMS,LAMBDA)
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
;   coord_pix1[i,*]=red_map(XMM,YMM,red_lam[i])
;
; HISTORY:
;   2012-11 (K Croxall): Written
;   2019-01-15 (R. Pogge): Added MODS2 Red Grating Transforms
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;=========================================================================
;   coef_red -- return the coefficient from the cubic fit
;               Each coefficient (CL##, CXnm, and CYnm) are 
;               functions of wavelength (l): 
;               C...(l) = c0 + c1*l + c2*l^2 + .. + c_deg*l^deg,  
;               where deg = degree of polynomial
;=========================================================================
function coef_red,coef_name,lambda
  coef=coef_name[0] + coef_name[1]*lambda + $
     coef_name[2]*lambda^2 + coef_name[3]*lambda^3 + $
     coef_name[4]*lambda^4
  return,coef
end

;=========================================================================
;   red_map -- map maks positions from milimeters on mask 
;              to pixels on chip 
;              coeff for MODS1 RED GRATING
;=========================================================================
function red_map,instrume,x_mms,y_mms,lambda
   if not keyword_set(verbose) then !Quiet=1

if instrume eq 'MODS1R' then begin
	; Cubic Fit Coefficients
	; Red Grating - 2012 Nov 6 [rwp]
	;id     deg  c0            c1            c2            c3            c4
	CL11 = [-5608.06,       1.38367,      -2.343280E-05,  8.242579E-10,  6.537931E-15]
	CL21 = [ 1466.73,       1.536392E-02, -8.659779E-07,  3.966486E-11, 0.] 
	CL12 = [-15.6864,       7.517714E-04, -6.370730E-08,  2.276665E-12, -6.799757E-17]
	CL22 = [-0.220149,      5.171989E-05, -6.201528E-09,  2.286377E-13, 0.] 
	CL13 = [-4.380066E-03, -3.944823E-05,  5.427036E-09, -2.297828E-13, 0.] 
	CL23 = [ 14.1854,      -1.486915E-04,  8.050403E-09,  1.624253E-13, 0.] 
	CX00 = [ 18.5459,      -4.092292E-03,  3.010620E-07, -1.270393E-11, 0.] 
	CY00 = [ 2.17849,      -5.945472E-04,  3.637781E-08,  0.          , 0.] 
	CX10 = [ 1.489887E-02, -3.328264E-06,  4.140425E-10, -2.091303E-14, 0.] 
	CY10 = [ 1.079413E-03, -1.696765E-07,  0.          ,  0.          , 0.] 
	CX20 = [-4.416710E-03,  9.656396E-07, -7.049647E-11,  2.869618E-15, 0.] 
	CY20 = [-2.214831E-04,  5.833243E-08, -3.441173E-12,  0.          , 0.] 
	CX30 = [-8.013121E-07, -6.482562E-11,  3.524057E-14, -2.562849E-18, 0.] 
	CY30 = [-1.403150E-07,  1.691218E-11,  0.          ,  0.          , 0.] 
	CX01 = [-2.749564E-05,  0.          ,  0.          ,  0.          , 0.] 
	CY01 = [-1.232022E-02,  1.619966E-06, -6.901868E-11,  0.          , 0.] 
	CX11 = [ 7.252693E-04, -2.140472E-07,  1.810156E-11, -3.656839E-16, 0.] 
	CY11 = [ 2.321225E-03, -3.251117E-07,  3.040293E-12,  0.          , 0.] 
	CX21 = [ 7.180726E-07, -9.757197E-11,  0.          ,  0.          , 0.] 
	CY21 = [ 1.478757E-06, -2.538724E-10,  1.454514E-14,  0.          , 0.] 
	CX02 = [-9.358168E-04,  1.387448E-07,  0.          ,  0.          , 0.] 
	CY02 = [-6.867777E-04,  2.159598E-07, -1.966841E-11,  5.237666E-16, 0.] 
	CX12 = [-9.609731E-07,  6.129878E-11, -3.482034E-15,  0.          , 0.]  
	CY12 = [-9.799198E-07,  1.352260E-10,  0.          ,  0.          , 0.] 
	CX03 = [ 2.045765E-07, -3.194970E-11,  0.          ,  0.          , 0.] 
	CY03 = [ 7.616682E-07, -9.273407E-11,  6.445709E-15,  0.          , 0.] 
endif else if instrume eq 'MODS2R' then begin
	CL11 = [    -5408.1481,       1.3857633,  -2.5123214e-05,   1.0264197e-09,   0]
	CL12 = [    -15.359014,   0.00061073032,  -3.8272371e-08,   1.0241438e-13,   0]
	CL13 = [   0.038206724,   -2.109187e-05,   3.4756591e-09,  -1.6037695e-13,   0]
	CL21 = [     1479.0191,     0.015653695,  -7.6107598e-07,   3.0880804e-11,   0]
	CL22 = [  -0.081296966,   3.6875878e-05,  -4.6594491e-09,   1.8911839e-13,   0]
	CL23 = [     14.208232,  -0.00016481286,   1.0486678e-08,   6.8531116e-14,   0]
	CX00 = [     17.825428,   -0.0039671641,   2.9297353e-07,  -1.2471933e-11,   0]
	CX10 = [   0.034605948,  -6.1102746e-06,   2.4236849e-10,   0,               0]
	CX20 = [ -0.0038492159,   7.6780663e-07,  -4.4766074e-11,   1.7483305e-15,   0]
	CX30 = [ 1.3171303e-07,  -1.6925205e-10,    1.424859e-14,   0,               0]
	CX01 = [   0.042007972,  -1.0991443e-05,   6.8911118e-10,   0,               0]
	CX11 = [ 0.00047705571,  -1.2770262e-07,    7.850898e-12,   0,               0]
	CX21 = [ 4.4934767e-07,  -4.0705655e-11,   0,               0,               0]
	CX02 = [-0.00093114483,   1.3746708e-07,   0,               0,               0]
	CX12 = [-3.8948911e-07,  -1.5098315e-11,   0,               0,               0]
	CX03 = [-1.2031598e-06,   3.1036055e-10,  -1.6853548e-14,   0,               0]
	CY00 = [     2.1860562,  -0.00051023016,   5.5355384e-09,   2.1955639e-12,   0]
	CY10 = [  -0.023975992,    7.109443e-06,  -4.8889618e-10,   0,               0]
	CY20 = [-0.00025039299,     7.29838e-08,  -4.6062358e-12,   0,               0]
	CY30 = [-2.7345037e-07,   -5.811821e-11,   9.2369725e-15,   0,               0]
	CY01 = [  -0.017025844,   1.7426525e-06,   0,               0,               0]
	CY11 = [  0.0022362704,  -3.0494924e-07,   1.2793253e-12,   0,               0]
	CY21 = [ 1.3619796e-06,  -2.5070315e-10,   1.5589193e-14,   0,               0]
	CY02 = [-0.00072962547,   2.4376286e-07,  -2.2535025e-11,   6.1540353e-16,   0]
	CY12 = [-6.7776444e-07,  -1.2252571e-10,     2.05408e-14,   0,               0]
	CY03 = [ 5.3061378e-07,  -1.9922516e-11,   0,               0,               0]
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
