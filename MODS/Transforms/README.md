# MODS Sudoku Transforms

R. Pogge & R. Stoll
2012 Sept 13

The purpose of the Sudoku mas is to measure how wavelengths map onto the MODS science detectors from a given (x,y)
mask location.  The apertures are 240um round pinholes which project to about 3-pixels on the science CCDs.  Pinholes 
are used because to first order any given spectral line will map to a "star" readily detected using source-finding 
algorithms, for example, SExtractor.

## The Sudoku Mask

The Sudoku Mask is a 2D mask meant to do a full empirical mapping of the mask (x,y) onto CCD (x,y) pixels at any given 
wavelength, lambda. For this mask, the pinhole placement is subject to the constraint that dispersed spectra from the 
pinholes must not overlap, and we want to sample the space in the X direction uniformly. If we impose a simple
echelon type pattern to avoid the vertical overlap we will introduce structured aliases into the eventual transformations.

OSU graduate student Rebecca Stoll noticed that the combined positional constraints were very similar to the rules 
of the Japanese combinatorial number-placement puzzle Sudoku.

A Sudoku matrix is a 9x9 matrix with each element containing an integer from 1 to 9 following these rules:
<ol>
<li>each integer can only appear once along a given row
<li>integer can only appear once down any given column
<li>each integer can only appear once within each of 9 3x3 subgrids that make up the overall grid
</ol>
The third rule of no repeats within the 3x3 subgrids is what makes the family of Sudoku matrices a subset of 9x9 Latin 
squares (which are created following the first two rules). Mathematicians have shown that there are ~6.7x10^21 
possible Sudoku matrixes (Felgenhauer & Jarvis 2006). 

Shannon entropy quantifies how random an ensemble of numbers is, developed by Claude Shannon, one of the founders of
information theory (see the seminal paper: Shannon, C., 1948, "A Mathematical Theory of Communication", Bell System 
Technical Journal, 27, 3).  Compared to generating a 9x9 grid of random numbers (of which there are 9^81
possible realizations), the Sudoku matrix has a higher Shannon entropy by about 3 orders of magnitude (see Newton & DeSalvo
2010, Proc Roy Soc A, DOI: 10.1098.rspa.2009.0052).  Greater Shannon entropy means greater intrinsic randomness
(in information theory, Shannon entropy estimates the number of data bits required to encode a string of
symbols calcluated based on the frequency of occurence of those symbols in the string).

We therefore created a 2D pinhole mask by finding a standard 9x9 Sudoku puzzle, and then used that number pattern to 
generate the pinhole positions following a simple algorithm for how to convert the number in the square to the 
location of the test pinhole along the diagonal within the cell. A couple of simple realizations showed that
this generates a mask that when used with the grating or prism disperser and a comparison lamp will give us a 
good, randomized coverage of the available science detector area to measure the transformations.

## Sudoku Transforms

We start by considering a pattern of spots from a single wavelength, lambda.  The 2D
transformation between the slit mask (x,y) pinhole positions in mm and CCD (X_ccd,Y_ccd) pixels
for these monochromatic spots is a 3rd-order coordinate transformation with linear and
polynomial components:
<pre>
Linear:
   X_ccd,lin = CL11 + CL12*x + CL13*y
   Y_ccd,lin = CL21 + CL22*x + CL23*y

Polynomial:
   X_ccd,pol = CX00 + CX10*x   + CX20*x^2   + CY30*x^3
                    + CX01*y   + CX11*x*y   + CX21*x^2*y
                    + CX02*y^2 + CX12*x*y^2 + CX03*y^3
   Y_ccd,pol = CY00 + CY10*x   + CY20*x^2   + CY30*x^3
                    + CY01*y   + CY11*x*y   + CY21*x^2*y
                    + CY02*y^2 + CY12*x*y^2 + CY03*y^3
</pre>
for a total of 6 linear and 20 polynomial coefficients.  The fitting is done
using the IRAF geomap task.

We then repeat this fitting procedure for a number of wavelengths convering the spectral range 
of the MODS channel mode (e.g., MODS2 Red Grating).  Each wavelength line has a different set
of linear and polynomial coordinate transformation coefficients.

We then fit this ensemble of coefficients as polynomials in wavelength, for example:
<pre>
   CX00(lambda) = c0 + c1*lambda + c2*lambda^2 + ... + cn*lambda^n
</pre>
where typically n=2 or n=3 is the most highest polynomial order that best fits the data.

These together give us a set of Sudoku Transforms of the form:
<pre>
    X_ccd = X_ccd,lin(lambda,x,y) + X_ccd,pol(lambda,x,y0
    Y_ccd = Y_ccd,lin(lambda,x,y) + Y_ccd,pol(lambda,x,y0
</pre>
That are then used to map a given location on the slit mask into CCD pixels as a function of wavelength.

The mapping is very robust in practice. Modulo small shifts because of a combination of
how accurately the masks are clamped in their mounting frames and small residual flexture in the spectrograph
channels not corrected out by the MODS active flexure compensation system, these transforms predictively map
where spectra will land on the CCD for a given slit mask position.

The application to MOS reduction is clear. While we did OK with more or less open-ended schemes to perform MOS mapping
or reductions, they lacked robustness.  More sophisticated algorithms (e.g., Carnegie COSMOS) attempt to improve this 
by doing forward optical modeling, but their internal models all assumed the optics are on-axis and therefore
optical distortions are radially symmetric.  MODS uses off-axis optics, so we have non axially-symmetric distortions, m
eaning the models have higher order asymmetric stretches (we cannot quantify the polynomials only in terms of 
even powers of x and y, but must include odd powers and xy cross terms). The xidl package that we were
evaluating at the time was much less sophisticated about its approach (with the penalty that it performed much worse
unmodified), but its simplicity allowed us greater latitude to break into the code and add a modeling approach based 
on the Sudoku mask transforms that, in the end, made it very robust and quite a bit faster as the guesses for
trace extraction and wavelength calibration were much closer to the mark on the first iteration, reducing the number of 
subsequent iterations and time lost to runaway divergences.

## Transform Scripts
<pre>
   blue_map.pro       - Blue Channel Grating mode Sudoku transforms
   red_map.pro        - Red Channel Grating mode Sudoku transforms
   blue_map_prism.pro - Blue Channel Prism mode Sudoku transforms
   red_map_prism.pro  - Red Channel Prism mode Sudoku transforms
</pre>
Each script contains, hardcoded, the transform coefficients for each instrument instance (MODS1 and MODS2).
