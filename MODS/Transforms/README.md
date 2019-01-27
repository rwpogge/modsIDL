# MODS Sudoku Transforms

R. Pogge & R. Stoll
2012 Sept 13

The purpose of the Sudoku and Pinhole Slit masks is to measure how
wavelengths map onto the MODS science detectors from a given (x,y)
mask location.  The apertures are 240um round pinholes which project
to about 3-pixels on the science CCDs.  Pinholes are used because to
first order any given spectral line will map to a "star" readily
detected using source-finding algorithms, for example, SExtractor.

## The Pinhole Slit

The Pinhole Slit is used to make a quick-look transform for the
facility long slit masks, and consists of an array of 17 pinholes
spaced every 12mm about (0,0) along the x=0 line of a typical centered
long slit.

## The Sudoku Mask

The Sudoku Mask is a 2D mask meant to do a full empirical mapping of
the mask (x,y) onto CCD (x,y) pixels at any given wavelength, lambda.
For this mask, the pinhole placement is subject to the constraint that
dispersed spectra from the pinholes must not overlap, and we want to
sample the space in the X direction uniformly. If we impose a simple
echelon type pattern to avoid the vertical overlap we will introduce
structured aliases into the eventual transformations.

Rebecca Stoll noticed that the combined positional constraints were
very similar to the rules of the Japanese combinatorial number-placement 
puzzle Sudoku.

A Sudoku matrix is a 9x9 matrix with each element containing an
integer from 1 to 9 following these rules:
<ol>
<li>each integer can only appear once along a given row
<li>integer can only appear once down any given column
<li>each integer can only appear once within each of 9 3x3 subgrids that make up the overall grid
</ol>
The third rule of no repeats within the 3x3 subgrids is what makes the family
of Sudoku matrices a subset of 9x9 Latin squares (which are created following
the first two rules). Mathematicians have shown that there are ~6.7x10^21 
possible Sudoku matrixes (Felgenhauer & Jarvis 2006). 

Shannon entropy quantifies how random an ensemble of numbers is,
developed by Claude Shannon, one of the founders of information theory
(see the seminal paper: Shannon, C., 1948, "A Mathematical Theory of
Communication", Bell System Technical Journal, 27, 3).  Compared to
generating a 9x9 grid of random numbers (of which there are 9^81
possible realizations), the Sudoku matrix has a higher Shannon entropy
by about 3 orders of magnitude (see Newton & DeSalvo 2010, Proc Roy
Soc A, DOI: 10.1098.rspa.2009.0052).  Greater Shannon entropy means
greater intrinsic randomness (in information theory, Shannon entropy
estimates the number of data bits required to encode a string of
symbols calcluated based on the frequency of occurence of those
symbols in the string).

We therefore created a 2D pinhole mask by finding a standard 9x9
Sudoku puzzle, and then used that number pattern to generate the
pinhole positions following a simple algorithm for how to convert the
number in the square to the location of the test pinhole along the
diagonal within the cell. A couple of simple realizations showed that
this generates a mask that when used with the grating or prism
disperser and a comparison lamp will give us a good, randomized
coverage of the available science detector area to measure the
transformations.

The ultimate transformations is of the form:
<pre>
   (x_CCD,y_CCD) = S(x_mask, y_mask, lambda)
</pre>
Given (x,y) coordinates of an aperture in the mask plane in
millimeters and a wavelength, the "Sudoku Transform", S(x,y,l) will
predict where that wavelength will land the CCD in (x,y) pixel
coordinates.

The application to MOS reduction is clear. While we did OK with more
or less open-ended schemes to perform MOS mapping for reductions, they
lacked robustness.  More sophisticated algorithms (e.g., Carnegie
COSMOS) attempt to improve this by doing forward optical modeling, but
their internal models all assumed the optics are on-axis and therefore
optical distortions are radially symmetric.  MODS uses off-axis
optics, so we have non axially-symmetric distortions, meaning the
models have higher order asymmetric stretches (we cannot quantify the
polynomials only in terms of even powers of x and y, but must include
odd powers and xy cross terms). The xidl package that we were
evaluating at the time was much less sophisticated about its approach
(with the penalty that it performed much worse unmodified), but its
simplicity allowed us greater latitude to break into the code and add
a modeling approach based on the Sudoku mask transforms that, in the
end, made it very robust and quite a bit faster as the guesses for
trace extraction and wavelength calibration were much closer to the
mark on the first iteration, reducing the number of subsequent
iterations and time lost to runaway divergences.

It was while wrestling with these issues in trying to improve our MOS
reductions that we realized we needed different calibration data to
better map the MODS system for MOS mode.  These masks are our new way
to do this.
