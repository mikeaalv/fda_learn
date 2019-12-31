# fda_learn
This is a list of relevant scripts and functions produced while learning FDA (functional data analysis).

## Contents
smooth_derivative.m: wrap up of the script for smoothing spline by GCV in:

Ramsay, James & Hooker, Giles & Graves, Spencer. (2009). Functional data analysis with R and MATLAB. 10.1007/978-0-387-98185-7.

smooth_derivative_test.m: example and test script for smooth_derivative.m

## Usage

1. Download and add the FDA matlab script (http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/)
2. Try run smooth_derivative_test.m, especially the first example.
3. Run similar script on your own data. Each column of Ymat is one curve, step size in xvec need to be close to 1, and choose nDer according to your need. You might also need to try different range of loglambda_vec.

For resources on FDA, refer to:

   http://www.psych.mcgill.ca/misc/fda/
   
   https://cran.r-project.org/web/packages/fda/index.html
   
   http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/
