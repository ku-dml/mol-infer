This directory contains a perl script by which
you can construct hyperplanes for data sets
that are feature-vectorized by our FV generator.

1. Open exp.pl and adjust $RESULTS and $DATADIR
   appropriately.
   - $RESULTS: the directory name where files are output;
   - $DATADIR: the directory name where FV data sets are placed;
   $LASSO_SCRIPT specifies the python script for Lasso learning,
   and you shouldn't change this one.
		    
2. Check @Mono in exp.pl, which is the list of properties
   to be learned.

3. Check %Lambda that are used as lambda in the objective
   of Lasso learning.

4. Run exp.pl by:
   $ perl exp.pl
   Then exp.sh is generated.

5. Run exp.sh as follows.
   $ chmod 755 exp.sh
   $ ./exp.sh

Then the results are output in $RESULTS directory. 
    
