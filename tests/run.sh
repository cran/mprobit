#!/bin/sh

# Tests after compilations from cmakefile
# Type programname to find command line usage

make -f cmakefile posdef
make -f cmakefile mprobit.exch
make -f cmakefile mprobit.ar
make -f cmakefile mprobit.unstr
make -f cmakefile ordprobit.exch
make -f cmakefile ordprobit.ar
make -f cmakefile ordprobit.unstr

# tests of C code
PATH=$PATH:.
export PATH

echo "test posdef...."
posdef < posdef.in > posdef.cout
echo "done."

echo "test mprobit.exch...."
mprobit.exch binaryex.dat 800 210 1 .3 .2 .6 > binaryex.cout
echo "done."

echo "test mprobit.ar...."
mprobit.ar binaryar.dat 800 210 1 -.3 .4 .7 > binaryar.cout
echo "done."

echo "test mprobit.unstr...."
mprobit.unstr binaryex.dat 800 1 4 .3 .2 .6 .6 .6 .6 .6 .6  > binunstr.cout
echo "done."

echo "test ordprobit.exch...."
ordprobit.exch ordinalex.dat 800 210 1 3 -.5 .3 .4 0.6 > ordprobitex.cout
echo "done."

echo "test ordprobit.ar...."
ordprobit.ar ordinalex.dat 800 210 1 3 -.5 .3 .4 0.7 > ordprobitar.cout
echo "done."

echo "test ordprobit.unstr...."
ordprobit.unstr ordinalex.dat 800 1 3 4 -.5 .3 .4 0.6 .6 .6 .6 .6 .6 > ordprobunstr.cout
echo "done."

#============================================================
# tests of R code

echo "test mvnder.test.R...."
R --no-save --slave -q < mvnder.test.R > mvnder.test.out 
echo "done."

echo "test exchmvn.test.R...."
R --no-save --slave -q < exchmvn.test.R > exchmvn.test.out 
echo "done."

echo "test bintest.R...."
R --no-save --slave -q < bintest.R > bintest.out 
echo "done."

echo "test ordtest.R...."
R --no-save --slave -q < ordtest.R > ordtest.out 
echo "done."

echo "test ordtest2.R...."
R --no-save --slave -q < ordtest2.R > ordtest2.out 
echo "done."
