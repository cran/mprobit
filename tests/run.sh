#!/bin/sh

# Tests after compilations from ../src/cmakefile
# Type programname to find command line usage

# In ../src

make -f cmakefile posdef
make -f cmakefile mprobit.exch
make -f cmakefile mprobit.ar
make -f cmakefile mprobit.unstr
make -f cmakefile ordprobit.exch
make -f cmakefile ordprobit.ar
make -f cmakefile ordprobit.unstr

# tests of C code
export PATH=$PATH:../src

posdef < posdef.in > posdef.cout

mprobit.exch binaryex.dat 800 210 1 .3 .2 .6 > binaryex.cout
mprobit.ar binaryar.dat 800 210 1 -.3 .4 .7 > binaryar.cout

mprobit.unstr binaryex.dat 800 1 4 .3 .2 .6 .6 .6 .6 .6 .6  > binunstr.cout

ordprobit.exch ordinalex.dat 800 210 1 3 -.5 .3 .4 0.6 > ordprobitex.cout
ordprobit.ar ordinalex.dat 800 210 1 3 -.5 .3 .4 0.7 > ordprobitar.cout

ordprobit.unstr ordinalex.dat 800 1 3 4 -.5 .3 .4 0.6 .6 .6 .6 .6 .6 > ordprobunstr.cout

#============================================================
# tests of R code

R --no-save --slave -q < mvnder.test.R > mvnder.test.out &
R --no-save --slave -q < exchmvn.test.R > exchmvn.test.out &
R --no-save --slave -q < bintest.R > bintest.out &
R --no-save --slave -q < ordtest.R > ordtest.out &
R --no-save --slave -q < ordtest2.R > ordtest2.out &
