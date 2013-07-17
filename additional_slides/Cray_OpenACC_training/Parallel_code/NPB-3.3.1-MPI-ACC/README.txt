This directory contains an OpenACC port of the NAS Parallel Benchmarks MG
code, based on the MPI version of this code.

The port is not optimised, but is included here to test the CUDA_PROXY
facility that allows multiple ranks on a node to share a single gpu.

To set up the compilation environment initially:

source ../../Tools/XK_setup.bash

To build the code to run on N MPI ranks:

make veryclean
make MG NPROCS=N

(N must be chosen at compile time). The default executable created is

./bin/mg.B.N

To build the code without OpenACC support (to run on the CPUs):

make MG NPROCS=2 ACC=no

The same executable name is created.

To submit the code to run with N ranks, with P ranks per node:

bash submit.bash mg B N P

If P>1, the CRAY_CUDA_PROXY environment variable is set in the jobscript.

EOF
