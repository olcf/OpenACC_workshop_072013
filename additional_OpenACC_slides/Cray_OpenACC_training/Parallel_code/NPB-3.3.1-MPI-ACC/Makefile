SHELL=/bin/sh
CLASS=B
NPROCS=1
SUBTYPE=
VERSION=
ACC=yes # controls whether we use OpenACC accelerator directives
SFILE=config/suite.def
NTHREADS=128 # Number of threads per block for GPU kernels

default: header
	@ sys/print_instructions

MG: mg
mg: header
	mkdir -p bin
	cd MG; $(MAKE) NPROCS=$(NPROCS) CLASS=$(CLASS) ACC=$(ACC) NTHREADS=$(NTHREADS)

# Awk script courtesy cmg@cray.com, modified by Haoqiang Jin
suite:
	@ awk -f sys/suite.awk SMAKE=$(MAKE) $(SFILE) | $(SHELL)

# It would be nice to make clean in each subdirectory (the targets
# are defined) but on a really clean system this will won't work
# because those makefiles need config/make.def
clean:
	- rm -f core 
	- rm -f *~ */core */*~ */*.o */npbparams.h */*.obj */*.exe
	- rm -f */*.cub */*.ptx */*.lst */*.cg */*.opt
	- rm -f MPI_dummy/test MPI_dummy/libmpi.a
	- rm -f sys/setparams sys/makesuite sys/setparams.h
	- rm -f btio.*.out*

veryclean: clean
#	- rm -f config/make.def config/suite.def 
	- rm -f bin/mg.*

header:
	@ sys/print_header

