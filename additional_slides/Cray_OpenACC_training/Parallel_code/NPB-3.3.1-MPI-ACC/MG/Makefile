# Makefile for MG code
# Revision: $Revision: 337 $
# Date:     $Date: 2011-11-21 11:06:28 -0600 (Mon, 21 Nov 2011) $

SHELL=/bin/sh
BENCHMARK=mg
BENCHMARKU=MG
ACC=yes
NTHREADS=128

include ../config/make.def

OBJS = mg.o $(COMMON)/print_results.o  \
       $(COMMON)/$(RAND).o $(COMMON)/timers.o

include ../sys/make.common

ifeq ($(ACC),no)
 FCOMPILE += -x acc
endif

FCOMPILE += -hlist=a

$(PROGRAM): config $(OBJS)
	$(FLINK) $(FLINKFLAGS) -o $(PROGRAM) $(OBJS) $(FMPI_LIB)

mg.o:		mg.F  globals.h mpinpb.h npbparams.h
	@ echo "Use OpenACC accelerator directives: $(ACC)"
	$(FCOMPILE) mg.F

clean:
	- rm -f *.o *~ *.ptx *.cub *.lst *.opt *.cg
	- rm -f npbparams.h core



