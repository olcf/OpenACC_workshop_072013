FC=ftn
#FCFLAGS=-O3 -rm -eF -hwp -hpl=vhone.aid -homp_analyze -DOMP_ANALYZE
#FCFLAGS=-O3 -rm -eF -hwp -hpl=vhone.aid -hprofile_generate
#FCFLAGS=-O3 -rm -eF -hwp -hpl=vhone.aid -DPAT_REGION
FCFLAGS=-rm -eF 

SRCS=vh1mods.f90 \
parabola.f90 \
boundary.f90 \
evolve.f90 \
flatten.f90 \
forces.f90 \
ppmlr.f90 \
remap.f90 \
riemann.f90\
states.f90 \
volume.f90 \
zonemod.f90 \
sweepx1.f90 \
sweepx2.f90 \
sweepy.f90 \
dtcon.f90 \
dump.f90 \
images.f90 \
init.f90 \
prin.f90 \
sweepz.f90


OBJS=$(SRCS:.f90=.o)

all: vhone vhone+pat

vhone+pat: vhone
	pat_build -w -g mpi $<

clean:
	rm -f *.o *.lst vhone vhone+pat
	rm -r vhone.aid

vhone: $(OBJS) vhone.f90
	$(FC) $(FCFLAGS)  vhone.f90 $(OBJS) -o $@

%.o: %.f90
	ftn $(FCFLAGS) -c $<

