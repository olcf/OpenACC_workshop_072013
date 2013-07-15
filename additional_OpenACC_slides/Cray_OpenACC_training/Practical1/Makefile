# You can override these with command line arguments when calling
# the Makefile, e.g. "make ACC=yes"
VERSION=00
ACC=yes

# You should not need to alter anything below this line
# The Cray wrappers for the different compilers
FC=ftn
CC=cc

FFLAGS=
CFLAGS=

# Options for the Cray compiler
# Write a compiler commentary in separate files with extension .lst
FFLAGS += -hlist=a
CFLAGS += -hlist=a

ifneq ($(ACC),yes)
# Turn off recognition of OpenACC directives
  FFLAGS += -xacc
  CFLAGS += -hnopragma=acc
endif

ifeq ($(VERSION),00)
# There is no Cdynamic version 00 (see README.txt for details)
all: Fstatic Fdynamic Cstatic
else
all: Fstatic Fdynamic Cstatic Cdynamic
endif

Fstatic: first_example_Fstatic_v$(VERSION).f90
	@echo "Executing: $(FC) $(FFLAGS) -o first_example_Fstatic_v$(VERSION).x $^"
	$(FC) $(FFLAGS) -o first_example_Fstatic_v$(VERSION).x $^

Fdynamic: first_example_Fdynamic_v$(VERSION).f90
	@echo "Executing: $(FC) $(FFLAGS) -o first_example_Fdynamic_v$(VERSION).x $^"
	$(FC) $(FFLAGS) -o first_example_Fdynamic_v$(VERSION).x $^

Cstatic: first_example_Cstatic_v$(VERSION).c
	@echo "Executing: $(CC) $(CFLAGS) -o first_example_Cstatic_v$(VERSION).x $^"
	$(CC) $(CFLAGS) -o first_example_Cstatic_v$(VERSION).x $^

Cdynamic: first_example_Cdynamic_v$(VERSION).c
	@echo "Executing: $(CC) $(CFLAGS) -o first_example_Cdynamic_v$(VERSION).x $^"
	$(CC) $(CFLAGS) -o first_example_Cdynamic_v$(VERSION).x $^

clean:
	@$(RM) *.cub *.ptx *.lst *.gpu *.mod *.x *~
