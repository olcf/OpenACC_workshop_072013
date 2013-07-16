#!/bin/bash
# USAGE: Add between aprun options and executable
# For Example: aprun -n 16 -N 1 ./foo arg1 arg2
# Becomes: aprun -n 16 -N 1 ./nvprof.sh ./foo arg1 arg2
 
# Give each *node* a separate file
LOG=profile_$(hostname).nvp
 
# Stripe each profile file by 1 to share the load on large runs
lfs setstripe -c 1 $LOG
 
# Execute the provided command.
exec nvprof .o $LOG $*
