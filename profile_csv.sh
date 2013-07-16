#!/bin/bash
# USAGE: Add between aprun options and executable
# For Example: aprun -n 16 -N 1 ./foo arg1 arg2
# Becomes: aprun -n 16 -N 1 ./profile.sh ./foo arg1 arg2
 
# Enable command-line profiler
export COMPUTE_PROFILE=1
 
# Set output to CSV (optional)
export COMPUTE_PROFILE_CSV=1
 
# Give each *node* a separate file
export COMPUTE_PROFILE_LOG=cuda_profile_$(hostname).csv
 
# Stripe each profile file by 1 to share the load on large runs
lfs setstripe -c 1 $COMPUTE_PROFILE_LOG
 
# Execute the provided command.
exec $*
