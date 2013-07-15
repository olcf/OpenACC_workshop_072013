#!/bin/bash

# Use this script to write a submission script and then submit the job

# Make sure there are the right amount of arguments given
case $# in
    1)
	echo "Defaulting to NTHREADS=1"
	NTHREADS=1;;
    2)
	NTHREADS=$2;;
    *)
	echo "Usage: $0 EXECUTABLE [NTHREADS]"
	exit
esac

EXECUTABLE=$1
app=$(basename $EXECUTABLE)
app_name=$(echo $app | sed 's/\.x//;s/\+pat//')

# Set the environment correctly (probably not strictly necessary)
MYPE=cray  # the Programming Environment to use
TOOLSDIR=../../Tools
source $TOOLSDIR/XK_setup.bash $MYPE

# A directory on the lustre in which to run the code
case $HOST in
    todi*)
	LUSDIR=/scratch/todi/$USER;;
    *)
	LUSDIR=/lus/scratch/$USER;;
esac
LUSBASE=${LUSDIR}/Cray_OpenACC_training/Practical3
LUS=${LUSBASE}/C/cray_${app}_$(date +%y%m%d_%H%M%S)
# Tell me where it is, for future reference
echo "Running in directory:"
echo $LUS
# Make the directory (and parents) if it isn't already there
mkdir -p $LUS
# Copy the application to the working directory
cp $EXECUTABLE $LUS
# Move the setup script there as well
cp $TOOLSDIR/XK_setup.bash $LUS
# Put the .lst file there (so we know what was compiled)
cp MG/*.lst $LUS

# Write a WLM submission script
#   We will use one node (mppwidth), with one thread per node (mppdepth) and
#   one MPI process per node (mppnppn).
jobscript=submit.wlm

# Write a personalised header first
sed "s/APPLICATION/$app_name/g;
     s/mppdepth=1/mppdepth=$NTHREADS/;
     s/cpus-per-task=1/cpus-per-task=$NTHREADS/" $TOOLSDIR/$WLM_HEADER > $LUS/$jobscript

# Move to the working directory
cd $LUS

# Now write the remainder of the jobscript
cat <<EOF >> $jobscript

# Make sure we are in the correct directory
cd $LUS

# Load appropriate programming environment
#   So appropriate runtime libraries are found
source XK_setup.bash $MYPE

# Uncomment this line for runtime commentary with CCE
#export CRAY_ACC_DEBUG=2
# Uncomment this line for basic profiling
#export COMPUTE_PROFILE=1

# Execute the code, inheriting the job parameters from PBS
export OMP_NUM_THREADS=$NTHREADS
aprun -B $app

EOF

# Submit the job (WLM_SUBMIT defined in XK_setup.bash)
echo "Submitting job: $WLM_SUBMIT $jobscript"
$WLM_SUBMIT $jobscript

#EOF
