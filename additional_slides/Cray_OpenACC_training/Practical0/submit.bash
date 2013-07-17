#!/bin/bash

# Use this script to write a submission script and then submit the job

# Set up the programming environment
MYPE=cray
TOOLSDIR=$(dirname $(pwd -P))/Tools
source $TOOLSDIR/XK_setup.bash $MYPE

# A directory on the lustre in which to run the code
case $HOST in
    todi*)
	LUSDIR=/scratch/todi/$USER;;
    *)
	LUSDIR=/lus/scratch/$USER;;
esac
LUSBASE=${LUSDIR}/Cray_OpenACC_training/Practical0
LUS=${LUSBASE}/$(date +%y%m%d_%H%M%S)
# Tell me where it is, for future reference
echo "Running in directory:"
echo $LUS
# Make the directory (and parents) if it isn't already there
mkdir -p $LUS
# Move to the working directory
cd $LUS

# Write a WLM submission script
#   We will use one node (mppwidth), with one thread per node (mppdepth) and
#   one MPI process per node (mppnppn).
jobscript=submit.wlm

# Write a personalised header first
app_name="practical0"
sed "s/APPLICATION/$app_name/g" $TOOLSDIR/$WLM_HEADER > $jobscript

# Now write the remainder of the jobscript
cat <<EOF >> $jobscript

# Make sure we are in the correct directory
cd $LUS

# If you are using OpenMP, set the number of threads here
#export OMP_NUM_THREADS=<whatever>

# Now run the job with aprun
#   aprun -B uses the same parameters (-N -n -d) as specified in WLM

# As an example, we'll run two utilities to provide information
# about the GPU(s): 

# nvidia-smi from Nvidia
# Don't worry about the specific options used here
echo -e "\nRunning nvidia-smi:"
CRAY_ROOTFS=INITRAMFS aprun -B -b /bin/sh -c "nvidia-smi -q"

# pgaccelinfo from PGI
# To access this, we need to load PrgEnv-pgi:
source \${MODULESHOME}/init/sh
module swap PrgEnv-cray PrgEnv-pgi

# Run the command:
echo -e "\nRunning pgaccelinfo:"
aprun -B pgaccelinfo

EOF

# Submit the job (WLM_SUBMIT defined in XK_setup.bash)
echo "Submitting job: $WLM_SUBMIT $jobscript"
$WLM_SUBMIT $jobscript

#EOF

