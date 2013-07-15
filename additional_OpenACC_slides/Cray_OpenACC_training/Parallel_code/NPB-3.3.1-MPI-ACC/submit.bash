#!/bin/bash

if [ $# -ne 4 ]
then
  echo "Usage: $0 BENCHMARK-NAME CLASS NRANKS NPPN"
  echo "    where BENCHMARK-NAME is 'mg'"
  echo "          CLASS is compiled problem size (e.g. 'B')"
  echo "          NRANKS is compiled value"
  echo "          NPPN is number of ranks per node"
  exit
fi

bmk=$1
class=$2
nranks=$3
nppn=$4

nodes=$(((nranks+nppn-1)/nppn))
if [ $nppn -gt $nranks ]
then
    nppn=$nranks
fi

EXECUTABLE=bin/${bmk}.${class}.${nranks}
app=$(basename $EXECUTABLE)
app_name=$app

# Set the environment
MYPE=cray
TOOLSDIR=../../Tools
source $TOOLSDIR/XK_setup.bash $MYPE

# A directory on the lustre in which to run the code
case $HOST in
    todi*)
	LUSDIR=/scratch/todi/$USER;;
    *)
	LUSDIR=/lus/scratch/$USER;;
esac
LUSBASE=${LUSDIR}/Cray_OpenACC_training/Parallel_code
LUS=${LUSBASE}/cray_${app}_${nppn}_$(date +%y%m%d_%H%M%S)
# Tell me where it is, for future reference
echo "Running in directory:"
echo $LUS
# Make the directory (and parents) if it isn't already there
mkdir -p $LUS
# Copy the application to the working directory
cp $EXECUTABLE $LUS
# Move the setup script there as well
cp $TOOLSDIR/XK_setup.bash $LUS
# Move the compiler feedback file there for book-keeping
cp MG/*.lst $LUS

# Write a WLM submission script
#   We will use one node (mppwidth), with one thread per node (mppdepth) and
#   one MPI process per node (mppnppn).
jobscript=submit.wlm

# Write a personalised header first
sed "s/APPLICATION/$app_name/g;
     s/mppwidth=1/mppwidth=$nranks/;
     s/mppnppn=1/mppnppn=$nppn/;
     s/nodes=1/nodes=$nodes/;
     s/ntasks=1/ntasks=$nranks/;
     s/ntasks-per-node=1/ntasks-per-node=$nppn/" $TOOLSDIR/$WLM_HEADER > $LUS/$jobscript

# Move to the working directory
cd $LUS

# Now write the remainder of the jobscript
cat <<EOF >> $jobscript

# Make sure we are in the correct directory
cd $LUS

# Load appropriate programming environment
#   So appropriate runtime libraries are found
. XK_setup.bash $MYPE

# Uncomment this line for runtime commentary
#export CRAY_ACC_DEBUG=2
# Uncomment this line for basic profiling
#export COMPUTE_PROFILE=1

# We'll need the CUDA_PROXY with multiple ranks
nppn=$nppn
if [ \$nppn -gt 1 ]
then
    export CRAY_CUDA_PROXY=1
fi

# Execute the code, inheriting the job parameters from PBS
aprun -B $app

EOF

# Submit the job (WLM_SUBMIT defined in XK_setup.bash)
echo "Submitting job: $WLM_SUBMIT $jobscript"
$WLM_SUBMIT $jobscript

#EOF
