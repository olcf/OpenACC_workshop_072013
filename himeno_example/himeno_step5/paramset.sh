#!/bin/sh
#
ndx=$2
ndy=$3
ndz=$4
#
case "$1" in
    ssmall | XS )
       mx0=65
       my0=33
       mz0=33 ;;
    small | S )
       mx0=129
       my0=65
       mz0=65 ;;
    midium | M )
       mx0=257
       my0=129
       mz0=129 ;;
    large| L )
       mx0=513
       my0=257
       mz0=257 ;;
    elarge| XL )
       mx0=1025
       my0=513
       mz0=513 ;;
    * )
       echo ' Invalid argument'
       echo ' Usage:: % program <Grid size> <ID> <JD> <KD>'
       echo '         Grid size= XS (64x32x32)'
       echo '                    S  (128x64x64)'
       echo '                    M  (256x128x128)'
       echo '                    L  (512x256x256)'
       echo '                    XL (1024x512x512)'
       echo ' '
       echo ' <ID> <JD> <KD> is partition size'
       echo '        <ID> is the number of partition for I-dimensional'
       echo '        <JD> is the number of partition for J-dimensional'
       echo '        <KD> is the number of partition for K-dimensional'
       echo ' '
       echo ' The number of PE is fixed partition size'
       echo '        Number of PE= <ID> x <JD> x <KD>'
       exit ;;
esac
#
if [ -f param.h ]
then
  rm param.h
fi
#
echo 'C' >> param.h
echo '      parameter(mx0='$mx0',my0='$my0',mz0='$mz0')' >> param.h
#
if [ $ndx -eq 1 ]
then
    itmp=$mx0
elif [ $ndx -ne 1 ]
then
    iib=`expr $mx0 / $ndx`
    itmp=`expr $iib + 3`
fi
#
if [ $ndy -eq 1 ]
then
    jtmp=$my0
elif [ $ndy -ne 1 ]
then
    iib=`expr $my0 / $ndy`
    jtmp=`expr $iib + 3`
fi
#
if [ $ndz -eq 1 ]
then
    ktmp=$mz0
elif [ $ndz -ne 1 ]
then
    iib=`expr $mz0 / $ndz`
    ktmp=`expr $iib + 3`
fi

echo '      parameter(mimax='$itmp',mjmax='$jtmp',mkmax='$ktmp')' >> param.h
echo '      parameter(ndx='$ndx',ndy='$ndy',ndz='$ndz',ndims=3)' >> param.h
echo 'C' >> param.h
unset mx0
unset my0
unset mz0
unset itmp
unset jtmp
unset ktmp
unset nxd
unset nyd
unset nzd
unset iib
#
echo '      dimension  iop(ndims)' >> param.h
echo '      dimension  npx(2),npy(2),npz(2)' >> param.h
echo 'CC Array' >> param.h
echo '      dimension  p(mimax,mjmax,mkmax)' >> param.h
echo '      dimension  a(mimax,mjmax,mkmax,4),' >> param.h
echo '     >           b(mimax,mjmax,mkmax,3),c(mimax,mjmax,mkmax,3)' >> param.h
echo '      dimension  bnd(mimax,mjmax,mkmax)' >> param.h
echo '      dimension  wrk1(mimax,mjmax,mkmax),wrk2(mimax,mjmax,mkmax)' >> param.h
echo 'C' >> param.h
echo 'C Communication parameter' >> param.h
echo '      common /icart/  iop,mpi_comm_cart' >> param.h
echo '      common /idrec/  npx,npy,npz' >> param.h
echo '      common /multi/  id,npe' >> param.h
echo '      common /nvect/  ijvec,ikvec,jkvec' >> param.h
echo 'C Other constants' >> param.h
echo '      common /indx/   imax,jmax,kmax' >> param.h
echo '      common /other/  omega' >> param.h
echo 'CC Array' >> param.h
echo '      common /pres/   p' >> param.h
echo '      common /mtrx/   a,b,c' >> param.h
echo '      common /bound/  bnd' >> param.h
echo '      common /work/   wrk1,wrk2' >> param.h
