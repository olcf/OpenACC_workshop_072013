C
      parameter(mx0=1025,my0=513,mz0=513)
      parameter(mimax=259,mjmax=131,mkmax=131)
      parameter(ndx=4,ndy=4,ndz=4,ndims=3)
C
      dimension  iop(ndims)
      dimension  npx(2),npy(2),npz(2)
CC Array
      dimension  p(mimax,mjmax,mkmax)
      dimension  a(mimax,mjmax,mkmax,4),
     >           b(mimax,mjmax,mkmax,3),c(mimax,mjmax,mkmax,3)
      dimension  bnd(mimax,mjmax,mkmax)
      dimension  wrk1(mimax,mjmax,mkmax),wrk2(mimax,mjmax,mkmax)
C
C Communication parameter
      common /icart/  iop,mpi_comm_cart
      common /idrec/  npx,npy,npz
      common /multi/  id,npe
      common /nvect/  ijvec,ikvec,jkvec
C Other constants
      common /indx/   imax,jmax,kmax
      common /other/  omega
CC Array
      common /pres/   p
      common /mtrx/   a,b,c
      common /bound/  bnd
      common /work/   wrk1,wrk2
