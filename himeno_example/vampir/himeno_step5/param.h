C
      parameter(mx0=257,my0=129,mz0=129)
      parameter(mimax=131,mjmax=67,mkmax=129)
      parameter(ndx=2,ndy=2,ndz=1,ndims=3)
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
