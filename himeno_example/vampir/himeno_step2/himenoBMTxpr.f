C*********************************************************************
C
C This benchmark test program is measuring a cpu performance
C of floating point operation by a Poisson equation solver.
CC
C If you have any question, please ask me via email.
C written by Ryutaro HIMENO, November 26, 2001.
C Version 3.0
C ----------------------------------------------
C Ryutaro Himeno, Dr. of Eng.
C Head of Computer Information Division,
C RIKEN (The Institute of Pysical and Chemical Research)
C Email : himeno@postman.riken.go.jp
C -----------------------------------------------------------
C You can adjust the size of this benchmark code to fit your target
C computer. In that case, please chose following sets of
C (mimax,mjmax,mkmax):
C small : 65,33,33
C small : 129,65,65
C midium: 257,129,129
C large : 513,257,257
C ext.large: 1025,513,513
C This program is to measure a computer performance in MFLOPS
C by using a kernel which appears in a linear solver of pressure
C Poisson eq. which appears in an incompressible Navier-Stokes solver.
C A point-Jacobi method is employed in this solver as this method can 
C be easyly vectrized and be parallelized.
C ------------------
C Finite-difference method, curvilinear coodinate system
C Vectorizable and parallelizable on each grid point
C No. of grid points : imax x jmax x kmax including boundaries
C ------------------
C A,B,C:coefficient matrix, wrk1: source term of Poisson equation
C wrk2 : working area, OMEGA : relaxation parameter
C BND:control variable for boundaries and objects ( = 0 or 1)
C P: pressure
C -------------------
      PROGRAM HIMENOBMTXP
C
      IMPLICIT REAL*4(a-h,o-z)
C
      include 'mpif.h'
      include 'param.h'
C
C     ttarget specifys the measuring period in sec
      PARAMETER (ttarget=60.0)
C
      real*8  cpu,cpu0,cpu1,xmflops2,flop
C
      omega=0.8
      mx= mx0-1
      my= my0-1
      mz= mz0-1
C
CC Initializing communicator
      call initcomm
C
CC Initializaing computational index
      call initmax(mx,my,mz,it)
C
CC Initializing matrixes
      call initmt(mz,it)
      if(id .eq. 0) then
         write(*,*) 'Sequential version array size'
         write(*,*) ' mimax=',mx0,' mjmax=',my0,' mkmax=',mz0
         write(*,*) 'Parallel version  array size'
         write(*,*) ' mimax=',mimax,' mjmax=',mjmax,' mkmax=',mkmax
         write(*,*) ' imax=',imax,' jmax=',jmax,' kmax=',kmax
         write(*,*) ' I-decomp= ',ndx,' J-decomp= ',ndy,
     >              ' K-decomp= ',ndz
         write(*,*)
      end if
C
CC Start measuring
C
      nn=3
      if(id .eq. 0) then
         write(*,*) ' Start rehearsal measurement process.'
         write(*,*) ' Measure the performance in 3 times.'
      end if
C
      gosa= 0.0
      cpu= 0.0
      call mpi_barrier(mpi_comm_world,ierr)
      cpu0= mpi_wtime()
C Jacobi iteration
      call jacobi(nn,gosa)
      cpu1= mpi_wtime() - cpu0
C
      call mpi_allreduce(cpu1,
     >                   cpu,
     >                   1,
     >                   mpi_real8,
     >                   mpi_max,
     >                   mpi_comm_world,
     >                   ierr)
C
      flop=real(mx-2)*real(my-2)*real(mz-2)*34.0
      if(cpu .ne. 0.0) xmflops2=flop/cpu*1.0d-6*real(nn)
      if(id .eq. 0) then
         write(*,*) '  MFLOPS:',xmflops2,'  time(s):',cpu,gosa
      end if
      nn= int(ttarget/(cpu/3.0))
C
C     end the test loop
      if(id .eq. 0) then
         write(*,*) 'Now, start the actual measurement process.'
         write(*,*) 'The loop will be excuted in',nn,' times.'
         write(*,*) 'This will take about one minute.'
         write(*,*) 'Wait for a while.'
      end if
C
      gosa= 0.0
      cpu= 0.0
      call mpi_barrier(mpi_comm_world,ierr)
      cpu0= mpi_wtime()
C Jacobi iteration
      call jacobi(nn,gosa)
      cpu1= mpi_wtime() - cpu0
C
      call mpi_reduce(cpu1,
     >                cpu,
     >                1,
     >                mpi_real8,
     >                mpi_max,
     >                0,
     >                mpi_comm_world,
     >                ierr)
C
      if(id .eq. 0) then
         if(cpu .ne. 0.0)  xmflops2=flop*1.0d-6/cpu*real(nn)
C
         write(*,*) ' Loop executed for ',nn,' times'
         write(*,*) ' Gosa :',gosa
         write(*,*) ' MFLOPS:',xmflops2, '  time(s):',cpu
         score=xmflops2/82.84
         write(*,*) ' Score based on Pentium III 600MHz :',score
      end if
      call mpi_finalize(ierr)
C
      stop
      END
C
C
C**************************************************************
      subroutine initmt(mz,it)
C**************************************************************
      IMPLICIT REAL*4(a-h,o-z)
C
      include 'param.h'
C
! Directive inserted by Cray Reveal.  May be incomplete.
!$OMP  parallel do default(none)                      
!$OMP&   private (i,j,k)                             
!$OMP&   shared  (a,b,bnd,c,p,wrk1,wrk2)
      do k=1,mkmax
         do j=1,mjmax
            do i=1,mimax
               a(i,j,k,1)=0.0
               a(i,j,k,2)=0.0
               a(i,j,k,3)=0.0
               a(i,j,k,4)=0.0
               b(i,j,k,1)=0.0
               b(i,j,k,2)=0.0
               b(i,j,k,3)=0.0
               c(i,j,k,1)=0.0
               c(i,j,k,2)=0.0
               c(i,j,k,3)=0.0
               p(i,j,k)=0.0
               wrk1(i,j,k)=0.0   
               wrk2(i,j,k)=0.0   
               bnd(i,j,k)=0.0 
            enddo
         enddo
      enddo
C
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
               a(i,j,k,1)=1.0
               a(i,j,k,2)=1.0
               a(i,j,k,3)=1.0
               a(i,j,k,4)=1.0/6.0
               b(i,j,k,1)=0.0
               b(i,j,k,2)=0.0
               b(i,j,k,3)=0.0
               c(i,j,k,1)=1.0
               c(i,j,k,2)=1.0
               c(i,j,k,3)=1.0
               p(i,j,k)=float((k-1+it)*(k-1+it))
     >                       /float((mz-1)*(mz-1))
               wrk1(i,j,k)=0.0   
               wrk2(i,j,k)=0.0   
               bnd(i,j,k)=1.0
            enddo
         enddo
      enddo
C
      return
      end
C
C*************************************************************
      subroutine jacobi(nn,gosa)
C*************************************************************
      IMPLICIT REAL*4(a-h,o-z)
C
      include 'mpif.h'
      include 'param.h'
C
      DO loop=1,nn
         gosa=0.0
         wgosa=0.0
! Directive inserted by Cray Reveal.  May be incomplete.
!$OMP  parallel do default(none)                                   
!$OMP&   private (i,j,k,s0,ss)                                    
!$OMP&   shared  (a,b,bnd,c,imax,jmax,kmax,omega,p,wrk1,wrk2)    
!$OMP&   reduction (+:wgosa)
         DO K=2,kmax-1
            DO J=2,jmax-1
               DO I=2,imax-1
                  S0=a(I,J,K,1)*p(I+1,J,K)+a(I,J,K,2)*p(I,J+1,K)
     1                 +a(I,J,K,3)*p(I,J,K+1)
     2                 +b(I,J,K,1)*(p(I+1,J+1,K)-p(I+1,J-1,K)
     3                 -p(I-1,J+1,K)+p(I-1,J-1,K))
     4                 +b(I,J,K,2)*(p(I,J+1,K+1)-p(I,J-1,K+1)
     5                 -p(I,J+1,K-1)+p(I,J-1,K-1))
     6                 +b(I,J,K,3)*(p(I+1,J,K+1)-p(I-1,J,K+1)
     7                 -p(I+1,J,K-1)+p(I-1,J,K-1))
     8                 +c(I,J,K,1)*p(I-1,J,K)+c(I,J,K,2)*p(I,J-1,K)
     9                 +c(I,J,K,3)*p(I,J,K-1)+wrk1(I,J,K)
                  SS=(S0*a(I,J,K,4)-p(I,J,K))*bnd(I,J,K)
                  WGOSA=WGOSA+SS*SS
                  wrk2(I,J,K)=p(I,J,K)+OMEGA *SS
               enddo
            enddo
         enddo
C     
!$OMP  parallel do 
!$OMP&   private (i,j,k)                                    
         DO K=2,kmax-1
            DO J=2,jmax-1
               DO I=2,imax-1
                  p(I,J,K)=wrk2(I,J,K)
               enddo
            enddo
         enddo
C
         call sendp(ndx,ndy,ndz)
C
         call mpi_allreduce(wgosa,
     >                      gosa,
     >                      1,
     >                      mpi_real4,
     >                      mpi_sum,
     >                      mpi_comm_world,
     >                      ierr)
C
      enddo
CC End of iteration
      return
      end
c
c
c
      subroutine initcomm
c
      IMPLICIT REAL*4(a-h,o-z)
c
      include 'mpif.h'
      include 'param.h'
c
      logical    ipd(3),ir
      dimension  idm(3)
c
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world,npe,ierr)
      call mpi_comm_rank(mpi_comm_world,id,ierr)
C
      if(ndx*ndy*ndz .ne. npe) then
         if(id .eq. 0) then
            write(*,*) 'Invalid number of PE'
            write(*,*) 'Please check partitioning pattern'
            write(*,*) '                 or number of  PE'
         end if
         call mpi_finalize(ierr)
         stop
      end if
C
      icomm= mpi_comm_world
c
      idm(1)= ndx
      idm(2)= ndy
      idm(3)= ndz
C
      ipd(1)= .false.
      ipd(2)= .false.
      ipd(3)= .false.
      ir= .false.
C
      call mpi_cart_create(icomm,
     >                     ndims,
     >                     idm,
     >                     ipd,
     >                     ir,
     >                     mpi_comm_cart,
     >                     ierr)
      call mpi_cart_get(mpi_comm_cart,
     >                  ndims,
     >                  idm,
     >                  ipd,
     >                  iop,
     >                  ierr)
c
c
      if(ndz .gt. 1) then
         call mpi_cart_shift(mpi_comm_cart,
     >                       2,
     >                       1,
     >                       npz(1),
     >                       npz(2),
     >                       ierr)
      end if
c
      if(ndy .gt. 1) then
         call mpi_cart_shift(mpi_comm_cart,
     >                       1,
     >                       1,
     >                       npy(1),
     >                       npy(2),
     >                       ierr)
      end if
c
      if(ndx .gt. 1) then
         call mpi_cart_shift(mpi_comm_cart,
     >                       0,
     >                       1,
     >                       npx(1),
     >                       npx(2),
     >                       ierr)
      end if
c
      return
      end
c
c
c
      subroutine initmax(mx,my,mz,ks)
c
      IMPLICIT REAL*4(a-h,o-z)
c
      include 'param.h'
      include 'mpif.h'
C
      integer  itmp,ks
      integer  mx1(0:ndx),my1(0:ndy),mz1(0:ndz)
      integer  mx2(0:ndx),my2(0:ndy),mz2(0:ndz)
C
CC    define imax, communication direction
      itmp= mx/ndx
      mx1(0)= 0
      do  i=1,ndx
         if(i .le. mod(mx,ndx)) then
            mx1(i)= mx1(i-1) + itmp + 1
         else
            mx1(i)= mx1(i-1) + itmp
         end if
      end do
      do i=0,ndx-1
         mx2(i)= mx1(i+1) - mx1(i)
         if(i .ne. 0)     mx2(i)= mx2(i) + 1
         if(i .ne. ndx-1) mx2(i)= mx2(i) + 1
      end do
c
      itmp= my/ndy
      my1(0)= 0
      do  i=1,ndy
         if(i .le. mod(my,ndy)) then
            my1(i)= my1(i-1) + itmp + 1
         else
            my1(i)= my1(i-1) + itmp
         end if
      end do
      do i=0,ndy-1
         my2(i)= my1(i+1) - my1(i)
         if(i .ne. 0)      my2(i)= my2(i) + 1
         if(i .ne. ndy-1)  my2(i)= my2(i) + 1
      end do
c
      itmp= mz/ndz
      mz1(0)= 0
      do  i=1,ndz
         if(i .le. mod(mz,ndz)) then
            mz1(i)= mz1(i-1) + itmp + 1
         else
            mz1(i)= mz1(i-1) + itmp
         end if
      end do
      do i=0,ndz-1
         mz2(i)= mz1(i+1) - mz1(i)
         if(i .ne. 0)      mz2(i)= mz2(i) + 1
         if(i .ne. ndz-1)  mz2(i)= mz2(i) + 1
      end do
c
      imax= mx2(iop(1))
      jmax= my2(iop(2))
      kmax= mz2(iop(3))
c
      if(iop(3) .eq. 0) then
         ks= mz1(iop(3))
      else
         ks= mz1(iop(3)) - 1
      end if
c
c     j-k plane  divied by i-direction
      if(ndx .gt. 1) then
         call mpi_type_vector(jmax*kmax,
     >                        1,
     >                        mimax,
     >                        mpi_real4,
     >                        jkvec,
     >                        ierr)
         call mpi_type_commit(jkvec,
     >                        ierr)
      end if
c
c     i-k plane  divied by j-direction
      if(ndy .gt. 1) then
         call mpi_type_vector(kmax,
     >                        imax,
     >                        mimax*mjmax,
     >                        mpi_real4,
     >                        ikvec,
     >                        ierr)
         call mpi_type_commit(ikvec,
     >                        ierr)
      end if
c
c     new vector k-direction
      if(ndz .gt. 1) then
         call mpi_type_vector(jmax,
     >                        imax,
     >                        mimax,
     >                        mpi_real4,
     >                        ijvec,
     >                        ierr)
         call mpi_type_commit(ijvec,
     >                        ierr)
      end if
c
      return
      end
c
c
c
      subroutine sendp(ndx,ndy,ndz)
c
      IMPLICIT REAL*4(a-h,o-z)
C
      if(ndz .gt. 1) then
         call sendp3()
      end if
c
      if(ndy .gt. 1) then
         call sendp2()
      end if
c
      if(ndx .gt. 1) then
         call sendp1()
      end if
c
      return
      end
c
c
c
      subroutine sendp3()
c
      IMPLICIT REAL*4(a-h,o-z)
c
      include 'mpif.h'
      include 'param.h'
c
      dimension ist(mpi_status_size,0:3),ireq(0:3)
      data ireq /4*mpi_request_null/
c
      call mpi_irecv(p(1,1,kmax),
     >               1,
     >               ijvec,
     >               npz(2),
     >               1,
     >               mpi_comm_cart,
     >               ireq(3),
     >               ierr)
c
      call mpi_irecv(p(1,1,1),
     >               1,
     >               ijvec,
     >               npz(1),
     >               2,
     >               mpi_comm_cart,
     >               ireq(2),
     >               ierr)
c
      call mpi_isend(p(1,1,2),
     >               1,
     >               ijvec,
     >               npz(1),
     >               1,
     >               mpi_comm_cart,
     >               ireq(0),
     >               ierr)
c
      call mpi_isend(p(1,1,kmax-1),
     >               1,
     >               ijvec,
     >               npz(2),
     >               2,
     >               mpi_comm_cart,
     >               ireq(1),
     >               ierr)
c
      call mpi_waitall(4,
     >                 ireq,
     >                 ist,
     >                 ierr)
c
      return
      end
c
c
c
      subroutine sendp2()
c
      IMPLICIT REAL*4(a-h,o-z)
c
      include 'mpif.h'
      include 'param.h'
c
      dimension ist(mpi_status_size,0:3),ireq(0:3)
      data ireq /4*mpi_request_null/
c
      call mpi_irecv(p(1,1,1),
     >               1,
     >               ikvec,
     >               npy(1),
     >               2,
     >               mpi_comm_cart,
     >               ireq(3),
     >               ierr)
c
      call mpi_irecv(p(1,jmax,1),
     >               1,
     >               ikvec,
     >               npy(2),
     >               1,
     >               mpi_comm_cart,
     >               ireq(2),
     >               ierr)
c
      call mpi_isend(p(1,2,1),
     >               1,
     >               ikvec,
     >               npy(1),
     >               1,
     >               mpi_comm_cart,
     >               ireq(0),
     >               ierr)
c
      call mpi_isend(p(1,jmax-1,1),
     >               1,
     >               ikvec,
     >               npy(2),
     >               2,
     >               mpi_comm_cart,
     >               ireq(1),
     >               ierr)
c
      call mpi_waitall(4,
     >                 ireq,
     >                 ist,
     >                 ierr)
c
      return
      end
c
c
c
      subroutine sendp1()
c
      IMPLICIT REAL*4(a-h,o-z)
c
      include 'mpif.h'
      include 'param.h'
c
      dimension ist(mpi_status_size,0:3),ireq(0:3)
      data ireq /4*mpi_request_null/
c
      call mpi_irecv(p(1,1,1),
     >               1,
     >               jkvec,
     >               npx(1),
     >               2,
     >               mpi_comm_cart,
     >               ireq(3),
     >               ierr)
c
      call mpi_irecv(p(imax,1,1),
     >               1,
     >               jkvec,
     >               npx(2),
     >               1,
     >               mpi_comm_cart,
     >               ireq(2),
     >               ierr)
c
      call mpi_isend(p(2,1,1),
     >               1,
     >               jkvec,
     >               npx(1),
     >               1,
     >               mpi_comm_cart,
     >               ireq(0),
     >               ierr)
c
      call mpi_isend(p(imax-1,1,1),
     >               1,
     >               jkvec,
     >               npx(2),
     >               2,
     >               mpi_comm_cart,
     >               ireq(1),
     >               ierr)
c
      call mpi_waitall(4,
     >                 ireq,
     >                 ist,
     >                 ierr)
c
      return
      end
