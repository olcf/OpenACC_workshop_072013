subroutine sweepx1

! This subroutine performs 1D hydro sweeps in the X direction, 
! looping over js by ks rows.  At the end, data is packed into
! a SEND buffer for data transpose in sweepy.
!-----------------------------------------------------------------------

! GLOBALS
use zone
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k, n, m,ierr

!-----------------------------------------------------------------------

sweep  = 'x'
ngeom  = ngeomx
nleft  = nleftx
nright = nrightx
nmin   = 7
nmax   = imax + 6

! Now Loop over each row...
#ifdef OMP_ANALYZE
!dir$ omp_analyze_loop
#endif
#ifdef PAT_REGION
call pat_region_begin(1,'sweepx1',ierr)
#endif
!$omp parallel do private(i,j,k,n)
do k = 1, ks
 do j = 1, js

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   do i = 1,imax
     n = i + 6
     r  (n) = zro(i,j,k)
     p  (n) = zpr(i,j,k)
     u  (n) = zux(i,j,k)
     v  (n) = zuy(i,j,k)
     w  (n) = zuz(i,j,k)
     f  (n) = zfl(i,j,k)

     xa0(n) = zxa(i)
     dx0(n) = zdx(i)
     xa (n) = zxa(i)
     dx (n) = zdx(i)
     p  (n) = max(smallp,p(n))
     e  (n) = p(n)/(r(n)*gamm)+0.5*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! Do 1D hydro update using PPMLR
   call ppmlr 
#ifdef DEBUGX
  print *,'In sweepx1',svel
#endif

   ! Put updated values into send buffer for sweepy, dropping ghost zones
   do i = 1, imax
     n = i + 6
     send1(1,k,j,i) = r(n)
     send1(2,k,j,i) = p(n)
     send1(3,k,j,i) = u(n)
     send1(4,k,j,i) = v(n)
     send1(5,k,j,i) = w(n)
     send1(6,k,j,i) = f(n)
   enddo

 enddo
enddo
#ifdef PAT_REGION
call pat_region_end(1,ierr)
#endif

return
end

