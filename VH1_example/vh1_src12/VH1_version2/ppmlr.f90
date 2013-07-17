subroutine ppmlr

! Using the 1D arrays of rho, u, and P, perform a 1D lagrangian hydrodynamics evolution:
!   - set boundary conditions by filling ghost zones
!   - obtain parabolic interpolations of rho, u, P
!   - compute input states from these interpolations for Riemann problem
!   - call the Riemann solver to find the time averages umid, pmid
!   - evolve the finite difference equations to get updated values of rho, u, and E
! Then perform a conservative remap back to the Eulerian grid
!-----------------------------------------------------------------------

! GLOBALS
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: n
REAL :: xwag
REAL, DIMENSION(maxsweep) :: dr, du, dp, r6, u6, p6, rl, ul, pl
REAL, DIMENSION(maxsweep) :: rrgh, urgh, prgh, rlft, ulft, plft, umid, pmid
REAL, DIMENSION(maxsweep) :: xaf, dxf

!-----------------------------------------------------------------------
if (ngeom<3) radius = 1.0

! Apply boundary conditions by filling ghost zones
call boundary
#ifdef DEBUG
print *,'after boundary',sum(r(nmin:nmax)),sum(u(nmin:nmax)),sum(v(nmin:nmax)),sum(w(nmin:nmax)),sum(p(nmin:nmax)),sum(e(nmin:nmax)),sum(f(nmin:nmax))
#endif

! Calculate flattening coefficients for smoothing near shocks
call flatten
#ifdef DEBUG
print *,'after flatten',sum(flat(nmin:nmax))
#endif

! Compute parabolic coefficients 
call paraset( nmin-4, nmax+5, para, dx, xa )
#ifdef DEBUG
print *,'after paraset',sum(para(nmin:nmax,:)),sum(dx(nmin:nmax)),sum(xa(nmin:nmax))
#endif

! Interpolate parabolae for fluid variables 
call parabola(nmin-4, nmax+4, para, p, dp, p6, pl, flat)
#ifdef DEBUG
print *,'after parabola',sum(dp(nmin:nmax)),sum(p6(nmin:nmax)),sum(pl(nmin:nmax))
#endif
call parabola(nmin-4, nmax+4, para, r, dr, r6, rl, flat)
#ifdef DEBUG
print *,'after parabola',sum(dr(nmin:nmax)),sum(r6(nmin:nmax)),sum(rl(nmin:nmax))
#endif
call parabola(nmin-4, nmax+4, para, u, du, u6, ul, flat)
#ifdef DEBUG
print *,'after parabola',sum(du(nmin:nmax)),sum(u6(nmin:nmax)),sum(ul(nmin:nmax))
#endif

! Integrate parabolae to get input states for Riemann problem
call states( pl, ul, rl, p6, u6, r6, dp, du, dr, plft, ulft, rlft, prgh, urgh, rrgh )
#ifdef DEBUG
print *,'after state',sum(prgh(nmin:nmax)),sum(urgh(nmin:nmax)),sum(rrgh(nmin:nmax))
#endif
 
! Call the Riemann solver to obtain the zone face averages, umid and pmid
call riemann( nmin-3, nmax+4, gam, prgh, urgh, rrgh, plft, ulft, rlft, pmid, umid )
#ifdef DEBUG
print *,'after riemann',sum(umid(nmin:nmax)),sum(pmid(nmin:nmax))
#endif

! do lagrangian update using umid and pmid
call evolve( umid, pmid )
#ifdef DEBUG
print *,'after evolve',sum(r(nmin:nmax)),sum(u(nmin:nmax)),sum(q(nmin:nmax))
#endif

!#########################################################################
! EXTRA DISSIPATION TO REDUCE CARBUNCLE NOISE
xwag = sum(flat)
if (xwag*xwig /= 0.0) then ! wiggle grid, remap, then remap back to Eulerian grid

 ! save Eulerian coordinates for second remap
 xaf = xa0
 dxf = dx0

 ! wiggle grid where there is a shock, leaving edges unchanged
 do n = nmin+1, nmax
  if (max(flat(n-1),flat(n)) > 0.0) xa0(n) = xa0(n) + xwig*dx0(n)
  dx0(n-1) = xa0(n) - xa0(n-1)
 enddo
 dx0(nmax) = xa0(nmax+1) - xa0(nmax)

 call remap

 ! put wiggled grid into xa(n), but keep ghost cells untouched
 do n = nmin, nmax+1
  xa(n)   = xa0(n)
  dx(n-1) = xa(n) - xa(n-1)
 enddo
 call volume(nmin, nmax, ngeom, radius, xa, dx, dvol)

 ! put Eulerian grid back into xa0
 xa0 = xaf
 dx0 = dxf

endif
!#########################################################################

! remap onto original Eulerian grid
call remap

return
end

