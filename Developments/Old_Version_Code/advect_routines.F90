! #################################################################### !
! ## this file contains various routines used from advect           ## !
! #################################################################### !
! --------------------------------------------------------------------------------
!> Subroutine to update velocities according to momenta change   
!!
!! @authors  Kai Wuennemann, MfN; Dirk Elbeshausen, MfN
!!
!!     Description                                     Programmer    Date
!!     ------------------------------------------------------------------
!!     Original version (5.7).............................KAI  2005/03/22
!!     Conversion into Fortran90..........................DE   2007/02/24
!!     Modified to include HIS algorithm..................GSC  2016/09/09
!!
!<--------------------------------------------------------------------------------
subroutine update_velocities(vmom_dim,vmomp)
  use mod_isale
  implicit none
  integer, intent(in) :: vmom_dim !> Dimension of 2nd vmom array co-ordinate
  real*8, intent(in) :: vmomp(X:Y,1:vmom_dim,igs:ige,jgs:jge) !< Momenta change
  integer i,j
  real*8 :: qmomv(X:Y)       ! Nodal momenta changes for SALE algorithm
  real*8 :: qmomv_b(X:Y,1:4) ! Nodal momenta changes for HIS algorithm

  if (vmom_dim .eq. 4) then
     do j=js,je
        do i=is,ie
           ! The HIS algorithm advects the change in each nodal momentum component separately
           ! Here the advected momentum change of the relevant node is used to update the
           ! nodal velocity. The factor of 1/4 is because each nodal momentum component is
           ! advected four times (from each of the neigbouring cells) and then averaged. . .
           qmomv_b(X:Y,1:4) = 0.25D0*vmomp(X:Y,1:4,i,j)
           V(X:Y,i+1,j)   = V(X:Y,i+1,j)   + qmomv_b(X:Y,1)*rmv(i+1,j)
           V(X:Y,i+1,j+1) = V(X:Y,i+1,j+1) + qmomv_b(X:Y,2)*rmv(i+1,j+1)
           V(X:Y,i,j+1)   = V(X:Y,i,j+1)   + qmomv_b(X:Y,3)*rmv(i,j+1)
           V(X:Y,i,j)     = V(X:Y,i,j)     + qmomv_b(X:Y,4)*rmv(i,j)           
        enddo
     enddo
  else
     do j=js,je
        do i=is,ie
           ! The SALE algorithm advects the change in cell-centered momenta. Here the
           ! momenta changes are partitioned equally among the cell nodes. . .
           qmomv(X:Y) = 0.25D0*vmomp(X:Y,1,i,j)
           V(X:Y,i+1,j)   = V(X:Y,i+1,j)   + qmomv(X:Y)*rmv(i+1,j)
           V(X:Y,i+1,j+1) = V(X:Y,i+1,j+1) + qmomv(X:Y)*rmv(i+1,j+1)
           V(X:Y,i,j+1)   = V(X:Y,i,j+1)   + qmomv(X:Y)*rmv(i,j+1)
           V(X:Y,i,j)     = V(X:Y,i,j)     + qmomv(X:Y)*rmv(i,j)           
        enddo
     enddo
  end if

end subroutine update_velocities
! --------------------------------------------------------------------------------

! --------------------------------------------------------------------------------
!> Subroutine to cap spuriously high vertex velocities   
!!
!! @authors  Kai Wuennemann, MfN; Dirk Elbeshausen, MfN; Gareth Collins, ICL
!!
!!     Description                                     Programmer    Date
!!     ------------------------------------------------------------------
!!     Original version (5.7).............................KAI  2005/03/22
!!     Conversion into Fortran90..........................DE   2007/02/24
!!     Cropping now based on velocity relative to
!!     maximum velocity of full-nodes ....................GSC  2008/05/20
!!
!<--------------------------------------------------------------------------------
subroutine crop_velocity()
  use mod_isale
  implicit none
  real*8 :: sfac       !Absolute maximum vertex velocity allowed
  real*8 :: vmod       !magnitude of velocity vector at node
  real*8 :: vmax       !Maximum node velocity
  real*8 :: mommax     !Maximum momentum (magnitude) of nodes surrounding vertex
  real*8 :: mom        !momentum (magnitude) of node adjacent to local node
  integer :: i,j       !Global node number
  integer :: iloc,jloc !Local node number
  real*8, parameter :: VEL_MAXFAC = 100.D0
  real*8, parameter :: CMV_CUTOFF = 0.D0 !1.D-3
  real*8, parameter :: VEL_MAX = 10.D0  ! times impact velocity

  if (VEL_CUTOFF .le. 0.D0) return

  ! Limit all velocities to VEL_CUTOFF or VEL_MIN
  ! (depending on which value is larger).
  sfac=max(VEL_CUTOFF,VEL_MIN)
  do j=js,jep
     do i=is,iep
        vmod=sqrt(V(X,i,j)**2+V(Y,i,j)**2)
        if (vmod.gt.sfac) then
           V(X:Y,i,j)=V(X:Y,i,j)*sfac/vmod
        endif
     enddo
  enddo

end subroutine crop_velocity
! --------------------------------------------------------------------------------

