! #################################################################### !
! ## this file contains various routines used from advect           ## !
! #################################################################### !
! --------------------------------------------------------------------------------
!> Subroutine to update velocities according to momenta change   
!!
!! @authors  Kai Wuennemann, MfN; Dirk Elbeshausen, MfN; Xianzheng Li, ICL
!! @GitHub: x1ng4me(XL)
!!     Description                                     Programmer    Date
!!     ------------------------------------------------------------------
!!     Original version (5.7).............................KAI  2005/03/22
!!     Conversion into Fortran90..........................DE   2007/02/24
!!     Modified to include HIS algorithm..................GSC  2016/09/09
!!     Remove race condition for parallelization..........XL   2019/08/28
!<--------------------------------------------------------------------------------
subroutine update_velocities(vmom_dim,vmomp)
  use mod_isale
  implicit none
  integer, intent(in) :: vmom_dim !> Dimension of 2nd vmom array co-ordinate
  real*8, intent(in) :: vmomp(X:Y,1:vmom_dim,igs:ige,jgs:jge) !< Momenta change
  integer i,j
  real*8 :: qmomv(X:Y)       ! Nodal momenta changes for SALE algorithm
  real*8 :: qmomv_b(X:Y,1:4) ! Nodal momenta changes for HIS algorithm
  if (vmom_dim.eq.4) then
     ! Refactored the velocities update loop, change loops over nodes to cells
     ! Make sure the code after adding OpenMP will not have conflicts or overwritting
     ! Consider boundary condition and use 'if statement' separate the code.
     !$omp parallel private(i,j)
     !$omp do
     do j=js,jep
        do i=is,iep    
           ! The HIS algorithm advects the change in each nodal momentum component separately
           ! Here the advected momentum change of the relevant node is used to update the
           ! nodal velocity. The factor of 1/4 is because each nodal momentum component is
           ! advected four times (from each of the neigbouring cells) and then averaged. . .
           if (j.eq.1) then
              if (i.eq.1) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,4,i,j)*rmv(i,j)
              else if (i.eq.iep) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i-1,j)*rmv(i,j)
              else
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,1,i-1,j) + vmomp(X:Y,4,i,j))
              end if
           else if (j.eq.jep) then
              if (i.eq.1) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,3,i,j-1)*rmv(i,j)
              else if (i.eq.iep) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,2,i-1,j-1)*rmv(i,j)
              else
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,2,i-1,j-1) + vmomp(X:Y,3,i,j-1))
              end if
           else if (i.eq.1) then
              if (j.gt.1 .AND. j.lt.jep) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,4,i,j) + vmomp(X:Y,3,i,j-1))
              end if
           else if (i.eq.iep) then
              if (j.gt.1 .AND. j.lt.jep) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,1,i-1,j) + vmomp(X:Y,2,i-1,j-1))
              end if
           else
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,4,i,j)*rmv(i,j)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,3,i,j-1)*rmv(i,j)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i-1,j)*rmv(i,j)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,2,i-1,j-1)*rmv(i,j)
           end if
        enddo
     enddo
     !$omp end do
     !$omp end parallel
  else
     do j=js,jep
        do i=is,iep
           ! The SALE algorithm advects the change in cell-centered momenta. Here the
           ! momenta changes are partitioned equally among the cell nodes. . .
           if (j.eq.1) then
              if (i.eq.1) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i,j)*rmv(i,j)
              else if (i.eq.iep) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i-1,j)*rmv(i,j)
              else
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,1,i-1,j) + vmomp(X:Y,1,i,j))
              end if
           else if (j.eq.jep) then
              if (i.eq.1) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i,j-1)*rmv(i,j)
              else if (i.eq.iep) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i-1,j-1)*rmv(i,j)
              else
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,1,i-1,j-1) + vmomp(X:Y,1,i,j-1))
              end if
           else if (i.eq.1) then
              if (j.gt.1 .AND. j.lt.jep) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,1,i,j) + vmomp(X:Y,1,i,j-1))
              end if
           else if (i.eq.iep) then
              if (j.gt.1 .AND. j.lt.jep) then
                 V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,1,i-1,j) + vmomp(X:Y,1,i-1,j-1))
              end if
           else
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i,j)*rmv(i,j)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i,j-1)*rmv(i,j)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i-1,j)*rmv(i,j)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i-1,j-1)*rmv(i,j)
           end if
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

