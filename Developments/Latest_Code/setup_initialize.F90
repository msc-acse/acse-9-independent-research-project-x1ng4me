! --------------------------------------------------------------------------------
!> Subroutine to initialize most array or matrix fields.
!!
!! @author Dirk Elbeshausen, MfN; Xianzheng Li, ICL
!! @GitHub: x1ng4me(XL)
!!     Description                                     Programmer    Date
!!     ------------------------------------------------------------------
!!     Original version (7.0).............................DE   2007/02/25
!!     Move arrays into one loop and add parallelization..XL   2019/08/28
!<--------------------------------------------------------------------------------
SUBROUTINE SETUP_INITIALIZE
  use mod_isale
  use ptool_interface
  use mod_store
  use mod_gradient
  IMPLICIT NONE
  integer :: i,j

  call ptool_text(io_oin,"","",2)
  call ptool_text(io_oin,"Allocating memory for fields","",1)
  
  ! Fundamental variables (advected; needed for restart)
  ! store coordinates globally to minimize communication
  ALLOCATE(C(1:2,1:nxp,1:nyp))        ! COORDS IN X,Y (store them global)
  ALLOCATE(V(1:2,IGS:IGE,JGS:JGE))        ! VELOCITY COMPONENTS
  ALLOCATE(MC(IGS:IGE,JGS:JGE))           ! MASS IN CELL
  ALLOCATE(CMC(NMAT,IGS:IGE,JGS:JGE))     ! MAT. CONCENTRATION IN CELL
  ALLOCATE(RO(NMAT,IGS:IGE,JGS:JGE))      ! DENSITY
  ALLOCATE(SIE(NMAT,IGS:IGE,JGS:JGE))     ! SPECIFIC INTERNAL ENERGY
  ALLOCATE(ALPHA(NMAT,IGS:IGE,JGS:JGE))   ! DISTENSION 1./(1.-POROSITY)
  ALLOCATE(VOLSTRAIN(IGS:IGE,JGS:JGE))    ! VOLUME STRAIN
  if (dummyfields) then
     ALLOCATE(DUMMYF(1:ndummy,IGS:IGE,JGS:JGE))       ! DUMMY FIELD #1
  end if
  ALLOCATE(INITCOORD(X:Y,IGS:IGE,JGS:JGE))  ! INITAL COORD FIELD
  ALLOCATE(DAMAGE(IGS:IGE,JGS:JGE))       ! SCALAR DAMAGE (needed in hydro)

  if (field_stored("Vis",varname_save,nvar_save)) then
     allocate(visco_eff(igs:ige,jgs:jge))
  endif

  ! Additional fundamental variables for stress
  if (.not. hydro) then
     ALLOCATE(STRESDEV(1:numstrcomp,IGS:IGE,JGS:JGE)) ! COMP. OF STRESS DEV. (pixx,piyy,pixy,pith)
     ALLOCATE(EPSTRAIN(IGS:IGE,JGS:JGE))  ! EQUIVALENT PLASTIC STRAIN
     ALLOCATE(VELO(IGS:IGE,JGS:JGE))      ! VIBRATIONAL VELOCITY (ACOUSTIC FLUIDIZATION)
     ALLOCATE(SCATTER(IGS:IGE,JGS:JGE))   ! VIBRATIONAL SCATTER ENERGY
     ALLOCATE(PLW(IGS:IGE,JGS:JGE))       ! PLASTIC WORK DONE
  endif

  if (field_stored("exx",varname_save,nvar_save)) then
     allocate(etam_tens(1:numstrcomp,IGS:IGE,JGS:JGE)) ! viscoelastic comp. of stress tensor
  end if
  if (field_stored("etm",varname_save,nvar_save)) then
     allocate(etam_field(IGS:IGE,JGS:JGE))
  end if
  ! Derived variables (from fundamental varibles)
  ALLOCATE(R(NXP,NYP))                    ! RADIAL COORDS (store them global)
  ALLOCATE(P(IGS:IGE,JGS:JGE))            ! PRESSURE
  ALLOCATE(MV(IGS:IGE,JGS:JGE))           ! MASS AT VERTICES
  ALLOCATE(RMV(IGS:IGE,JGS:JGE))          ! = 1/MV
  ALLOCATE(CMV(NMAT,IGS:IGE,JGS:JGE))     ! MAT. CONCENTRATION AT VERTICES
  ALLOCATE(VOIDCONC(IGS:IGE,JGS:JGE))     ! CONCENTRATION OF VOID
  ALLOCATE(TMPR(IGS:IGE,JGS:JGE))         ! TEMPERATURE
  ALLOCATE(TMELT(IGS:IGE,JGS:JGE))        ! MELTING TEMPERATURE
  ALLOCATE(CSOUND(IGS:IGE,JGS:JGE))       ! SPEED OF SOUND
  ALLOCATE(ENTROPY(NMAT,IGS:IGE,JGS:JGE)) ! SPECIFIC ENTROPY (J/K/kg)
  ALLOCATE(Q(IGS:IGE,JGS:JGE))            ! ARTIFICIAL VISCOSITY
  ALLOCATE(VOL(IGS:IGE,JGS:JGE))          ! CELL VOLUME
  ALLOCATE(PMAX(IGS:IGE,JGS:JGE))         ! MAXIMUM PRESSURE
  ALLOCATE(A_S(IGS:IGE,JGS:JGE))          ! BULK MODULUS
  ALLOCATE(MU(IGS:IGE,JGS:JGE))           ! SHEAR MODULUS
  ALLOCATE(RRSUM(IGS:IGE,JGS:JGE))        ! 1./(SUM OVER LOCAL VERTICES OF R)

  if (.not.hydro) then
     ALLOCATE(YIELD(IGS:IGE,JGS:JGE))        ! YIELD STRENGTH
     ALLOCATE(YAC(IGS:IGE,JGS:JGE))          ! ACOUSTIC FLUIDIZATION INTENSITY
     ALLOCATE(ERATE(IGS:IGE,JGS:JGE))        ! STRAIN RATE
     ALLOCATE(PVIBR(IGS:IGE,JGS:JGE))        ! VIBRATIONAL PRESSURE
     ALLOCATE(FAILURE(IGS:IGE,JGS:JGE))      ! FAILURE MODE
     failure = FAIL_NONE
     if (field_stored("FaP",varname_save,nvar_save)) then
        ALLOCATE(failure_first(IGS:IGE,JGS:JGE))
     endif

     if (      field_stored("Exx",varname_save,nvar_save) &
          .or. field_stored("Eyy",varname_save,nvar_save) &
          .or. field_stored("Exy",varname_save,nvar_save) &
          .or. field_stored("Eth",varname_save,nvar_save)) then
        ALLOCATE(erate_tens(1:numstrcomp,igs:ige,jgs:jge))
     endif
  end if

  ! Temporary variables
  ALLOCATE(VL(1:2,IGS:IGE,JGS:JGE))       ! VELOCITY COMP. (TEMP. STORED)
  ALLOCATE(VREL(1:2,IGS:IGE,JGS:JGE))     ! VELOCITY COMP. RELATIVE TO GRID
  ALLOCATE(ROL(NMAT,IGS:IGE,JGS:JGE))     ! DENSITY (TEMP. STORED)
  ALLOCATE(PL(IGS:IGE,JGS:JGE))           ! PRESSURE (TEMP. STORED)

  ! Allocate memory for gravity vector field
  if (gradient_type==GRADTYPE_CENTRAL .or. gradient_type==GRADTYPE_SELF) then
     ALLOCATE(GRAVITY(1:2,IGS:IGE,JGS:JGE))
  end if
  ! To add OpenMP to the initialization, separate the allocation and initialization part
  ! Then add OpenMP to initialize all temporary arrays in iSALE-2D.
  !$omp parallel private(i,j)
  !$omp do 
  do j = 1,nyp
     do i = 1,nxp
        C(:,i,j) = 0.D0
        V(:,i,j) = 0.D0
        mc(i,j) = 0.D0
        cmc(:,i,j) = 0.D0
        ro(:,i,j) = 0.D0
        sie(:,i,j) = 0.D0
        alpha(:,i,j) = 0.D0
        volstrain(i,j) = 0.D0
        if (dummyfields) then
           dummyf(:,i,j) = 0.D0
        end if
        initcoord(:,i,j) = 0.D0
        damage(i,j)=0.D0 
        if (field_stored("Vis",varname_save,nvar_save)) then
           visco_eff(i,j) = 0.D0
        end if
        if (.not. hydro) then
           stresdev(:,i,j)=0.D0
           epstrain(i,j)=0.D0
           velo(i,j) = 0.D0
           scatter(i,j)=0.D0 
           plw(i,j)=0.D0
        end if
        if (field_stored("exx",varname_save,nvar_save)) then
           etam_tens(:,i,j) = 0.D0
        end if
        if (field_stored("etm",varname_save,nvar_save)) then
           etam_field(i,j) = 0.D0
        end if
        r(i,j)=1.D0
        p(i,j)=0.D0
        mv(i,j)=0.D0
        rmv(i,j) = 0.D0
        cmv(:,i,j)=0.D0
        voidconc(i,j)=1.D0
        tmpr(i,j)=0.D0
        tmelt(i,j)=0.D0
        csound(i,j)=0.D0
        entropy(:,i,j)=0.D0
        q(i,j)=0.D0
        vol(i,j)=0.D0
        pmax(i,j)=0.D0
        A_s(i,j)=0.D0
        mu(i,j)=0.D0
        rrsum(i,j)=0.D0
        if (.not.hydro) then
           yield(i,j)=0.D0
           YAc(i,j)=0.D0
           erate(i,j)=0.D0
           pvibr(i,j)=0.D0
           if (field_stored("FaP",varname_save,nvar_save)) then
              failure_first(i,j)=FAIL_NONE
           end if
           if (      field_stored("Exx",varname_save,nvar_save) &
                .or. field_stored("Eyy",varname_save,nvar_save) &
                .or. field_stored("Exy",varname_save,nvar_save) &
                .or. field_stored("Eth",varname_save,nvar_save)) then
              erate_tens(:,i,j) = 0.D0
           end if
        end if
        VL(:,i,j)=0.D0
        VREL(:,i,j)=0.D0
        ROL(:,i,j)=0.D0
        PL(i,j)=0.D0
        if (gradient_type==GRADTYPE_CENTRAL .or. gradient_type==GRADTYPE_SELF) then
           GRAVITY(:,i,j) = 0.d0
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
END SUBROUTINE SETUP_INITIALIZE

subroutine destroymesh()
  use mod_isale
  use ptool_interface
  IMPLICIT NONE

  DEALLOCATE(C,V)
  DEALLOCATE(MC)           ! MASS IN CELL
  DEALLOCATE(CMC)          ! MAT. CONCENTRATION IN CELL
  DEALLOCATE(RO)           ! DENSITY
  DEALLOCATE(SIE)          ! SPECIFIC INTERNAL ENERGY
  DEALLOCATE(ALPHA)        ! DISTENSION 1./(1.-POROSITY)
  DEALLOCATE(VOLSTRAIN)    ! VOLUME STRAIN
  if (dummyfields) DEALLOCATE(DUMMYF)
  DEALLOCATE(INITCOORD)
  DEALLOCATE(DAMAGE)       ! SCALAR DAMAGE

  ! Additional fundamental variables for stress
  if (.not. hydro) then
     DEALLOCATE(STRESDEV)  ! COMP. OF STRESS DEV. (pixx,piyy,pixy,pith)
     DEALLOCATE(EPSTRAIN)  ! EQUIVALENT PLASTIC STRAIN
     DEALLOCATE(VELO)      ! VIBRATIONAL VELOCITY (ACOUSTIC FLUIDIZATION)
     DEALLOCATE(SCATTER)   ! SCATTERED VIBRATIONAL ENERGY
     DEALLOCATE(PLW)       ! PLASTIC WORK DONE
  endif

  ! Derived variables (from fundamental varibles)
  DEALLOCATE(R)            ! RADIAL COORDS (store them global)
  DEALLOCATE(P)            ! PRESSURE
  DEALLOCATE(MV)           ! MASS AT VERTICES
  DEALLOCATE(RMV)          ! = 1/MV
  DEALLOCATE(CMV)          ! MAT. CONCENTRATION AT VERTICES
  DEALLOCATE(VOIDCONC)     ! CONCENTRATION OF VOID
  DEALLOCATE(TMPR)         ! TEMPERATURE
  DEALLOCATE(TMELT)        ! MELTING TEMPERATURE
  DEALLOCATE(CSOUND)       ! SPEED OF SOUND
  DEALLOCATE(ENTROPY)      ! SPECIFIC ENTROPY
  DEALLOCATE(Q)            ! ARTIFICIAL VISCOSITY
  DEALLOCATE(VOL)          ! CELL VOLUME
  DEALLOCATE(PMAX)         ! MAXIMUM PRESSURE
  DEALLOCATE(A_S)          ! BULK MODULUS
  DEALLOCATE(MU)           ! SHEAR MODULUS
  DEALLOCATE(RRSUM)        ! 1./(SUM OVER LOCAL VERTICES OF R)

  if (.not.hydro) then
     DEALLOCATE(YIELD,YAC,ERATE,PVIBR,FAILURE)
  end if

  ! Temporary variables
  DEALLOCATE(VL,VREL,ROL,PL)

  ! Gravity
  if (allocated(gravity)) deallocate(gravity)

  call ptool_text(io_out,"Finished deallocating memory for fields","")

end subroutine destroymesh
