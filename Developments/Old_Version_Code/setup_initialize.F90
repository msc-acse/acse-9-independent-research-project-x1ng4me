! --------------------------------------------------------------------------------
!> Subroutine to initialize most array or matrix fields.
!!
!! @author Dirk Elbeshausen, MfN
!!
!!     Description                                     Programmer    Date
!!     ------------------------------------------------------------------
!!     Original version (7.0).............................DE   2007/02/25
!!
!<--------------------------------------------------------------------------------
SUBROUTINE SETUP_INITIALIZE
  use mod_isale
  use ptool_interface
  use mod_store
  use mod_gradient
  IMPLICIT NONE

  call ptool_text(io_oin,"","",2)
  call ptool_text(io_oin,"Allocating memory for fields","",1)

  ! Fundamental variables (advected; needed for restart)
  ! store coordinates globally to minimize communication
  ALLOCATE(C(1:2,1:nxp,1:nyp))        ! COORDS IN X,Y (store them global)
  C=0.D0
  ! does this need to be initialized?
  ALLOCATE(V(1:2,IGS:IGE,JGS:JGE))        ! VELOCITY COMPONENTS
  V=0.D0
  ALLOCATE(MC(IGS:IGE,JGS:JGE))           ! MASS IN CELL
  mc=0.D0
  ALLOCATE(CMC(NMAT,IGS:IGE,JGS:JGE))     ! MAT. CONCENTRATION IN CELL
  cmc=0.D0
  ALLOCATE(RO(NMAT,IGS:IGE,JGS:JGE))      ! DENSITY
  ro=0.D0
  ALLOCATE(SIE(NMAT,IGS:IGE,JGS:JGE))     ! SPECIFIC INTERNAL ENERGY
  sie=0.D0
  ALLOCATE(ALPHA(NMAT,IGS:IGE,JGS:JGE))   ! DISTENSION 1./(1.-POROSITY)
  alpha=1.D0
  ALLOCATE(VOLSTRAIN(IGS:IGE,JGS:JGE))    ! VOLUME STRAIN
  volstrain=0.D0
  if (dummyfields) then
     ALLOCATE(DUMMYF(1:ndummy,IGS:IGE,JGS:JGE))       ! DUMMY FIELD #1
     dummyf=0.D0
  end if
  ALLOCATE(INITCOORD(X:Y,IGS:IGE,JGS:JGE))  ! INITAL COORD FIELD
  initcoord=0.D0
  ALLOCATE(DAMAGE(IGS:IGE,JGS:JGE))       ! SCALAR DAMAGE (needed in hydro)
  damage=0.D0

  if (field_stored("Vis",varname_save,nvar_save)) then
     allocate(visco_eff(igs:ige,jgs:jge))
     visco_eff = 0.D0
  endif

  ! Additional fundamental variables for stress
  if (.not. hydro) then
     ALLOCATE(STRESDEV(1:numstrcomp,IGS:IGE,JGS:JGE)) ! COMP. OF STRESS DEV. (pixx,piyy,pixy,pith)
     stresdev=0.D0
     ALLOCATE(EPSTRAIN(IGS:IGE,JGS:JGE))  ! EQUIVALENT PLASTIC STRAIN
     epstrain=0.D0
     ALLOCATE(VELO(IGS:IGE,JGS:JGE))      ! VIBRATIONAL VELOCITY (ACOUSTIC FLUIDIZATION)
     velo(:,:)=0.D0
     ALLOCATE(SCATTER(IGS:IGE,JGS:JGE))   ! VIBRATIONAL SCATTER ENERGY
     scatter=0.D0
     ALLOCATE(PLW(IGS:IGE,JGS:JGE))       ! PLASTIC WORK DONE
     plw=0.D0
  endif

  if (field_stored("exx",varname_save,nvar_save)) then
     allocate(etam_tens(1:numstrcomp,IGS:IGE,JGS:JGE)) ! viscoelastic comp. of stress tensor
     etam_tens = 0.D0
  endif

  if (field_stored("etm",varname_save,nvar_save)) then
     allocate(etam_field(IGS:IGE,JGS:JGE))
     etam_field = 0.D0
  endif

  ! Derived variables (from fundamental varibles)
  ALLOCATE(R(NXP,NYP))                    ! RADIAL COORDS (store them global)
  r=1.D0
  ALLOCATE(P(IGS:IGE,JGS:JGE))            ! PRESSURE
  p=0.D0
  ALLOCATE(MV(IGS:IGE,JGS:JGE))           ! MASS AT VERTICES
  mv=0.D0
  ALLOCATE(RMV(IGS:IGE,JGS:JGE))          ! = 1/MV
  rmv=0.D0
  ALLOCATE(CMV(NMAT,IGS:IGE,JGS:JGE))     ! MAT. CONCENTRATION AT VERTICES
  cmv=0.D0
  ALLOCATE(VOIDCONC(IGS:IGE,JGS:JGE))     ! CONCENTRATION OF VOID
  voidconc=1.D0
  ALLOCATE(TMPR(IGS:IGE,JGS:JGE))         ! TEMPERATURE
  tmpr=0.D0
  ALLOCATE(TMELT(IGS:IGE,JGS:JGE))        ! MELTING TEMPERATURE
  tmelt=0.D0
  ALLOCATE(CSOUND(IGS:IGE,JGS:JGE))       ! SPEED OF SOUND
  csound=0.D0
  ALLOCATE(ENTROPY(NMAT,IGS:IGE,JGS:JGE)) ! SPECIFIC ENTROPY (J/K/kg)
  entropy=0.D0
  ALLOCATE(Q(IGS:IGE,JGS:JGE))            ! ARTIFICIAL VISCOSITY
  q=0.D0
  ALLOCATE(VOL(IGS:IGE,JGS:JGE))          ! CELL VOLUME
  vol=0.D0
  ALLOCATE(PMAX(IGS:IGE,JGS:JGE))         ! MAXIMUM PRESSURE
  pmax=0.D0
  ALLOCATE(A_S(IGS:IGE,JGS:JGE))          ! BULK MODULUS
  A_s=0.D0
  ALLOCATE(MU(IGS:IGE,JGS:JGE))           ! SHEAR MODULUS
  mu=0.D0
  ALLOCATE(RRSUM(IGS:IGE,JGS:JGE))        ! 1./(SUM OVER LOCAL VERTICES OF R)
  rrsum=0.D0

  if (.not.hydro) then
     ALLOCATE(YIELD(IGS:IGE,JGS:JGE))        ! YIELD STRENGTH
     yield=0.D0
     ALLOCATE(YAC(IGS:IGE,JGS:JGE))          ! ACOUSTIC FLUIDIZATION INTENSITY
     YAc=0.D0
     ALLOCATE(ERATE(IGS:IGE,JGS:JGE))        ! STRAIN RATE
     erate=0.D0
     ALLOCATE(PVIBR(IGS:IGE,JGS:JGE))        ! VIBRATIONAL PRESSURE
     pvibr=0.D0
     ALLOCATE(FAILURE(IGS:IGE,JGS:JGE))      ! FAILURE MODE
     failure = FAIL_NONE
     if (field_stored("FaP",varname_save,nvar_save)) then
        ALLOCATE(failure_first(IGS:IGE,JGS:JGE))
        failure_first=FAIL_NONE
     endif

     if (      field_stored("Exx",varname_save,nvar_save) &
          .or. field_stored("Eyy",varname_save,nvar_save) &
          .or. field_stored("Exy",varname_save,nvar_save) &
          .or. field_stored("Eth",varname_save,nvar_save)) then
        ALLOCATE(erate_tens(1:numstrcomp,igs:ige,jgs:jge))
        erate_tens = 0.D0
     endif
     
  end if

  ! Temporary variables
  ALLOCATE(VL(1:2,IGS:IGE,JGS:JGE))       ! VELOCITY COMP. (TEMP. STORED)
  VL=0.D0
  ALLOCATE(VREL(1:2,IGS:IGE,JGS:JGE))     ! VELOCITY COMP. RELATIVE TO GRID
  VREL=0.D0
  ALLOCATE(ROL(NMAT,IGS:IGE,JGS:JGE))     ! DENSITY (TEMP. STORED)
  ROL=0.D0
  ALLOCATE(PL(IGS:IGE,JGS:JGE))           ! PRESSURE (TEMP. STORED)
  PL=0.D0

  ! Allocate memory for gravity vector field
  if (gradient_type==GRADTYPE_CENTRAL .or. gradient_type==GRADTYPE_SELF) then
     ALLOCATE(GRAVITY(1:2,IGS:IGE,JGS:JGE))
     GRAVITY = 0.d0
  end if

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
