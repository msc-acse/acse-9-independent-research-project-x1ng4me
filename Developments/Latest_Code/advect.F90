!--------------------------------------------------------------------------------
!> Fluxing of cell centered quantities mass, energy,
!! and momentum. Then convert momenta to velocities.
!! Phase boundaries are included into this version:
!! Speration between matter1, matter2 and vacuum.
!!
!! This version fluxes by either volume or mass.
!! Stresses are advected if necessary.
!!
!! Uses full-donor cell advection.
!!
!! @authors Kai Wuennemann, MfN; Dirk Elbeshausen, MfN; Gareth Collins, ICL; Xianzheng Li, ICL
!! @Github Alias : x1ng4me(XL)
!!     Description                                     Programmer    Date
!!     ------------------------------------------------------------------
!!     Original version (5.7)                             KAI  2005/03/22
!!     Incorporate porosity (6.0)                         KAI  2005/03/30  
!!     Multiple interface tracking from Gareth (6.1)      KAI  2006/05/25
!!     Conversion into Fortran90                          DE   2007/02/23
!!     Merged volume and mass advection                   DE   2007/05/20
!!     Move out subroutines and add parallelization       XL   2019/08/28
!<--------------------------------------------------------------------------------
subroutine advect
  use mod_isale
  use mod_advect
  use mod_identify_mat
  use omp_lib
  implicit none
  integer i,j,m,l,non(nmat)
  integer cc11,cc12,cc21,cc22,sumcc(nmat)
  integer :: filling(nmat)
  integer :: nomic,nom(nmat),mmax
  integer :: matincell(nmat),countmat,m1,m2
  integer ck11(nmat),ck12(nmat),ck21(nmat),ck22(nmat)
  real*8  rx1,rx2,zx1,zx2,c11,c21,c12,c22,qm
  real*8  r1,r2,r3,r4,tvmat,tvmat_plus,conc,cmax
  real*8  sfac,cc(nmat)
  real*8  s(nmat,9),sumout,sumin,matflux
  real*8  masse(nmat),vmat(nmat)
  real*8  fr,ar,vmod,sfac1
  real*8  zufz,ranval
  real*8  rx1k(nmat),rx2k(nmat),zx1k(nmat),zx2k(nmat)
  real*8, allocatable ::  VolFlux(:,:,:)
  real*8, allocatable ::  initcoordp(:,:,:)
  real*8, allocatable ::  vmom(:,:,:,:,:),vmomp(:,:,:,:)
  real*8, allocatable ::  mp(:,:,:),mvp(:,:)
  real*8, allocatable ::  siep(:,:,:),plwp(:,:)
  real*8, allocatable ::  damagep(:,:),epstrp(:,:)
  real*8, allocatable ::  volstrp(:,:),alphap(:,:,:)
  real*8, allocatable :: stresdevp(:,:,:)
  real*8, allocatable :: velp(:,:)
  real*8, allocatable :: dummyfp(:,:,:)
  character*1, allocatable :: vstat(:,:)
  real*8, parameter :: NO_OUTFLUX = 0.0d0
  real*8, external :: total_flux, mirror_flux
  real*8 fluxpart

  integer, external :: int_bound,priority_material
  integer  nompri !Material number of the priority material for interface construction

  integer, parameter :: f_right=2,f_top=4,f_left=6,f_bottom=8,FL_OUT=9
  integer, parameter :: ADVMOM_HIS=1
  integer :: vmom_dim ! Dimension of vmom and vmomp arrays (2nd index)

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++ Allocate and initialise memory for temporary arrays +++
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  allocate(VolFlux(1:8,igs:ige,jgs:jge))
  allocate(vstat(igs:ige,jgs:jge))
  allocate(voltran_r(nmat,igs:ige,jgs:jge), voltran_t(nmat,igs:ige,jgs:jge))
  allocate(voltran_l(nmat,jgs:jge), voltran_b(nmat,igs:ige))
  allocate(damagep(igs:ige,jgs:jge))
  allocate(initcoordp(X:Y,igs:ige,jgs:jge))
  if (dummyfields) then
     allocate(dummyfp(1:ndummy,igs:ige,jgs:jge))
  end if
   if (ADVMOM .eq. ADVMOM_HIS) then
     vmom_dim = 4 ! HIS method advects all four nodal momenta
  else
     vmom_dim = 1 ! SALE method advects one cell-centred set of momenta
  end if 
  allocate(vmom(X:Y,1:vmom_dim,nmat,igs:ige,jgs:jge), vmomp(X:Y,1:vmom_dim,igs:ige,jgs:jge)) 
  allocate(mp(nmat,igs:ige,jgs:jge), mvp(igs:ige,jgs:jge), siep(nmat,igs:ige,jgs:jge)) 
  allocate(volstrp(igs:ige,jgs:jge), alphap(nmat,igs:ige,jgs:jge)) 
  if (.not. hydro) then
     allocate(stresdevp(1:4,igs:ige,jgs:jge))
     allocate(plwp(igs:ige,jgs:jge), epstrp(igs:ige,jgs:jge), velp(igs:ige,jgs:jge))
  end if
  ! If Eulerian calculation (not ALE), set the velocity of flow relative to
  ! grid, VREL = VL - VG, equal to flow velocity as grid velocity is zero
  if (ALE_MODE == EULERIAN) VREL = VL
  
  ! +++ Put initialization into first loop, add OpenMP in the first loop.
  ! +++ Initialise memory for temporary arrays.
  ! +++ Calculate and store volume transfer between cells 
  !$omp parallel private(i,j,fr,ar)
  !$omp do 
  do j=jgs,jge
     do i=igs,ige
        VolFlux(:,i,j) = 0.D0
        voltran_r(:,i,j) = 0.D0
        voltran_t(:,i,j) = 0.D0
        voltran_l(i,j) = 0.D0
        voltran_b(i,j) = 0.D0
        damagep(i,j) = 0.D0
        mp(:,i,j) = 0.D0
        mvp(i,j) = 0.D0
        siep(:,i,j) = 0.D0
        volstrp(i,j) = 0.D0
        alphap(:,i,j) = 0.D0
        if (.not. hydro) then
           stresdevp(:,i,j) = 0.D0
           plwp(i,j) = 0.D0
           epstrp(i,j) = 0.D0
           velp(i,j) = 0.D0
        end if
        vmom(:,:,:,i,j) = 0.D0
        vmomp(:,:,i,j) = 0.D0
        if (dummyfields) then
           dummyfp(:,i,j) = 0.D0
        end if
        vstat(i,j) = 'E'
        initcoordp(:,i,j) = 0.D0
        fr = 0.D0
        ar = 0.D0
        ! Left side of cell
        if (.not. i.eq.ige .OR. j.eq.jge) then
           if (i.eq.1) then ! Left boundary
              if (BCTYPE(LEFT) == 2) then ! Cont. outflow through left boundary
                 ! Compute flux (fr,ar) for right side of cell to mirror
                 VolFlux(FL_LEFT,i,j) = total_flux(i,j,FL_RIGHT,fr,ar)
              else ! All other boundary conditions
                 ! Zero flux to mirror on boundary
                 VolFlux(FL_LEFT,i,j) = 0.D0
              end if
           else
              VolFlux(FL_LEFT,i,j) = total_flux(i-1,j,FL_RIGHT,fr,ar)
              VolFlux(FL_LEFT,i,j) = mirror_flux(fr, ar) 
           end if

           ! Bottom side of cell
           if (j.eq.1) then ! Bottom boundary
              if (BCTYPE(BOTTOM) == 2) then ! Cont. outflow through bottom boundary
                 ! Compute flux (ft,at) for top side of cell to mirror
                 VolFlux(FL_BOTTOM,i,j) = total_flux(i,j,FL_TOP,fr,ar)
              else ! All other boundary conditions
                 ! Zero flux to mirror on boundary
                 VolFlux(FL_BOTTOM,i,j) = 0.D0
              end if
           else
              VolFlux(FL_BOTTOM,i,j) = total_flux(i,j-1,FL_TOP,fr,ar)
              VolFlux(FL_BOTTOM,i,j) = mirror_flux(fr, ar)
           end if

           VolFlux(FL_RIGHT,i,j) =total_flux(i,j,FL_RIGHT,fr,ar)
           VolFlux(FL_TOP,i,j) =total_flux(i,j,FL_TOP,fr,ar)
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel



  ! **************************************************************************
  !                      declaration of some quantities
  ! **************************************************************************
  !     i,j,ij,ipj,ijp,ipjp,imj,ijm: pointer on the main array
  !     non: =1 for simple cell (full or empty), =0 if internal boundary in cell
  !     cc11,cc12,cc21,cc22: =1 if vertex is part of the matter filled volume
  !     sumcc: sum of cc11,cc12,cc21,cc22, if =4 than no internal boundary
  !     mat: to distinguish different cases of material filling 
  !          if mat(n)=1 type of matter available, else =0
  !     rx1,rx2,zx1,zx2: rel. position of int.Bound. crossing cell sides
  !     c11,c12,c21,c22: concentration of matter at vertices (=cmv(i,j))
  !     qm: 1/4 of mass in cell
  !     qmomu,qmomv: 1/4 of momentum
  !     sfac, sbuf, ssbuf: buffers, used at different places
  !                        for buffering and fact.
  !     ddx,ddy: spatial increment
  !     s: dimension(9,3), 9=4+4+1 (4-sides with in- and ouflux, 
  !        1 for total values), 3=different materials, 
  !        1=vacuum,2=mat_1,3=mat_2, ... n=mat_n
  !     VolFlux: dimension(8), Volume outflux through the sides (1,3,5,7)
  !     sumout,sumin,matflux: sum of in- and outflux
  !
  !     ucenter,vcenter: centeres values of velocity components (cell center)
  !     donor : array (1-8) to store the koord. i,j of right,left,top,bottom side
  !
  ! --- Array variables with dimension of nmat: 1=total content or vacuum
  !                                             2=material 1
  !                                             3=material 2
  !                                             ...
  !                                             n=material n
  !
  !     vmat: volume of matter with 2=mat_1, 3=mat_2, ...n=mat_n,
  !           1=volume of void
  !     tvmat: total volume of matter (mat_1 - mat_n)
  !     masse: total mass of each material (1=sum of mat_1+mat_2+...mat_n)
  !     vmom: momentum with 1=total, 2=mat_1, 3=mat_2, ...n=mat_n
  !
  !
  ! --- Array variables with dimension of ([nmat],nxp,nyp):
  !     voltran_r,voltran_t: uesed to store volume transport for each material
  !     siep(nmat,i,j) : to store internal energy for each material during
  !                     the computation of advection (1=energy of total matter)
  !     vmomp: to store momentum during the computation of advection
  !     stresdevp: to store stress components during the computation of advection.
  !     velp: to store vibrations during the computation of advection
  !

  ! +++++++++++++++++++++++++++++++++++++++++++++++
  ! +++ Calculate volume transfer for each cell +++
  ! +++++++++++++++++++++++++++++++++++++++++++++++
  ! Initialise fluxes
  ! +++ Construct material boundaries in mixed cells and calculate individual volume fluxes
  ! +++ for each material based on material boundary position
  ! +++ Apply OpenMP for the large loop, using omp do which equals to omp_scheduling(static)

  !$omp parallel private(i,j,fr,ar,r1,r2,r3,r4,rx1,rx2,zx1,zx2,sumcc,nompri,cc,nomic,nom,filling),&
  !$omp& private(cmax,mmax,l,m,conc,c11,c21,c12,c22,cc11,cc12,cc21,cc22,ck11,ck12,ck21,ck22,rx1k,rx2k,zx1k,zx2k),&
  !$omp& private(tvmat,matflux,vmat,non,s,ranval)
  !$omp do
  do j=js,je
     fr=0.D0
     ar=0.D0
     do i=is,ie

        ! Set radial coordinates
        r1  =r(i+1,j)    
        r2  =r(i+1,j+1)
        r3  =r(i,j+1)
        r4  =r(i,j)

        ! Initialise interface positions
        rx1=0.D0
        rx2=0.D0
        zx1=0.D0
        zx2=0.D0            

        ! Initialize other variables
        sumcc = 0
        nompri= 0   !Priority material number

        ! Store material Concentration in new variable
        cc(1:nmat) = cmc(1:nmat,i,j) 

        ! Distinguish different cases of material filling in a cell
        nomic = cell(i,j,1)%nm
        nom = cell(i,j,1)%matid(1:nmat) !DErelic: note: cell%matid(1:MAXMAT) --> change later...
        call advect_filling(cell(i,j,1),filling)

        ! Find maximum of material concentration
        cmax=0.D0
        mmax=0
        cc(VOID) = 1.D0-cmc(TOTAL,i,j) ! Redefine cc(1) as concentration of void
        do l=1,nomic
           m=nom(l)
           if (cc(m).gt.cmax) then
              mmax=m
              cmax=cc(m)
           endif
        enddo

        !-----------Simple Cell, only one material or vacuum
        if (nomic.eq.1) then 
           if (nom(1).ne.VOID) then
              vstat(i+1,j) ='M'      ! Mark material vertices
              vstat(i+1,j+1)='M'
              vstat(i,j+1) ='M'
              vstat(i,j)  ='M'
           endif
           m=nom(1)
           non(m)=1  ! No internal boundary

           !-----------Mixed Cell; multiple materials
        elseif (nomic.gt.1 .and. nomic.le.nmat) then

           ! Let's try something new here... Why not calculate the
           ! internal boundary location for each material and
           ! account for this in the fluxing?  Previously we only
           ! determined one boundary: that between the dominant material
           ! component and the rest...

           ! In simple vac+mat cell, set priority material as non vacuum
           !if (nomic.eq.2.and.nom(1).eq.1) nompri=nom(2)

           ranval=zufz()! Define ran here so that in case of multiple isolated 
           ! matter the interfaces are consistent...

           ! So, first, let's loop over the materials present...
           do l=1,nomic

              ! Recall the material number for the species present
              m=nom(l)

              if (m<1) then
                 if (nomic==1) write(IOERR,*) "YES NOMIC = 1",nomic
                 write(IOERR,*) "m < 1 = ",m
                 write(IOERR,*) "CELL  = ",i,j,ncyc
                 write(IOERR,*) "nomic = ",nomic
                 write(IOERR,*) "matid = ",filling
                 write(IOERR,*) "NOM   = ",nom(1:nomic)
                 write(IOERR,*) "CONC  = ",cc(1:nmat)
              endif

              ! define the inputs (vertex mass concentrations) 
              ! for the internal-boundary finding routine.
              !
              ! Note change: in the case that vacuum exists in the cell, the 
              ! routine will use the void fraction to compute one interface.
              ! Previously, the vacuum interface was calculated implicitly,
              ! by calculating the interface for the total material, which is
              ! 1.-(void fraction).
              !
              ! In all other cases, the routine uses the volume fraction
              ! of the material species to define the boundary
              conc=cc(m)
              if (m.eq.VOID) then
                 c11 =1.D0-cmv(m,i,j) ! Set concentration of matter at vertices
                 c21 =1.D0-cmv(m,i+1,j) ! to new variables
                 c12 =1.D0-cmv(m,i,j+1)
                 c22 =1.D0-cmv(m,i+1,j+1)
              else
                 c11 =cmv(m,i,j) ! Set concentration of matter at vertices
                 c21 =cmv(m,i+1,j) ! to new variables
                 c12 =cmv(m,i,j+1)
                 c22 =cmv(m,i+1,j+1)
              end if

              cc11=0        ! Initialize variables if vertex is included
              cc12=0        ! in material filling of the cell
              cc21=0
              cc22=0

              ! Calculate the position (rx1/2, zx1/2) of internal
              ! boundary between different materials.
              ! ccxx are set to 1 if vertex is included in
              ! material filling of the cell; 0 if it is not
              ! non(m) == 0 for successful interface construction
              ! non(m) \= 0 for unsuccessful interface construction
              ! non(m) == 1 for isolated matter or cell full
              ! non(m) == 2 for split matter
              non(m)=int_bound(i,j,rx1,rx2,zx1,zx2,&
                   cc11,cc12,cc21,cc22,c11,c12,c21,c22,conc)

              ! Check whether "maximum material" has an interface
              ! defined; if it does this should be the priority material
              !if (nompri.eq.0.and.m.eq.mmax.and.non(m).eq.0) nompri=mmax

              ! Now account for "isolated" cells: cells that
              ! should contain a boundary but for which none
              ! can easily be defined...
              sumcc(m)=cc11+cc12+cc21+cc22
              if (isolated_treatment>0) then
                 if (sumcc(m).eq.0.or.non(m).eq.2) then
                    call advect_isolate(i,j,m,nomic,filling,nom,non(m),&
                         rx1,rx2,zx1,zx2,cc11,cc12,cc21,cc22,ranval)
                    ! If cell contents is deleted, isolate returns
                    ! with nomic==1; in this case, there is no longer 
                    ! an interface, so skip interface construction
                    if (nomic.eq.1) then
                       mmax = nom(1)
                       goto 20
                    end if
                 endif
              endif

              ! Store the parameters needed to compute fluxes
              ck11(m)=cc11
              ck12(m)=cc12
              ck21(m)=cc21
              ck22(m)=cc22
              rx1k(m)=rx1
              rx2k(m)=rx2
              zx1k(m)=zx1
              zx2k(m)=zx2

              ! For cells with some vacuum present record the
              ! vertices that have matter attached
              if (m.eq.VOID) then !Note that here cc## == 0 implies material vertex
                 if (cc11.eq.0) vstat(i,j)  ='M'
                 if (cc21.eq.0) vstat(i+1,j) ='M'
                 if (cc12.eq.0) vstat(i,j+1) ='M'
                 if (cc22.eq.0) vstat(i+1,j+1)='M'
              endif

           enddo

           ! If material is not vacuum then all the
           ! vertices should be considered
           ! for the momentum calculation. . .
           if (nom(1).ne.VOID) then
              vstat(i  ,j)  ='M'
              vstat(i+1,j)  ='M'
              vstat(i  ,j+1)='M'
              vstat(i+1,j+1)='M'
           endif

        else
           write(IOERR,*) 'Wrong number of mat in cell: ',nomic
           call mystop("ERROR IN ADVECT: Wrong number of materials in cell (nomic)")
        endif

20      continue

        ! +++ Start calculation flux +++++++++++++++++++++++++++++++
        !

        ! Define the priority material (0 if no prioritization),
        ! depending on the user option "priority_option"
        nompri = priority_material(nomic,nom,mmax,non)


        ! For each cell side, determine the total volume outflux through the
        ! side (total_flux/mirror_flux), then modify this total to account for 
        ! the position of the material interfaces (cases_flux)...

        call cases_flux(FL_LEFT,s,VolFlux(:,i,j),ck11,ck12,zx1k,r4,r4,nom,nomic,nompri,cc)

        call cases_flux(FL_BOTTOM,s,VolFlux(:,i,j),ck11,ck21,rx1k,r4,r1,nom,nomic,nompri,cc)

        ! Compute right-side outflux and store flux (fr,ar) for left-side of right neighbour 
        call cases_flux(FL_RIGHT,s,VolFlux(:,i,j),ck21,ck22,zx2k,r4,r4,nom,nomic,nompri,cc)

        ! Compute top-side outflux and store flux (ft,at) for bottom-side of top neighbour 
        call cases_flux(FL_TOP,s,VolFlux(:,i,j),ck12,ck22,rx2k,r3,r2,nom,nomic,nompri,cc)

        ! Total outflux of volume of matter/vacuum +++++++++++++++++++++
        tvmat=vol(i,j)*cmc(TOTAL,i,j) ! Calc total volume of matter

        matflux=0.D0
        do m=1,nmat
           s(m,FL_OUT)=s(m,FL_RIGHT)+s(m,FL_TOP)+s(m,FL_LEFT)+s(m,FL_BOTTOM) !total flux for each mat.
           matflux=matflux+s(m,FL_OUT)
           if (m.eq.VOID) then 
              vmat(VOID)=vol(i,j)-tvmat ! old volume of void 
           else
              vmat(m)=vol(i,j)*cmc(m,i,j) ! old volume of matter
           endif
        enddo

        ! Check if donor cell has enough matter/vacuum (only cells with int.Bound.)
        if (nomic.ge.2) then            
           if (nmat.ge.5) then 
              call advect_adjust_flux_simple(s,nomic,nom,vmat)
           else
              call advect_adjust_flux(s,nomic,nom,vmat)
           end if
        endif

        ! Store fluxes if they present outflux
        do m=1,nmat ! loop for each material
           if (s(m,FL_RIGHT).lt.0.D0) voltran_r(m,i,j)=s(m,FL_RIGHT)
           if (s(m,FL_TOP).lt.0.D0)   voltran_t(m,i,j)=s(m,FL_TOP)
           if (i.gt.1.and.s(m,FL_LEFT).lt.0.D0)   voltran_r(m,i-1,j)=-s(m,FL_LEFT)
           if (j.gt.1.and.s(m,FL_BOTTOM).lt.0.D0) voltran_t(m,i,j-1)=-s(m,FL_BOTTOM)
           if (i.eq.1.and.s(m,FL_LEFT).lt.0.D0)   voltran_l(m,j)=s(m,FL_LEFT)
           if (j.eq.1.and.s(m,FL_BOTTOM).lt.0.D0) voltran_b(m,i)=s(m,FL_BOTTOM)

        enddo

        ! Calculate the cell-centered momentum
        if (ADVMOM .eq. ADVMOM_HIS) then
           call advect_momentum_his(i, j, nom, nomic, vmom)
        else
           call advect_momentum(i, j, sumcc, nom, nomic, vstat, vmom)
        end if

     enddo
  enddo
  !$omp end do
  !$omp end parallel
  deallocate(VolFlux)
  ! +++ End first array loop

  ! +++++++++++++++++++++++++++++++++++++
  ! +++      Calculate advection      +++
  ! +++++++++++++++++++++++++++++++++++++
  ! +++ Calculate advection: account for outflow and influx of fields across cell boundaries;
  ! +++ update material densities and volume fractions.
  ! +++ Use same parallel scheduling method, list all private parameters to avoid conflicts.
  
  !$omp parallel private(i,j,sumout,s,sumin,tvmat,vmat,tvmat_plus,masse,countmat,matincell),&
  !$omp& private(m,m1,m2,sfac,fluxpart)
  !$omp do
  do j=js,je
     do i=is,ie
        ! restore in- and outflux of volume
        call advect_restore_fluxes(i, j, s) 

        ! calculate total amount of in- and outfluxes...
        s(:,FL_OUT)=s(:,FL_RIGHT)+s(:,FL_TOP)+s(:,FL_LEFT)+s(:,FL_BOTTOM)

        sumout =0.D0
        sumin  =0.D0
        tvmat  =0.D0
        do m=firstmat,nmat
           sumout =sumout+s(m,FL_OUT)
           sumin  =sumin +s(m,F_RIGHT)+s(m,F_TOP)+s(m,F_LEFT)+s(m,F_BOTTOM) !total influx of matter
           vmat(m)=vol(i,j)*cmc(m,i,j)                 ! amount of each material
           tvmat=tvmat+vmat(m)                       !total outflux of matter
        enddo
       
        vmat(VOID)=vol(i,j)-tvmat

        ! --------- 1. compute effects of outflux only                        
        call advect_outflux(nomic, nom, vmom_dim, vmom, siep, alphap, masse, fluxpart, volstrp, i, j, sumout, s, tvmat, vmat, velp, plwp, dummyfp, stresdevp, damagep, epstrp, initcoordp, vmomp)
        ! --------- 2. compute effects of influx
        if (i.gt.1 ) call advect_influx(f_left  ,i-1,j  ,nomic, nom, vmom_dim, vmom, siep, alphap, masse, fluxpart, volstrp, i, j, sumout, s,&
          &tvmat, vmat, velp, plwp, dummyfp, stresdevp, damagep, epstrp, initcoordp, vmomp)
        if (i.lt.ie) call advect_influx(f_right ,i+1,j  ,nomic, nom, vmom_dim, vmom, siep, alphap, masse, fluxpart, volstrp, i, j, sumout, s,&
          &tvmat, vmat, velp, plwp, dummyfp, stresdevp, damagep, epstrp, initcoordp, vmomp)
        if (j.gt.1 ) call advect_influx(f_bottom,i  ,j-1,nomic, nom, vmom_dim, vmom, siep, alphap, masse, fluxpart, volstrp, i, j, sumout, s,&
          &tvmat, vmat, velp, plwp, dummyfp, stresdevp, damagep, epstrp, initcoordp, vmomp)
        if (j.lt.je) call advect_influx(f_top   ,i  ,j+1,nomic, nom, vmom_dim, vmom, siep, alphap, masse, fluxpart, volstrp, i, j, sumout, s,&
          &tvmat, vmat, velp, plwp, dummyfp, stresdevp, damagep, epstrp, initcoordp, vmomp)
        
        ! Retain the total advected volume (before accounting for compression/expansion)
        ! needed to recover advected variables when using transport volumes. . .
        tvmat_plus = tvmat

        ! Capping and totals of mass, volume and energy. . .
        vmat(VOID)     =min(vol(i,j),max(0.D0,vmat(VOID)))
        tvmat          =0.D0
        masse(TOTAL)   =0.D0
        siep(TOTAL,i,j)=0.D0
        countmat = 0
        matincell = 0
        do m=firstmat,nmat
           masse(m)       =max(0.D0,masse(m))
           siep(m,i,j)    =max(0.D0,siep(m,i,j))
           vmat(m)        =min(vol(i,j),max(0.D0,vmat(m)))
           tvmat          =tvmat+vmat(m)
           masse(TOTAL)   =masse(TOTAL)+masse(m)
           siep(TOTAL,i,j)=siep(TOTAL,i,j)+siep(m,i,j)
           ! Identify the material numbers in cell to help account for
           ! volume change. . . 
           if (vmat(m) .gt. eps_min*vol(i,j)) then 
              countmat = countmat + 1
              matincell(countmat) = m
           end if
        enddo

        ! --------- Account for compression / expansion of non-void material -----------
        ! If the material volume exceeds the cell volume (compression) or if there is
        ! insufficient material to fill the cell and no void (expansion) then
        ! adjust the volumes of each component to fill the cell with matter. . .
        if ( tvmat.gt.oeps_min*vol(i,j) .or. &
             (tvmat.gt.0.D0 .and. vmat(VOID).lt.eps_min*vol(i,j))) then
           if (balance_pressure .eq. 1 .and. countmat .eq. 2) then
              ! If there are only two materials in the cell, we can compress/expand
              ! the materials based on their compressibility to achieve an approximately
              ! equal pressure change for the two components. . .
              m1 = matincell(1)
              m2 = matincell(2)
              vmat(m1) = min((mat(m1)%bulkmod0-mat(m2)%bulkmod0*(1.D0-vol(i,j)/vmat(m2)))/ &
                         (mat(m1)%bulkmod0/vmat(m1)+mat(m2)%bulkmod0/vmat(m2)),vol(i,j))
              vmat(m2) = max(0.D0,vol(i,j)-vmat(m1))
              ! Another approach, if one material is much more compressible than the other
              ! is to take up all the volume change in the more compressible component. . .
              !vmat(m1)=min(vol(i,j),vmat(m1))      ! Volume of asteroid smaller than cell.
              !vmat(m2)=max(0.D0,vol(i,j)-vmat(m1)) ! Volume of air is the remainder.
           else
              ! Compress or expand each material proportionally with
              ! volume fraction of the material.  In essence, this assumes
              ! that the compressibility of each material is the same. . .
              sfac=vol(i,j)/tvmat
              do m=firstmat,nmat
                 vmat(m)=vmat(m)*sfac
              enddo
           end if
           vmat(VOID)=0.D0
           tvmat  =vol(i,j)

           ! If 'part_pres' is true, we treat vacuum consistently
           ! with other materials; i.e. so that the compressibility
           ! of vacuum is assumed to be the same as that of all
           ! other materials.  This has been found to help improve
           ! interface construction when the target has high strength
           ! when used in combination with computing pressure, and
           ! the effects of pressure gradients in partially full cells.
           !
           ! If there is vacuum in the cell
        else
           if (part_pres) then
              sfac=vol(i,j)/(tvmat+vmat(VOID))
              tvmat=0.D0
              do m=firstmat,nmat
                 vmat(m)=vmat(m)*sfac
                 tvmat=tvmat+vmat(m)
              enddo
              vmat(VOID)=min(vol(i,j),max(0.,vol(i,j)-tvmat))
           endif
        endif
        ! ---- End of accounting for compression / expansion of non-void material --------
        

        ! --------- compute new values due to in- and ouflux
        masse(TOTAL)   =0.D0
        siep(TOTAL,i,j)=0.D0
        do m=firstmat,nmat
           masse(TOTAL)   =masse(TOTAL)+masse(m)
           siep(TOTAL,i,j)=siep(TOTAL,i,j)+siep(m,i,j)
        enddo

        do m=1,nmat
           masse(m)=max(0.D0,masse(m))
           mp(m,i,j) =masse(m)
           if (m.eq.TOTAL) then
              cmc(m,i,j)=tvmat/vol(i,j) ! only for total matter
           else
              cmc(m,i,j)=vmat(m)/vol(i,j)
           endif
        enddo

        ! --------- fill or empty cells ----------------------- 
        ! --------- with respect to the amount of matter (mass)
        call advect_cutoff(i,j,vmat,tvmat,siep,mp)

        ! Determine all other parameters that are advected 
        ! accordingly to volume transport
        if (mp(TOTAL,i,j).gt.0.D0 .and. tvmat.gt.0.D0) then                              
           ! Calc new distension according to mass flux
           do m=firstmat,nmat
              if (mp(m,i,j).gt.0.D0) then
                 alphap(m,i,j) =max(1.D0,alphap(m,i,j) /mp(m,i,j))
              else
                 alphap(m,i,j) =1.D0
              endif
           enddo

           select case (ADVTYPE)
           case (ADVTYPE_VOL)    ! according to volume flux
              fluxpart = tvmat_plus ! This is the total volume of material after advection
           case (ADVTYPE_MASS)   ! according to mass flux
              fluxpart = mp(TOTAL,i,j)
           case default
              write(IOERR,*) "ADVTYPE UNKNOWN IN ADVECT : ",advtype
              call mystop("ERROR IN ADVECT: ADVTYPE UNKNOWN")
           end select

           ! Calc new other parameters
           ! according to transport mass OR volume
           ! (fluxpart contains the scaling value)
           volstrp(i,j)       = volstrp(i,j)       / fluxpart
           damagep(i,j)       = damagep(i,j)       / fluxpart
           if (dummyfields) &
              dummyfp(1:ndummy,i,j)       = dummyfp(1:ndummy,i,j)/ fluxpart
           initcoordp(X:Y,i,j)=initcoordp(X:Y,i,j) / fluxpart

           ! stresses
           if (.not. HYDRO) then 
              plwp(i,j)          = plwp(i,j)          / mp(TOTAL,i,j)
              epstrp(i,j)        = epstrp(i,j)        / fluxpart
              epstrp(i,j)        = min(epstrp(i,j),  10.D0)
              stresdevp(1:4,i,j) = stresdevp(1:4,i,j) / fluxpart
              velp(i,j)          = velp(i,j)          / mp(TOTAL,i,j)
           end if
           
           ! Limit shear and volume strain to a reasonable range
           ! to avoid spurious results
           volstrp(i,j)=max(volstrp(i,j),-10.D0)
           volstrp(i,j)=min(volstrp(i,j), 10.D0)

        else               

           alphap(firstmat:nmat,i,j) = 1.D0
           volstrp(i,j) = 0.d0
           damagep(i,j) = 0.d0
           if (dummyfields) dummyfp(1:ndummy,i,j) = 0.d0
           initcoordp(X,i,j) =initcoord(X,i,j)
           initcoordp(Y,i,j) =initcoord(Y,i,j)

           ! stresses
           if (.not. hydro) then 
              stresdevp(1:4,i,j) = 0.d0
              epstrp(i,j)        = 0.D0
              velp(i,j)          = 0.d0               
           end if

        endif

        voidconc(i,j)=vmat(VOID) ! Store amount of void !

     enddo
  enddo
  !$omp end do
  !$omp end parallel

  ! +++ End secound array loop

  ! +++ Update global fields from temporary fields; 
  ! +++ Calculate bulk fields for material specific quantities;
  ! +++ Identify materials in cell based on volume fractions.
  ! +++ Applied OpenMP
  call advect_finalize()
     
  ! +++ Calculate nodal volume fractions from cell volume fractions (REFACTORED);
  ! +++ Compute nodal mass (REFACTORED) and reciprocal mass, based on cell mass.
  ! +++ Prepare velocity for update;
  ! +++ Update nodal velocities based on cell-entered momenta fluxes (REFACTORED);
  ! +++ Cap nodal velocity to eliminate spurious velocities;
  ! +++ Remove calling advect_routine.F90 file to reduce communication time (REFACTORED)
  !$omp parallel private(i,j,sfac1,vmod)
  sfac1=max(VEL_CUTOFF,VEL_MIN)
  if (VEL_CUTOFF .le. 0.D0) then
     sfac1 = 1D32
  endif
  !$omp do
  do j=js,jep
     do i=is,iep
        cmv(:,i,j) = 0.D0
        if (j.eq.1) then
           if (i.eq.1) then
              mvp(i,j) = mvp(i,j) + 0.25d0*mc(i,j)
              cmv(:,i,j) = cmv(:,i,j) + cmc(:,i,j)
              call cal_cmv(i,j,mvp)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,4,i,j)*rmv(i,j)
           else if (i.eq.iep) then
              mvp(i,j) = mvp(i,j) + 0.25d0*mc(i-1,j)
              cmv(:,i,j) = cmv(:,i,j) + cmc(:,i-1,j)
              call cal_cmv(i,j,mvp)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i-1,j)*rmv(i,j)
           else
              mvp(i,j) = mvp(i,j) + 0.25d0*(mc(i-1,j) + mc(i,j))
              cmv(:,i,j) = cmv(:,i,j) +0.5D0*(cmc(:,i-1,j)+cmc(:,i,j))
              call cal_cmv(i,j,mvp)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,1,i-1,j) + vmomp(X:Y,4,i,j))
           end if
        else if (j.eq.jep) then
           if (i.eq.1) then
              mvp(i,j) = mvp(i,j) + 0.25d0*mc(i,j-1)
              cmv(:,i,j) = cmv(:,i,j) + cmc(:,i,j-1)
              call cal_cmv(i,j,mvp)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,3,i,j-1)*rmv(i,j)
           else if (i.eq.iep) then
              mvp(i,j) = mvp(i,j) + 0.25d0*mc(i-1,j-1)
              cmv(:,i,j) = cmv(:,i,j) + cmc(:,i-1,j-1)
              call cal_cmv(i,j,mvp)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,2,i-1,j-1)*rmv(i,j)
           else
              mvp(i,j) = mvp(i,j) + 0.25d0*(mc(i-1,j-1) + mc(i,j-1))
              cmv(:,i,j) = cmv(:,i,j) +0.5D0*(cmc(:,i-1,j-1)+cmc(:,i,j-1))
              call cal_cmv(i,j,mvp)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,2,i-1,j-1) + vmomp(X:Y,3,i,j-1))
           end if
        else if (i.eq.1) then
           if (j.gt.1 .AND. j.lt.jep) then
              mvp(i,j) = mvp(i,j) + 0.25d0*(mc(i,j) + mc(i,j-1))
              cmv(:,i,j) = cmv(:,i,j) +0.5D0*(cmc(:,i,j)+cmc(:,i,j-1))
              call cal_cmv(i,j,mvp)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,4,i,j) + vmomp(X:Y,3,i,j-1))
           end if
        else if (i.eq.iep) then
           if (j.gt.1 .AND. j.lt.jep) then
              mvp(i,j) = mvp(i,j) + 0.25d0*(mc(i-1,j) + mc(i-1,j-1))
              cmv(:,i,j) = cmv(:,i,j) +0.5D0*(cmc(:,i-1,j)+cmc(:,i-1,j-1))
              call cal_cmv(i,j,mvp)
              V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*rmv(i,j)*(vmomp(X:Y,1,i-1,j) + vmomp(X:Y,2,i-1,j-1))
           end if
        else
           mvp(i,j) = mvp(i,j) + 0.25d0*mc(i,j)
           mvp(i,j) = mvp(i,j) + 0.25d0*mc(i-1,j)
           mvp(i,j) = mvp(i,j) + 0.25d0*mc(i,j-1)
           mvp(i,j) = mvp(i,j) + 0.25d0*mc(i-1,j-1)
           cmv(:,i,j) = cmv(:,i,j) +0.25d0*cmc(:,i,j)
           cmv(:,i,j) = cmv(:,i,j) +0.25d0*cmc(:,i-1,j)
           cmv(:,i,j) = cmv(:,i,j) +0.25d0*cmc(:,i,j-1)
           cmv(:,i,j) = cmv(:,i,j) +0.25d0*cmc(:,i-1,j-1)
           call cal_cmv(i,j,mvp)
           V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,4,i,j)*rmv(i,j)
           V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,3,i,j-1)*rmv(i,j)
           V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,1,i-1,j)*rmv(i,j)
           V(X:Y,i,j) = V(X:Y,i,j) + 0.25D0*vmomp(X:Y,2,i-1,j-1)*rmv(i,j)
        end if

        vmod=sqrt(V(X,i,j)**2+V(Y,i,j)**2)
        if (vmod.gt.sfac1) then
           V(X:Y,i,j)=V(X:Y,i,j)*sfac1/vmod
        endif

     enddo
  enddo
  !$omp end do
  !$omp end parallel

  ! +++ Set Boundary
  call update_boundary(v)

  ! +++ Move tracers here if moving according to material fluxes
  !     Need to call now, before voltran arrays are deallocated
  if (tracer_motion.eq.TR_MAT) call movetracer
  deallocate(vstat, voltran_r, voltran_t, voltran_l, voltran_b)
  deallocate(vmom)
  deallocate(vmomp)
  deallocate(mp)
  deallocate(mvp)
  deallocate(siep)
  deallocate(volstrp, alphap, initcoordp, damagep)
  if (dummyfields) deallocate(dummyfp)
  if (.not. hydro) deallocate(stresdevp, plwp, epstrp, velp)


contains

  ! ------------------------------------------------------------------------------
  !> This subroutine overwrites the original fields with the advected
  !! quantities (dummy arrays) and calculates the new vertex mass
  !!
  !! @authors Kai Wuennemann, MfN; Dirk Elbeshausen, MfN; Xianzheng Li, ICL
  !! @GitHub: x1ng4me(XL)
  !!     Description                                     Programmer    Date
  !!     ------------------------------------------------------------------
  !!     Original version (5.7).............................KAI  2005/03/22
  !!     Conversion into Fortran90..........................DE   2007/02/24
  !!     Add Parallelization................................XL   2019/08/28
  !< ------------------------------------------------------------------------------
  !The reason for reamining this loop inside the advect is it has a entire loop.
  !Applied OpenMP and refactored part of the subroutine advect_finalize()
  subroutine advect_finalize()
    use mod_cell, only : cell
    implicit none
    integer :: nump
    real*8, external :: calc_vertex_mass

    !$omp parallel private(i,j,m)
    !$omp do
    do j=js,je
       do i=is,ie
          ! Put Cell mass Calculation inside the loop to enable OpenMP for this part
          mc(i,j) = mp(TOTAL,i,j)
          ! Distension and specific internal energy
          alpha(TOTAL,i,j)=0.D0
          do m=firstmat,nmat
             sie(m,i,j)=siep(m,i,j) 
             alpha(m,i,j) =max(1.D0,min(mat(m)%alphamax,alphap(m,i,j)))
             if (cmc(m,i,j).gt.eps_min) &
                  alpha(TOTAL,i,j)=alpha(TOTAL,i,j) + & 
                  cmc(m,i,j)/(alpha(m,i,j)*cmc(TOTAL,i,j))  
          enddo
          if (alpha(TOTAL,i,j).gt.0.D0) alpha(TOTAL,i,j)=1.D0/alpha(TOTAL,i,j)
          sie(TOTAL,i,j)  =siep(TOTAL,i,j) 

          ! Non-stress-related fields
          volstrain(i,j)=volstrp(i,j)
          if (dummyfields) dummyf(1:ndummy,i,j)   =dummyfp(1:ndummy,i,j)
          initcoord(X:Y,i,j) =initcoordp(X:Y,i,j)
          damage(i,j)   =max(0.D0,min(1.D0,damagep(i,j)))

          ! Stress-related fields
          if (.not. hydro) then 
             stresdev(1:4,i,j) = stresdevp(1:4,i,j)
             velo(i,j)         = dsqrt(max(0.D0,velp(i,j)))
             plw(i,j)          = plwp(i,j)
             epstrain(i,j)     = epstrp(i,j)
          end if
          call identify_mat_cell(eps_min,oeps_min,cell(i,j,1),i,j,1)

       enddo
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine advect_finalize

end subroutine advect

  ! ------------------------------------------------------------------------------
  !> Subroutine to restore in- and outfluxes for current cell
  !!
  !! @authors Kai Wuennemann, MfN; Dirk Elbeshausen, MfN; Xianzheng Li, ICL
  !! @GitHub: x1ng4me(XL)
  !!     Description                                     Programmer    Date
  !!     ------------------------------------------------------------------
  !!     Original version (5.7).............................KAI  2005/03/22
  !!     Conversion into Fortran90..........................DE   2007/02/24
  !!     Move subroutine outside advect.....................XL   2019/08/28
  !< ------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------
  ! Move outside the subroutine advect for parallelization
  subroutine advect_restore_fluxes(i, j, s)
    use mod_isale
    use mod_advect

    implicit none
    
    integer :: m
    integer, intent(in) :: i, j
    real*8, intent(out) :: s(nmat,9)
    !DE when using parallel version, restore fluxes based on STORE_FLUX
    do m=1,nmat

       !              ----- Left side -----
       if (i.eq.1) then ! left boundary
          s(m,5)=min(0.D0,voltran_l(m,j)) ! outflux
          s(m,6)=max(0.D0,voltran_l(m,j)) ! influx
       else
          s(m,5)=min(0.D0,-voltran_r(m,i-1,j)) ! outflux
          s(m,6)=max(0.D0,-voltran_r(m,i-1,j)) ! influx
       endif
       !              ----- Bottom side -----
       if (j.eq.1) then ! bottom boundary
          s(m,7)=min(0.D0,voltran_b(m,i)) ! outflux
          s(m,8)=max(0.D0,voltran_b(m,i)) ! influx
       else
          s(m,7)=min(0.D0,-voltran_t(m,i,j-1)) ! outflux
          s(m,8)=max(0.D0,-voltran_t(m,i,j-1)) ! influx
       endif
       !              ----- Right side -----
       s(m,1)=min(0.D0,voltran_r(m,i,j))
       s(m,2)=max(0.D0,voltran_r(m,i,j))
       !              ----- Top side -----
       s(m,3)=min(0.D0,voltran_t(m,i,j))
       s(m,4)=max(0.D0,voltran_t(m,i,j))

    enddo

  end subroutine advect_restore_fluxes
  

  ! ------------------------------------------------------------------------------
  !> Subroutine to calculate momentum (cell centered) for cell i,j
  !! and material k
  !!
  !! @authors Kai Wuennemann, MfN; Dirk Elbeshausen, MfN; Xianzheng Li, ICL
  !! @GitHub: x1ng4me(XL)
  !!     Description                                     Programmer    Date
  !!     ------------------------------------------------------------------
  !!     Original version (5.7).............................KAI  2005/03/22
  !!     Conversion into Fortran90..........................DE   2007/02/24
  !!     Move subroutine outside advect.....................XL   2019/08/28
  !< ------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------  
  ! Move outside the subroutine advect for parallelization
  subroutine advect_momentum(i, j, sumcc, nom, nomic, vstat, vmom)
    use mod_isale
    use mod_advect
    use mod_identify_mat

    implicit none
    integer :: nvm
    integer :: l, m, sumcc(nmat)
    integer :: nom(nmat), nomic
    character*1 :: vstat(igs:ige,jgs:jge)
    ! character*1, intent(in) :: vstat(:,:)
    integer :: i, j
    real*8 :: vmom(X:Y,1,nmat,igs:ige,jgs:jge)

    if (cmc(TOTAL,i,j).eq.0.D0) then

       ! if the cell is empty, there are no cell momenta
       vmom(X:Y,1,firstmat:nmat,i,j) = 0.d0

    elseif (sumcc(VOID).eq.4) then

       ! If the cell has material, but all four vertices are flagged 
       ! as void use all vertex velocities to calculate cell momenta
       !do k=nom(1),nom(nomic)
       do l=1,nomic
          m=nom(l)
          vmom(X:Y,1,m,i,j)=0.25*rol(m,i,j)* &
               (VL(X:Y,i+1,j)+VL(X:Y,i+1,j+1)+VL(X:Y,i,j+1)+VL(X:Y,i,j))
       enddo

    else

       ! Otherwise, look at each vertex--if it is flagged as having 
       ! material, include that vertex in the calculation of the 
       ! cell-centered momenta
       do l=1,nomic
          m=nom(l)
          nvm=0
          vmom(X:Y,1,m,i,j)=0.D0
          if (vstat(i,j).eq.'M') then
             vmom(X:Y,1,m,i,j)=vmom(X:Y,1,m,i,j)+VL(X:Y,i,j)
             nvm=nvm+1
          endif
          if (vstat(i+1,j).eq.'M') then
             vmom(X:Y,1,m,i,j)=vmom(X:Y,1,m,i,j)+VL(X:Y,i+1,j)
             nvm=nvm+1
          endif
          if (vstat(i,j+1).eq.'M') then
             vmom(X:Y,1,m,i,j)=vmom(X:Y,1,m,i,j)+VL(X:Y,i,j+1)
             nvm=nvm+1
          endif
          if (vstat(i+1,j+1).eq.'M') then
             vmom(X:Y,1,m,i,j)=vmom(X:Y,1,m,i,j)+VL(X:Y,i+1,j+1)
             nvm=nvm+1
          endif
          vmom(X:Y,1,m,i,j)=rol(m,i,j)*vmom(X:Y,1,m,i,j)/dfloat(nvm)
       enddo

    endif

  end subroutine advect_momentum
  ! ------------------------------------------------------------------------------
  !> Subroutine to calculate momentum for cell i,j and material k for use in the
  !! Half-Index Shift algorithm, where each nodal momenta are advected
  !! separately. See Benson (1992, JCP, 100: 143-162) for details.
  !!
  !! @authors Gareth S. Collins, Imperial College, Xianzheng Li, ICL
  !! @Github: x1ng4me(XL)
  !!     Description                                     Programmer    Date
  !!     ------------------------------------------------------------------
  !!     Original version (5.7).............................GSC  2016/09/09
  !!     Move subroutine outside advect.....................XL   2019/08/28
  !< ------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------
  ! Move out from subroutine advect for parallelization
  subroutine advect_momentum_his(i, j, nom, nomic, vmom)
    use mod_isale
    use mod_advect
    use mod_identify_mat
 
    implicit none
    integer :: l, m
    ! character*1, intent(in) :: vstat(:,:)
    integer :: i, j
    integer :: nom(nmat),nomic
    real*8 :: vmom(X:Y,4,nmat,igs:ige,jgs:jge)
 
    ! Only cells with material have momenta to be advected
    vmom(:,:,:,i,j)=0.D0
    if (cmc(TOTAL,i,j).gt.0.D0) then
       do l=1,nomic
          m=nom(l)
          vmom(X:Y,1,m,i,j)=rol(m,i,j)*VL(X:Y,i+1,j)
          vmom(X:Y,3,m,i,j)=rol(m,i,j)*VL(X:Y,i,  j+1)
          vmom(X:Y,2,m,i,j)=rol(m,i,j)*VL(X:Y,i+1,j+1)
          vmom(X:Y,4,m,i,j)=rol(m,i,j)*VL(X:Y,i,  j)
       end do
    end if
 
  end subroutine advect_momentum_his

  ! ------------------------------------------------------------------------------
  !> Subroutine for performing advection of outfluxes for current cell
  !!
  !! @authors Kai Wuennemann, MfN; Dirk Elbeshausen, MfN; Xianzheng Li, ICL
  !! @GitHub: x1ng4me(XL)
  !!     Description                                     Programmer    Date
  !!     ------------------------------------------------------------------
  !!     Original version (5.7).............................KAI  2005/03/22
  !!     Conversion into Fortran90..........................DE   2007/02/24
  !!     Move subroutine outside advect.....................XL   2019/08/28
  !< ------------------------------------------------------------------------------
  ! Move out from subroutine advect for parallelization
  subroutine advect_outflux(nomic, nom, vmom_dim, vmom, siep, alphap, masse, fluxpart, volstrp, i, j, sumout, s, tvmat, vmat, velp, plwp, dummyfp, stresdevp, damagep, epstrp, initcoordp, vmomp)
    use mod_isale
    use mod_advect
    use mod_identify_mat

    implicit none
    integer :: m
    integer :: nomic, nom(nmat), vmom_dim
    real*8 :: vmom(X:Y,1:vmom_dim,nmat,igs:ige,jgs:jge)
    ! character*1, intent(in) :: vstat(:,:)
    real*8 :: siep(nmat,igs:ige,jgs:jge), alphap(nmat,igs:ige,jgs:jge)
    real*8 :: masse(nmat)
    real*8 :: fluxpart
    real*8 :: volstrp(igs:ige,jgs:jge)
    integer :: i, j
    real*8 :: sumout, s(nmat,9)
    real*8 :: tvmat, vmat(nmat)
    real*8 :: vmomp(X:Y,1:vmom_dim,igs:ige,jgs:jge)
    real*8 :: velp(igs:ige,jgs:jge), plwp(igs:ige,jgs:jge)
    real*8 :: dummyfp(1:ndummy,igs:ige,jgs:jge)
    real*8 :: stresdevp(1:4,igs:ige,jgs:jge)
    real*8 :: damagep(igs:ige,jgs:jge),epstrp(igs:ige,jgs:jge)
    real*8 :: initcoordp(X:Y,igs:ige,jgs:jge)
        
    ! The total volume and mass of material in cell after outflux
    tvmat   =tvmat+sumout
    masse(TOTAL)=max(0.D0,tvmat*rol(TOTAL,i,j))
    vmat(VOID) =vmat(VOID) +s(VOID,9)

    ! Volume, mass, energy and distension in cell after outflux
    do m=firstmat,nmat
       vmat(m) =vmat(m) +s(m,9)
       masse(m)=max(0.D0,vmat(m)*rol(m,i,j))
       siep(m,i,j)=masse(m)*sie(m,i,j)
       alphap(m,i,j) =masse(m)*alpha(m,i,j)
    enddo

    ! Other fields are advected by mass or volume
    select case (ADVTYPE)
    case (ADVTYPE_VOL)    ! according to volume flux
       fluxpart = tvmat
    case (ADVTYPE_MASS)   ! according to mass flux
       fluxpart = masse(TOTAL)
    case default
       write(IOERR,*) "ADVTYPE UNKNOWN INADVECT_ OUTFLUX : ",advtype
       call mystop("ERROR IN ADVECT_OUTFLUX: ADVTYPE UNKNOWN")
    end select

    ! outflux according to volume OR mass
    ! (fluxing amount is stored in fluxpart)
    volstrp(i,j)=fluxpart*volstrain(i,j)
    if(dummyfields) dummyfp(1:ndummy,i,j)=fluxpart*dummyf(1:ndummy,i,j)
    initcoordp(X:Y,i,j)=fluxpart*initcoord(X:Y,i,j)
    damagep(i,j)=fluxpart*damage(i,j)

    ! stresses after outflux
    if (.not. hydro) then
       stresdevp(1:4,i,j)=fluxpart*stresdev(1:4,i,j)
       velp(i,j)         =masse(TOTAL)*velo(i,j)*velo(i,j)
       epstrp(i,j)       =fluxpart*epstrain(i,j)
       plwp(i,j)         =masse(TOTAL)*plw(i,j)
    end if

    ! momentum after outflux (according to volume)
    vmomp(X:Y,:,i,j)=0.D0
    do m=firstmat,nmat
       vmomp(X:Y,:,i,j)=vmomp(X:Y,:,i,j)+s(m,9)*vmom(X:Y,:,m,i,j)
    end do

  end subroutine advect_outflux

  ! ------------------------------------------------------------------------------
  !> Subroutine for performing advection of fluxes from donor cell
  !! (idon,jdon) to current cell (i,j) through given face (side)...
  !!
  !! @authors Kai Wuennemann, MfN; Dirk Elbeshausen, MfN; Gareth Collins, ICL; Xianzheng Li, ICL
  !! @Github: x1ng4me(XL)
  !!     Description                                     Programmer    Date
  !!     ------------------------------------------------------------------
  !!     Original version (5.7).............................KAI  2005/03/22
  !!     Conversion into Fortran90..........................DE   2007/02/24
  !!     Modified to improve efficiency.....................GSC  2008/05/22
  !!     Move subroutine outside advect.....................XL   2019/08/28
  !< ------------------------------------------------------------------------------  
  !Move out from subroutine advect for parallelization
  subroutine advect_influx(side,idon,jdon,nomic, nom, vmom_dim, vmom, siep, alphap, masse, fluxpart, volstrp, i, j, sumout, s,&
          &tvmat, vmat, velp, plwp, dummyfp, stresdevp, damagep, epstrp, initcoordp, vmomp)
    use mod_isale
    use mod_advect
    use mod_identify_mat

    implicit none
    integer, intent(in) :: side !< Integer flag for cell face
    integer, intent(in) :: idon !< i-cell index for donor cell
    integer, intent(in) :: jdon !< j-cell index for donor cell
    ! Local memory
    real*8 :: tflux,volmatin

    integer :: m
    integer :: nomic, nom(nmat), vmom_dim
    real*8 :: vmom(X:Y,1:vmom_dim,nmat,igs:ige,jgs:jge)
    ! character*1, intent(in) :: vstat(:,:)
    real*8 :: siep(nmat,igs:ige,jgs:jge), alphap(nmat,igs:ige,jgs:jge)
    real*8 :: masse(nmat)
    real*8 :: fluxpart
    real*8 :: volstrp(igs:ige,jgs:jge)
    integer :: i, j
    real*8 :: sumout, s(nmat,9)
    real*8 :: tvmat, vmat(nmat)
    real*8 :: vmomp(X:Y,1:vmom_dim,igs:ige,jgs:jge)
    real*8 :: velp(igs:ige,jgs:jge), plwp(igs:ige,jgs:jge)
    real*8 :: dummyfp(1:ndummy,igs:ige,jgs:jge)
    real*8 :: stresdevp(1:4,igs:ige,jgs:jge)
    real*8 :: damagep(igs:ige,jgs:jge),epstrp(igs:ige,jgs:jge)
    real*8 :: initcoordp(X:Y,igs:ige,jgs:jge)
    ! Add influx of vacuum to void volume
    vmat(VOID) =vmat(VOID)+s(VOID,side)

    ! Account for influx by the mass of each material
    tflux = 0.D0        ! Total influxing mass (or volume)
    do m=firstmat,nmat

       ! flux according to mass
       volmatin = s(m,side)
       fluxpart = volmatin*rol(m,idon,jdon)
       if (s(m,side).lt.0.D0) then
          write(IOERR,*) 'S<0: ',i,j,volmatin
       endif

       ! Flux each material mass, energy, distension and momenta
       if (fluxpart.gt.0.D0) then
          masse(m)      = masse(m)      + fluxpart
          siep(m,i,j)   = siep(m,i,j)   + fluxpart*sie(m,idon,jdon)
          alphap(m,i,j) = alphap(m,i,j) + fluxpart*alpha(m,idon,jdon)

          ! momentum (advection according to volume)
          vmomp(X:Y,:,i,j) = vmomp(X:Y,:,i,j) + volmatin*vmom(X:Y,:,m,idon,jdon)

          ! Plastic work and velp are always advected by mass
          if (.not.hydro) then
             plwp(i,j) = plwp(i,j) + fluxpart*plw(idon,jdon)
             velp(i,j) = velp(i,j) + fluxpart*velo(idon,jdon)**2.D0
          end if
       end if

       ! Record material volume, and sum total volume
       vmat(m)  = vmat(m) + volmatin
       tvmat    = tvmat   + volmatin

       ! switch now to advection according to volume if desired
       if (advtype==ADVTYPE_VOL) fluxpart = volmatin
       ! Sum total mass (or volume) fluxes over materials
       tflux    = tflux   + fluxpart

    end do

    ! Flux the remaining variables by the total mass (or volume) flux
    if (tflux.gt.0.D0) then

       volstrp(i,j)  = volstrp(i,j)  + tflux*volstrain(idon,jdon)
       if (dummyfields) &
          dummyfp(1:ndummy,i,j)  = dummyfp(1:ndummy,i,j)  + tflux*dummyf(1:ndummy,idon,jdon)
       initcoordp(X:Y,i,j)=initcoordp(X:Y,i,j) +tflux*initcoord(X:Y,idon,jdon)
       damagep(i,j)  = damagep(i,j)  + tflux*damage(idon,jdon)

       ! shear stresses and strains
       if (.not. hydro) then
          stresdevp(1:4,i,j) = stresdevp(1:4,i,j) + tflux*stresdev(1:4,idon,jdon)
          epstrp(i,j)   = epstrp(i,j)   + tflux*epstrain(idon,jdon)
       end if

    end if

  end subroutine advect_influx
  ! ------------------------------------------------------------------------------
  !! Use for initialize mass and volume fractions
  !! needed to be written as a single subroutine 
  !! because of the order-specific problem in the last loop.
  !! @authors Xianzheng Li, ICL
  !! @Github: x1ng4me(XL)
  !!     Description                                     Programmer    Date
  !!     ------------------------------------------------------------------
  !!     Original version (5.7).............................XL      2019/08/28
  !< ------------------------------------------------------------------------------  
  subroutine cal_cmv(i,j,mvp)
    use mod_isale
    use mod_advect
    use mod_identify_mat

    implicit none
    integer, intent(in) :: i,j
    real*8 :: mvp(igs:ige,jgs:jge)
    !real*8 :: rmv(igs:ige,jgs:jge)
    !real*8 :: mv(igs:ige,jgs:jge)
    !real*8 :: V
    if (mvp(i,j).gt.0.D0) then
        rmv(i,j)=1.D0/mvp(i,j)
        V(X:Y,i,j)  =VL(X:Y,i,j)*mv(i,j)*rmv(i,j)
        mv(i,j)=mvp(i,j)
    else
        rmv(i,j)  = 0.D0
        V(:,i,j)  = 0.D0
        mv(i,j)   = 0.D0
        VL(:,i,j) = 0.D0
        cmv(:,i,j) = 0.D0
    endif
  end subroutine cal_cmv
! --------------------------------------------------------------------------------
!> This function determines the priority material number
!! based on the number of materials in the cell (nomic), 
!! the materials in the cell nom(i), and the material in the cell
!! that occupies the largest volume fraction (mmax).
!!
!! ibflag is an integer that identifies the state of any internal
!! boundary in the cell:
!! ibflag = 0 => Internal boundary found
!!        = 1 => No internal boundary found (isolated material)
!!        = 2 => No internal boundary found (split material)
!!
!! Different prioritisation strategies are possible, chosen by
!! priority_option
!!
!! @authors Gareth Collins, ICL
!!
!!     Description                                     Programmer    Date
!!     ------------------------------------------------------------------
!!     Original version ..................................GSC  2009/01/15
!!
!< --------------------------------------------------------------------------------
integer function priority_material(nomic,nom,mmax,ibflag)
  use mod_isale
  implicit none
  integer, intent(in) :: nomic,nom(nmat),ibflag(nmat),mmax
  integer, parameter :: INTERFACE_EXISTS=0
  integer :: l,m

  ! Default is no prioritisation
  priority_material = 0

  ! No prioritisation; old advection scheme or one-material cell
  if (priority_option == 0 .or. nomic == 1 ) return

  ! Prioritise only VOID and do not use interface construction for other materials.
  if (priority_option == -1) then
     priority_material = -1
     return
  end if
  
  ! The priority-material schemes
  if (priority_option == 1) then

     ! Strategy (1): Prioritise maximum material, apart from cases:
     !               (a) two-material incl. vacuum
     !               (b) no interface for maximum material and interface
     !                   for one of other materials in cell is defined
     if (nomic.eq.2.and.nom(1).eq.VOID) then ! Two materials, one==vacuum

        ! Simple mat+vac case, prioritise material
        priority_material = nom(2)
        
     elseif (ibflag(mmax)==INTERFACE_EXISTS) then 
        
        ! If maximum material has an interface, use this
        priority_material = mmax
        
     else
        
        ! If all else fails, use the "maximum material"
        priority_material=mmax
        
        ! If vacuum is present and interface was found set this as priority
        if (nom(1).eq.VOID.and.ibflag(VOID).eq.INTERFACE_EXISTS) then 

           priority_material=VOID
           
           ! Find any material where the interface is found
        else
           do l=1,nomic
              m=nom(l)
              if (ibflag(m).eq.INTERFACE_EXISTS) then 
                 priority_material=m
                 return
              end if
           end do
        end if
        
     end if

  elseif (priority_option == 2) then

     ! Strategy (2): Always prioritise maximum material
     priority_material = mmax

  elseif (priority_option == 3) then

     ! Strategy (3): Prioritise maxium material, apart from case:
     !               (a) two-material incl. vacuum
     if (nomic.eq.2.and.nom(1).eq.VOID) then ! Two materials, one==vacuum
        ! Simple mat+vac case, prioritise material
        priority_material = nom(2)
     else
        priority_material = mmax
     end if

  elseif (priority_option == 4) then

     ! Strategy (4): Prioritise upper target material
     if (nomic.eq.2.and.nom(1).eq.VOID) then
        ! Simple mat+vac case, prioritise material
        priority_material = nom(2)
     else
        ! Prioritise target surface material, if in cell
        do l=1,nomic
           if (nom(l)==obj(obj_num+lay_num)%mat) priority_material = nom(l)
        end do
        ! Otherwise, go with maximum
        if (priority_material == 0) priority_material = mmax
     end if

  elseif (priority_option == 5) then

     ! Strategy (5): Prioritise upper target material (alt)
     if (nomic.eq.2.and.nom(1).eq.VOID) then
        ! Simple mat+vac case, prioritise maximum
        priority_material = mmax
     else
        ! Prioritise target surface material, if in cell
        do l=1,nomic
           if (nom(l)==obj(obj_num+lay_num)%mat) priority_material = nom(l)
        end do
        ! Otherwise, go with maximum
        if (priority_material == 0) priority_material = mmax
     end if

  elseif (priority_option == 6) then

     ! Strategy (6): Prioritise projectile material
     if (nomic.eq.2.and.nom(1).eq.VOID) then
        ! Simple mat+vac case, prioritise material
        priority_material = nom(2)
     else
        ! Prioritise target surface material, if in cell
        do l=1,nomic
           if (nom(l)==obj(1)%mat) priority_material = nom(l)
        end do
        ! Otherwise, go with maximum
        if (priority_material == 0) priority_material = mmax
     end if

  elseif (priority_option == 7) then

     ! Strategy (7): Prioritise projectile material (alt)
     if (nomic.eq.2.and.nom(1).eq.VOID) then
        ! Simple mat+vac case, prioritise maximum
        priority_material = mmax
     else
        ! Prioritise target surface material, if in cell
        do l=1,nomic
           if (nom(l)==obj(1)%mat) priority_material = nom(l)
        end do
        ! Otherwise, go with maximum
        if (priority_material == 0) priority_material = mmax
     end if

  else

     write(IOERR,*) 'Priority option not recognised: ',priority_option
     call mystop("ERROR IN PRIORITY_MATERIAL: Priority option not recognized")

  end if

end function priority_material
