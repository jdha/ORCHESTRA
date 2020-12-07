MODULE obs_oper
   !!======================================================================
   !!                       ***  MODULE obs_oper  ***
   !! Observation diagnostics: Observation operators for various observation
   !!                          types
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   obs_prof_opt :    Compute the model counterpart of profile data
   !!   obs_surf_opt :    Compute the model counterpart of surface data
   !!   obs_pro_sco_opt: Compute the model counterpart of temperature and 
   !!                    salinity observations from profiles in generalised 
   !!                    vertical coordinates 
   !!----------------------------------------------------------------------

   !! * Modules used
   USE par_kind, ONLY : &         ! Precision variables
      & wp
   USE in_out_manager             ! I/O manager
   USE obs_inter_sup              ! Interpolation support
   USE obs_inter_h2d, ONLY : &    ! Horizontal interpolation to the obs pt
      & obs_int_h2d, &
      & obs_int_h2d_init
   USE obs_inter_z1d, ONLY : &    ! Vertical interpolation to the obs pt
      & obs_int_z1d,    &
      & obs_int_z1d_spl
   USE obs_const,  ONLY :     &
      & obfillflt		  ! Fillvalue   
   USE dom_oce,       ONLY : &
      & glamt, glamu, glamv, &
      & gphit, gphiu, gphiv, & 
      & gdept_n, gdept_0 
   USE lib_mpp,       ONLY : &
      & ctl_warn, ctl_stop
   USE obs_grid,      ONLY : & 
      & obs_level_search     
   USE sbcdcy,        ONLY : &    ! For calculation of where it is night-time
      & sbc_dcy, nday_qsr

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   PUBLIC obs_prof_opt, &  ! Compute the model counterpart of profile obs
      &   obs_pro_sco_opt, &  ! Compute the model counterpart of profile observations 
      &   obs_surf_opt     ! Compute the model counterpart of surface obs

   INTEGER, PARAMETER, PUBLIC :: &
      & imaxavtypes = 20   ! Max number of daily avgd obs types

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_oper.F90 7646 2017-02-06 09:25:03Z timgraham $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE obs_prof_opt( prodatqc, kt, kpi, kpj, kpk,          &
      &                     kit000, kdaystp,                      &
      &                     pvar1, pvar2, pgdept, pmask1, pmask2, &
      &                     plam1, plam2, pphi1, pphi2,           &
      &                     k1dint, k2dint, kdailyavtypes )

      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_pro_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of profiles
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    First, a vertical profile of horizontally interpolated model
      !!    now values is computed at the obs (lon, lat) point.
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!
      !!    Next, the vertical profile is interpolated to the
      !!    data depth points. Two vertical interpolation schemes are
      !!    available:
      !!        - linear       (k1dint = 0)
      !!        - Cubic spline (k1dint = 1)
      !!
      !!    For the cubic spline the 2nd derivative of the interpolating 
      !!    polynomial is computed before entering the vertical interpolation 
      !!    routine.
      !!
      !!    If the logical is switched on, the model equivalent is
      !!    a daily mean model temperature field. So, we first compute
      !!    the mean, then interpolate only at the end of the day.
      !!
      !!    Note: in situ temperature observations must be converted
      !!    to potential temperature (the model variable) prior to
      !!    assimilation. 
      !!
      !! ** Action  :
      !!
      !! History :
      !!      ! 97-11 (A. Weaver, S. Ricci, N. Daget)
      !!      ! 06-03 (G. Smith) NEMOVAR migration
      !!      ! 06-10 (A. Weaver) Cleanup
      !!      ! 07-01 (K. Mogensen) Merge of temperature and salinity
      !!      ! 07-03 (K. Mogensen) General handling of profiles
      !!      ! 15-02 (M. Martin) Combined routine for all profile types
      !!-----------------------------------------------------------------------

      !! * Modules used
      USE obs_profiles_def ! Definition of storage space for profile obs.

      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_prof), INTENT(INOUT) :: &
         & prodatqc                  ! Subset of profile data passing QC
      INTEGER, INTENT(IN) :: kt      ! Time step
      INTEGER, INTENT(IN) :: kpi     ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kpk
      INTEGER, INTENT(IN) :: kit000  ! Number of the first time step
                                     !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: k1dint  ! Vertical interpolation type (see header)
      INTEGER, INTENT(IN) :: k2dint  ! Horizontal interpolation type (see header)
      INTEGER, INTENT(IN) :: kdaystp ! Number of time steps per day
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj,kpk) :: &
         & pvar1,    &               ! Model field 1
         & pvar2,    &               ! Model field 2
         & pmask1,   &               ! Land-sea mask 1
         & pmask2                    ! Land-sea mask 2
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj) :: &
         & plam1,    &               ! Model longitudes for variable 1
         & plam2,    &               ! Model longitudes for variable 2
         & pphi1,    &               ! Model latitudes for variable 1
         & pphi2                     ! Model latitudes for variable 2
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpk) :: &
         & pgdept                    ! Model array of depth levels
      INTEGER, DIMENSION(imaxavtypes), OPTIONAL :: &
         & kdailyavtypes             ! Types for daily averages

      !! * Local declarations
      INTEGER ::   ji
      INTEGER ::   jj
      INTEGER ::   jk
      INTEGER ::   jobs
      INTEGER ::   inrc
      INTEGER ::   ipro
      INTEGER ::   idayend
      INTEGER ::   ista
      INTEGER ::   iend
      INTEGER ::   iobs
      INTEGER, DIMENSION(imaxavtypes) :: &
         & idailyavtypes
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi1, &
         & igrdi2, &
         & igrdj1, &
         & igrdj2
      REAL(KIND=wp) :: zlam
      REAL(KIND=wp) :: zphi
      REAL(KIND=wp) :: zdaystp
      REAL(KIND=wp), DIMENSION(kpk) :: &
         & zobsmask1, &
         & zobsmask2, &
         & zobsk,    &
         & zobs2k
      REAL(KIND=wp), DIMENSION(2,2,kpk) :: &
         & zweig1, &
         & zweig2
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
         & zmask1, &
         & zmask2, &
         & zint1, &
         & zint2, &
         & zinm1, &
         & zinm2
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zglam1, &
         & zglam2, &
         & zgphi1, &
         & zgphi2
      LOGICAL :: ld_dailyav

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! Record and data counters
      inrc = kt - kit000 + 2
      ipro = prodatqc%npstp(inrc)

      ! Daily average types
      ld_dailyav = .FALSE.
      IF ( PRESENT(kdailyavtypes) ) THEN
         idailyavtypes(:) = kdailyavtypes(:)
         IF ( ANY (idailyavtypes(:) /= -1) ) ld_dailyav = .TRUE.
      ELSE
         idailyavtypes(:) = -1
      ENDIF

      ! Daily means are calculated for values over timesteps:
      !  [1 <= kt <= kdaystp], [kdaystp+1 <= kt <= 2*kdaystp], ...
      idayend = MOD( kt - kit000 + 1, kdaystp )

      IF ( ld_dailyav ) THEN

         ! Initialize daily mean for first timestep of the day
         IF ( idayend == 1 .OR. kt == 0 ) THEN
            DO jk = 1, jpk
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     prodatqc%vdmean(ji,jj,jk,1) = 0.0
                     prodatqc%vdmean(ji,jj,jk,2) = 0.0
                  END DO
               END DO
            END DO
         ENDIF

         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! Increment field 1 for computing daily mean
                  prodatqc%vdmean(ji,jj,jk,1) = prodatqc%vdmean(ji,jj,jk,1) &
                     &                        + pvar1(ji,jj,jk)
                  ! Increment field 2 for computing daily mean
                  prodatqc%vdmean(ji,jj,jk,2) = prodatqc%vdmean(ji,jj,jk,2) &
                     &                        + pvar2(ji,jj,jk)
               END DO
            END DO
         END DO

         ! Compute the daily mean at the end of day
         zdaystp = 1.0 / REAL( kdaystp )
         IF ( idayend == 0 ) THEN
            IF (lwp) WRITE(numout,*) 'Calculating prodatqc%vdmean on time-step: ',kt
            CALL FLUSH(numout)
            DO jk = 1, jpk
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     prodatqc%vdmean(ji,jj,jk,1) = prodatqc%vdmean(ji,jj,jk,1) &
                        &                        * zdaystp
                     prodatqc%vdmean(ji,jj,jk,2) = prodatqc%vdmean(ji,jj,jk,2) &
                        &                        * zdaystp
                  END DO
               END DO
            END DO
         ENDIF

      ENDIF

      ! Get the data for interpolation
      ALLOCATE( &
         & igrdi1(2,2,ipro),      &
         & igrdi2(2,2,ipro),      &
         & igrdj1(2,2,ipro),      &
         & igrdj2(2,2,ipro),      &
         & zglam1(2,2,ipro),      &
         & zglam2(2,2,ipro),      &
         & zgphi1(2,2,ipro),      &
         & zgphi2(2,2,ipro),      &
         & zmask1(2,2,kpk,ipro),  &
         & zmask2(2,2,kpk,ipro),  &
         & zint1(2,2,kpk,ipro),  &
         & zint2(2,2,kpk,ipro)   &
         & )

      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro
         iobs = jobs - prodatqc%nprofup
         igrdi1(1,1,iobs) = prodatqc%mi(jobs,1)-1
         igrdj1(1,1,iobs) = prodatqc%mj(jobs,1)-1
         igrdi1(1,2,iobs) = prodatqc%mi(jobs,1)-1
         igrdj1(1,2,iobs) = prodatqc%mj(jobs,1)
         igrdi1(2,1,iobs) = prodatqc%mi(jobs,1)
         igrdj1(2,1,iobs) = prodatqc%mj(jobs,1)-1
         igrdi1(2,2,iobs) = prodatqc%mi(jobs,1)
         igrdj1(2,2,iobs) = prodatqc%mj(jobs,1)
         igrdi2(1,1,iobs) = prodatqc%mi(jobs,2)-1
         igrdj2(1,1,iobs) = prodatqc%mj(jobs,2)-1
         igrdi2(1,2,iobs) = prodatqc%mi(jobs,2)-1
         igrdj2(1,2,iobs) = prodatqc%mj(jobs,2)
         igrdi2(2,1,iobs) = prodatqc%mi(jobs,2)
         igrdj2(2,1,iobs) = prodatqc%mj(jobs,2)-1
         igrdi2(2,2,iobs) = prodatqc%mi(jobs,2)
         igrdj2(2,2,iobs) = prodatqc%mj(jobs,2)
      END DO

      CALL obs_int_comm_2d( 2, 2, ipro, kpi, kpj, igrdi1, igrdj1, plam1, zglam1 )
      CALL obs_int_comm_2d( 2, 2, ipro, kpi, kpj, igrdi1, igrdj1, pphi1, zgphi1 )
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, pmask1, zmask1 )
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, pvar1,   zint1 )
      
      CALL obs_int_comm_2d( 2, 2, ipro, kpi, kpj, igrdi2, igrdj2, plam2, zglam2 )
      CALL obs_int_comm_2d( 2, 2, ipro, kpi, kpj, igrdi2, igrdj2, pphi2, zgphi2 )
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi2, igrdj2, pmask2, zmask2 )
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi2, igrdj2, pvar2,   zint2 )

      ! At the end of the day also get interpolated means
      IF ( ld_dailyav .AND. idayend == 0 ) THEN

         ALLOCATE( &
            & zinm1(2,2,kpk,ipro),  &
            & zinm2(2,2,kpk,ipro)   &
            & )

         CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, &
            &                  prodatqc%vdmean(:,:,:,1), zinm1 )
         CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi2, igrdj2, &
            &                  prodatqc%vdmean(:,:,:,2), zinm2 )

      ENDIF

      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro

         iobs = jobs - prodatqc%nprofup

         IF ( kt /= prodatqc%mstp(jobs) ) THEN

            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                    &
                  &            ' kt      = ', kt,                      &
                  &            ' mstp    = ', prodatqc%mstp(jobs), &
                  &            ' ntyp    = ', prodatqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_pro_opt', 'Inconsistent time' )
         ENDIF

         zlam = prodatqc%rlam(jobs)
         zphi = prodatqc%rphi(jobs)

         ! Horizontal weights and vertical mask

         IF ( prodatqc%npvend(jobs,1) > 0 ) THEN

            CALL obs_int_h2d_init( kpk, kpk, k2dint, zlam, zphi,     &
               &                   zglam1(:,:,iobs), zgphi1(:,:,iobs), &
               &                   zmask1(:,:,:,iobs), zweig1, zobsmask1 )

         ENDIF

         IF ( prodatqc%npvend(jobs,2) > 0 ) THEN

            CALL obs_int_h2d_init( kpk, kpk, k2dint, zlam, zphi,     &
               &                   zglam2(:,:,iobs), zgphi2(:,:,iobs), &
               &                   zmask2(:,:,:,iobs), zweig2, zobsmask2 )
 
         ENDIF

         IF ( prodatqc%npvend(jobs,1) > 0 ) THEN

            zobsk(:) = obfillflt

            IF ( ANY (idailyavtypes(:) == prodatqc%ntyp(jobs)) ) THEN

               IF ( idayend == 0 )  THEN
                  ! Daily averaged data
                  CALL obs_int_h2d( kpk, kpk,      &
                     &              zweig1, zinm1(:,:,:,iobs), zobsk )

               ENDIF

            ELSE 

               ! Point data
               CALL obs_int_h2d( kpk, kpk,      &
                  &              zweig1, zint1(:,:,:,iobs), zobsk )

            ENDIF

            !-------------------------------------------------------------
            ! Compute vertical second-derivative of the interpolating 
            ! polynomial at obs points
            !-------------------------------------------------------------

            IF ( k1dint == 1 ) THEN
               CALL obs_int_z1d_spl( kpk, zobsk, zobs2k,   &
                  &                  pgdept, zobsmask1 )
            ENDIF

            !-----------------------------------------------------------------
            !  Vertical interpolation to the observation point
            !-----------------------------------------------------------------
            ista = prodatqc%npvsta(jobs,1)
            iend = prodatqc%npvend(jobs,1)
            CALL obs_int_z1d( kpk,                &
               & prodatqc%var(1)%mvk(ista:iend),  &
               & k1dint, iend - ista + 1,         &
               & prodatqc%var(1)%vdep(ista:iend), &
               & zobsk, zobs2k,                   &
               & prodatqc%var(1)%vmod(ista:iend), &
               & pgdept, zobsmask1 )

         ENDIF

         IF ( prodatqc%npvend(jobs,2) > 0 ) THEN

            zobsk(:) = obfillflt

            IF ( ANY (idailyavtypes(:) == prodatqc%ntyp(jobs)) ) THEN

               IF ( idayend == 0 )  THEN

                  ! Daily averaged data
                  CALL obs_int_h2d( kpk, kpk,      &
                     &              zweig2, zinm2(:,:,:,iobs), zobsk )

               ENDIF

            ELSE

               ! Point data
               CALL obs_int_h2d( kpk, kpk,      &
                  &              zweig2, zint2(:,:,:,iobs), zobsk )

            ENDIF


            !-------------------------------------------------------------
            ! Compute vertical second-derivative of the interpolating 
            ! polynomial at obs points
            !-------------------------------------------------------------

            IF ( k1dint == 1 ) THEN
               CALL obs_int_z1d_spl( kpk, zobsk, zobs2k, &
                  &                  pgdept, zobsmask2 )
            ENDIF

            !----------------------------------------------------------------
            !  Vertical interpolation to the observation point
            !----------------------------------------------------------------
            ista = prodatqc%npvsta(jobs,2)
            iend = prodatqc%npvend(jobs,2)
            CALL obs_int_z1d( kpk, &
               & prodatqc%var(2)%mvk(ista:iend),&
               & k1dint, iend - ista + 1, &
               & prodatqc%var(2)%vdep(ista:iend),&
               & zobsk, zobs2k, &
               & prodatqc%var(2)%vmod(ista:iend),&
               & pgdept, zobsmask2 )

         ENDIF

      END DO

      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdi1, &
         & igrdi2, &
         & igrdj1, &
         & igrdj2, &
         & zglam1, &
         & zglam2, &
         & zgphi1, &
         & zgphi2, &
         & zmask1, &
         & zmask2, &
         & zint1,  &
         & zint2   &
         & )

      ! At the end of the day also get interpolated means
      IF ( ld_dailyav .AND. idayend == 0 ) THEN
         DEALLOCATE( &
            & zinm1,  &
            & zinm2   &
            & )
      ENDIF

      prodatqc%nprofup = prodatqc%nprofup + ipro 

   END SUBROUTINE obs_prof_opt

   SUBROUTINE obs_pro_sco_opt( prodatqc, kt, kpi, kpj, kpk, kit000, kdaystp, & 
      &                    ptn, psn, pgdept, pgdepw, ptmask, k1dint, k2dint, & 
      &                    kdailyavtypes ) 
      !!----------------------------------------------------------------------- 
      !! 
      !!                     ***  ROUTINE obs_pro_opt  *** 
      !! 
      !! ** Purpose : Compute the model counterpart of profiles 
      !!              data by interpolating from the model grid to the  
      !!              observation point. Generalised vertical coordinate version 
      !! 
      !! ** Method  : Linearly interpolate to each observation point using  
      !!              the model values at the corners of the surrounding grid box. 
      !! 
      !!          First, model values on the model grid are interpolated vertically to the 
      !!          Depths of the profile observations.  Two vertical interpolation schemes are 
      !!          available: 
      !!          - linear       (k1dint = 0) 
      !!          - Cubic spline (k1dint = 1)    
      !! 
      !! 
      !!         Secondly the interpolated values are interpolated horizontally to the  
      !!         obs (lon, lat) point. 
      !!         Several horizontal interpolation schemes are available: 
      !!        - distance-weighted (great circle) (k2dint = 0) 
      !!        - distance-weighted (small angle)  (k2dint = 1) 
      !!        - bilinear (geographical grid)     (k2dint = 2) 
      !!        - bilinear (quadrilateral grid)    (k2dint = 3) 
      !!        - polynomial (quadrilateral grid)  (k2dint = 4) 
      !! 
      !!    For the cubic spline the 2nd derivative of the interpolating  
      !!    polynomial is computed before entering the vertical interpolation  
      !!    routine. 
      !! 
      !!    For ENACT moored buoy data (e.g., TAO), the model equivalent is 
      !!    a daily mean model temperature field. So, we first compute 
      !!    the mean, then interpolate only at the end of the day. 
      !! 
      !!    This is the procedure to be used with generalised vertical model  
      !!    coordinates (ie s-coordinates. It is ~4x slower than the equivalent 
      !!    horizontal then vertical interpolation algorithm, but can deal with situations 
      !!    where the model levels are not flat. 
      !!    ONLY PERFORMED if ln_sco=.TRUE.  
      !!       
      !!    Note: the in situ temperature observations must be converted 
      !!    to potential temperature (the model variable) prior to 
      !!    assimilation.  
      !!?????????????????????????????????????????????????????????????? 
      !!    INCLUDE POTENTIAL TEMP -> IN SITU TEMP IN OBS OPERATOR??? 
      !!?????????????????????????????????????????????????????????????? 
      !! 
      !! ** Action  : 
      !! 
      !! History : 
      !!      ! 2014-08 (J. While) Adapted from obs_pro_opt to handel generalised 
      !!                           vertical coordinates
      !!----------------------------------------------------------------------- 
   
      !! * Modules used 
      USE obs_profiles_def   ! Definition of storage space for profile obs. 
       
      IMPLICIT NONE 
 
      !! * Arguments 
      TYPE(obs_prof), INTENT(INOUT) :: prodatqc   ! Subset of profile data not failing screening 
      INTEGER, INTENT(IN) :: kt        ! Time step 
      INTEGER, INTENT(IN) :: kpi       ! Model grid parameters 
      INTEGER, INTENT(IN) :: kpj 
      INTEGER, INTENT(IN) :: kpk 
      INTEGER, INTENT(IN) :: kit000    ! Number of the first time step  
                                       !   (kit000-1 = restart time) 
      INTEGER, INTENT(IN) :: k1dint    ! Vertical interpolation type (see header) 
      INTEGER, INTENT(IN) :: k2dint    ! Horizontal interpolation type (see header) 
      INTEGER, INTENT(IN) :: kdaystp   ! Number of time steps per day                     
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj,kpk) :: & 
         & ptn,    &    ! Model temperature field 
         & psn,    &    ! Model salinity field 
         & ptmask       ! Land-sea mask 
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj,kpk) :: & 
         & pgdept,  &    ! Model array of depth T levels    
         & pgdepw       ! Model array of depth W levels 
      INTEGER, DIMENSION(imaxavtypes), OPTIONAL :: & 
         & kdailyavtypes   ! Types for daily averages 
      
      !! * Local declarations 
      INTEGER ::   ji 
      INTEGER ::   jj 
      INTEGER ::   jk 
      INTEGER ::   iico, ijco 
      INTEGER ::   jobs 
      INTEGER ::   inrc 
      INTEGER ::   ipro 
      INTEGER ::   idayend 
      INTEGER ::   ista 
      INTEGER ::   iend 
      INTEGER ::   iobs 
      INTEGER ::   iin, ijn, ikn, ik   ! looping indices over interpolation nodes 
      INTEGER, DIMENSION(imaxavtypes) :: & 
         & idailyavtypes 
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: & 
         & igrdi, & 
         & igrdj 
      INTEGER :: & 
         & inum_obs
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iv_indic    
      REAL(KIND=wp) :: zlam 
      REAL(KIND=wp) :: zphi 
      REAL(KIND=wp) :: zdaystp 
      REAL(KIND=wp), DIMENSION(kpk) :: & 
         & zobsmask, & 
         & zobsk,    & 
         & zobs2k 
      REAL(KIND=wp), DIMENSION(2,2,1) :: & 
         & zweig, & 
         & l_zweig 
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: & 
         & zmask, & 
         & zintt, & 
         & zints, & 
         & zinmt, & 
         & zgdept,& 
         & zgdepw,& 
         & zinms 
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: & 
         & zglam, & 
         & zgphi    
      REAL(KIND=wp), DIMENSION(1) :: zmsk_1       
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: interp_corner       
 
      !------------------------------------------------------------------------ 
      ! Local initialization  
      !------------------------------------------------------------------------ 
      ! ... Record and data counters 
      inrc = kt - kit000 + 2 
      ipro = prodatqc%npstp(inrc) 
  
      ! Daily average types 
      IF ( PRESENT(kdailyavtypes) ) THEN 
         idailyavtypes(:) = kdailyavtypes(:) 
      ELSE 
         idailyavtypes(:) = -1 
      ENDIF 
 
      ! Initialize daily mean for first time-step 
      idayend = MOD( kt - kit000 + 1, kdaystp ) 
 
      ! Added kt == 0 test to catch restart case  
      IF ( idayend == 1 .OR. kt == 0) THEN 
          
         IF (lwp) WRITE(numout,*) 'Reset prodatqc%vdmean on time-step: ',kt 
         DO jk = 1, jpk 
            DO jj = 1, jpj 
               DO ji = 1, jpi 
                  prodatqc%vdmean(ji,jj,jk,1) = 0.0 
                  prodatqc%vdmean(ji,jj,jk,2) = 0.0 
               END DO 
            END DO 
         END DO 
       
      ENDIF 
       
      DO jk = 1, jpk 
         DO jj = 1, jpj 
            DO ji = 1, jpi 
               ! Increment the temperature field for computing daily mean 
               prodatqc%vdmean(ji,jj,jk,1) = prodatqc%vdmean(ji,jj,jk,1) & 
               &                        + ptn(ji,jj,jk) 
               ! Increment the salinity field for computing daily mean 
               prodatqc%vdmean(ji,jj,jk,2) = prodatqc%vdmean(ji,jj,jk,2) & 
               &                        + psn(ji,jj,jk) 
            END DO 
         END DO 
      END DO 
    
      ! Compute the daily mean at the end of day 
      zdaystp = 1.0 / REAL( kdaystp ) 
      IF ( idayend == 0 ) THEN 
         DO jk = 1, jpk 
            DO jj = 1, jpj 
               DO ji = 1, jpi 
                  prodatqc%vdmean(ji,jj,jk,1) = prodatqc%vdmean(ji,jj,jk,1) & 
                  &                        * zdaystp 
                  prodatqc%vdmean(ji,jj,jk,2) = prodatqc%vdmean(ji,jj,jk,2) & 
                  &                           * zdaystp 
               END DO 
            END DO 
         END DO 
      ENDIF 
 
      ! Get the data for interpolation 
      ALLOCATE( & 
         & igrdi(2,2,ipro),      & 
         & igrdj(2,2,ipro),      & 
         & zglam(2,2,ipro),      & 
         & zgphi(2,2,ipro),      & 
         & zmask(2,2,kpk,ipro),  & 
         & zintt(2,2,kpk,ipro),  & 
         & zints(2,2,kpk,ipro),  & 
         & zgdept(2,2,kpk,ipro), & 
         & zgdepw(2,2,kpk,ipro)  & 
         & ) 
 
      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro 
         iobs = jobs - prodatqc%nprofup 
         igrdi(1,1,iobs) = prodatqc%mi(jobs,1)-1 
         igrdj(1,1,iobs) = prodatqc%mj(jobs,1)-1 
         igrdi(1,2,iobs) = prodatqc%mi(jobs,1)-1 
         igrdj(1,2,iobs) = prodatqc%mj(jobs,1) 
         igrdi(2,1,iobs) = prodatqc%mi(jobs,1) 
         igrdj(2,1,iobs) = prodatqc%mj(jobs,1)-1 
         igrdi(2,2,iobs) = prodatqc%mi(jobs,1) 
         igrdj(2,2,iobs) = prodatqc%mj(jobs,1) 
      END DO 
     
      ! Initialise depth arrays
      zgdept = 0.0
      zgdepw = 0.0
 
      CALL obs_int_comm_2d( 2, 2, ipro, kpi, kpj, igrdi, igrdj, glamt, zglam ) 
      CALL obs_int_comm_2d( 2, 2, ipro, kpi, kpj, igrdi, igrdj, gphit, zgphi ) 
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi, igrdj, ptmask,zmask ) 
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi, igrdj, ptn,   zintt ) 
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi, igrdj, psn,   zints ) 
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi, igrdj, pgdept(:,:,:), & 
        &                     zgdept ) 
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi, igrdj, pgdepw(:,:,:), & 
        &                     zgdepw ) 
 
      ! At the end of the day also get interpolated means 
      IF ( idayend == 0 ) THEN 
 
         ALLOCATE( & 
            & zinmt(2,2,kpk,ipro),  & 
            & zinms(2,2,kpk,ipro)   & 
            & ) 
 
         CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi, igrdj, & 
            &                  prodatqc%vdmean(:,:,:,1), zinmt ) 
         CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi, igrdj, & 
            &                  prodatqc%vdmean(:,:,:,2), zinms ) 
 
      ENDIF 
       
      ! Return if no observations to process 
      ! Has to be done after comm commands to ensure processors 
      ! stay in sync 
      IF ( ipro == 0 ) RETURN 
 
      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro 
    
         iobs = jobs - prodatqc%nprofup 
    
         IF ( kt /= prodatqc%mstp(jobs) ) THEN 
             
            IF(lwp) THEN 
               WRITE(numout,*) 
               WRITE(numout,*) ' E R R O R : Observation',              & 
                  &            ' time step is not consistent with the', & 
                  &            ' model time step' 
               WRITE(numout,*) ' =========' 
               WRITE(numout,*) 
               WRITE(numout,*) ' Record  = ', jobs,                    & 
                  &            ' kt      = ', kt,                      & 
                  &            ' mstp    = ', prodatqc%mstp(jobs), & 
                  &            ' ntyp    = ', prodatqc%ntyp(jobs) 
            ENDIF 
            CALL ctl_stop( 'obs_pro_opt', 'Inconsistent time' ) 
         ENDIF 
          
         zlam = prodatqc%rlam(jobs) 
         zphi = prodatqc%rphi(jobs) 
          
         ! Horizontal weights 
         ! Only calculated once, for both T and S. 
         ! Masked values are calculated later.  
 
         IF ( ( prodatqc%npvend(jobs,1) > 0 ) .OR. & 
            & ( prodatqc%npvend(jobs,2) > 0 ) ) THEN 
 
            CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,     & 
               &                   zglam(:,:,iobs), zgphi(:,:,iobs), & 
               &                   zmask(:,:,1,iobs), zweig, zmsk_1 ) 
 
         ENDIF 
         
         ! IF zmsk_1 = 0; then ob is on land 
         IF (zmsk_1(1) < 0.1) THEN 
            WRITE(numout,*) 'WARNING (obs_oper) :- profile found within landmask' 
   
         ELSE  
             
            ! Temperature 
             
            IF ( prodatqc%npvend(jobs,1) > 0 ) THEN  
    
               zobsk(:) = obfillflt 
    
               IF ( ANY (idailyavtypes(:) == prodatqc%ntyp(jobs)) ) THEN 
    
                  IF ( idayend == 0 )  THEN 
                   
                     ! Daily averaged moored buoy (MRB) data 
                   
                     ! vertically interpolate all 4 corners 
                     ista = prodatqc%npvsta(jobs,1) 
                     iend = prodatqc%npvend(jobs,1) 
                     inum_obs = iend - ista + 1 
                     ALLOCATE(interp_corner(2,2,inum_obs),iv_indic(inum_obs)) 
      
                     DO iin=1,2 
                        DO ijn=1,2 
                                       
                                       
           
                           IF ( k1dint == 1 ) THEN 
                              CALL obs_int_z1d_spl( kpk, & 
                                 &     zinmt(iin,ijn,:,iobs), & 
                                 &     zobs2k, zgdept(iin,ijn,:,iobs), & 
                                 &     zmask(iin,ijn,:,iobs)) 
                           ENDIF 
       
                           CALL obs_level_search(kpk, & 
                              &    zgdept(iin,ijn,:,iobs), & 
                              &    inum_obs, prodatqc%var(1)%vdep(ista:iend), & 
                              &    iv_indic) 
                           CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs, & 
                              &    prodatqc%var(1)%vdep(ista:iend), & 
                              &    zinmt(iin,ijn,:,iobs), & 
                              &    zobs2k, interp_corner(iin,ijn,:), & 
                              &    zgdept(iin,ijn,:,iobs), & 
                              &    zmask(iin,ijn,:,iobs)) 
       
                        ENDDO 
                     ENDDO 
                   
                   
                  ELSE 
                
                     CALL ctl_stop( ' A nonzero' //     & 
                        &           ' number of profile T BUOY data should' // & 
                        &           ' only occur at the end of a given day' ) 
    
                  ENDIF 
         
               ELSE  
                
                  ! Point data 
     
                  ! vertically interpolate all 4 corners 
                  ista = prodatqc%npvsta(jobs,1) 
                  iend = prodatqc%npvend(jobs,1) 
                  inum_obs = iend - ista + 1 
                  ALLOCATE(interp_corner(2,2,inum_obs), iv_indic(inum_obs)) 
                  DO iin=1,2  
                     DO ijn=1,2 
                                    
                                    
                        IF ( k1dint == 1 ) THEN 
                           CALL obs_int_z1d_spl( kpk, & 
                              &    zintt(iin,ijn,:,iobs),& 
                              &    zobs2k, zgdept(iin,ijn,:,iobs), & 
                              &    zmask(iin,ijn,:,iobs)) 
  
                        ENDIF 
       
                        CALL obs_level_search(kpk, & 
                            &        zgdept(iin,ijn,:,iobs),& 
                            &        inum_obs, prodatqc%var(1)%vdep(ista:iend), & 
                            &         iv_indic) 
                        CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs,     & 
                            &          prodatqc%var(1)%vdep(ista:iend),     & 
                            &          zintt(iin,ijn,:,iobs),            & 
                            &          zobs2k,interp_corner(iin,ijn,:), & 
                            &          zgdept(iin,ijn,:,iobs),         & 
                            &          zmask(iin,ijn,:,iobs) )      
         
                     ENDDO 
                  ENDDO 
             
               ENDIF 
       
               !------------------------------------------------------------- 
               ! Compute the horizontal interpolation for every profile level 
               !------------------------------------------------------------- 
             
               DO ikn=1,inum_obs 
                  iend=ista+ikn-1 

                  l_zweig(:,:,1) = 0._wp 

                  ! This code forces the horizontal weights to be  
                  ! zero IF the observation is below the bottom of the  
                  ! corners of the interpolation nodes, Or if it is in  
                  ! the mask. This is important for observations are near  
                  ! steep bathymetry 
                  DO iin=1,2 
                     DO ijn=1,2 
     
                        depth_loop1: DO ik=kpk,2,-1 
                           IF(zmask(iin,ijn,ik-1,iobs ) > 0.9 )THEN   
                            
                              l_zweig(iin,ijn,1) = &  
                                 & zweig(iin,ijn,1) * & 
                                 & MAX( SIGN(1._wp,(zgdepw(iin,ijn,ik,iobs) ) & 
                                 &  - prodatqc%var(1)%vdep(iend)),0._wp) 
                            
                              EXIT depth_loop1 
                           ENDIF 
                        ENDDO depth_loop1 
     
                     ENDDO 
                  ENDDO 
   
                  CALL obs_int_h2d( 1, 1, l_zweig, interp_corner(:,:,ikn), & 
                  &          prodatqc%var(1)%vmod(iend:iend) ) 
 
               ENDDO 
 
 
               DEALLOCATE(interp_corner,iv_indic) 
          
            ENDIF 
       
 
            ! Salinity  
          
            IF ( prodatqc%npvend(jobs,2) > 0 ) THEN  
    
               zobsk(:) = obfillflt 
    
               IF ( ANY (idailyavtypes(:) == prodatqc%ntyp(jobs)) ) THEN 
    
                  IF ( idayend == 0 )  THEN 
                   
                     ! Daily averaged moored buoy (MRB) data 
                   
                     ! vertically interpolate all 4 corners 
                     ista = prodatqc%npvsta(iobs,2) 
                     iend = prodatqc%npvend(iobs,2) 
                     inum_obs = iend - ista + 1 
                     ALLOCATE(interp_corner(2,2,inum_obs),iv_indic(inum_obs)) 
      
                     DO iin=1,2 
                        DO ijn=1,2 
                                       
                                       
           
                           IF ( k1dint == 1 ) THEN 
                              CALL obs_int_z1d_spl( kpk, & 
                                 &     zinms(iin,ijn,:,iobs), & 
                                 &     zobs2k, zgdept(iin,ijn,:,iobs), & 
                                 &     zmask(iin,ijn,:,iobs)) 
                           ENDIF 
       
                           CALL obs_level_search(kpk, & 
                              &    zgdept(iin,ijn,:,iobs), & 
                              &    inum_obs, prodatqc%var(2)%vdep(ista:iend), & 
                              &    iv_indic) 
                           CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs, & 
                              &    prodatqc%var(2)%vdep(ista:iend), & 
                              &    zinms(iin,ijn,:,iobs), & 
                              &    zobs2k, interp_corner(iin,ijn,:), & 
                              &    zgdept(iin,ijn,:,iobs), & 
                              &    zmask(iin,ijn,:,iobs)) 
       
                        ENDDO 
                     ENDDO 
                   
                   
                  ELSE 
                
                     CALL ctl_stop( ' A nonzero' //     & 
                        &           ' number of profile T BUOY data should' // & 
                        &           ' only occur at the end of a given day' ) 
    
                  ENDIF 
         
               ELSE  
                
                  ! Point data 
     
                  ! vertically interpolate all 4 corners 
                  ista = prodatqc%npvsta(jobs,2) 
                  iend = prodatqc%npvend(jobs,2) 
                  inum_obs = iend - ista + 1 
                  ALLOCATE(interp_corner(2,2,inum_obs), iv_indic(inum_obs)) 
                   
                  DO iin=1,2     
                     DO ijn=1,2  
                                 
                                 
                        IF ( k1dint == 1 ) THEN 
                           CALL obs_int_z1d_spl( kpk, & 
                              &    zints(iin,ijn,:,iobs),& 
                              &    zobs2k, zgdept(iin,ijn,:,iobs), & 
                              &    zmask(iin,ijn,:,iobs)) 
  
                        ENDIF 
       
                        CALL obs_level_search(kpk, & 
                           &        zgdept(iin,ijn,:,iobs),& 
                           &        inum_obs, prodatqc%var(2)%vdep(ista:iend), & 
                           &         iv_indic) 
                        CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs,  & 
                           &          prodatqc%var(2)%vdep(ista:iend),     & 
                           &          zints(iin,ijn,:,iobs),               & 
                           &          zobs2k,interp_corner(iin,ijn,:),     & 
                           &          zgdept(iin,ijn,:,iobs),              & 
                           &          zmask(iin,ijn,:,iobs) )      
         
                     ENDDO 
                  ENDDO 
             
               ENDIF 
       
               !------------------------------------------------------------- 
               ! Compute the horizontal interpolation for every profile level 
               !------------------------------------------------------------- 
             
               DO ikn=1,inum_obs 
                  iend=ista+ikn-1 

                  l_zweig(:,:,1) = 0._wp
   
                  ! This code forces the horizontal weights to be  
                  ! zero IF the observation is below the bottom of the  
                  ! corners of the interpolation nodes, Or if it is in  
                  ! the mask. This is important for observations are near  
                  ! steep bathymetry 
                  DO iin=1,2 
                     DO ijn=1,2 
     
                        depth_loop2: DO ik=kpk,2,-1 
                           IF(zmask(iin,ijn,ik-1,iobs ) > 0.9 )THEN   
                            
                              l_zweig(iin,ijn,1) = &  
                                 &  zweig(iin,ijn,1) * & 
                                 &  MAX( SIGN(1._wp,(zgdepw(iin,ijn,ik,iobs) ) & 
                                 &  - prodatqc%var(2)%vdep(iend)),0._wp) 
                            
                              EXIT depth_loop2 
                           ENDIF 
                        ENDDO depth_loop2 
     
                     ENDDO 
                  ENDDO 
   
                  CALL obs_int_h2d( 1, 1, l_zweig, interp_corner(:,:,ikn), & 
                  &          prodatqc%var(2)%vmod(iend:iend) ) 
 
               ENDDO 
 
 
               DEALLOCATE(interp_corner,iv_indic) 
          
            ENDIF 
          
         ENDIF 
       
      END DO 
     
      ! Deallocate the data for interpolation 
      DEALLOCATE( & 
         & igrdi, & 
         & igrdj, & 
         & zglam, & 
         & zgphi, & 
         & zmask, & 
         & zintt, & 
         & zints, & 
         & zgdept,&
         & zgdepw &
         & ) 
      ! At the end of the day also get interpolated means 
      IF ( idayend == 0 ) THEN 
         DEALLOCATE( & 
            & zinmt,  & 
            & zinms   & 
            & ) 
      ENDIF 
    
      prodatqc%nprofup = prodatqc%nprofup + ipro  
       
   END SUBROUTINE obs_pro_sco_opt 
 
   SUBROUTINE obs_surf_opt( surfdataqc, kt, kpi, kpj,         &
      &                    kit000, kdaystp, psurf, psurfmask, &
      &                    k2dint, ldnightav )

      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_surf_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of surface
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    The new model value is first computed at the obs (lon, lat) point.
      !!
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!
      !!
      !! ** Action  :
      !!
      !! History :
      !!      ! 07-03 (A. Weaver)
      !!      ! 15-02 (M. Martin) Combined routine for surface types
      !!-----------------------------------------------------------------------

      !! * Modules used
      USE obs_surf_def  ! Definition of storage space for surface observations

      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: &
         & surfdataqc                  ! Subset of surface data passing QC
      INTEGER, INTENT(IN) :: kt        ! Time step
      INTEGER, INTENT(IN) :: kpi       ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kit000    ! Number of the first time step 
                                       !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: kdaystp   ! Number of time steps per day
      INTEGER, INTENT(IN) :: k2dint    ! Horizontal interpolation type (see header)
      REAL(wp), INTENT(IN), DIMENSION(kpi,kpj) :: &
         & psurf,  &                   ! Model surface field
         & psurfmask                   ! Land-sea mask
      LOGICAL, INTENT(IN) :: ldnightav ! Logical for averaging night-time data

      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jobs
      INTEGER :: inrc
      INTEGER :: isurf
      INTEGER :: iobs
      INTEGER :: idayend
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi, &
         & igrdj
      INTEGER, DIMENSION(:,:), SAVE, ALLOCATABLE :: &
         & icount_night,      &
         & imask_night
      REAL(wp) :: zlam
      REAL(wp) :: zphi
      REAL(wp), DIMENSION(1) :: zext, zobsmask
      REAL(wp) :: zdaystp
      REAL(wp), DIMENSION(2,2,1) :: &
         & zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmask,  &
         & zsurf,  &
         & zsurfm, &
         & zglam,  &
         & zgphi
      REAL(wp), DIMENSION(:,:), SAVE, ALLOCATABLE :: &
         & zintmp,  &
         & zouttmp, &
         & zmeanday    ! to compute model sst in region of 24h daylight (pole)

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! Record and data counters
      inrc = kt - kit000 + 2
      isurf = surfdataqc%nsstp(inrc)

      IF ( ldnightav ) THEN

      ! Initialize array for night mean
         IF ( kt == 0 ) THEN
            ALLOCATE ( icount_night(kpi,kpj) )
            ALLOCATE ( imask_night(kpi,kpj) )
            ALLOCATE ( zintmp(kpi,kpj) )
            ALLOCATE ( zouttmp(kpi,kpj) )
            ALLOCATE ( zmeanday(kpi,kpj) )
            nday_qsr = -1   ! initialisation flag for nbc_dcy
         ENDIF

         ! Night-time means are calculated for night-time values over timesteps:
         !  [1 <= kt <= kdaystp], [kdaystp+1 <= kt <= 2*kdaystp], .....
         idayend = MOD( kt - kit000 + 1, kdaystp )

         ! Initialize night-time mean for first timestep of the day
         IF ( idayend == 1 .OR. kt == 0 ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  surfdataqc%vdmean(ji,jj) = 0.0
                  zmeanday(ji,jj) = 0.0
                  icount_night(ji,jj) = 0
               END DO
            END DO
         ENDIF

         zintmp(:,:) = 0.0
         zouttmp(:,:) = sbc_dcy( zintmp(:,:), .TRUE. )
         imask_night(:,:) = INT( zouttmp(:,:) )

         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Increment the temperature field for computing night mean and counter
               surfdataqc%vdmean(ji,jj) = surfdataqc%vdmean(ji,jj)  &
                      &                    + psurf(ji,jj) * REAL( imask_night(ji,jj) )
               zmeanday(ji,jj)          = zmeanday(ji,jj) + psurf(ji,jj)
               icount_night(ji,jj)      = icount_night(ji,jj) + imask_night(ji,jj)
            END DO
         END DO

         ! Compute the night-time mean at the end of the day
         zdaystp = 1.0 / REAL( kdaystp )
         IF ( idayend == 0 ) THEN
            IF (lwp) WRITE(numout,*) 'Calculating surfdataqc%vdmean on time-step: ',kt
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! Test if "no night" point
                  IF ( icount_night(ji,jj) > 0 ) THEN
                     surfdataqc%vdmean(ji,jj) = surfdataqc%vdmean(ji,jj) &
                       &                        / REAL( icount_night(ji,jj) )
                  ELSE
                     !At locations where there is no night (e.g. poles),
                     ! calculate daily mean instead of night-time mean.
                     surfdataqc%vdmean(ji,jj) = zmeanday(ji,jj) * zdaystp
                  ENDIF
               END DO
            END DO
         ENDIF

      ENDIF

      ! Get the data for interpolation

      ALLOCATE( &
         & igrdi(2,2,isurf), &
         & igrdj(2,2,isurf), &
         & zglam(2,2,isurf), &
         & zgphi(2,2,isurf), &
         & zmask(2,2,isurf), &
         & zsurf(2,2,isurf)  &
         & )

      DO jobs = surfdataqc%nsurfup + 1, surfdataqc%nsurfup + isurf
         iobs = jobs - surfdataqc%nsurfup
         igrdi(1,1,iobs) = surfdataqc%mi(jobs)-1
         igrdj(1,1,iobs) = surfdataqc%mj(jobs)-1
         igrdi(1,2,iobs) = surfdataqc%mi(jobs)-1
         igrdj(1,2,iobs) = surfdataqc%mj(jobs)
         igrdi(2,1,iobs) = surfdataqc%mi(jobs)
         igrdj(2,1,iobs) = surfdataqc%mj(jobs)-1
         igrdi(2,2,iobs) = surfdataqc%mi(jobs)
         igrdj(2,2,iobs) = surfdataqc%mj(jobs)
      END DO

      CALL obs_int_comm_2d( 2, 2, isurf, kpi, kpj, &
         &                  igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( 2, 2, isurf, kpi, kpj, &
         &                  igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_2d( 2, 2, isurf, kpi, kpj, &
         &                  igrdi, igrdj, psurfmask, zmask )
      CALL obs_int_comm_2d( 2, 2, isurf, kpi, kpj, &
         &                  igrdi, igrdj, psurf, zsurf )

      ! At the end of the day get interpolated means
      IF (ldnightav ) THEN
         IF ( idayend == 0 ) THEN

            ALLOCATE( &
               & zsurfm(2,2,isurf)  &
               & )

            CALL obs_int_comm_2d( 2, 2, isurf, kpi, kpj, igrdi, igrdj, &
               &               surfdataqc%vdmean(:,:), zsurfm )

         ENDIF
      ENDIF

      ! Loop over observations
      DO jobs = surfdataqc%nsurfup + 1, surfdataqc%nsurfup + isurf

         iobs = jobs - surfdataqc%nsurfup

         IF ( kt /= surfdataqc%mstp(jobs) ) THEN

            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                &
                  &            ' kt      = ', kt,                  &
                  &            ' mstp    = ', surfdataqc%mstp(jobs), &
                  &            ' ntyp    = ', surfdataqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_surf_opt', 'Inconsistent time' )

         ENDIF

         zlam = surfdataqc%rlam(jobs)
         zphi = surfdataqc%rphi(jobs)

         ! Get weights to interpolate the model value to the observation point
         CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
            &                   zglam(:,:,iobs), zgphi(:,:,iobs), &
            &                   zmask(:,:,iobs), zweig, zobsmask )

         ! Interpolate the model field to the observation point
         IF ( ldnightav .AND. idayend == 0 ) THEN
            ! Night-time averaged data
            CALL obs_int_h2d( 1, 1, zweig, zsurfm(:,:,iobs), zext )
         ELSE
            CALL obs_int_h2d( 1, 1, zweig, zsurf(:,:,iobs),  zext )
         ENDIF

         IF ( TRIM(surfdataqc%cvars(1)) == 'SLA' .AND. surfdataqc%nextra == 2 ) THEN
            ! ... Remove the MDT from the SSH at the observation point to get the SLA
            surfdataqc%rext(jobs,1) = zext(1)
            surfdataqc%rmod(jobs,1) = surfdataqc%rext(jobs,1) - surfdataqc%rext(jobs,2)
         ELSE
            surfdataqc%rmod(jobs,1) = zext(1)
         ENDIF

      END DO

      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdi, &
         & igrdj, &
         & zglam, &
         & zgphi, &
         & zmask, &
         & zsurf  &
         & )

      ! At the end of the day also deallocate night-time mean array
      IF ( ldnightav ) THEN
         IF ( idayend == 0 ) THEN
            DEALLOCATE( &
               & zsurfm  &
               & )
         ENDIF
      ENDIF

      surfdataqc%nsurfup = surfdataqc%nsurfup + isurf

   END SUBROUTINE obs_surf_opt

END MODULE obs_oper
