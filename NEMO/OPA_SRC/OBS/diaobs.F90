MODULE diaobs
   !!======================================================================
   !!                       ***  MODULE diaobs  ***
   !! Observation diagnostics: Computation of the misfit between data and
   !!                          their model equivalent 
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   dia_obs_init : Reading and prepare observations
   !!   dia_obs      : Compute model equivalent to observations
   !!   dia_obs_wri  : Write observational diagnostics
   !!   calc_date    : Compute the date of timestep in YYYYMMDD.HHMMSS format
   !!   ini_date     : Compute the initial date YYYYMMDD.HHMMSS
   !!   fin_date     : Compute the final date YYYYMMDD.HHMMSS
   !!----------------------------------------------------------------------
   !! * Modules used
   USE wrk_nemo                 ! Memory Allocation
   USE par_kind                 ! Precision variables
   USE in_out_manager           ! I/O manager
   USE par_oce
   USE dom_oce                  ! Ocean space and time domain variables
   USE obs_read_prof            ! Reading and allocation of profile obs
   USE obs_read_surf            ! Reading and allocation of surface obs
   USE obs_sstbias              ! Bias correction routine for SST 
   USE obs_readmdt              ! Reading and allocation of MDT for SLA.
   USE obs_prep                 ! Preparation of obs. (grid search etc).
   USE obs_oper                 ! Observation operators
   USE obs_write                ! Writing of observation related diagnostics
   USE obs_grid                 ! Grid searching
   USE obs_read_altbias         ! Bias treatment for altimeter
   USE obs_profiles_def         ! Profile data definitions
   USE obs_surf_def             ! Surface data definitions
   USE obs_types                ! Definitions for observation types
   USE mpp_map                  ! MPP mapping
   USE lib_mpp                  ! For ctl_warn/stop

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE
   PUBLIC dia_obs_init, &  ! Initialize and read observations
      &   dia_obs,      &  ! Compute model equivalent to observations
      &   dia_obs_wri,  &  ! Write model equivalent to observations
      &   dia_obs_dealloc, &  ! Deallocate dia_obs data
      &   calc_date           ! Compute the date of a timestep

   !! * Module variables
   LOGICAL, PUBLIC :: ln_diaobs   !: Logical switch for the obs operator
   LOGICAL :: ln_sstnight         !: Logical switch for night mean SST obs
   
   INTEGER :: nn_1dint       !: Vertical interpolation method
   INTEGER :: nn_2dint       !: Horizontal interpolation method
   INTEGER, DIMENSION(imaxavtypes) :: &
      & nn_profdavtypes      !: Profile data types representing a daily average
   INTEGER :: nproftypes     !: Number of profile obs types
   INTEGER :: nsurftypes     !: Number of surface obs types
   INTEGER, DIMENSION(:), ALLOCATABLE :: &
      & nvarsprof, &         !: Number of profile variables
      & nvarssurf            !: Number of surface variables
   INTEGER, DIMENSION(:), ALLOCATABLE :: &
      & nextrprof, &         !: Number of profile extra variables
      & nextrsurf            !: Number of surface extra variables
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: sstbias_type !SST bias type    
   TYPE(obs_surf), PUBLIC, POINTER, DIMENSION(:) :: &
      & surfdata, &          !: Initial surface data
      & surfdataqc           !: Surface data after quality control
   TYPE(obs_prof), PUBLIC, POINTER, DIMENSION(:) :: &
      & profdata, &          !: Initial profile data
      & profdataqc           !: Profile data after quality control

   CHARACTER(len=6), PUBLIC, DIMENSION(:), ALLOCATABLE :: &
      & cobstypesprof, &     !: Profile obs types
      & cobstypessurf        !: Surface obs types

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: diaobs.F90 6140 2015-12-21 11:35:23Z timgraham $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dia_obs_init
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs_init  ***
      !!          
      !! ** Purpose : Initialize and read observations
      !!
      !! ** Method  : Read the namelist and call reading routines
      !!
      !! ** Action  : Read the namelist and call reading routines
      !!
      !! History :
      !!        !  06-03  (K. Mogensen) Original code
      !!        !  06-05  (A. Weaver) Reformatted
      !!        !  06-10  (A. Weaver) Cleaning and add controls
      !!        !  07-03  (K. Mogensen) General handling of profiles
      !!        !  14-08  (J.While) Incorporated SST bias correction  
      !!        !  15-02  (M. Martin) Simplification of namelist and code
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      !! * Local declarations
      INTEGER, PARAMETER :: &
         & jpmaxnfiles = 1000    ! Maximum number of files for each obs type
      INTEGER, DIMENSION(:), ALLOCATABLE :: &
         & ifilesprof, &         ! Number of profile files
         & ifilessurf            ! Number of surface files
      INTEGER :: ios             ! Local integer output status for namelist read
      INTEGER :: jtype           ! Counter for obs types
      INTEGER :: jvar            ! Counter for variables
      INTEGER :: jfile           ! Counter for files

      CHARACTER(len=128), DIMENSION(jpmaxnfiles) :: &
         & cn_profbfiles, &      ! T/S profile input filenames
         & cn_sstfbfiles, &      ! Sea surface temperature input filenames
         & cn_slafbfiles, &      ! Sea level anomaly input filenames
         & cn_sicfbfiles, &      ! Seaice concentration input filenames
         & cn_velfbfiles, &      ! Velocity profile input filenames
         & cn_sstbias_files      ! SST bias input filenames
      CHARACTER(LEN=128) :: &
         & cn_altbiasfile        ! Altimeter bias input filename
      CHARACTER(len=128), DIMENSION(:,:), ALLOCATABLE :: &
         & clproffiles, &        ! Profile filenames
         & clsurffiles           ! Surface filenames

      LOGICAL :: ln_t3d          ! Logical switch for temperature profiles
      LOGICAL :: ln_s3d          ! Logical switch for salinity profiles
      LOGICAL :: ln_sla          ! Logical switch for sea level anomalies 
      LOGICAL :: ln_sst          ! Logical switch for sea surface temperature
      LOGICAL :: ln_sic          ! Logical switch for sea ice concentration
      LOGICAL :: ln_vel3d        ! Logical switch for velocity (u,v) obs
      LOGICAL :: ln_nea          ! Logical switch to remove obs near land
      LOGICAL :: ln_altbias      ! Logical switch for altimeter bias
      LOGICAL :: ln_sstbias     !: Logical switch for bias corection of SST 
      LOGICAL :: ln_ignmis       ! Logical switch for ignoring missing files
      LOGICAL :: ln_s_at_t       ! Logical switch to compute model S at T obs
      LOGICAL :: llvar1          ! Logical for profile variable 1
      LOGICAL :: llvar2          ! Logical for profile variable 1
      LOGICAL :: llnightav       ! Logical for calculating night-time averages
      LOGICAL, DIMENSION(jpmaxnfiles) :: lmask ! Used for finding number of sstbias files

      REAL(dp) :: rn_dobsini     ! Obs window start date YYYYMMDD.HHMMSS
      REAL(dp) :: rn_dobsend     ! Obs window end date   YYYYMMDD.HHMMSS
      REAL(wp), POINTER, DIMENSION(:,:) :: &
         & zglam1, &             ! Model longitudes for profile variable 1
         & zglam2                ! Model longitudes for profile variable 2
      REAL(wp), POINTER, DIMENSION(:,:) :: &
         & zgphi1, &             ! Model latitudes for profile variable 1
         & zgphi2                ! Model latitudes for profile variable 2
      REAL(wp), POINTER, DIMENSION(:,:,:) :: &
         & zmask1, &             ! Model land/sea mask associated with variable 1
         & zmask2                ! Model land/sea mask associated with variable 2

      NAMELIST/namobs/ln_diaobs, ln_t3d, ln_s3d, ln_sla,              &
         &            ln_sst, ln_sic, ln_vel3d,                       &
         &            ln_altbias, ln_nea, ln_grid_global,             &
         &            ln_grid_search_lookup,                          &
         &            ln_ignmis, ln_s_at_t, ln_sstnight,              &
         &            cn_profbfiles, cn_slafbfiles,                   &
         &            cn_sstfbfiles, cn_sicfbfiles,                   &
         &            cn_velfbfiles, cn_altbiasfile,                  &
         &            cn_gridsearchfile, rn_gridsearchres,            &
         &            rn_dobsini, rn_dobsend, nn_1dint, nn_2dint,     &
         &            nn_msshc, rn_mdtcorr, rn_mdtcutoff,             &
         &            nn_profdavtypes, ln_sstbias, cn_sstbias_files

      INTEGER :: jnumsstbias
      CALL wrk_alloc( jpi, jpj, zglam1 )
      CALL wrk_alloc( jpi, jpj, zglam2 )
      CALL wrk_alloc( jpi, jpj, zgphi1 )
      CALL wrk_alloc( jpi, jpj, zgphi2 )
      CALL wrk_alloc( jpi, jpj, jpk, zmask1 )
      CALL wrk_alloc( jpi, jpj, jpk, zmask2 )

      !-----------------------------------------------------------------------
      ! Read namelist parameters
      !-----------------------------------------------------------------------
      
      !Initalise all values in namelist arrays
      ALLOCATE(sstbias_type(jpmaxnfiles))
      ! Some namelist arrays need initialising
      cn_profbfiles(:) = ''
      cn_slafbfiles(:) = ''
      cn_sstfbfiles(:) = ''
      cn_sicfbfiles(:) = ''
      cn_velfbfiles(:) = ''
      cn_sstbias_files(:) = ''
      nn_profdavtypes(:) = -1

      CALL ini_date( rn_dobsini )
      CALL fin_date( rn_dobsend )

      ! Read namelist namobs : control observation diagnostics
      REWIND( numnam_ref )   ! Namelist namobs in reference namelist
      READ  ( numnam_ref, namobs, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namobs in reference namelist', lwp )

      REWIND( numnam_cfg )   ! Namelist namobs in configuration namelist
      READ  ( numnam_cfg, namobs, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namobs in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namobs )

      IF ( .NOT. ln_diaobs ) THEN
         IF(lwp) WRITE(numout,cform_war)
         IF(lwp) WRITE(numout,*)' ln_diaobs is set to false so not calling dia_obs'
         RETURN
      ENDIF
      
      !-----------------------------------------------------------------------
      ! Set up list of observation types to be used
      ! and the files associated with each type
      !-----------------------------------------------------------------------

      nproftypes = COUNT( (/ln_t3d .OR. ln_s3d, ln_vel3d /) )
      nsurftypes = COUNT( (/ln_sla, ln_sst, ln_sic /) )

      IF (ln_sstbias) THEN 
         lmask(:) = .FALSE. 
         WHERE (cn_sstbias_files(:) /= '') lmask(:) = .TRUE. 
         jnumsstbias = COUNT(lmask) 
         lmask(:) = .FALSE. 
      ENDIF      

      IF ( nproftypes == 0 .AND. nsurftypes == 0 ) THEN
         IF(lwp) WRITE(numout,cform_war)
         IF(lwp) WRITE(numout,*) ' ln_diaobs is set to true, but all obs operator logical flags', &
            &                    ' ln_t3d, ln_s3d, ln_sla, ln_sst, ln_sic, ln_vel3d', &
            &                    ' are set to .FALSE. so turning off calls to dia_obs'
         nwarn = nwarn + 1
         ln_diaobs = .FALSE.
         RETURN
      ENDIF

      IF ( nproftypes > 0 ) THEN

         ALLOCATE( cobstypesprof(nproftypes) )
         ALLOCATE( ifilesprof(nproftypes) )
         ALLOCATE( clproffiles(nproftypes,jpmaxnfiles) )

         jtype = 0
         IF (ln_t3d .OR. ln_s3d) THEN
            jtype = jtype + 1
            clproffiles(jtype,:) = cn_profbfiles(:)
            cobstypesprof(jtype) = 'prof  '
            ifilesprof(jtype) = 0
            DO jfile = 1, jpmaxnfiles
               IF ( trim(clproffiles(jtype,jfile)) /= '' ) &
                  ifilesprof(jtype) = ifilesprof(jtype) + 1
            END DO
         ENDIF
         IF (ln_vel3d) THEN
            jtype = jtype + 1
            clproffiles(jtype,:) = cn_velfbfiles(:)
            cobstypesprof(jtype) = 'vel   '
            ifilesprof(jtype) = 0
            DO jfile = 1, jpmaxnfiles
               IF ( trim(clproffiles(jtype,jfile)) /= '' ) &
                  ifilesprof(jtype) = ifilesprof(jtype) + 1
            END DO
         ENDIF

      ENDIF

      IF ( nsurftypes > 0 ) THEN

         ALLOCATE( cobstypessurf(nsurftypes) )
         ALLOCATE( ifilessurf(nsurftypes) )
         ALLOCATE( clsurffiles(nsurftypes, jpmaxnfiles) )

         jtype = 0
         IF (ln_sla) THEN
            jtype = jtype + 1
            clsurffiles(jtype,:) = cn_slafbfiles(:)
            cobstypessurf(jtype) = 'sla   '
            ifilessurf(jtype) = 0
            DO jfile = 1, jpmaxnfiles
               IF ( trim(clsurffiles(jtype,jfile)) /= '' ) &
                  ifilessurf(jtype) = ifilessurf(jtype) + 1
            END DO
         ENDIF
         IF (ln_sst) THEN
            jtype = jtype + 1
            clsurffiles(jtype,:) = cn_sstfbfiles(:)
            cobstypessurf(jtype) = 'sst   '
            ifilessurf(jtype) = 0
            DO jfile = 1, jpmaxnfiles
               IF ( trim(clsurffiles(jtype,jfile)) /= '' ) &
                  ifilessurf(jtype) = ifilessurf(jtype) + 1
            END DO
         ENDIF
#if defined key_lim2 || defined key_lim3
         IF (ln_sic) THEN
            jtype = jtype + 1
            clsurffiles(jtype,:) = cn_sicfbfiles(:)
            cobstypessurf(jtype) = 'sic   '
            ifilessurf(jtype) = 0
            DO jfile = 1, jpmaxnfiles
               IF ( trim(clsurffiles(jtype,jfile)) /= '' ) &
                  ifilessurf(jtype) = ifilessurf(jtype) + 1
            END DO
         ENDIF
#endif

      ENDIF

      !Write namelist settings to stdout
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_obs_init : Observation diagnostic initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '          Namelist namobs : set observation diagnostic parameters' 
         WRITE(numout,*) '             Logical switch for T profile observations                ln_t3d = ', ln_t3d
         WRITE(numout,*) '             Logical switch for S profile observations                ln_s3d = ', ln_s3d
         WRITE(numout,*) '             Logical switch for SLA observations                      ln_sla = ', ln_sla
         WRITE(numout,*) '             Logical switch for SST observations                      ln_sst = ', ln_sst
         WRITE(numout,*) '             Logical switch for Sea Ice observations                  ln_sic = ', ln_sic
         WRITE(numout,*) '             Logical switch for velocity observations               ln_vel3d = ', ln_vel3d
         WRITE(numout,*) '             Global distribution of observations              ln_grid_global = ',ln_grid_global
         WRITE(numout,*) '             Logical switch for SST bias correction         ln_sstbias = ', ln_sstbias 
         WRITE(numout,*) '             Logical switch for obs grid search lookup ln_grid_search_lookup = ',ln_grid_search_lookup
         IF (ln_grid_search_lookup) &
            WRITE(numout,*) '             Grid search lookup file header                cn_gridsearchfile = ', cn_gridsearchfile
         WRITE(numout,*) '             Initial date in window YYYYMMDD.HHMMSS               rn_dobsini = ', rn_dobsini
         WRITE(numout,*) '             Final date in window YYYYMMDD.HHMMSS                 rn_dobsend = ', rn_dobsend
         WRITE(numout,*) '             Type of vertical interpolation method                  nn_1dint = ', nn_1dint
         WRITE(numout,*) '             Type of horizontal interpolation method                nn_2dint = ', nn_2dint
         WRITE(numout,*) '             Rejection of observations near land switch               ln_nea = ', ln_nea
         WRITE(numout,*) '             MSSH correction scheme                                 nn_msshc = ', nn_msshc
         WRITE(numout,*) '             MDT  correction                                      rn_mdtcorr = ', rn_mdtcorr
         WRITE(numout,*) '             MDT cutoff for computed correction                 rn_mdtcutoff = ', rn_mdtcutoff
         WRITE(numout,*) '             Logical switch for alt bias                          ln_altbias = ', ln_altbias
         WRITE(numout,*) '             Logical switch for ignoring missing files             ln_ignmis = ', ln_ignmis
         WRITE(numout,*) '             Daily average types                             nn_profdavtypes = ', nn_profdavtypes
         WRITE(numout,*) '             Logical switch for night-time SST obs               ln_sstnight = ', ln_sstnight
         WRITE(numout,*) '          Number of profile obs types: ',nproftypes

         IF ( nproftypes > 0 ) THEN
            DO jtype = 1, nproftypes
               DO jfile = 1, ifilesprof(jtype)
                  WRITE(numout,'(1X,2A)') '             '//cobstypesprof(jtype)//' input observation file names  = ', &
                     TRIM(clproffiles(jtype,jfile))
               END DO
            END DO
         ENDIF

         WRITE(numout,*)'          Number of surface obs types: ',nsurftypes
         IF ( nsurftypes > 0 ) THEN
            DO jtype = 1, nsurftypes
               DO jfile = 1, ifilessurf(jtype)
                  WRITE(numout,'(1X,2A)') '             '//cobstypessurf(jtype)//' input observation file names  = ', &
                     TRIM(clsurffiles(jtype,jfile))
               END DO
            END DO
         ENDIF
         WRITE(numout,*) '~~~~~~~~~~~~'

      ENDIF

      !-----------------------------------------------------------------------
      ! Obs operator parameter checking and initialisations
      !-----------------------------------------------------------------------

      IF ( ln_vel3d .AND. ( .NOT. ln_grid_global ) ) THEN
         CALL ctl_stop( 'Velocity data only works with ln_grid_global=.true.' )
         RETURN
      ENDIF

      IF ( ln_grid_global ) THEN
         CALL ctl_warn( 'ln_grid_global=T may cause memory issues when used with a large number of processors' )
      ENDIF

      IF ( ( nn_1dint < 0 ) .OR. ( nn_1dint > 1 ) ) THEN
         CALL ctl_stop(' Choice of vertical (1D) interpolation method', &
            &                    ' is not available')
      ENDIF

      IF ( ( nn_2dint < 0 ) .OR. ( nn_2dint > 4 ) ) THEN
         CALL ctl_stop(' Choice of horizontal (2D) interpolation method', &
            &                    ' is not available')
      ENDIF

      CALL obs_typ_init
      IF(ln_grid_global) THEN
         CALL mppmap_init
      ENDIF

      CALL obs_grid_setup( )

      !-----------------------------------------------------------------------
      ! Depending on switches read the various observation types
      !-----------------------------------------------------------------------

      IF ( nproftypes > 0 ) THEN

         ALLOCATE(profdata(nproftypes))
         ALLOCATE(profdataqc(nproftypes))
         ALLOCATE(nvarsprof(nproftypes))
         ALLOCATE(nextrprof(nproftypes))

         DO jtype = 1, nproftypes

            nvarsprof(jtype) = 2
            IF ( TRIM(cobstypesprof(jtype)) == 'prof' ) THEN
               nextrprof(jtype) = 1
               llvar1 = ln_t3d
               llvar2 = ln_s3d
               zglam1 = glamt
               zgphi1 = gphit
               zmask1 = tmask
               zglam2 = glamt
               zgphi2 = gphit
               zmask2 = tmask
            ENDIF
            IF ( TRIM(cobstypesprof(jtype)) == 'vel' )  THEN
               nextrprof(jtype) = 2
               llvar1 = ln_vel3d
               llvar2 = ln_vel3d
               zglam1 = glamu
               zgphi1 = gphiu
               zmask1 = umask
               zglam2 = glamv
               zgphi2 = gphiv
               zmask2 = vmask
            ENDIF

            !Read in profile or profile obs types
            CALL obs_rea_prof( profdata(jtype), ifilesprof(jtype),       &
               &               clproffiles(jtype,1:ifilesprof(jtype)), &
               &               nvarsprof(jtype), nextrprof(jtype), nitend-nit000+2, &
               &               rn_dobsini, rn_dobsend, llvar1, llvar2, &
               &               ln_ignmis, ln_s_at_t, .FALSE., &
               &               kdailyavtypes = nn_profdavtypes )

            DO jvar = 1, nvarsprof(jtype)
               CALL obs_prof_staend( profdata(jtype), jvar )
            END DO

            CALL obs_pre_prof( profdata(jtype), profdataqc(jtype), &
               &               llvar1, llvar2, &
               &               jpi, jpj, jpk, &
               &               zmask1, zglam1, zgphi1, zmask2, zglam2, zgphi2,  &
               &               ln_nea, kdailyavtypes = nn_profdavtypes )

         END DO

         DEALLOCATE( ifilesprof, clproffiles )

      ENDIF

      IF ( nsurftypes > 0 ) THEN

         ALLOCATE(surfdata(nsurftypes))
         ALLOCATE(surfdataqc(nsurftypes))
         ALLOCATE(nvarssurf(nsurftypes))
         ALLOCATE(nextrsurf(nsurftypes))

         DO jtype = 1, nsurftypes

            nvarssurf(jtype) = 1
            nextrsurf(jtype) = 0
            llnightav = .FALSE.
            IF ( TRIM(cobstypessurf(jtype)) == 'sla' ) nextrsurf(jtype) = 2
            IF ( TRIM(cobstypessurf(jtype)) == 'sst' ) llnightav = ln_sstnight

            !Read in surface obs types
            CALL obs_rea_surf( surfdata(jtype), ifilessurf(jtype), &
               &               clsurffiles(jtype,1:ifilessurf(jtype)), &
               &               nvarssurf(jtype), nextrsurf(jtype), nitend-nit000+2, &
               &               rn_dobsini, rn_dobsend, ln_ignmis, .FALSE., llnightav )
         
         
            CALL obs_pre_surf( surfdata(jtype), surfdataqc(jtype), ln_nea )

            IF ( TRIM(cobstypessurf(jtype)) == 'sla' ) THEN
               CALL obs_rea_mdt( surfdataqc(jtype), nn_2dint )
               IF ( ln_altbias ) CALL obs_rea_altbias ( surfdataqc(jtype), nn_2dint, cn_altbiasfile )
            ENDIF
            IF ( TRIM(cobstypessurf(jtype)) == 'sst' .AND. ln_sstbias ) THEN
               !Read in bias field and correct SST.
               IF ( jnumsstbias == 0 ) CALL ctl_stop("ln_sstbias set,"// &
                                                     "  but no bias"// &
                                                     " files to read in")   
                  CALL obs_app_sstbias( surfdataqc(jtype), nn_2dint, &
                                        jnumsstbias, cn_sstbias_files(1:jnumsstbias) )
            ENDIF
         END DO

         DEALLOCATE( ifilessurf, clsurffiles )

      ENDIF

      CALL wrk_dealloc( jpi, jpj, zglam1 )
      CALL wrk_dealloc( jpi, jpj, zglam2 )
      CALL wrk_dealloc( jpi, jpj, zgphi1 )
      CALL wrk_dealloc( jpi, jpj, zgphi2 )
      CALL wrk_dealloc( jpi, jpj, jpk, zmask1 )
      CALL wrk_dealloc( jpi, jpj, jpk, zmask2 )

   END SUBROUTINE dia_obs_init

   SUBROUTINE dia_obs( kstp )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs  ***
      !!          
      !! ** Purpose : Call the observation operators on each time step
      !!
      !! ** Method  : Call the observation operators on each time step to
      !!              compute the model equivalent of the following data:
      !!               - Profile data, currently T/S or U/V
      !!               - Surface data, currently SST, SLA or sea-ice concentration.
      !!
      !! ** Action  :
      !!
      !! History :
      !!        !  06-03  (K. Mogensen) Original code
      !!        !  06-05  (K. Mogensen) Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  07-03  (K. Mogensen) General handling of profiles
      !!        !  07-04  (G. Smith) Generalized surface operators
      !!        !  08-10  (M. Valdivieso) obs operator for velocity profiles
      !!        !  14-08  (J. While) observation operator for profiles in 
      !!                             generalised vertical coordinates
      !!        !  15-08  (M. Martin) Combined surface/profile routines.
      !!----------------------------------------------------------------------
      !! * Modules used
      USE dom_oce, ONLY : &             ! Ocean space and time domain variables
         & gdept_n,       &      
         & gdept_1d      
      USE phycst, ONLY : &              ! Physical constants
         & rday                         
      USE oce, ONLY : &                 ! Ocean dynamics and tracers variables
         & tsn,  &             
         & un, vn, &
         & sshn  
      USE phycst, ONLY : &         ! Physical constants
         & rday
#if defined  key_lim3
      USE ice, ONLY : &            ! LIM3 Ice model variables
         & frld
#endif
#if defined key_lim2
      USE ice_2, ONLY : &          ! LIM2 Ice model variables
         & frld
#endif
      IMPLICIT NONE

      !! * Arguments
      INTEGER, INTENT(IN) :: kstp  ! Current timestep
      !! * Local declarations
      INTEGER :: idaystp           ! Number of timesteps per day
      INTEGER :: jtype             ! Data loop variable
      INTEGER :: jvar              ! Variable number
      INTEGER :: ji, jj            ! Loop counters
      REAL(wp), POINTER, DIMENSION(:,:,:) :: &
         & zprofvar1, &            ! Model values for 1st variable in a prof ob
         & zprofvar2               ! Model values for 2nd variable in a prof ob
      REAL(wp), POINTER, DIMENSION(:,:,:) :: &
         & zprofmask1, &           ! Mask associated with zprofvar1
         & zprofmask2              ! Mask associated with zprofvar2
      REAL(wp), POINTER, DIMENSION(:,:) :: &
         & zsurfvar                ! Model values equivalent to surface ob.
      REAL(wp), POINTER, DIMENSION(:,:) :: &
         & zglam1,    &            ! Model longitudes for prof variable 1
         & zglam2,    &            ! Model longitudes for prof variable 2
         & zgphi1,    &            ! Model latitudes for prof variable 1
         & zgphi2                  ! Model latitudes for prof variable 2
#if ! defined key_lim2 && ! defined key_lim3
      REAL(wp), POINTER, DIMENSION(:,:) :: frld
#endif
      LOGICAL :: llnightav        ! Logical for calculating night-time average

      !Allocate local work arrays
      CALL wrk_alloc( jpi, jpj, jpk, zprofvar1 )
      CALL wrk_alloc( jpi, jpj, jpk, zprofvar2 )
      CALL wrk_alloc( jpi, jpj, jpk, zprofmask1 )
      CALL wrk_alloc( jpi, jpj, jpk, zprofmask2 )
      CALL wrk_alloc( jpi, jpj, zsurfvar )
      CALL wrk_alloc( jpi, jpj, zglam1 )
      CALL wrk_alloc( jpi, jpj, zglam2 )
      CALL wrk_alloc( jpi, jpj, zgphi1 )
      CALL wrk_alloc( jpi, jpj, zgphi2 )
#if ! defined key_lim2 && ! defined key_lim3
      CALL wrk_alloc(jpi,jpj,frld) 
#endif

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_obs : Call the observation operators', kstp
         WRITE(numout,*) '~~~~~~~'
      ENDIF

      idaystp = NINT( rday / rdt )

      !-----------------------------------------------------------------------
      ! No LIM => frld == 0.0_wp
      !-----------------------------------------------------------------------
#if ! defined key_lim2 && ! defined key_lim3
      frld(:,:) = 0.0_wp
#endif
      !-----------------------------------------------------------------------
      ! Call the profile and surface observation operators
      !-----------------------------------------------------------------------

      IF ( nproftypes > 0 ) THEN

         DO jtype = 1, nproftypes

            SELECT CASE ( TRIM(cobstypesprof(jtype)) )
            CASE('prof')
               zprofvar1(:,:,:) = tsn(:,:,:,jp_tem)
               zprofvar2(:,:,:) = tsn(:,:,:,jp_sal)
               zprofmask1(:,:,:) = tmask(:,:,:)
               zprofmask2(:,:,:) = tmask(:,:,:)
               zglam1(:,:) = glamt(:,:)
               zglam2(:,:) = glamt(:,:)
               zgphi1(:,:) = gphit(:,:)
               zgphi2(:,:) = gphit(:,:)
            CASE('vel')
               zprofvar1(:,:,:) = un(:,:,:)
               zprofvar2(:,:,:) = vn(:,:,:)
               zprofmask1(:,:,:) = umask(:,:,:)
               zprofmask2(:,:,:) = vmask(:,:,:)
               zglam1(:,:) = glamu(:,:)
               zglam2(:,:) = glamv(:,:)
               zgphi1(:,:) = gphiu(:,:)
               zgphi2(:,:) = gphiv(:,:)
            END SELECT

            IF( ln_zco .OR. ln_zps ) THEN 
               CALL obs_prof_opt( profdataqc(jtype), kstp, jpi, jpj, jpk,  &
                  &               nit000, idaystp,                         &
                  &               zprofvar1, zprofvar2,                    &
                  &               gdept_1d, zprofmask1, zprofmask2,        &
                  &               zglam1, zglam2, zgphi1, zgphi2,          &
                  &               nn_1dint, nn_2dint,                      &
                  &               kdailyavtypes = nn_profdavtypes )
            ELSE IF(TRIM(cobstypesprof(jtype)) == 'prof') THEN
               !TG - THIS NEEDS MODIFICATION TO MATCH SIMPLIFICATION
               CALL obs_pro_sco_opt( profdataqc(jtype),                    & 
                  &              kstp, jpi, jpj, jpk, nit000, idaystp,   & 
                  &              zprofvar1, zprofvar2,                   & 
                  &              gdept_n(:,:,:), gdepw_n(:,:,:),           &
                  &              tmask, nn_1dint, nn_2dint,              & 
                  &              kdailyavtypes = nn_profdavtypes ) 
            ELSE
               CALL ctl_stop('DIA_OBS: Generalised vertical interpolation not'// &
                         'yet working for velocity data (turn off velocity observations')
            ENDIF

         END DO

      ENDIF

      IF ( nsurftypes > 0 ) THEN

         DO jtype = 1, nsurftypes

            SELECT CASE ( TRIM(cobstypessurf(jtype)) )
            CASE('sst')
               zsurfvar(:,:) = tsn(:,:,1,jp_tem)
               llnightav = ln_sstnight
            CASE('sla')
               zsurfvar(:,:) = sshn(:,:)
               llnightav = .FALSE.
#if defined key_lim2 || defined key_lim3
            CASE('sic')
               IF ( kstp == 0 ) THEN
                  IF ( lwp .AND. surfdataqc(jtype)%nsstpmpp(1) > 0 ) THEN
                     CALL ctl_warn( 'Sea-ice not initialised on zeroth '// &
                        &           'time-step but some obs are valid then.' )
                     WRITE(numout,*)surfdataqc(jtype)%nsstpmpp(1), &
                        &           ' sea-ice obs will be missed'
                  ENDIF
                  surfdataqc(jtype)%nsurfup = surfdataqc(jtype)%nsurfup + &
                     &                        surfdataqc(jtype)%nsstp(1)
                  CYCLE
               ELSE
                  zsurfvar(:,:) = 1._wp - frld(:,:)
               ENDIF

               llnightav = .FALSE.
#endif
            END SELECT

            CALL obs_surf_opt( surfdataqc(jtype), kstp, jpi, jpj,       &
               &               nit000, idaystp, zsurfvar, tmask(:,:,1), &
               &               nn_2dint, llnightav )

         END DO

      ENDIF

      CALL wrk_dealloc( jpi, jpj, jpk, zprofvar1 )
      CALL wrk_dealloc( jpi, jpj, jpk, zprofvar2 )
      CALL wrk_dealloc( jpi, jpj, jpk, zprofmask1 )
      CALL wrk_dealloc( jpi, jpj, jpk, zprofmask2 )
      CALL wrk_dealloc( jpi, jpj, zsurfvar )
      CALL wrk_dealloc( jpi, jpj, zglam1 )
      CALL wrk_dealloc( jpi, jpj, zglam2 )
      CALL wrk_dealloc( jpi, jpj, zgphi1 )
      CALL wrk_dealloc( jpi, jpj, zgphi2 )
#if ! defined key_lim2 && ! defined key_lim3
      CALL wrk_dealloc(jpi,jpj,frld)
#endif

   END SUBROUTINE dia_obs

   SUBROUTINE dia_obs_wri
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs_wri  ***
      !!          
      !! ** Purpose : Call observation diagnostic output routines
      !!
      !! ** Method  : Call observation diagnostic output routines
      !!
      !! ** Action  :
      !!
      !! History :
      !!        !  06-03  (K. Mogensen) Original code
      !!        !  06-05  (K. Mogensen) Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  07-03  (K. Mogensen) General handling of profiles
      !!        !  08-09  (M. Valdivieso) Velocity component (U,V) profiles
      !!        !  15-08  (M. Martin) Combined writing for prof and surf types
      !!----------------------------------------------------------------------
      !! * Modules used
      USE obs_rot_vel          ! Rotation of velocities

      IMPLICIT NONE

      !! * Local declarations
      INTEGER :: jtype                    ! Data set loop variable
      INTEGER :: jo, jvar, jk
      REAL(wp), DIMENSION(:), ALLOCATABLE :: &
         & zu, &
         & zv

      !-----------------------------------------------------------------------
      ! Depending on switches call various observation output routines
      !-----------------------------------------------------------------------

      IF ( nproftypes > 0 ) THEN

         DO jtype = 1, nproftypes

            IF ( TRIM(cobstypesprof(jtype)) == 'vel' ) THEN

               ! For velocity data, rotate the model velocities to N/S, E/W
               ! using the compressed data structure.
               ALLOCATE( &
                  & zu(profdataqc(jtype)%nvprot(1)), &
                  & zv(profdataqc(jtype)%nvprot(2))  &
                  & )

               CALL obs_rotvel( profdataqc(jtype), nn_2dint, zu, zv )

               DO jo = 1, profdataqc(jtype)%nprof
                  DO jvar = 1, 2
                     DO jk = profdataqc(jtype)%npvsta(jo,jvar), profdataqc(jtype)%npvend(jo,jvar)

                        IF ( jvar == 1 ) THEN
                           profdataqc(jtype)%var(jvar)%vmod(jk) = zu(jk)
                        ELSE
                           profdataqc(jtype)%var(jvar)%vmod(jk) = zv(jk)
                        ENDIF

                     END DO
                  END DO
               END DO

               DEALLOCATE( zu )
               DEALLOCATE( zv )

            END IF

            CALL obs_prof_decompress( profdataqc(jtype), &
               &                      profdata(jtype), .TRUE., numout )

            CALL obs_wri_prof( profdata(jtype) )

         END DO

      ENDIF

      IF ( nsurftypes > 0 ) THEN

         DO jtype = 1, nsurftypes

            CALL obs_surf_decompress( surfdataqc(jtype), &
               &                      surfdata(jtype), .TRUE., numout )

            CALL obs_wri_surf( surfdata(jtype) )

         END DO

      ENDIF

   END SUBROUTINE dia_obs_wri

   SUBROUTINE dia_obs_dealloc
      IMPLICIT NONE
      !!----------------------------------------------------------------------
      !!                    *** ROUTINE dia_obs_dealloc ***
      !!
      !!  ** Purpose : To deallocate data to enable the obs_oper online loop.
      !!               Specifically: dia_obs_init --> dia_obs --> dia_obs_wri
      !!
      !!  ** Method : Clean up various arrays left behind by the obs_oper.
      !!
      !!  ** Action :
      !!
      !!----------------------------------------------------------------------
      ! obs_grid deallocation
      CALL obs_grid_deallocate

      ! diaobs deallocation
      IF ( nproftypes > 0 ) &
         &   DEALLOCATE( cobstypesprof, profdata, profdataqc, nvarsprof, nextrprof )

      IF ( nsurftypes > 0 ) &
         &   DEALLOCATE( cobstypessurf, surfdata, surfdataqc, nvarssurf, nextrsurf )

   END SUBROUTINE dia_obs_dealloc

   SUBROUTINE calc_date( kstp, ddobs )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE calc_date  ***
      !!          
      !! ** Purpose : Get date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Method  : Get date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Action  : Get date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Action  : Get initial date in double precision YYYYMMDD.HHMMSS format
      !!
      !! History :
      !!        !  06-03  (K. Mogensen)  Original code
      !!        !  06-05  (K. Mogensen)  Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  06-10  (G. Smith) Calculates initial date the same as method for final date
      !!        !  10-05  (D. Lea) Update to month length calculation for NEMO vn3.2
      !!        !  2014-09  (D. Lea) New generic routine now deals with arbitrary initial time of day
      !!----------------------------------------------------------------------
      USE phycst, ONLY : &            ! Physical constants
         & rday
      USE dom_oce, ONLY : &           ! Ocean space and time domain variables
         & rdt

      IMPLICIT NONE

      !! * Arguments
      REAL(KIND=dp), INTENT(OUT) :: ddobs                        ! Date in YYYYMMDD.HHMMSS
      INTEGER :: kstp

      !! * Local declarations
      INTEGER :: iyea        ! date - (year, month, day, hour, minute)
      INTEGER :: imon
      INTEGER :: iday
      INTEGER :: ihou
      INTEGER :: imin
      INTEGER :: imday       ! Number of days in month.
      REAL(wp) :: zdayfrc    ! Fraction of day

      INTEGER, DIMENSION(12) ::   imonth_len    !: length in days of the months of the current year

      !!----------------------------------------------------------------------
      !! Initial date initialization (year, month, day, hour, minute)
      !!----------------------------------------------------------------------
      iyea =   ndate0 / 10000
      imon = ( ndate0 - iyea * 10000 ) / 100
      iday =   ndate0 - iyea * 10000 - imon * 100
      ihou =   nn_time0 / 100
      imin = ( nn_time0 - ihou * 100 ) 

      !!----------------------------------------------------------------------
      !! Compute number of days + number of hours + min since initial time
      !!----------------------------------------------------------------------
      zdayfrc = kstp * rdt / rday
      zdayfrc = zdayfrc - aint(zdayfrc)
      imin = imin + int( zdayfrc * 24 * 60 ) 
      DO WHILE (imin >= 60) 
        imin=imin-60
        ihou=ihou+1
      END DO
      DO WHILE (ihou >= 24)
        ihou=ihou-24
        iday=iday+1
      END DO 
      iday = iday + kstp * rdt / rday 

      !-----------------------------------------------------------------------
      ! Convert number of days (iday) into a real date
      !----------------------------------------------------------------------

      CALL calc_month_len( iyea, imonth_len )

      DO WHILE ( iday > imonth_len(imon) )
         iday = iday - imonth_len(imon)
         imon = imon + 1 
         IF ( imon > 12 ) THEN
            imon = 1
            iyea = iyea + 1
            CALL calc_month_len( iyea, imonth_len )  ! update month lengths
         ENDIF
      END DO

      !----------------------------------------------------------------------
      ! Convert it into YYYYMMDD.HHMMSS format.
      !----------------------------------------------------------------------
      ddobs = iyea * 10000_dp + imon * 100_dp + &
         &    iday + ihou * 0.01_dp + imin * 0.0001_dp

   END SUBROUTINE calc_date

   SUBROUTINE ini_date( ddobsini )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ini_date  ***
      !!          
      !! ** Purpose : Get initial date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  06-03  (K. Mogensen)  Original code
      !!        !  06-05  (K. Mogensen)  Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  10-05  (D. Lea) Update to month length calculation for NEMO vn3.2
      !!        !  2014-09  (D. Lea) Change to call generic routine calc_date
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      REAL(KIND=dp), INTENT(OUT) :: ddobsini                   ! Initial date in YYYYMMDD.HHMMSS

      CALL calc_date( nit000 - 1, ddobsini )

   END SUBROUTINE ini_date

   SUBROUTINE fin_date( ddobsfin )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE fin_date  ***
      !!          
      !! ** Purpose : Get final date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  06-03  (K. Mogensen)  Original code
      !!        !  06-05  (K. Mogensen)  Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  10-05  (D. Lea) Update to month length calculation for NEMO vn3.2
      !!        !  2014-09  (D. Lea) Change to call generic routine calc_date
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      REAL(dp), INTENT(OUT) :: ddobsfin ! Final date in YYYYMMDD.HHMMSS

      CALL calc_date( nitend, ddobsfin )

   END SUBROUTINE fin_date
   
END MODULE diaobs
