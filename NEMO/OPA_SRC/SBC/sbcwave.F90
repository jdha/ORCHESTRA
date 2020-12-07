MODULE sbcwave
   !!======================================================================
   !!                       ***  MODULE  sbcwave  ***
   !! Wave module 
   !!======================================================================
   !! History :  3.3  !  2011-09  (M. Adani)  Original code: Drag Coefficient 
   !!         :  3.4  !  2012-10  (M. Adani)  Stokes Drift 
   !!            3.6  !  2014-09  (E. Clementi,P. Oddo) New Stokes Drift Computation
   !!             -   !  2016-12  (G. Madec, E. Clementi) update Stoke drift computation
   !!                                                    + add sbc_wave_ini routine
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_stokes    : calculate 3D Stokes-drift velocities
   !!   sbc_wave      : wave data from wave model in netcdf files 
   !!   sbc_wave_init : initialisation fo surface waves 
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants 
   USE oce            ! ocean variables
   USE sbc_oce	       ! Surface boundary condition: ocean fields
   USE zdf_oce,  ONLY : ln_zdfqiao
   USE bdy_oce        ! open boundary condition variables
   USE domvvl         ! domain: variable volume layers
   !
   USE iom            ! I/O manager library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE fldread	       ! read input fields
   USE wrk_nemo       !

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_stokes      ! routine called in sbccpl
   PUBLIC   sbc_wave        ! routine called in sbcmod
   PUBLIC   sbc_wave_init   ! routine called in sbcmod
   
   ! Variables checking if the wave parameters are coupled (if not, they are read from file)
   LOGICAL, PUBLIC ::   cpl_hsig   = .FALSE.
   LOGICAL, PUBLIC ::   cpl_phioc  = .FALSE.
   LOGICAL, PUBLIC ::   cpl_sdrftx = .FALSE.
   LOGICAL, PUBLIC ::   cpl_sdrfty = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wper   = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wnum   = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wstrf  = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wdrag  = .FALSE.

   INTEGER ::   jpfld    ! number of files to read for stokes drift
   INTEGER ::   jp_usd   ! index of stokes drift  (i-component) (m/s)    at T-point
   INTEGER ::   jp_vsd   ! index of stokes drift  (j-component) (m/s)    at T-point
   INTEGER ::   jp_hsw   ! index of significant wave hight      (m)      at T-point
   INTEGER ::   jp_wmp   ! index of mean wave period            (s)      at T-point

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_cd      ! structure of input fields (file informations, fields read) Drag Coefficient
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sd      ! structure of input fields (file informations, fields read) Stokes Drift
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_wn      ! structure of input fields (file informations, fields read) wave number for Qiao
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tauoc   ! structure of input fields (file informations, fields read) normalized wave stress into the ocean
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   cdn_wave            !:
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   hsw, wmp, wnum      !: 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tauoc_wave          !:  
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tsd2d               !: 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   div_sd              !: barotropic stokes drift divergence
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   ut0sd, vt0sd        !: surface Stokes drift velocities at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   usd  , vsd  , wsd   !: Stokes drift velocities at u-, v- & w-points, resp.

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014) 
   !! $Id: sbcwave.F90 7864 2017-03-31 15:40:41Z emanuelaclementi $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_stokes( )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_stokes  ***
      !!
      !! ** Purpose :   compute the 3d Stokes Drift according to Breivik et al.,
      !!                2014 (DOI: 10.1175/JPO-D-14-0020.1)
      !!
      !! ** Method  : - Calculate Stokes transport speed 
      !!              - Calculate horizontal divergence 
      !!              - Integrate the horizontal divergenze from the bottom 
      !! ** action  
      !!---------------------------------------------------------------------
      INTEGER  ::   jj, ji, jk   ! dummy loop argument
      INTEGER  ::   ik           ! local integer 
      REAL(wp) ::  ztransp, zfac, ztemp, zsp0
      REAL(wp) ::  zdep_u, zdep_v, zkh_u, zkh_v, zda_u, zda_v
      REAL(wp), DIMENSION(:,:)  , POINTER ::   zk_t, zk_u, zk_v, zu0_sd, zv0_sd   ! 2D workspace
      REAL(wp), DIMENSION(:,:,:), POINTER ::   ze3divh                            ! 3D workspace
      !!---------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi,jpj,jpk,   ze3divh )
      CALL wrk_alloc( jpi,jpj,       zk_t, zk_u, zk_v, zu0_sd, zv0_sd )
      !
      !
      zfac =  2.0_wp * rpi / 16.0_wp
      DO jj = 1, jpj                ! exp. wave number at t-point    (Eq. (19) in Breivick et al. (2014) )
         DO ji = 1, jpi
               ! Stokes drift velocity estimated from Hs and Tmean
               ztransp = zfac * hsw(ji,jj)*hsw(ji,jj) / MAX( wmp(ji,jj) , 0.0000001_wp )
               ! Stokes surface speed
               zsp0 = SQRT( ut0sd(ji,jj)*ut0sd(ji,jj) + vt0sd(ji,jj)*vt0sd(ji,jj) )
               tsd2d(ji,jj) = zsp0
               ! Wavenumber scale
               zk_t(ji,jj) = ABS( zsp0 ) / MAX( ABS( 5.97_wp*ztransp ) , 0.0000001_wp )
         END DO
      END DO      
      DO jj = 1, jpjm1              ! exp. wave number & Stokes drift velocity at u- & v-points
         DO ji = 1, jpim1
            zk_u(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji+1,jj) )
            zk_v(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji,jj+1) )
            !
            zu0_sd(ji,jj) = 0.5_wp * ( ut0sd(ji,jj) + ut0sd(ji+1,jj) )
            zv0_sd(ji,jj) = 0.5_wp * ( vt0sd(ji,jj) + vt0sd(ji,jj+1) )
         END DO
      END DO
      !
      !                       !==  horizontal Stokes Drift 3D velocity  ==!
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zdep_u = 0.5_wp * ( gdept_n(ji,jj,jk) + gdept_n(ji+1,jj,jk) )
               zdep_v = 0.5_wp * ( gdept_n(ji,jj,jk) + gdept_n(ji,jj+1,jk) )
               !                          
               zkh_u = zk_u(ji,jj) * zdep_u     ! k * depth
               zkh_v = zk_v(ji,jj) * zdep_v
               !                                ! Depth attenuation
               zda_u = EXP( -2.0_wp*zkh_u ) / ( 1.0_wp + 8.0_wp*zkh_u )
               zda_v = EXP( -2.0_wp*zkh_v ) / ( 1.0_wp + 8.0_wp*zkh_v )
               !
               usd(ji,jj,jk) = zda_u * zu0_sd(ji,jj) * umask(ji,jj,jk)
               vsd(ji,jj,jk) = zda_v * zv0_sd(ji,jj) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO   
      CALL lbc_lnk( usd(:,:,:), 'U', vsd(:,:,:), 'V', -1. )
      !
      !                       !==  vertical Stokes Drift 3D velocity  ==!
      !
      DO jk = 1, jpkm1               ! Horizontal e3*divergence
         DO jj = 2, jpj
            DO ji = fs_2, jpi
               ze3divh(ji,jj,jk) = (  e2u(ji  ,jj) * e3u_n(ji  ,jj,jk) * usd(ji  ,jj,jk)    &
                  &                 - e2u(ji-1,jj) * e3u_n(ji-1,jj,jk) * usd(ji-1,jj,jk)    &
                  &                 + e1v(ji,jj  ) * e3v_n(ji,jj  ,jk) * vsd(ji,jj  ,jk)    &
                  &                 - e1v(ji,jj-1) * e3v_n(ji,jj-1,jk) * vsd(ji,jj-1,jk)  ) * r1_e1e2t(ji,jj)
            END DO
         END DO
      END DO
      !
      IF( .NOT. AGRIF_Root() ) THEN
         IF( nbondi ==  1 .OR. nbondi == 2 )   ze3divh(nlci-1,   :  ,:) = 0._wp      ! east
         IF( nbondi == -1 .OR. nbondi == 2 )   ze3divh(  2   ,   :  ,:) = 0._wp      ! west
         IF( nbondj ==  1 .OR. nbondj == 2 )   ze3divh(  :   ,nlcj-1,:) = 0._wp      ! north
         IF( nbondj == -1 .OR. nbondj == 2 )   ze3divh(  :   ,  2   ,:) = 0._wp      ! south
      ENDIF
      !
      CALL lbc_lnk( ze3divh, 'T', 1. )
      !
      IF( ln_linssh ) THEN   ;   ik = 1   ! none zero velocity through the sea surface
      ELSE                   ;   ik = 2   ! w=0 at the surface (set one for all in sbc_wave_init)
      ENDIF
      DO jk = jpkm1, ik, -1          ! integrate from the bottom the hor. divergence (NB: at k=jpk w is always zero)
         wsd(:,:,jk) = wsd(:,:,jk+1) - ze3divh(:,:,jk)
      END DO
      !
      IF( ln_bdy ) THEN
         DO jk = 1, jpkm1
            wsd(:,:,jk) = wsd(:,:,jk) * bdytmask(:,:)
         END DO
      ENDIF
      !                       !==  Horizontal divergence of barotropic Stokes transport  ==!
      div_sd(:,:) = 0._wp
      DO jk = 1, jpkm1                                 ! 
        div_sd(:,:) = div_sd(:,:) + ze3divh(:,:,jk)
      END DO
      !
      CALL iom_put( "ustokes",  usd  )
      CALL iom_put( "vstokes",  vsd  )
      CALL iom_put( "wstokes",  wsd  )
      !
      CALL wrk_dealloc( jpi,jpj,jpk,   ze3divh )
      CALL wrk_dealloc( jpi,jpj,       zk_t, zk_u, zk_v, zu0_sd, zv0_sd )
      !
   END SUBROUTINE sbc_stokes


   SUBROUTINE sbc_wave( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wave  ***
      !!
      !! ** Purpose :   read wave parameters from wave model  in netcdf files.
      !!
      !! ** Method  : - Read namelist namsbc_wave
      !!              - Read Cd_n10 fields in netcdf files 
      !!              - Read stokes drift 2d in netcdf files 
      !!              - Read wave number in netcdf files 
      !!              - Compute 3d stokes drift using Breivik et al.,2014
      !!                formulation
      !! ** action  
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      !!---------------------------------------------------------------------
      !
      IF( ln_cdgw .AND. .NOT. cpl_wdrag ) THEN     !==  Neutral drag coefficient  ==!
         CALL fld_read( kt, nn_fsbc, sf_cd )             ! read from external forcing
         cdn_wave(:,:) = sf_cd(1)%fnow(:,:,1)
      ENDIF

      IF( ln_tauoc .AND. .NOT. cpl_wstrf ) THEN    !==  Wave induced stress  ==!
         CALL fld_read( kt, nn_fsbc, sf_tauoc )          ! read wave norm stress from external forcing
         tauoc_wave(:,:) = sf_tauoc(1)%fnow(:,:,1)
      ENDIF

      IF( ln_sdw )  THEN                           !==  Computation of the 3d Stokes Drift  ==! 
         !
         IF( jpfld > 0 ) THEN                            ! Read from file only if the field is not coupled
            CALL fld_read( kt, nn_fsbc, sf_sd )          ! read wave parameters from external forcing
            IF( jp_hsw > 0 )   hsw  (:,:) = sf_sd(jp_hsw)%fnow(:,:,1)   ! significant wave height
            IF( jp_wmp > 0 )   wmp  (:,:) = sf_sd(jp_wmp)%fnow(:,:,1)   ! wave mean period
            IF( jp_usd > 0 )   ut0sd(:,:) = sf_sd(jp_usd)%fnow(:,:,1)   ! 2D zonal Stokes Drift at T point
            IF( jp_vsd > 0 )   vt0sd(:,:) = sf_sd(jp_vsd)%fnow(:,:,1)   ! 2D meridional Stokes Drift at T point
         ENDIF
         !
         ! Read also wave number if needed, so that it is available in coupling routines
         IF( ln_zdfqiao .AND. .NOT.cpl_wnum ) THEN
            CALL fld_read( kt, nn_fsbc, sf_wn )          ! read wave parameters from external forcing
            wnum(:,:) = sf_wn(1)%fnow(:,:,1)
         ENDIF
           
         !                                         !==  Computation of the 3d Stokes Drift  ==! 
         !
         IF( jpfld == 4 )   CALL sbc_stokes()            ! Calculate only if required fields are read
         !                                               ! In coupled wave model-NEMO case the call is done after coupling
         !
      ENDIF
      !
   END SUBROUTINE sbc_wave


   SUBROUTINE sbc_wave_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wave_init  ***
      !!
      !! ** Purpose :   read wave parameters from wave model  in netcdf files.
      !!
      !! ** Method  : - Read namelist namsbc_wave
      !!              - Read Cd_n10 fields in netcdf files 
      !!              - Read stokes drift 2d in netcdf files 
      !!              - Read wave number in netcdf files 
      !!              - Compute 3d stokes drift using Breivik et al.,2014
      !!                formulation
      !! ** action  
      !!---------------------------------------------------------------------
      INTEGER ::   ierror, ios   ! local integer
      INTEGER ::   ifpr
      !!
      CHARACTER(len=100)     ::  cn_dir                          ! Root directory for location of drag coefficient files
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::   slf_i     ! array of namelist informations on the fields to read
      TYPE(FLD_N)            ::  sn_cdg, sn_usd, sn_vsd,  &
                             &   sn_hsw, sn_wmp, sn_wnum, sn_tauoc      ! informations about the fields to be read
      !
      NAMELIST/namsbc_wave/  sn_cdg, cn_dir, sn_usd, sn_vsd, sn_hsw, sn_wmp, sn_wnum, sn_tauoc
      !!---------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namsbc_wave in reference namelist : File for drag coeff. from wave model
      READ  ( numnam_ref, namsbc_wave, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_wave in reference namelist', lwp )
         
      REWIND( numnam_cfg )              ! Namelist namsbc_wave in configuration namelist : File for drag coeff. from wave model
      READ  ( numnam_cfg, namsbc_wave, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_wave in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namsbc_wave )
      !
      IF( ln_cdgw ) THEN
         IF( .NOT. cpl_wdrag ) THEN
            ALLOCATE( sf_cd(1), STAT=ierror )           !* allocate and fill sf_wave with sn_cdg
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
            !
                                   ALLOCATE( sf_cd(1)%fnow(jpi,jpj,1)   )
            IF( sn_cdg%ln_tint )   ALLOCATE( sf_cd(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_cd, (/ sn_cdg /), cn_dir, 'sbc_wave_init', 'Wave module ', 'namsbc_wave' )
         ENDIF
         ALLOCATE( cdn_wave(jpi,jpj) )
      ENDIF

      IF( ln_tauoc ) THEN
         IF( .NOT. cpl_wstrf ) THEN
            ALLOCATE( sf_tauoc(1), STAT=ierror )           !* allocate and fill sf_wave with sn_tauoc
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
            !
                                    ALLOCATE( sf_tauoc(1)%fnow(jpi,jpj,1)   )
            IF( sn_tauoc%ln_tint )  ALLOCATE( sf_tauoc(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_tauoc, (/ sn_tauoc /), cn_dir, 'sbc_wave_init', 'Wave module', 'namsbc_wave' )
         ENDIF
         ALLOCATE( tauoc_wave(jpi,jpj) )
      ENDIF

      IF( ln_sdw ) THEN   ! Find out how many fields have to be read from file if not coupled
         jpfld=0
         jp_usd=0   ;   jp_vsd=0   ;   jp_hsw=0   ;   jp_wmp=0
         IF( .NOT. cpl_sdrftx ) THEN
            jpfld  = jpfld + 1
            jp_usd = jpfld
         ENDIF
         IF( .NOT. cpl_sdrfty ) THEN
            jpfld  = jpfld + 1
            jp_vsd = jpfld
         ENDIF
         IF( .NOT. cpl_hsig ) THEN
            jpfld  = jpfld + 1
            jp_hsw = jpfld
         ENDIF
         IF( .NOT. cpl_wper ) THEN
            jpfld  = jpfld + 1
            jp_wmp = jpfld
         ENDIF

         ! Read from file only the non-coupled fields 
         IF( jpfld > 0 ) THEN
            ALLOCATE( slf_i(jpfld) )
            IF( jp_usd > 0 )   slf_i(jp_usd) = sn_usd
            IF( jp_vsd > 0 )   slf_i(jp_vsd) = sn_vsd
            IF( jp_hsw > 0 )   slf_i(jp_hsw) = sn_hsw
            IF( jp_wmp > 0 )   slf_i(jp_wmp) = sn_wmp
            ALLOCATE( sf_sd(jpfld), STAT=ierror )   !* allocate and fill sf_sd with stokes drift
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
            !
            DO ifpr= 1, jpfld
               ALLOCATE( sf_sd(ifpr)%fnow(jpi,jpj,1) )
               IF( slf_i(ifpr)%ln_tint )   ALLOCATE( sf_sd(ifpr)%fdta(jpi,jpj,1,2) )
            END DO
            !
            CALL fld_fill( sf_sd, slf_i, cn_dir, 'sbc_wave_init', 'Wave module ', 'namsbc_wave' )
         ENDIF
         ALLOCATE( usd  (jpi,jpj,jpk), vsd  (jpi,jpj,jpk), wsd(jpi,jpj,jpk) )
         ALLOCATE( hsw  (jpi,jpj)    , wmp  (jpi,jpj)     )
         ALLOCATE( ut0sd(jpi,jpj)    , vt0sd(jpi,jpj)     )
         ALLOCATE( div_sd(jpi,jpj) )
         ALLOCATE( tsd2d (jpi,jpj) )
         usd(:,:,:) = 0._wp
         vsd(:,:,:) = 0._wp
         wsd(:,:,:) = 0._wp
         ! Wave number needed only if ln_zdfqiao=T
         IF( .NOT. cpl_wnum ) THEN
            ALLOCATE( sf_wn(1), STAT=ierror )           !* allocate and fill sf_wave with sn_wnum
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable toallocate sf_wave structure' )
                                   ALLOCATE( sf_wn(1)%fnow(jpi,jpj,1)   )
            IF( sn_wnum%ln_tint )  ALLOCATE( sf_wn(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_wn, (/ sn_wnum /), cn_dir, 'sbc_wave', 'Wave module', 'namsbc_wave' )
         ENDIF
         ALLOCATE( wnum(jpi,jpj) )
      ENDIF
      !
   END SUBROUTINE sbc_wave_init

   !!======================================================================
END MODULE sbcwave
