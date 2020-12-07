MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  PISCES  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_sink_alloc

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking, sinking2  !: POC sinking fluxes 
   !                                                          !  (different meanings depending on the parameterization)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkingn, sinking2n  !: POC sinking fluxes 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkingp, sinking2p  !: POC sinking fluxes 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkcal, sinksil   !: CaCO3 and BSi sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer            !: Small BFe sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer2           !: Big iron sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfep      !: Fep sinking fluxes

   INTEGER  :: ik100

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zsink.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!   'standard sinking parameterisation'                  ???
   !!----------------------------------------------------------------------

   SUBROUTINE p4z_sink ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER  ::   ji, jj, jk, jit
      INTEGER  ::   iiter1, iiter2
      REAL(wp) ::   zagg1, zagg2, zagg3, zagg4
      REAL(wp) ::   zagg , zaggfe, zaggdoc, zaggdoc2, zaggdoc3
      REAL(wp) ::   zfact, zwsmax, zmax
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zw3d
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')


      ! Initialization of some global variables
      ! ---------------------------------------
      prodpoc(:,:,:) = 0.
      conspoc(:,:,:) = 0.
      prodgoc(:,:,:) = 0.
      consgoc(:,:,:) = 0.

      !
      !    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      !    -----------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zmax  = MAX( heup_01(ji,jj), hmld(ji,jj) )
               zfact = MAX( 0., gdepw_n(ji,jj,jk+1) - zmax ) / wsbio2scale
               wsbio4(ji,jj,jk) = wsbio2 + MAX(0., ( wsbio2max - wsbio2 )) * zfact
            END DO
         END DO
      END DO

      ! limit the values of the sinking speeds to avoid numerical instabilities  
      wsbio3(:,:,:) = wsbio

      !
      ! OA This is (I hope) a temporary solution for the problem that may 
      ! OA arise in specific situation where the CFL criterion is broken 
      ! OA for vertical sedimentation of particles. To avoid this, a time
      ! OA splitting algorithm has been coded. A specific maximum
      ! OA iteration number is provided and may be specified in the namelist 
      ! OA This is to avoid very large iteration number when explicit free
      ! OA surface is used (for instance). When niter?max is set to 1, 
      ! OA this computation is skipped. The crude old threshold method is 
      ! OA then applied. This also happens when niter exceeds nitermax.
      IF( MAX( niter1max, niter2max ) == 1 ) THEN
        iiter1 = 1
        iiter2 = 1
      ELSE
        iiter1 = 1
        iiter2 = 1
        DO jk = 1, jpkm1
          DO jj = 1, jpj
             DO ji = 1, jpi
                IF( tmask(ji,jj,jk) == 1) THEN
                   zwsmax =  0.5 * e3t_n(ji,jj,jk) / xstep
                   iiter1 =  MAX( iiter1, INT( wsbio3(ji,jj,jk) / zwsmax ) )
                   iiter2 =  MAX( iiter2, INT( wsbio4(ji,jj,jk) / zwsmax ) )
                ENDIF
             END DO
          END DO
        END DO
        IF( lk_mpp ) THEN
           CALL mpp_max( iiter1 )
           CALL mpp_max( iiter2 )
        ENDIF
        iiter1 = MIN( iiter1, niter1max )
        iiter2 = MIN( iiter2, niter2max )
      ENDIF

      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,jk) == 1 ) THEN
                 zwsmax = 0.5 * e3t_n(ji,jj,jk) / xstep
                 wsbio3(ji,jj,jk) = MIN( wsbio3(ji,jj,jk), zwsmax * REAL( iiter1, wp ) )
                 wsbio4(ji,jj,jk) = MIN( wsbio4(ji,jj,jk), zwsmax * REAL( iiter2, wp ) )
               ENDIF
            END DO
         END DO
      END DO

      wscal (:,:,:) = wsbio4(:,:,:)

      !  Initializa to zero all the sinking arrays 
      !   -----------------------------------------
      sinking (:,:,:) = 0.e0
      sinking2(:,:,:) = 0.e0
      sinkcal (:,:,:) = 0.e0
      sinkfer (:,:,:) = 0.e0
      sinksil (:,:,:) = 0.e0
      sinkfer2(:,:,:) = 0.e0

      !   Compute the sedimentation term using p4zsink2 for all the sinking particles
      !   -----------------------------------------------------
      DO jit = 1, iiter1
        CALL p4z_sink2( wsbio3, sinking , jppoc, iiter1 )
        CALL p4z_sink2( wsbio3, sinkfer , jpsfe, iiter1 )
      END DO

      DO jit = 1, iiter2
        CALL p4z_sink2( wsbio4, sinking2, jpgoc, iiter2 )
        CALL p4z_sink2( wsbio4, sinkfer2, jpbfe, iiter2 )
        CALL p4z_sink2( wsbio4, sinksil , jpgsi, iiter2 )
        CALL p4z_sink2( wscal , sinkcal , jpcal, iiter2 )
      END DO

      IF( ln_p5z ) THEN
         sinkingn (:,:,:) = 0.e0
         sinking2n(:,:,:) = 0.e0
         sinkingp (:,:,:) = 0.e0
         sinking2p(:,:,:) = 0.e0

         !   Compute the sedimentation term using p4zsink2 for all the sinking particles
         !   -----------------------------------------------------
         DO jit = 1, iiter1
           CALL p4z_sink2( wsbio3, sinkingn , jppon, iiter1 )
           CALL p4z_sink2( wsbio3, sinkingp , jppop, iiter1 )
         END DO

         DO jit = 1, iiter2
           CALL p4z_sink2( wsbio4, sinking2n, jpgon, iiter2 )
           CALL p4z_sink2( wsbio4, sinking2p, jpgop, iiter2 )
         END DO
      ENDIF

      IF( ln_ligand ) THEN
         wsfep (:,:,:) = wfep
         DO jk = 1,jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( tmask(ji,jj,jk) == 1 ) THEN
                    zwsmax = 0.5 * e3t_n(ji,jj,jk) / xstep
                    wsfep(ji,jj,jk) = MIN( wsfep(ji,jj,jk), zwsmax * REAL( iiter1, wp ) )
                  ENDIF
               END DO
            END DO
         END DO
         !
         sinkfep(:,:,:) = 0.e0
         DO jit = 1, iiter1
           CALL p4z_sink2( wsfep, sinkfep , jpfep, iiter1 )
         END DO
      ENDIF

     ! Total carbon export per year
     IF( iom_use( "tcexp" ) .OR. ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  )  &
        &   t_oce_co2_exp = glob_sum( ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * e1e2t(:,:) * tmask(:,:,1) )
     !
     IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          CALL wrk_alloc( jpi, jpj,      zw2d )
          CALL wrk_alloc( jpi, jpj, jpk, zw3d )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "EPC100" ) )  THEN
              zw2d(:,:) = ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of carbon at 100m
              CALL iom_put( "EPC100"  , zw2d )
          ENDIF
          IF( iom_use( "EPFE100" ) )  THEN
              zw2d(:,:) = ( sinkfer(:,:,ik100) + sinkfer2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of iron at 100m
              CALL iom_put( "EPFE100"  , zw2d )
          ENDIF
          IF( iom_use( "EPCAL100" ) )  THEN
              zw2d(:,:) = sinkcal(:,:,ik100) * zfact * tmask(:,:,1) ! Export of calcite at 100m
              CALL iom_put( "EPCAL100"  , zw2d )
          ENDIF
          IF( iom_use( "EPSI100" ) )  THEN
              zw2d(:,:) =  sinksil(:,:,ik100) * zfact * tmask(:,:,1) ! Export of bigenic silica at 100m
              CALL iom_put( "EPSI100"  , zw2d )
          ENDIF
          IF( iom_use( "EXPC" ) )  THEN
              zw3d(:,:,:) = ( sinking(:,:,:) + sinking2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of carbon in the water column
              CALL iom_put( "EXPC"  , zw3d )
          ENDIF
          IF( iom_use( "EXPFE" ) )  THEN
              zw3d(:,:,:) = ( sinkfer(:,:,:) + sinkfer2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of iron 
              CALL iom_put( "EXPFE"  , zw3d )
          ENDIF
          IF( iom_use( "EXPCAL" ) )  THEN
              zw3d(:,:,:) = sinkcal(:,:,:) * zfact * tmask(:,:,:) ! Export of calcite 
              CALL iom_put( "EXPCAL"  , zw3d )
          ENDIF
          IF( iom_use( "EXPSI" ) )  THEN
              zw3d(:,:,:) = sinksil(:,:,:) * zfact * tmask(:,:,:) ! Export of bigenic silica
              CALL iom_put( "EXPSI"  , zw3d )
          ENDIF
          IF( iom_use( "tcexp" ) )  CALL iom_put( "tcexp" , t_oce_co2_exp * zfact )   ! molC/s
          ! 
          CALL wrk_dealloc( jpi, jpj,      zw2d )
          CALL wrk_dealloc( jpi, jpj, jpk, zw3d )
        ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink

   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!----------------------------------------------------------------------
      INTEGER :: jk

      ik100 = 10        !  last level where depth less than 100 m
      DO jk = jpkm1, 1, -1
         IF( gdept_1d(jk) > 100. )  ik100 = jk - 1
      END DO
      IF (lwp) WRITE(numout,*)
      IF (lwp) WRITE(numout,*) ' Level corresponding to 100m depth ',  ik100 + 1
      IF (lwp) WRITE(numout,*)
      !
      t_oce_co2_exp = 0._wp
      !
   END SUBROUTINE p4z_sink_init

   SUBROUTINE p4z_sink2( pwsink, psinkflx, jp_tra, kiter )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink2  ***
      !!
      !! ** Purpose :   Compute the sedimentation terms for the various sinking
      !!     particles. The scheme used to compute the trends is based
      !!     on MUSCL.
      !!
      !! ** Method  : - this ROUTINE compute not exactly the advection but the
      !!      transport term, i.e.  div(u*tra).
      !!---------------------------------------------------------------------
      !
      INTEGER , INTENT(in   )                         ::   jp_tra    ! tracer index index      
      INTEGER , INTENT(in   )                         ::   kiter     ! number of iterations for time-splitting 
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj,jpk) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   psinkflx  ! sinking fluxe
      !!
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zigma,zew,zign, zflx, zstep
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztraz, zakz, zwsink2, ztrb 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink2')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )

      zstep = rfact2 / REAL( kiter, wp ) / 2.

      ztraz(:,:,:) = 0.e0
      zakz (:,:,:) = 0.e0
      ztrb (:,:,:) = trb(:,:,:,jp_tra)

      DO jk = 1, jpkm1
         zwsink2(:,:,jk+1) = -pwsink(:,:,jk) / rday * tmask(:,:,jk+1) 
      END DO
      zwsink2(:,:,1) = 0.e0


      ! Vertical advective flux
      DO jn = 1, 2
         !  first guess of the slopes interior values
         DO jk = 2, jpkm1
            ztraz(:,:,jk) = ( trb(:,:,jk-1,jp_tra) - trb(:,:,jk,jp_tra) ) * tmask(:,:,jk)
         END DO
         ztraz(:,:,1  ) = 0.0
         ztraz(:,:,jpk) = 0.0

         ! slopes
         DO jk = 2, jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zign = 0.25 + SIGN( 0.25, ztraz(ji,jj,jk) * ztraz(ji,jj,jk+1) )
                  zakz(ji,jj,jk) = ( ztraz(ji,jj,jk) + ztraz(ji,jj,jk+1) ) * zign
               END DO
            END DO
         END DO
         
         ! Slopes limitation
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zakz(ji,jj,jk) = SIGN( 1., zakz(ji,jj,jk) ) *        &
                     &             MIN( ABS( zakz(ji,jj,jk) ), 2. * ABS(ztraz(ji,jj,jk+1)), 2. * ABS(ztraz(ji,jj,jk) ) )
               END DO
            END DO
         END DO
         
         ! vertical advective flux
         DO jk = 1, jpkm1
            DO jj = 1, jpj      
               DO ji = 1, jpi    
                  zigma = zwsink2(ji,jj,jk+1) * zstep / e3w_n(ji,jj,jk+1)
                  zew   = zwsink2(ji,jj,jk+1)
                  psinkflx(ji,jj,jk+1) = -zew * ( trb(ji,jj,jk,jp_tra) - 0.5 * ( 1 + zigma ) * zakz(ji,jj,jk) ) * zstep
               END DO
            END DO
         END DO
         !
         ! Boundary conditions
         psinkflx(:,:,1  ) = 0.e0
         psinkflx(:,:,jpk) = 0.e0
         
         DO jk=1,jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / e3t_n(ji,jj,jk)
                  trb(ji,jj,jk,jp_tra) = trb(ji,jj,jk,jp_tra) + zflx
               END DO
            END DO
         END DO

      ENDDO

      DO jk = 1,jpkm1
         DO jj = 1,jpj
            DO ji = 1, jpi
               zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / e3t_n(ji,jj,jk)
               ztrb(ji,jj,jk) = ztrb(ji,jj,jk) + 2. * zflx
            END DO
         END DO
      END DO

      trb(:,:,:,jp_tra) = ztrb(:,:,:)
      psinkflx(:,:,:)   = 2. * psinkflx(:,:,:)
      !
      CALL wrk_dealloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink2')
      !
   END SUBROUTINE p4z_sink2


   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(3)

      ierr(:) = 0
      !
      ALLOCATE( sinking(jpi,jpj,jpk) , sinking2(jpi,jpj,jpk)                    ,     &                
         &      sinkcal(jpi,jpj,jpk) , sinksil (jpi,jpj,jpk)                    ,     &                
         &      sinkfer2(jpi,jpj,jpk)                                           ,     &                
         &      sinkfer(jpi,jpj,jpk)                                            , STAT=ierr(1) )                
         !
      IF( ln_ligand ) ALLOCATE( sinkfep(jpi,jpj,jpk)                            , STAT=ierr(2) )  
         
      IF( ln_p5z    ) ALLOCATE( sinkingn(jpi,jpj,jpk), sinking2n(jpi,jpj,jpk)   ,     &
         &                      sinkingp(jpi,jpj,jpk), sinking2p(jpi,jpj,jpk)   , STAT=ierr(3) )
      !
      p4z_sink_alloc = MAXVAL( ierr )
      IF( p4z_sink_alloc /= 0 ) CALL ctl_warn('p4z_sink_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_sink_alloc
   
   !!======================================================================
END MODULE p4zsink
