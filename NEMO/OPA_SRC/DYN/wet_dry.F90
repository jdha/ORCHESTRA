MODULE wet_dry
   !!==============================================================================
   !!                       ***  MODULE  wet_dry  ***
   !! Wetting and drying includes initialisation routine and routines to
   !! compute and apply flux limiters and preserve water depth positivity
   !! only effects if wetting/drying is on (ln_wd == .true.)
   !!==============================================================================
   !! History :  3.6  ! 2014-09  ((H.Liu)  Original code
   !!                 ! will add the runoff and periodic BC case later
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   wad_lmt    : Compute the horizontal flux limiter and the limited velocity
   !!                when wetting and drying happens 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce, ONLY : ln_rnf   ! surface boundary condition: ocean
   USE sbcrnf          ! river runoff 
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   !!----------------------------------------------------------------------
   !! critical depths,filters, limiters,and masks for  Wetting and Drying
   !! ---------------------------------------------------------------------

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   wdmask         !: u- and v- limiter 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   ht_wd          !: wetting and drying t-pt depths
                                                                     !  (can include negative depths)

   LOGICAL,  PUBLIC  ::   ln_wd       !: Wetting/drying activation switch (T:on,F:off)
   REAL(wp), PUBLIC  ::   rn_wdmin1   !: minimum water depth on dried cells
   REAL(wp), PUBLIC  ::   rn_wdmin2   !: tolerrance of minimum water depth on dried cells
   REAL(wp), PUBLIC  ::   rn_wdld     !: land elevation below which wetting/drying 
                                           !: will be considered
   INTEGER , PUBLIC  ::   nn_wdit     !: maximum number of iteration for W/D limiter

   PUBLIC   wad_init                  ! initialisation routine called by step.F90
   PUBLIC   wad_lmt                   ! routine called by sshwzv.F90
   PUBLIC   wad_lmt_bt                ! routine called by dynspg_ts.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
CONTAINS

   SUBROUTINE wad_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE wad_init  ***
      !!                    
      !! ** Purpose :   read wetting and drying namelist and print the variables.
      !!
      !! ** input   : - namwad namelist
      !!----------------------------------------------------------------------
      NAMELIST/namwad/ ln_wd, rn_wdmin1, rn_wdmin2, rn_wdld, nn_wdit
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      INTEGER  ::   ierr                ! Local integer status array allocation 
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist namwad in reference namelist 
                                        ! : Parameters for Wetting/Drying
      READ  ( numnam_ref, namwad, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namwad in reference namelist', .TRUE.) 
      REWIND( numnam_cfg )              ! Namelist namwad in configuration namelist 
                                        ! : Parameters for Wetting/Drying
      READ  ( numnam_cfg, namwad, IOSTAT = ios, ERR = 906)
906   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namwad in configuration namelist', .TRUE. )
      IF(lwm) WRITE ( numond, namwad )

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'wad_init : Wetting and drying initialization through namelist read'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namwad'
         WRITE(numout,*) '      Logical activation                 ln_wd      = ', ln_wd
         WRITE(numout,*) '      Minimum wet depth on dried cells rn_wdmin1    = ', rn_wdmin1
         WRITE(numout,*) '      Tolerance of min wet depth     rn_wdmin2      = ', rn_wdmin2
         WRITE(numout,*) '      land elevation threshold         rn_wdld      = ', rn_wdld
         WRITE(numout,*) '      Max iteration for W/D limiter    nn_wdit      = ', nn_wdit
      ENDIF
      !
      IF(ln_wd) THEN
         ALLOCATE( wdmask(jpi,jpj), ht_wd(jpi,jpj),  STAT=ierr )
         IF( ierr /= 0 ) CALL ctl_stop('STOP', 'wad_init : Array allocation error')
      ENDIF
      !
   END SUBROUTINE wad_init


   SUBROUTINE wad_lmt( sshb1, sshemp, z2dt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_lmt  ***
      !!                    
      !! ** Purpose :   generate flux limiters for wetting/drying
      !!
      !! ** Method  : - Prevent negative depth occurring (Not ready for Agrif) 
      !!
      !! ** Action  : - calculate flux limiter and W/D flag
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   sshb1
      REAL(wp), DIMENSION(:,:), INTENT(in)    ::   sshemp
      REAL(wp), INTENT(in) :: z2dt
      !
      INTEGER  ::   ji, jj, jk, jk1     ! dummy loop indices
      INTEGER  ::   jflag               ! local scalar
      REAL(wp) ::   zcoef, zdep1, zdep2 ! local scalars
      REAL(wp) ::   zzflxp, zzflxn      ! local scalars
      REAL(wp) ::   zdepwd              ! local scalar, always wet cell depth
      REAL(wp) ::   ztmp                ! local scalars
      REAL(wp), POINTER,  DIMENSION(:,:) ::   zwdlmtu, zwdlmtv         !: W/D flux limiters
      REAL(wp), POINTER,  DIMENSION(:,:) ::   zflxp,  zflxn            ! local 2D workspace
      REAL(wp), POINTER,  DIMENSION(:,:) ::   zflxu,  zflxv            ! local 2D workspace
      REAL(wp), POINTER,  DIMENSION(:,:) ::   zflxu1, zflxv1           ! local 2D workspace
      !!----------------------------------------------------------------------
      !

      IF( nn_timing == 1 )  CALL timing_start('wad_lmt')

      IF(ln_wd) THEN

        CALL wrk_alloc( jpi, jpj, zflxp, zflxn, zflxu, zflxv, zflxu1, zflxv1 )
        CALL wrk_alloc( jpi, jpj, zwdlmtu, zwdlmtv)
        !
       
        !IF(lwp) WRITE(numout,*)
        !IF(lwp) WRITE(numout,*) 'wad_lmt : wetting/drying limiters and velocity limiting'
       
        jflag  = 0
        zdepwd = 50._wp   !maximum depth on which that W/D could possibly happen

       
        zflxp(:,:)   = 0._wp
        zflxn(:,:)   = 0._wp
        zflxu(:,:)   = 0._wp
        zflxv(:,:)   = 0._wp

        zwdlmtu(:,:)  = 1._wp
        zwdlmtv(:,:)  = 1._wp
       
        ! Horizontal Flux in u and v direction
        DO jk = 1, jpkm1  
           DO jj = 1, jpj
              DO ji = 1, jpi
                 zflxu(ji,jj) = zflxu(ji,jj) + e3u_n(ji,jj,jk) * un(ji,jj,jk) * umask(ji,jj,jk)
                 zflxv(ji,jj) = zflxv(ji,jj) + e3v_n(ji,jj,jk) * vn(ji,jj,jk) * vmask(ji,jj,jk)
              END DO  
           END DO  
        END DO
       
        zflxu(:,:) = zflxu(:,:) * e2u(:,:)
        zflxv(:,:) = zflxv(:,:) * e1v(:,:)
       
        wdmask(:,:) = 1
        DO jj = 2, jpj
           DO ji = 2, jpi 

             IF( tmask(ji, jj, 1) < 0.5_wp ) CYCLE   ! we don't care about land cells
             IF( ht_wd(ji,jj) > zdepwd )      CYCLE   ! and cells which are unlikely to dry

              zflxp(ji,jj) = max(zflxu(ji,jj), 0._wp) - min(zflxu(ji-1,jj),   0._wp) + &
                           & max(zflxv(ji,jj), 0._wp) - min(zflxv(ji,  jj-1), 0._wp) 
              zflxn(ji,jj) = min(zflxu(ji,jj), 0._wp) - max(zflxu(ji-1,jj),   0._wp) + &
                           & min(zflxv(ji,jj), 0._wp) - max(zflxv(ji,  jj-1), 0._wp) 

              zdep2 = ht_wd(ji,jj) + sshb1(ji,jj) - rn_wdmin1
              IF(zdep2 .le. 0._wp) THEN  !add more safty, but not necessary
                sshb1(ji,jj) = rn_wdmin1 - ht_wd(ji,jj)
                IF(zflxu(ji,  jj) > 0._wp) zwdlmtu(ji  ,jj) = 0._wp
                IF(zflxu(ji-1,jj) < 0._wp) zwdlmtu(ji-1,jj) = 0._wp
                IF(zflxv(ji,  jj) > 0._wp) zwdlmtv(ji  ,jj) = 0._wp
                IF(zflxv(ji,jj-1) < 0._wp) zwdlmtv(ji,jj-1) = 0._wp 
                wdmask(ji,jj) = 0._wp
              END IF
           ENDDO
        END DO

      
        !! start limiter iterations 
        DO jk1 = 1, nn_wdit + 1
       
          
           zflxu1(:,:) = zflxu(:,:) * zwdlmtu(:,:)
           zflxv1(:,:) = zflxv(:,:) * zwdlmtv(:,:)
           jflag = 0     ! flag indicating if any further iterations are needed
          
           DO jj = 2, jpj
              DO ji = 2, jpi 
        
                 IF( tmask(ji, jj, 1) < 0.5_wp ) CYCLE 
                 IF( ht_wd(ji,jj) > zdepwd )      CYCLE
        
                 ztmp = e1e2t(ji,jj)

                 zzflxp = max(zflxu1(ji,jj), 0._wp) - min(zflxu1(ji-1,jj),   0._wp) + &
                        & max(zflxv1(ji,jj), 0._wp) - min(zflxv1(ji,  jj-1), 0._wp) 
                 zzflxn = min(zflxu1(ji,jj), 0._wp) - max(zflxu1(ji-1,jj),   0._wp) + &
                        & min(zflxv1(ji,jj), 0._wp) - max(zflxv1(ji,  jj-1), 0._wp) 
          
                 zdep1 = (zzflxp + zzflxn) * z2dt / ztmp
                 zdep2 = ht_wd(ji,jj) + sshb1(ji,jj) - rn_wdmin1 - z2dt * sshemp(ji,jj)
          
                 IF( zdep1 > zdep2 ) THEN
                   wdmask(ji, jj) = 0
                   zcoef = ( ( zdep2 - rn_wdmin2 ) * ztmp - zzflxn * z2dt ) / ( zflxp(ji,jj) * z2dt )
                   !zcoef = ( ( zdep2 - rn_wdmin2 ) * ztmp - zzflxn * z2dt ) / ( zzflxp * z2dt )
                   ! flag if the limiter has been used but stop flagging if the only
                   ! changes have zeroed the coefficient since further iterations will
                   ! not change anything
                   IF( zcoef > 0._wp ) THEN
                      jflag = 1 
                   ELSE
                      zcoef = 0._wp
                   ENDIF
                   IF(jk1 > nn_wdit) zcoef = 0._wp
                   IF(zflxu1(ji,  jj) > 0._wp) zwdlmtu(ji  ,jj) = zcoef
                   IF(zflxu1(ji-1,jj) < 0._wp) zwdlmtu(ji-1,jj) = zcoef
                   IF(zflxv1(ji,  jj) > 0._wp) zwdlmtv(ji  ,jj) = zcoef
                   IF(zflxv1(ji,jj-1) < 0._wp) zwdlmtv(ji,jj-1) = zcoef
                 END IF
              END DO ! ji loop
           END DO  ! jj loop

           CALL lbc_lnk( zwdlmtu, 'U', 1. )
           CALL lbc_lnk( zwdlmtv, 'V', 1. )

           IF(lk_mpp) CALL mpp_max(jflag)   !max over the global domain

           IF(jflag == 0) EXIT
          
        END DO  ! jk1 loop
       
        DO jk = 1, jpkm1
          un(:,:,jk) = un(:,:,jk) * zwdlmtu(:, :) 
          vn(:,:,jk) = vn(:,:,jk) * zwdlmtv(:, :) 
        END DO

        CALL lbc_lnk( un, 'U', -1. )
        CALL lbc_lnk( vn, 'V', -1. )
      !
        un_b(:,:) = un_b(:,:) * zwdlmtu(:, :)
        vn_b(:,:) = vn_b(:,:) * zwdlmtv(:, :)
        CALL lbc_lnk( un_b, 'U', -1. )
        CALL lbc_lnk( vn_b, 'V', -1. )
       
        IF(jflag == 1 .AND. lwp) WRITE(numout,*) 'Need more iterations in wad_lmt!!!'
       
        !IF( ln_rnf      )   CALL sbc_rnf_div( hdivn )          ! runoffs (update hdivn field)
        !IF( nn_cla == 1 )   CALL cla_div    ( kt )             ! Cross Land Advection (update hdivn field)
        !
        !
        CALL wrk_dealloc( jpi, jpj, zflxp, zflxn, zflxu, zflxv, zflxu1, zflxv1 )
        CALL wrk_dealloc( jpi, jpj, zwdlmtu, zwdlmtv)
        !
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('wad_lmt')
      !
   END SUBROUTINE wad_lmt


   SUBROUTINE wad_lmt_bt( zflxu, zflxv, sshn_e, zssh_frc, rdtbt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_lmt  ***
      !!                    
      !! ** Purpose :   limiting flux in the barotropic stepping (dynspg_ts)
      !!
      !! ** Method  : - Prevent negative depth occurring (Not ready for Agrif) 
      !!
      !! ** Action  : - calculate flux limiter and W/D flag
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in)    ::   rdtbt    ! ocean time-step index
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   zflxu,  zflxv, sshn_e, zssh_frc  
      !
      INTEGER  ::   ji, jj, jk, jk1     ! dummy loop indices
      INTEGER  ::   jflag               ! local scalar
      REAL(wp) ::   z2dt
      REAL(wp) ::   zcoef, zdep1, zdep2 ! local scalars
      REAL(wp) ::   zzflxp, zzflxn      ! local scalars
      REAL(wp) ::   zdepwd              ! local scalar, always wet cell depth
      REAL(wp) ::   ztmp                ! local scalars
      REAL(wp), POINTER,  DIMENSION(:,:) ::   zwdlmtu, zwdlmtv         !: W/D flux limiters
      REAL(wp), POINTER,  DIMENSION(:,:) ::   zflxp,  zflxn            ! local 2D workspace
      REAL(wp), POINTER,  DIMENSION(:,:) ::   zflxu1, zflxv1           ! local 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('wad_lmt_bt')

      IF(ln_wd) THEN

        CALL wrk_alloc( jpi, jpj, zflxp, zflxn, zflxu1, zflxv1 )
        CALL wrk_alloc( jpi, jpj, zwdlmtu, zwdlmtv)
        !
       
        !IF(lwp) WRITE(numout,*)
        !IF(lwp) WRITE(numout,*) 'wad_lmt_bt : wetting/drying limiters and velocity limiting'
       
        jflag  = 0
        zdepwd = 50._wp   !maximum depth that ocean cells can have W/D processes

        z2dt = rdtbt   
       
        zflxp(:,:)   = 0._wp
        zflxn(:,:)   = 0._wp

        zwdlmtu(:,:)  = 1._wp
        zwdlmtv(:,:)  = 1._wp
       
        ! Horizontal Flux in u and v direction
       
        DO jj = 2, jpj
           DO ji = 2, jpi 

             IF( tmask(ji, jj, 1 ) < 0.5_wp) CYCLE   ! we don't care about land cells
             IF( ht_wd(ji,jj) > zdepwd )      CYCLE   ! and cells which are unlikely to dry

              zflxp(ji,jj) = max(zflxu(ji,jj), 0._wp) - min(zflxu(ji-1,jj),   0._wp) + &
                           & max(zflxv(ji,jj), 0._wp) - min(zflxv(ji,  jj-1), 0._wp) 
              zflxn(ji,jj) = min(zflxu(ji,jj), 0._wp) - max(zflxu(ji-1,jj),   0._wp) + &
                           & min(zflxv(ji,jj), 0._wp) - max(zflxv(ji,  jj-1), 0._wp) 

              zdep2 = ht_wd(ji,jj) + sshn_e(ji,jj) - rn_wdmin1
              IF(zdep2 .le. 0._wp) THEN  !add more safety, but not necessary
                sshn_e(ji,jj) = rn_wdmin1 - ht_wd(ji,jj)
                IF(zflxu(ji,  jj) > 0._wp) zwdlmtu(ji  ,jj) = 0._wp
                IF(zflxu(ji-1,jj) < 0._wp) zwdlmtu(ji-1,jj) = 0._wp
                IF(zflxv(ji,  jj) > 0._wp) zwdlmtv(ji  ,jj) = 0._wp
                IF(zflxv(ji,jj-1) < 0._wp) zwdlmtv(ji,jj-1) = 0._wp 
              END IF
           ENDDO
        END DO

      
        !! start limiter iterations 
        DO jk1 = 1, nn_wdit + 1
       
          
           zflxu1(:,:) = zflxu(:,:) * zwdlmtu(:,:)
           zflxv1(:,:) = zflxv(:,:) * zwdlmtv(:,:)
           jflag = 0     ! flag indicating if any further iterations are needed
          
           DO jj = 2, jpj
              DO ji = 2, jpi 
        
                 IF( tmask(ji, jj, 1 ) < 0.5_wp) CYCLE 
                 IF( ht_wd(ji,jj) > zdepwd )      CYCLE
        
                 ztmp = e1e2t(ji,jj)

                 zzflxp = max(zflxu1(ji,jj), 0._wp) - min(zflxu1(ji-1,jj),   0._wp) + &
                        & max(zflxv1(ji,jj), 0._wp) - min(zflxv1(ji,  jj-1), 0._wp) 
                 zzflxn = min(zflxu1(ji,jj), 0._wp) - max(zflxu1(ji-1,jj),   0._wp) + &
                        & min(zflxv1(ji,jj), 0._wp) - max(zflxv1(ji,  jj-1), 0._wp) 
          
                 zdep1 = (zzflxp + zzflxn) * z2dt / ztmp
                 zdep2 = ht_wd(ji,jj) + sshn_e(ji,jj) - rn_wdmin1 - z2dt * zssh_frc(ji,jj)
          
                 IF(zdep1 > zdep2) THEN
                   zcoef = ( ( zdep2 - rn_wdmin2 ) * ztmp - zzflxn * z2dt ) / ( zflxp(ji,jj) * z2dt )
                   !zcoef = ( ( zdep2 - rn_wdmin2 ) * ztmp - zzflxn * z2dt ) / ( zzflxp * z2dt )
                   ! flag if the limiter has been used but stop flagging if the only
                   ! changes have zeroed the coefficient since further iterations will
                   ! not change anything
                   IF( zcoef > 0._wp ) THEN
                      jflag = 1 
                   ELSE
                      zcoef = 0._wp
                   ENDIF
                   IF(jk1 > nn_wdit) zcoef = 0._wp
                   IF(zflxu1(ji,  jj) > 0._wp) zwdlmtu(ji  ,jj) = zcoef
                   IF(zflxu1(ji-1,jj) < 0._wp) zwdlmtu(ji-1,jj) = zcoef
                   IF(zflxv1(ji,  jj) > 0._wp) zwdlmtv(ji  ,jj) = zcoef
                   IF(zflxv1(ji,jj-1) < 0._wp) zwdlmtv(ji,jj-1) = zcoef
                 END IF
              END DO ! ji loop
           END DO  ! jj loop

           CALL lbc_lnk( zwdlmtu, 'U', 1. )
           CALL lbc_lnk( zwdlmtv, 'V', 1. )

           IF(lk_mpp) CALL mpp_max(jflag)   !max over the global domain

           IF(jflag == 0) EXIT
          
        END DO  ! jk1 loop
       
        zflxu(:,:) = zflxu(:,:) * zwdlmtu(:, :) 
        zflxv(:,:) = zflxv(:,:) * zwdlmtv(:, :) 

        CALL lbc_lnk( zflxu, 'U', -1. )
        CALL lbc_lnk( zflxv, 'V', -1. )
       
        IF(jflag == 1 .AND. lwp) WRITE(numout,*) 'Need more iterations in wad_lmt_bt!!!'
       
        !IF( ln_rnf      )   CALL sbc_rnf_div( hdivn )          ! runoffs (update hdivn field)
        !IF( nn_cla == 1 )   CALL cla_div    ( kt )             ! Cross Land Advection (update hdivn field)
        !
        !
        CALL wrk_dealloc( jpi, jpj, zflxp, zflxn, zflxu1, zflxv1 )
        CALL wrk_dealloc( jpi, jpj, zwdlmtu, zwdlmtv)
        !
      END IF

      IF( nn_timing == 1 )  CALL timing_stop('wad_lmt')
   END SUBROUTINE wad_lmt_bt

   !!==============================================================================
END MODULE wet_dry
