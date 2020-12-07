MODULE trcadv
   !!==============================================================================
   !!                       ***  MODULE  trcadv  ***
   !! Ocean passive tracers:  advection trend 
   !!==============================================================================
   !! History :  2.0  !  2005-11  (G. Madec)  Original code
   !!            3.0  !  2010-06  (C. Ethe)   Adapted to passive tracers
   !!            3.7  !  2014-05  (G. Madec, C. Ethe)  Add 2nd/4th order cases for CEN and FCT schemes 
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_adv       : compute ocean tracer advection trend
   !!   trc_adv_ini   : control the different options of advection scheme
   !!----------------------------------------------------------------------
   USE oce_trc        ! ocean dynamics and active tracers
   USE trc            ! ocean passive tracers variables
   USE traadv_cen     ! centered scheme           (tra_adv_cen  routine)
   USE traadv_fct     ! FCT      scheme           (tra_adv_fct  routine)
   USE traadv_mus     ! MUSCL    scheme           (tra_adv_mus  routine)
   USE traadv_ubs     ! UBS      scheme           (tra_adv_ubs  routine)
   USE traadv_qck     ! QUICKEST scheme           (tra_adv_qck  routine)
   USE traadv_mle     ! ML eddy induced velocity  (tra_adv_mle  routine)
   USE ldftra         ! lateral diffusion coefficient on tracers
   USE ldfslp         ! Lateral diffusion: slopes of neutral surfaces
   !
   USE prtctl_trc     ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_adv       
   PUBLIC   trc_adv_ini  

   !                            !!* Namelist namtrc_adv *
   LOGICAL ::   ln_trcadv_cen    ! centered scheme flag
   INTEGER ::      nn_cen_h, nn_cen_v   ! =2/4 : horizontal and vertical choices of the order of CEN scheme
   LOGICAL ::   ln_trcadv_fct    ! FCT scheme flag
   INTEGER ::      nn_fct_h, nn_fct_v   ! =2/4 : horizontal and vertical choices of the order of FCT scheme
   INTEGER ::      nn_fct_zts           ! >=1 : 2nd order FCT with vertical sub-timestepping
   LOGICAL ::   ln_trcadv_mus    ! MUSCL scheme flag
   LOGICAL ::      ln_mus_ups           ! use upstream scheme in vivcinity of river mouths
   LOGICAL ::   ln_trcadv_ubs    ! UBS scheme flag
   INTEGER ::      nn_ubs_v             ! =2/4 : vertical choice of the order of UBS scheme
   LOGICAL ::   ln_trcadv_qck    ! QUICKEST scheme flag

   !                                        ! choices of advection scheme:
   INTEGER, PARAMETER ::   np_NO_adv  = 0   ! no T-S advection
   INTEGER, PARAMETER ::   np_CEN     = 1   ! 2nd/4th order centered scheme
   INTEGER, PARAMETER ::   np_FCT     = 2   ! 2nd/4th order Flux Corrected Transport scheme
   INTEGER, PARAMETER ::   np_FCT_zts = 3   ! 2nd order FCT scheme with vertical sub-timestepping
   INTEGER, PARAMETER ::   np_MUS     = 4   ! MUSCL scheme
   INTEGER, PARAMETER ::   np_UBS     = 5   ! 3rd order Upstream Biased Scheme
   INTEGER, PARAMETER ::   np_QCK     = 6   ! QUICK scheme

   INTEGER ::              nadv             ! chosen advection scheme
   !
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.7 , NEMO Consortium (2015)
   !! $Id: trcadv.F90 7753 2017-03-03 11:46:59Z mocavero $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_adv( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_adv  ***
      !!
      !! ** Purpose :   compute the ocean tracer advection trend.
      !!
      !! ** Method  : - Update the tracer with the advection term following nadv
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER ::   jk   ! dummy loop index
      CHARACTER (len=22) ::   charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zun, zvn, zwn  ! effective velocity
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_adv')
      !
      CALL wrk_alloc( jpi,jpj,jpk,   zun, zvn, zwn )
      !                                               !==  effective transport  ==!
      IF( l_offline ) THEN
         zun(:,:,:) = un(:,:,:)     ! effective transport already in un/vn/wn
         zvn(:,:,:) = vn(:,:,:)
         zwn(:,:,:) = wn(:,:,:)
      ELSE
         !       
         DO jk = 1, jpkm1
            zun(:,:,jk) = e2u  (:,:) * e3u_n(:,:,jk) * un(:,:,jk)                   ! eulerian transport
            zvn(:,:,jk) = e1v  (:,:) * e3v_n(:,:,jk) * vn(:,:,jk)
            zwn(:,:,jk) = e1e2t(:,:)                 * wn(:,:,jk)
         END DO
         !
         IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN                                 ! add z-tilde and/or vvl corrections
            zun(:,:,:) = zun(:,:,:) + un_td(:,:,:)
            zvn(:,:,:) = zvn(:,:,:) + vn_td(:,:,:)
         ENDIF
         !
         IF( ln_ldfeiv .AND. .NOT. ln_traldf_triad )   & 
            &              CALL ldf_eiv_trp( kt, nittrc000, zun, zvn, zwn, 'TRC' )  ! add the eiv transport
         !
         IF( ln_mle    )   CALL tra_adv_mle( kt, nittrc000, zun, zvn, zwn, 'TRC' )  ! add the mle transport
         !
         zun(:,:,jpk) = 0._wp                                                       ! no transport trough the bottom
         zvn(:,:,jpk) = 0._wp
         zwn(:,:,jpk) = 0._wp
         !
      ENDIF
      !
      SELECT CASE ( nadv )                      !==  compute advection trend and add it to general trend  ==!
      !
      CASE ( np_CEN )                                    ! Centered : 2nd / 4th order
         CALL tra_adv_cen    ( kt, nittrc000,'TRC',       zun, zvn, zwn     , trn, tra, jptra, nn_cen_h, nn_cen_v )
      CASE ( np_FCT )                                    ! FCT      : 2nd / 4th order
         CALL tra_adv_fct    ( kt, nittrc000,'TRC', r2dttrc, zun, zvn, zwn, trb, trn, tra, jptra, nn_fct_h, nn_fct_v )
      CASE ( np_FCT_zts )                                ! 2nd order FCT with vertical time-splitting
         CALL tra_adv_fct_zts( kt, nittrc000,'TRC', r2dttrc, zun, zvn, zwn, trb, trn, tra, jptra        , nn_fct_zts )
      CASE ( np_MUS )                                    ! MUSCL
         CALL tra_adv_mus    ( kt, nittrc000,'TRC', r2dttrc, zun, zvn, zwn, trb,      tra, jptra        , ln_mus_ups ) 
      CASE ( np_UBS )                                    ! UBS
         CALL tra_adv_ubs    ( kt, nittrc000,'TRC', r2dttrc, zun, zvn, zwn, trb, trn, tra, jptra        , nn_ubs_v   )
      CASE ( np_QCK )                                    ! QUICKEST
         CALL tra_adv_qck    ( kt, nittrc000,'TRC', r2dttrc, zun, zvn, zwn, trb, trn, tra, jptra                     )
      !
      END SELECT
      !                  
      IF( ln_ctl )   THEN                             !== print mean trends (used for debugging)
         WRITE(charout, FMT="('adv ')")  ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      END IF
      !
      CALL wrk_dealloc( jpi,jpj,jpk,   zun, zvn, zwn )
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_adv')
      !
   END SUBROUTINE trc_adv


   SUBROUTINE trc_adv_ini
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trc_adv_ini  ***
      !!                
      !! ** Purpose : Control the consistency between namelist options for 
      !!              passive tracer advection schemes and set nadv
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio
      INTEGER ::  ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/namtrc_adv/ ln_trcadv_cen, nn_cen_h, nn_cen_v,               &   ! CEN
         &                 ln_trcadv_fct, nn_fct_h, nn_fct_v, nn_fct_zts,   &   ! FCT
         &                 ln_trcadv_mus,                     ln_mus_ups,   &   ! MUSCL
         &                 ln_trcadv_ubs,           nn_ubs_v,               &   ! UBS
         &                 ln_trcadv_qck                                        ! QCK
      !!----------------------------------------------------------------------
      !
      REWIND( numnat_ref )              !  namtrc_adv in reference namelist 
      READ  ( numnat_ref, namtrc_adv, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_adv in reference namelist', lwp )

      REWIND( numnat_cfg )              ! namtrc_adv in configuration namelist
      READ  ( numnat_cfg, namtrc_adv, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_adv in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_adv )

      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'trc_adv_ini : choice/control of the tracer advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtrc_adv : chose a advection scheme for tracers'
         WRITE(numout,*) '      centered scheme                           ln_trcadv_cen = ', ln_trcadv_cen
         WRITE(numout,*) '            horizontal 2nd/4th order               nn_cen_h   = ', nn_fct_h
         WRITE(numout,*) '            vertical   2nd/4th order               nn_cen_v   = ', nn_fct_v
         WRITE(numout,*) '      Flux Corrected Transport scheme           ln_trcadv_fct = ', ln_trcadv_fct
         WRITE(numout,*) '            horizontal 2nd/4th order               nn_fct_h   = ', nn_fct_h
         WRITE(numout,*) '            vertical   2nd/4th order               nn_fct_v   = ', nn_fct_v
         WRITE(numout,*) '            2nd order + vertical sub-timestepping  nn_fct_zts = ', nn_fct_zts
         WRITE(numout,*) '      MUSCL scheme                              ln_trcadv_mus = ', ln_trcadv_mus
         WRITE(numout,*) '            + upstream scheme near river mouths    ln_mus_ups = ', ln_mus_ups
         WRITE(numout,*) '      UBS scheme                                ln_trcadv_ubs = ', ln_trcadv_ubs
         WRITE(numout,*) '            vertical   2nd/4th order               nn_ubs_v   = ', nn_ubs_v
         WRITE(numout,*) '      QUICKEST scheme                           ln_trcadv_qck = ', ln_trcadv_qck
      ENDIF
      !

      ioptio = 0                       !==  Parameter control  ==!
      IF( ln_trcadv_cen )   ioptio = ioptio + 1
      IF( ln_trcadv_fct )   ioptio = ioptio + 1
      IF( ln_trcadv_mus )   ioptio = ioptio + 1
      IF( ln_trcadv_ubs )   ioptio = ioptio + 1
      IF( ln_trcadv_qck )   ioptio = ioptio + 1

      !
      IF( ioptio == 0 ) THEN
         nadv = np_NO_adv
         CALL ctl_warn( 'trc_adv_init: You are running without tracer advection.' )
      ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE advection scheme in namelist namtrc_adv' )
      !
      IF( ln_trcadv_cen .AND. ( nn_cen_h /= 2 .AND. nn_cen_h /= 4 )   &
                        .AND. ( nn_cen_v /= 2 .AND. nn_cen_v /= 4 )   ) THEN
        CALL ctl_stop( 'trc_adv_init: CEN scheme, choose 2nd or 4th order' )
      ENDIF
      IF( ln_trcadv_fct .AND. ( nn_fct_h /= 2 .AND. nn_fct_h /= 4 )   &
                        .AND. ( nn_fct_v /= 2 .AND. nn_fct_v /= 4 )   ) THEN
        CALL ctl_stop( 'trc_adv_init: FCT scheme, choose 2nd or 4th order' )
      ENDIF
      IF( ln_trcadv_fct .AND. nn_fct_zts > 0 ) THEN
         IF( nn_fct_h == 4 ) THEN
            nn_fct_h = 2
            CALL ctl_stop( 'trc_adv_init: force 2nd order FCT scheme, 4th order does not exist with sub-timestepping' )
         ENDIF
         IF( .NOT.ln_linssh ) THEN
            CALL ctl_stop( 'trc_adv_init: vertical sub-timestepping not allow in non-linear free surface' )
         ENDIF
         IF( nn_fct_zts == 1 )   CALL ctl_warn( 'trc_adv_init: FCT with ONE sub-timestep = FCT without sub-timestep' )
      ENDIF
      IF( ln_trcadv_ubs .AND. ( nn_ubs_v /= 2 .AND. nn_ubs_v /= 4 )   ) THEN
        CALL ctl_stop( 'trc_adv_init: UBS scheme, choose 2nd or 4th order' )
      ENDIF
      IF( ln_trcadv_ubs .AND. nn_ubs_v == 4 ) THEN
         CALL ctl_warn( 'trc_adv_init: UBS scheme, only 2nd FCT scheme available on the vertical. It will be used' )
      ENDIF
      IF( ln_isfcav ) THEN                                                       ! ice-shelf cavities
         IF(  ln_trcadv_cen .AND. nn_cen_v /= 4    .OR.   &                            ! NO 4th order with ISF
            & ln_trcadv_fct .AND. nn_fct_v /= 4   )   CALL ctl_stop( 'tra_adv_init: 4th order COMPACT scheme not allowed with ISF' )
      ENDIF
      !
      !                                !==  used advection scheme  ==!
      !                                      ! set nadv
      IF( ln_trcadv_cen                      )   nadv = np_CEN
      IF( ln_trcadv_fct                      )   nadv = np_FCT
      IF( ln_trcadv_fct .AND. nn_fct_zts > 0 )   nadv = np_FCT_zts
      IF( ln_trcadv_mus                      )   nadv = np_MUS
      IF( ln_trcadv_ubs                      )   nadv = np_UBS
      IF( ln_trcadv_qck                      )   nadv = np_QCK
      !
      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         IF( nadv == np_NO_adv  )   WRITE(numout,*) '         NO passive tracer advection'
         IF( nadv == np_CEN     )   WRITE(numout,*) '         CEN      scheme is used. Horizontal order: ', nn_cen_h,   &
            &                                                                        ' Vertical   order: ', nn_cen_v
         IF( nadv == np_FCT     )   WRITE(numout,*) '         FCT      scheme is used. Horizontal order: ', nn_fct_h,   &
            &                                                                       ' Vertical   order: ', nn_fct_v
         IF( nadv == np_FCT_zts )   WRITE(numout,*) '         use 2nd order FCT with ', nn_fct_zts,'vertical sub-timestepping'
         IF( nadv == np_MUS     )   WRITE(numout,*) '         MUSCL    scheme is used'
         IF( nadv == np_UBS     )   WRITE(numout,*) '         UBS      scheme is used'
         IF( nadv == np_QCK     )   WRITE(numout,*) '         QUICKEST scheme is used'
      ENDIF
      !
   END SUBROUTINE trc_adv_ini
   
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_adv( kt )
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_adv: You should not have seen this print! error?', kt
   END SUBROUTINE trc_adv
#endif

  !!======================================================================
END MODULE trcadv
