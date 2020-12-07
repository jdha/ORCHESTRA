MODULE traadv_fct
   !!==============================================================================
   !!                       ***  MODULE  traadv_fct  ***
   !! Ocean  tracers:  horizontal & vertical advective trend (2nd/4th order Flux Corrected Transport method)
   !!==============================================================================
   !! History :  3.7  !  2015-09  (L. Debreu, G. Madec)  original code (inspired from traadv_tvd.F90)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  tra_adv_fct    : update the tracer trend with a 3D advective trends using a 2nd or 4th order FCT scheme
   !!  tra_adv_fct_zts: update the tracer trend with a 3D advective trends using a 2nd order FCT scheme 
   !!                   with sub-time-stepping in the vertical direction
   !!  nonosc         : compute monotonic tracer fluxes by a non-oscillatory algorithm 
   !!  interp_4th_cpt : 4th order compact scheme for the vertical component of the advection
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE trc_oce        ! share passive tracers/Ocean variables
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! tracers trends
   USE diaptr         ! poleward transport diagnostics
   USE diaar5         ! AR5 diagnostics
   USE phycst, ONLY: rau0_rcp
   !
   USE in_out_manager ! I/O manager
   USE iom
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary condition (or mpp link) 
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_fct        ! routine called by traadv.F90
   PUBLIC   tra_adv_fct_zts    ! routine called by traadv.F90
   PUBLIC   interp_4th_cpt     ! routine called by traadv_cen.F90

   LOGICAL  ::   l_trd   ! flag to compute trends
   LOGICAL  ::   l_ptr   ! flag to compute poleward transport
   LOGICAL  ::   l_hst   ! flag to compute heat/salt transport
   REAL(wp) ::   r1_6 = 1._wp / 6._wp   ! =1/6

   !                                        ! tridiag solver associated indices:
   INTEGER, PARAMETER ::   np_NH   = 0   ! Neumann homogeneous boundary condition
   INTEGER, PARAMETER ::   np_CEN2 = 1   ! 2nd order centered  boundary condition

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: traadv_fct.F90 7753 2017-03-03 11:46:59Z mocavero $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_fct( kt, kit000, cdtype, p2dt, pun, pvn, pwn,       &
      &                                              ptb, ptn, pta, kjpt, kn_fct_h, kn_fct_v )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_fct  ***
      !! 
      !! **  Purpose :   Compute the now trend due to total advection of tracers
      !!               and add it to the general trend of tracer equations
      !!
      !! **  Method  : - 2nd or 4th FCT scheme on the horizontal direction
      !!               (choice through the value of kn_fct)
      !!               - on the vertical the 4th order is a compact scheme 
      !!               - corrected flux (monotonic correction) 
      !!
      !! ** Action : - update pta  with the now advective tracer trends
      !!             - send trends to trdtra module for further diagnostcs (l_trdtra=T)
      !!             - htr_adv, str_adv : poleward advective heat and salt transport (ln_diaptr=T)
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
      INTEGER                              , INTENT(in   ) ::   kn_fct_h        ! order of the FCT scheme (=2 or 4)
      INTEGER                              , INTENT(in   ) ::   kn_fct_v        ! order of the FCT scheme (=2 or 4)
      REAL(wp)                             , INTENT(in   ) ::   p2dt            ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb, ptn        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk, jn                           ! dummy loop indices  
      REAL(wp) ::   ztra                                     ! local scalar
      REAL(wp) ::   zfp_ui, zfp_vj, zfp_wk, zC2t_u, zC4t_u   !   -      -
      REAL(wp) ::   zfm_ui, zfm_vj, zfm_wk, zC2t_v, zC4t_v   !   -      -
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwi, zwx, zwy, zwz, ztu, ztv, zltu, zltv, ztw
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ztrdx, ztrdy, ztrdz, zptry
      REAL(wp), POINTER, DIMENSION(:,:)   :: z2d
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_adv_fct')
      !
      CALL wrk_alloc( jpi,jpj,jpk,   zwi, zwx, zwy, zwz, ztu, ztv, zltu, zltv, ztw )
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_fct : FCT advection scheme on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      l_trd = .FALSE.
      l_hst = .FALSE.
      l_ptr = .FALSE.
      IF( ( cdtype == 'TRA'   .AND. l_trdtra ) .OR. ( cdtype == 'TRC' .AND. l_trdtrc ) )     l_trd = .TRUE.
      IF(   cdtype == 'TRA'   .AND. ln_diaptr )                                              l_ptr = .TRUE. 
      IF(   cdtype == 'TRA' .AND. ( iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. &
         &                          iom_use("uadv_salttr") .OR. iom_use("vadv_salttr")  ) )  l_hst = .TRUE.
      !
      IF( l_trd .OR. l_hst )  THEN
         CALL wrk_alloc( jpi, jpj, jpk, ztrdx, ztrdy, ztrdz )
         ztrdx(:,:,:) = 0._wp   ;    ztrdy(:,:,:) = 0._wp   ;   ztrdz(:,:,:) = 0._wp
      ENDIF
      !
      IF( l_ptr ) THEN  
         CALL wrk_alloc( jpi, jpj, jpk, zptry )
         zptry(:,:,:) = 0._wp
      ENDIF
      !                          ! surface & bottom value : flux set to zero one for all
      zwz(:,:, 1 ) = 0._wp            
      zwx(:,:,jpk) = 0._wp   ;   zwy(:,:,jpk) = 0._wp    ;    zwz(:,:,jpk) = 0._wp
      !
      zwi(:,:,:) = 0._wp        
      !
      DO jn = 1, kjpt            !==  loop over the tracers  ==!
         !
         !        !==  upstream advection with initial mass fluxes & intermediate update  ==!
         !                    !* upstream tracer flux in the i and j direction 
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  ! upstream scheme
                  zfp_ui = pun(ji,jj,jk) + ABS( pun(ji,jj,jk) )
                  zfm_ui = pun(ji,jj,jk) - ABS( pun(ji,jj,jk) )
                  zfp_vj = pvn(ji,jj,jk) + ABS( pvn(ji,jj,jk) )
                  zfm_vj = pvn(ji,jj,jk) - ABS( pvn(ji,jj,jk) )
                  zwx(ji,jj,jk) = 0.5 * ( zfp_ui * ptb(ji,jj,jk,jn) + zfm_ui * ptb(ji+1,jj  ,jk,jn) )
                  zwy(ji,jj,jk) = 0.5 * ( zfp_vj * ptb(ji,jj,jk,jn) + zfm_vj * ptb(ji  ,jj+1,jk,jn) )
               END DO
            END DO
         END DO
         !                    !* upstream tracer flux in the k direction *!
         DO jk = 2, jpkm1        ! Interior value ( multiplied by wmask)
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zfp_wk = pwn(ji,jj,jk) + ABS( pwn(ji,jj,jk) )
                  zfm_wk = pwn(ji,jj,jk) - ABS( pwn(ji,jj,jk) )
                  zwz(ji,jj,jk) = 0.5 * ( zfp_wk * ptb(ji,jj,jk,jn) + zfm_wk * ptb(ji,jj,jk-1,jn) ) * wmask(ji,jj,jk)
               END DO
            END DO
         END DO
         IF( ln_linssh ) THEN    ! top ocean value (only in linear free surface as zwz has been w-masked)
            IF( ln_isfcav ) THEN             ! top of the ice-shelf cavities and at the ocean surface
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zwz(ji,jj, mikt(ji,jj) ) = pwn(ji,jj,mikt(ji,jj)) * ptb(ji,jj,mikt(ji,jj),jn)   ! linear free surface 
                  END DO
               END DO   
            ELSE                             ! no cavities: only at the ocean surface
               zwz(:,:,1) = pwn(:,:,1) * ptb(:,:,1,jn)
            ENDIF
         ENDIF
         !               
         DO jk = 1, jpkm1     !* trend and after field with monotonic scheme
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  !                             ! total intermediate advective trends
                  ztra = - (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                     &      + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )   &
                     &      + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1) ) * r1_e1e2t(ji,jj)
                  !                             ! update and guess with monotonic sheme
                  pta(ji,jj,jk,jn) =                     pta(ji,jj,jk,jn) +        ztra   / e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                  zwi(ji,jj,jk)    = ( e3t_b(ji,jj,jk) * ptb(ji,jj,jk,jn) + p2dt * ztra ) / e3t_a(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL lbc_lnk( zwi, 'T', 1. )  ! Lateral boundary conditions on zwi  (unchanged sign)
         !                
         IF( l_trd .OR. l_hst )  THEN             ! trend diagnostics (contribution of upstream fluxes)
            ztrdx(:,:,:) = zwx(:,:,:)   ;   ztrdy(:,:,:) = zwy(:,:,:)   ;   ztrdz(:,:,:) = zwz(:,:,:)
         END IF
         !                             ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( l_ptr )  zptry(:,:,:) = zwy(:,:,:) 
         !
         !        !==  anti-diffusive flux : high order minus low order  ==!
         !
         SELECT CASE( kn_fct_h )    !* horizontal anti-diffusive fluxes
         !
         CASE(  2  )                   !- 2nd order centered
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     zwx(ji,jj,jk) = 0.5_wp * pun(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji+1,jj,jk,jn) ) - zwx(ji,jj,jk)
                     zwy(ji,jj,jk) = 0.5_wp * pvn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji,jj+1,jk,jn) ) - zwy(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
         CASE(  4  )                   !- 4th order centered
            zltu(:,:,jpk) = 0._wp            ! Bottom value : flux set to zero
            zltv(:,:,jpk) = 0._wp
            DO jk = 1, jpkm1                 ! Laplacian
               DO jj = 1, jpjm1                    ! 1st derivative (gradient)
                  DO ji = 1, fs_jpim1   ! vector opt.
                     ztu(ji,jj,jk) = ( ptn(ji+1,jj  ,jk,jn) - ptn(ji,jj,jk,jn) ) * umask(ji,jj,jk)
                     ztv(ji,jj,jk) = ( ptn(ji  ,jj+1,jk,jn) - ptn(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
                  END DO
               END DO
               DO jj = 2, jpjm1                    ! 2nd derivative * 1/ 6
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     zltu(ji,jj,jk) = (  ztu(ji,jj,jk) + ztu(ji-1,jj,jk)  ) * r1_6
                     zltv(ji,jj,jk) = (  ztv(ji,jj,jk) + ztv(ji,jj-1,jk)  ) * r1_6
                  END DO
               END DO
            END DO
            CALL lbc_lnk( zltu, 'T', 1. )   ;    CALL lbc_lnk( zltv, 'T', 1. )   ! Lateral boundary cond. (unchanged sgn)
            !
            DO jk = 1, jpkm1                 ! Horizontal advective fluxes
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     zC2t_u = ptn(ji,jj,jk,jn) + ptn(ji+1,jj  ,jk,jn)   ! 2 x C2 interpolation of T at u- & v-points
                     zC2t_v = ptn(ji,jj,jk,jn) + ptn(ji  ,jj+1,jk,jn)
                     !                                                  ! C4 minus upstream advective fluxes 
                     zwx(ji,jj,jk) =  0.5_wp * pun(ji,jj,jk) * ( zC2t_u + zltu(ji,jj,jk) - zltu(ji+1,jj,jk) ) - zwx(ji,jj,jk)
                     zwy(ji,jj,jk) =  0.5_wp * pvn(ji,jj,jk) * ( zC2t_v + zltv(ji,jj,jk) - zltv(ji,jj+1,jk) ) - zwy(ji,jj,jk)
                  END DO
               END DO
            END DO         
            !
         CASE(  41 )                   !- 4th order centered       ==>>   !!gm coding attempt   need to be tested
            ztu(:,:,jpk) = 0._wp             ! Bottom value : flux set to zero
            ztv(:,:,jpk) = 0._wp
            DO jk = 1, jpkm1                 ! 1st derivative (gradient)
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     ztu(ji,jj,jk) = ( ptn(ji+1,jj  ,jk,jn) - ptn(ji,jj,jk,jn) ) * umask(ji,jj,jk)
                     ztv(ji,jj,jk) = ( ptn(ji  ,jj+1,jk,jn) - ptn(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            CALL lbc_lnk( ztu, 'U', -1. )   ;    CALL lbc_lnk( ztv, 'V', -1. )   ! Lateral boundary cond. (unchanged sgn)
            !
            DO jk = 1, jpkm1                 ! Horizontal advective fluxes
               DO jj = 2, jpjm1
                  DO ji = 2, fs_jpim1   ! vector opt.
                     zC2t_u = ptn(ji,jj,jk,jn) + ptn(ji+1,jj  ,jk,jn)   ! 2 x C2 interpolation of T at u- & v-points (x2)
                     zC2t_v = ptn(ji,jj,jk,jn) + ptn(ji  ,jj+1,jk,jn)
                     !                                                  ! C4 interpolation of T at u- & v-points (x2)
                     zC4t_u =  zC2t_u + r1_6 * ( ztu(ji-1,jj  ,jk) - ztu(ji+1,jj  ,jk) )
                     zC4t_v =  zC2t_v + r1_6 * ( ztv(ji  ,jj-1,jk) - ztv(ji  ,jj+1,jk) )
                     !                                                  ! C4 minus upstream advective fluxes 
                     zwx(ji,jj,jk) =  0.5_wp * pun(ji,jj,jk) * zC4t_u - zwx(ji,jj,jk)
                     zwy(ji,jj,jk) =  0.5_wp * pvn(ji,jj,jk) * zC4t_v - zwy(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
         END SELECT
         !                      
         SELECT CASE( kn_fct_v )    !* vertical anti-diffusive fluxes (w-masked interior values)
         !
         CASE(  2  )                   !- 2nd order centered
            DO jk = 2, jpkm1    
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zwz(ji,jj,jk) =  (  pwn(ji,jj,jk) * 0.5_wp * ( ptn(ji,jj,jk,jn) + ptn(ji,jj,jk-1,jn) )   &
                        &              - zwz(ji,jj,jk)  ) * wmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
         CASE(  4  )                   !- 4th order COMPACT
            CALL interp_4th_cpt( ptn(:,:,:,jn) , ztw )   ! zwt = COMPACT interpolation of T at w-point
            DO jk = 2, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zwz(ji,jj,jk) = ( pwn(ji,jj,jk) * ztw(ji,jj,jk) - zwz(ji,jj,jk) ) * wmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
         END SELECT
         IF( ln_linssh ) THEN    ! top ocean value: high order = upstream  ==>>  zwz=0
            zwz(:,:,1) = 0._wp   ! only ocean surface as interior zwz values have been w-masked
         ENDIF
         !
         CALL lbc_lnk( zwx, 'U', -1. )   ;   CALL lbc_lnk( zwy, 'V', -1. )         ! Lateral bondary conditions
         CALL lbc_lnk( zwz, 'W',  1. )
         !
         !        !==  monotonicity algorithm  ==!
         !
         CALL nonosc( ptb(:,:,:,jn), zwx, zwy, zwz, zwi, p2dt )
         !
         !        !==  final trend with corrected fluxes  ==!
         !
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.  
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) - (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                     &                                   + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )   &
                     &                                   + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1) ) &
                     &                                * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
               END DO
            END DO
         END DO
         !
         IF( l_trd .OR. l_hst ) THEN     ! trend diagnostics (contribution of upstream fluxes)
            ztrdx(:,:,:) = ztrdx(:,:,:) + zwx(:,:,:)  ! <<< Add to previously computed
            ztrdy(:,:,:) = ztrdy(:,:,:) + zwy(:,:,:)  ! <<< Add to previously computed
            ztrdz(:,:,:) = ztrdz(:,:,:) + zwz(:,:,:)  ! <<< Add to previously computed
         ENDIF
            !
         IF( l_trd ) THEN 
            CALL trd_tra( kt, cdtype, jn, jptra_xad, ztrdx, pun, ptn(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_yad, ztrdy, pvn, ptn(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_zad, ztrdz, pwn, ptn(:,:,:,jn) )
            !
         END IF
         !                                !  heat/salt transport
         IF( l_hst )  CALL dia_ar5_hst( jn, 'adv', ztrdx(:,:,:), ztrdy(:,:,:) )

         !                                ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( l_ptr ) THEN  
            zptry(:,:,:) = zptry(:,:,:) + zwy(:,:,:)  ! <<< Add to previously computed
            CALL dia_ptr_hst( jn, 'adv', zptry(:,:,:) )
         ENDIF
         !
      END DO                     ! end of tracer loop
      !
                              CALL wrk_dealloc( jpi,jpj,jpk,    zwi, zwx, zwy, zwz, ztu, ztv, zltu, zltv, ztw )
      IF( l_trd .OR. l_hst )  CALL wrk_dealloc( jpi, jpj, jpk, ztrdx, ztrdy, ztrdz )
      IF( l_ptr )             CALL wrk_dealloc( jpi, jpj, jpk, zptry )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_adv_fct')
      !
   END SUBROUTINE tra_adv_fct


   SUBROUTINE tra_adv_fct_zts( kt, kit000, cdtype, p2dt, pun, pvn, pwn,      &
      &                                                  ptb, ptn, pta, kjpt, kn_fct_zts )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_fct_zts  ***
      !! 
      !! **  Purpose :   Compute the now trend due to total advection of 
      !!       tracers and add it to the general trend of tracer equations
      !!
      !! **  Method  :   TVD ZTS scheme, i.e. 2nd order centered scheme with
      !!       corrected flux (monotonic correction). This version use sub-
      !!       timestepping for the vertical advection which increases stability
      !!       when vertical metrics are small.
      !!       note: - this advection scheme needs a leap-frog time scheme
      !!
      !! ** Action : - update (pta) with the now advective tracer trends
      !!             - save the trends 
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
      INTEGER                              , INTENT(in   ) ::   kn_fct_zts      ! number of number of vertical sub-timesteps
      REAL(wp)                             , INTENT(in   ) ::   p2dt            ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb, ptn        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
      !
      REAL(wp), DIMENSION( jpk )                           ::   zts             ! length of sub-timestep for vertical advection
      REAL(wp)                                             ::   zr_p2dt         ! reciprocal of tracer timestep
      INTEGER  ::   ji, jj, jk, jl, jn       ! dummy loop indices  
      INTEGER  ::   jtb, jtn, jta   ! sub timestep pointers for leap-frog/euler forward steps
      INTEGER  ::   jtaken          ! toggle for collecting appropriate fluxes from sub timesteps
      REAL(wp) ::   z_rzts          ! Fractional length of Euler forward sub-timestep for vertical advection
      REAL(wp) ::   ztra            ! local scalar
      REAL(wp) ::   zfp_ui, zfp_vj, zfp_wk   !   -      -
      REAL(wp) ::   zfm_ui, zfm_vj, zfm_wk   !   -      -
      REAL(wp), POINTER, DIMENSION(:,:  )   ::   zwx_sav , zwy_sav
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zwi, zwx, zwy, zwz, zhdiv, zwzts, zwz_sav
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   ztrdx, ztrdy, ztrdz
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zptry
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   ztrs
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_adv_fct_zts')
      !
      CALL wrk_alloc( jpi,jpj,             zwx_sav, zwy_sav )
      CALL wrk_alloc( jpi,jpj,jpk,         zwx, zwy, zwz, zwi, zhdiv, zwzts, zwz_sav )
      CALL wrk_alloc( jpi,jpj,jpk,kjpt+1,  ztrs )
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_fct_zts : 2nd order FCT scheme with ', kn_fct_zts, ' vertical sub-timestep on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      l_trd = .FALSE.
      l_hst = .FALSE.
      l_ptr = .FALSE.
      IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR. ( cdtype == 'TRC' .AND. l_trdtrc ) )      l_trd = .TRUE.
      IF(   cdtype == 'TRA' .AND. ln_diaptr )                                               l_ptr = .TRUE. 
      IF(   cdtype == 'TRA' .AND. ( iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. &
         &                          iom_use("uadv_salttr") .OR. iom_use("vadv_salttr")  ) ) l_hst = .TRUE.
      !
      IF( l_trd .OR. l_hst )  THEN
         CALL wrk_alloc( jpi,jpj,jpk,   ztrdx, ztrdy, ztrdz )
         ztrdx(:,:,:) = 0._wp  ;    ztrdy(:,:,:) = 0._wp  ;   ztrdz(:,:,:) = 0._wp
      ENDIF
      !
      IF( l_ptr ) THEN  
         CALL wrk_alloc( jpi, jpj,jpk, zptry )
         zptry(:,:,:) = 0._wp
      ENDIF
      zwi(:,:,:) = 0._wp
      z_rzts = 1._wp / REAL( kn_fct_zts, wp )
      zr_p2dt = 1._wp / p2dt
      !
      ! surface & Bottom value : flux set to zero for all tracers
      zwz(:,:, 1 ) = 0._wp
      zwx(:,:,jpk) = 0._wp   ;    zwz(:,:,jpk) = 0._wp
      zwy(:,:,jpk) = 0._wp   ;    zwi(:,:,jpk) = 0._wp
      !
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         !
         ! Upstream advection with initial mass fluxes & intermediate update
         DO jk = 1, jpkm1        ! upstream tracer flux in the i and j direction
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  ! upstream scheme
                  zfp_ui = pun(ji,jj,jk) + ABS( pun(ji,jj,jk) )
                  zfm_ui = pun(ji,jj,jk) - ABS( pun(ji,jj,jk) )
                  zfp_vj = pvn(ji,jj,jk) + ABS( pvn(ji,jj,jk) )
                  zfm_vj = pvn(ji,jj,jk) - ABS( pvn(ji,jj,jk) )
                  zwx(ji,jj,jk) = 0.5_wp * ( zfp_ui * ptb(ji,jj,jk,jn) + zfm_ui * ptb(ji+1,jj  ,jk,jn) )
                  zwy(ji,jj,jk) = 0.5_wp * ( zfp_vj * ptb(ji,jj,jk,jn) + zfm_vj * ptb(ji  ,jj+1,jk,jn) )
               END DO
            END DO
         END DO
         !                       ! upstream tracer flux in the k direction
         DO jk = 2, jpkm1              ! Interior value
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zfp_wk = pwn(ji,jj,jk) + ABS( pwn(ji,jj,jk) )
                  zfm_wk = pwn(ji,jj,jk) - ABS( pwn(ji,jj,jk) )
                  zwz(ji,jj,jk) = 0.5_wp * ( zfp_wk * ptb(ji,jj,jk,jn) + zfm_wk * ptb(ji,jj,jk-1,jn) ) * wmask(ji,jj,jk)
               END DO
            END DO
         END DO
         IF( ln_linssh ) THEN          ! top value : linear free surface case only (as zwz is multiplied by wmask)
            IF( ln_isfcav ) THEN             ! ice-shelf cavities: top value
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zwz(ji,jj, mikt(ji,jj) ) = pwn(ji,jj,mikt(ji,jj)) * ptb(ji,jj,mikt(ji,jj),jn) 
                  END DO
               END DO   
            ELSE                             ! no cavities, surface value
               zwz(:,:,1) = pwn(:,:,1) * ptb(:,:,1,jn)
            ENDIF
         ENDIF
         !
         DO jk = 1, jpkm1         ! total advective trend
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  !                             ! total intermediate advective trends
                  ztra = - (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                     &      + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )   &
                     &      + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1)   ) * r1_e1e2t(ji,jj)
                  !                             ! update and guess with monotonic sheme
                  pta(ji,jj,jk,jn) =                     pta(ji,jj,jk,jn) +        ztra   / e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                  zwi(ji,jj,jk)    = ( e3t_b(ji,jj,jk) * ptb(ji,jj,jk,jn) + p2dt * ztra ) / e3t_a(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !                           
         CALL lbc_lnk( zwi, 'T', 1. )     ! Lateral boundary conditions on zwi  (unchanged sign)
         !                
         IF( l_trd .OR. l_hst )  THEN                ! trend diagnostics (contribution of upstream fluxes)
            ztrdx(:,:,:) = zwx(:,:,:)   ;    ztrdy(:,:,:) = zwy(:,:,:)  ;   ztrdz(:,:,:) = zwz(:,:,:)
         END IF
         !                                ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( l_ptr )  zptry(:,:,:) = zwy(:,:,:)

         ! 3. anti-diffusive flux : high order minus low order
         ! ---------------------------------------------------

         DO jk = 1, jpkm1                    !* horizontal anti-diffusive fluxes
            !
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zwx_sav(ji,jj) = zwx(ji,jj,jk)
                  zwy_sav(ji,jj) = zwy(ji,jj,jk)
                  !
                  zwx(ji,jj,jk) = 0.5_wp * pun(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji+1,jj,jk,jn) )
                  zwy(ji,jj,jk) = 0.5_wp * pvn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji,jj+1,jk,jn) )
               END DO
            END DO
            !
            DO jj = 2, jpjm1                    ! partial horizontal divergence
               DO ji = fs_2, fs_jpim1
                  zhdiv(ji,jj,jk) = (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk)   &
                     &               + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk)  )
               END DO
            END DO
            !
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zwx(ji,jj,jk) = zwx(ji,jj,jk) - zwx_sav(ji,jj)
                  zwy(ji,jj,jk) = zwy(ji,jj,jk) - zwy_sav(ji,jj)
               END DO
            END DO
         END DO
         !
         !                                !* vertical anti-diffusive flux
         zwz_sav(:,:,:)   = zwz(:,:,:)
         ztrs   (:,:,:,1) = ptb(:,:,:,jn)
         ztrs   (:,:,1,2) = ptb(:,:,1,jn)
         ztrs   (:,:,1,3) = ptb(:,:,1,jn)
         zwzts  (:,:,:)   = 0._wp
         !
         DO jl = 1, kn_fct_zts                  ! Start of sub timestepping loop
            !
            IF( jl == 1 ) THEN                        ! Euler forward to kick things off
               jtb = 1   ;   jtn = 1   ;   jta = 2
               zts(:) = p2dt * z_rzts
               jtaken = MOD( kn_fct_zts + 1 , 2)            ! Toggle to collect every second flux
               !                                            ! starting at jl =1 if kn_fct_zts is odd; 
               !                                            ! starting at jl =2 otherwise
            ELSEIF( jl == 2 ) THEN                    ! First leapfrog step
               jtb = 1   ;   jtn = 2   ;   jta = 3
               zts(:) = 2._wp * p2dt * z_rzts
            ELSE                                      ! Shuffle pointers for subsequent leapfrog steps
               jtb = MOD(jtb,3) + 1
               jtn = MOD(jtn,3) + 1
               jta = MOD(jta,3) + 1
            ENDIF
            DO jk = 2, jpkm1                          ! interior value
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zwz(ji,jj,jk) = 0.5_wp * pwn(ji,jj,jk) * ( ztrs(ji,jj,jk,jtn) + ztrs(ji,jj,jk-1,jtn) ) * wmask(ji,jj,jk)
                     IF( jtaken == 0 )   zwzts(ji,jj,jk) = zwzts(ji,jj,jk) + zwz(ji,jj,jk) * zts(jk)    ! Accumulate time-weighted vertcal flux
                  END DO
               END DO
            END DO
            IF( ln_linssh ) THEN                    ! top value (only in linear free surface case)
               IF( ln_isfcav ) THEN                      ! ice-shelf cavities
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        zwz(ji,jj, mikt(ji,jj) ) = pwn(ji,jj,mikt(ji,jj)) * ptb(ji,jj,mikt(ji,jj),jn)   ! linear free surface 
                     END DO
                  END DO   
               ELSE                                      ! no ocean cavities
                  zwz(:,:,1) = pwn(:,:,1) * ptb(:,:,1,jn)
               ENDIF
            ENDIF
            !
            jtaken = MOD( jtaken + 1 , 2 )
            !
            DO jk = 2, jpkm1                             ! total advective trends
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ztrs(ji,jj,jk,jta) = ztrs(ji,jj,jk,jtb)                                                 &
                        &               - zts(jk) * (  zhdiv(ji,jj,jk) + zwz(ji,jj,jk) - zwz(ji,jj,jk+1) )   &
                        &                         * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
         END DO

         DO jk = 2, jpkm1          ! Anti-diffusive vertical flux using average flux from the sub-timestepping
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  zwz(ji,jj,jk) = ( zwzts(ji,jj,jk) * zr_p2dt - zwz_sav(ji,jj,jk) ) * wmask(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL lbc_lnk( zwx, 'U', -1. )   ;   CALL lbc_lnk( zwy, 'V', -1. )         ! Lateral bondary conditions
         CALL lbc_lnk( zwz, 'W',  1. )

         ! 4. monotonicity algorithm
         ! -------------------------
         CALL nonosc( ptb(:,:,:,jn), zwx, zwy, zwz, zwi, p2dt )


         ! 5. final trend with corrected fluxes
         ! ------------------------------------
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.  
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + (   zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )       &
                     &                                    + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1)   )   &
                     &                                * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
               END DO
            END DO
         END DO

        !
         IF( l_trd .OR. l_hst ) THEN     ! trend diagnostics (contribution of upstream fluxes)
            ztrdx(:,:,:) = ztrdx(:,:,:) + zwx(:,:,:)  ! <<< Add to previously computed
            ztrdy(:,:,:) = ztrdy(:,:,:) + zwy(:,:,:)  ! <<< Add to previously computed
            ztrdz(:,:,:) = ztrdz(:,:,:) + zwz(:,:,:)  ! <<< Add to previously computed
         ENDIF
            !
         IF( l_trd ) THEN 
            CALL trd_tra( kt, cdtype, jn, jptra_xad, ztrdx, pun, ptn(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_yad, ztrdy, pvn, ptn(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_zad, ztrdz, pwn, ptn(:,:,:,jn) )
            !
         END IF
         !                                             ! heat/salt transport
         IF( l_hst )  CALL dia_ar5_hst( jn, 'adv', ztrdx(:,:,:), ztrdy(:,:,:) )

         !                                            ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( l_ptr ) THEN  
            zptry(:,:,:) = zptry(:,:,:) + zwy(:,:,:)  ! <<< Add to previously computed
            CALL dia_ptr_hst( jn, 'adv', zptry(:,:,:) )
         ENDIF
         !
      END DO
      !
                              CALL wrk_alloc( jpi,jpj,             zwx_sav, zwy_sav )
                              CALL wrk_alloc( jpi,jpj, jpk,        zwx, zwy, zwz, zwi, zhdiv, zwzts, zwz_sav )
                              CALL wrk_alloc( jpi,jpj,jpk,kjpt+1,  ztrs )
      IF( l_trd .OR. l_hst )  CALL wrk_dealloc( jpi, jpj, jpk, ztrdx, ztrdy, ztrdz )
      IF( l_ptr )             CALL wrk_dealloc( jpi, jpj, jpk, zptry )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_adv_fct_zts')
      !
   END SUBROUTINE tra_adv_fct_zts


   SUBROUTINE nonosc( pbef, paa, pbb, pcc, paft, p2dt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc  ***
      !!     
      !! **  Purpose :   compute monotonic tracer fluxes from the upstream 
      !!       scheme and the before field by a nonoscillatory algorithm 
      !!
      !! **  Method  :   ... ???
      !!       warning : pbef and paft must be masked, but the boundaries
      !!       conditions on the fluxes are not necessary zalezak (1979)
      !!       drange (1995) multi-dimensional forward-in-time and upstream-
      !!       in-space based differencing for fluid
      !!----------------------------------------------------------------------
      REAL(wp)                         , INTENT(in   ) ::   p2dt            ! tracer time-step
      REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT(in   ) ::   pbef, paft      ! before & after field
      REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT(inout) ::   paa, pbb, pcc   ! monotonic fluxes in the 3 directions
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikm1         ! local integer
      REAL(wp) ::   zpos, zneg, zbt, za, zb, zc, zbig, zrtrn    ! local scalars
      REAL(wp) ::   zau, zbu, zcu, zav, zbv, zcv, zup, zdo            !   -      -
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zbetup, zbetdo, zbup, zbdo
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('nonosc')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zbetup, zbetdo, zbup, zbdo )
      !
      zbig  = 1.e+40_wp
      zrtrn = 1.e-15_wp
      zbetup(:,:,:) = 0._wp   ;   zbetdo(:,:,:) = 0._wp

      ! Search local extrema
      ! --------------------
      ! max/min of pbef & paft with large negative/positive value (-/+zbig) inside land
      zbup = MAX( pbef * tmask - zbig * ( 1._wp - tmask ),   &
         &        paft * tmask - zbig * ( 1._wp - tmask )  )
      zbdo = MIN( pbef * tmask + zbig * ( 1._wp - tmask ),   &
         &        paft * tmask + zbig * ( 1._wp - tmask )  )

      DO jk = 1, jpkm1
         ikm1 = MAX(jk-1,1)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.

               ! search maximum in neighbourhood
               zup = MAX(  zbup(ji  ,jj  ,jk  ),   &
                  &        zbup(ji-1,jj  ,jk  ), zbup(ji+1,jj  ,jk  ),   &
                  &        zbup(ji  ,jj-1,jk  ), zbup(ji  ,jj+1,jk  ),   &
                  &        zbup(ji  ,jj  ,ikm1), zbup(ji  ,jj  ,jk+1)  )

               ! search minimum in neighbourhood
               zdo = MIN(  zbdo(ji  ,jj  ,jk  ),   &
                  &        zbdo(ji-1,jj  ,jk  ), zbdo(ji+1,jj  ,jk  ),   &
                  &        zbdo(ji  ,jj-1,jk  ), zbdo(ji  ,jj+1,jk  ),   &
                  &        zbdo(ji  ,jj  ,ikm1), zbdo(ji  ,jj  ,jk+1)  )

               ! positive part of the flux
               zpos = MAX( 0., paa(ji-1,jj  ,jk  ) ) - MIN( 0., paa(ji  ,jj  ,jk  ) )   &
                  & + MAX( 0., pbb(ji  ,jj-1,jk  ) ) - MIN( 0., pbb(ji  ,jj  ,jk  ) )   &
                  & + MAX( 0., pcc(ji  ,jj  ,jk+1) ) - MIN( 0., pcc(ji  ,jj  ,jk  ) )

               ! negative part of the flux
               zneg = MAX( 0., paa(ji  ,jj  ,jk  ) ) - MIN( 0., paa(ji-1,jj  ,jk  ) )   &
                  & + MAX( 0., pbb(ji  ,jj  ,jk  ) ) - MIN( 0., pbb(ji  ,jj-1,jk  ) )   &
                  & + MAX( 0., pcc(ji  ,jj  ,jk  ) ) - MIN( 0., pcc(ji  ,jj  ,jk+1) )

               ! up & down beta terms
               zbt = e1e2t(ji,jj) * e3t_n(ji,jj,jk) / p2dt
               zbetup(ji,jj,jk) = ( zup            - paft(ji,jj,jk) ) / ( zpos + zrtrn ) * zbt
               zbetdo(ji,jj,jk) = ( paft(ji,jj,jk) - zdo            ) / ( zneg + zrtrn ) * zbt
            END DO
         END DO
      END DO
      CALL lbc_lnk( zbetup, 'T', 1. )   ;   CALL lbc_lnk( zbetdo, 'T', 1. )   ! lateral boundary cond. (unchanged sign)

      ! 3. monotonic flux in the i & j direction (paa & pbb)
      ! ----------------------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zau = MIN( 1._wp, zbetdo(ji,jj,jk), zbetup(ji+1,jj,jk) )
               zbu = MIN( 1._wp, zbetup(ji,jj,jk), zbetdo(ji+1,jj,jk) )
               zcu =       ( 0.5  + SIGN( 0.5 , paa(ji,jj,jk) ) )
               paa(ji,jj,jk) = paa(ji,jj,jk) * ( zcu * zau + ( 1._wp - zcu) * zbu )

               zav = MIN( 1._wp, zbetdo(ji,jj,jk), zbetup(ji,jj+1,jk) )
               zbv = MIN( 1._wp, zbetup(ji,jj,jk), zbetdo(ji,jj+1,jk) )
               zcv =       ( 0.5  + SIGN( 0.5 , pbb(ji,jj,jk) ) )
               pbb(ji,jj,jk) = pbb(ji,jj,jk) * ( zcv * zav + ( 1._wp - zcv) * zbv )

      ! monotonic flux in the k direction, i.e. pcc
      ! -------------------------------------------
               za = MIN( 1., zbetdo(ji,jj,jk+1), zbetup(ji,jj,jk) )
               zb = MIN( 1., zbetup(ji,jj,jk+1), zbetdo(ji,jj,jk) )
               zc =       ( 0.5  + SIGN( 0.5 , pcc(ji,jj,jk+1) ) )
               pcc(ji,jj,jk+1) = pcc(ji,jj,jk+1) * ( zc * za + ( 1._wp - zc) * zb )
            END DO
         END DO
      END DO
      CALL lbc_lnk( paa, 'U', -1. )   ;   CALL lbc_lnk( pbb, 'V', -1. )   ! lateral boundary condition (changed sign)
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zbetup, zbetdo, zbup, zbdo )
      !
      IF( nn_timing == 1 )  CALL timing_stop('nonosc')
      !
   END SUBROUTINE nonosc


   SUBROUTINE interp_4th_cpt_org( pt_in, pt_out )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interp_4th_cpt_org  ***
      !! 
      !! **  Purpose :   Compute the interpolation of tracer at w-point
      !!
      !! **  Method  :   4th order compact interpolation
      !!----------------------------------------------------------------------
      REAL(wp),DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pt_in    ! now tracer fields
      REAL(wp),DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   pt_out   ! now tracer field interpolated at w-pts
      !
      INTEGER :: ji, jj, jk   ! dummy loop integers
      REAL(wp),DIMENSION(jpi,jpj,jpk) :: zwd, zwi, zws, zwrm, zwt
      !!----------------------------------------------------------------------
      
      DO jk = 3, jpkm1        !==  build the three diagonal matrix  ==!
         DO jj = 1, jpj
            DO ji = 1, jpi
               zwd (ji,jj,jk) = 4._wp
               zwi (ji,jj,jk) = 1._wp
               zws (ji,jj,jk) = 1._wp
               zwrm(ji,jj,jk) = 3._wp * ( pt_in(ji,jj,jk-1) + pt_in(ji,jj,jk) )
               !
               IF( tmask(ji,jj,jk+1) == 0._wp) THEN   ! Switch to second order centered at bottom
                  zwd (ji,jj,jk) = 1._wp
                  zwi (ji,jj,jk) = 0._wp
                  zws (ji,jj,jk) = 0._wp
                  zwrm(ji,jj,jk) = 0.5 * ( pt_in(ji,jj,jk-1) + pt_in(ji,jj,jk) )    
               ENDIF
            END DO
         END DO
      END DO
      !
      jk = 2                                          ! Switch to second order centered at top
      DO jj = 1, jpj
         DO ji = 1, jpi
            zwd (ji,jj,jk) = 1._wp
            zwi (ji,jj,jk) = 0._wp
            zws (ji,jj,jk) = 0._wp
            zwrm(ji,jj,jk) = 0.5 * ( pt_in(ji,jj,jk-1) + pt_in(ji,jj,jk) )
         END DO
      END DO   
      !
      !                       !==  tridiagonal solve  ==!
      DO jj = 1, jpj                ! first recurrence
         DO ji = 1, jpi
            zwt(ji,jj,2) = zwd(ji,jj,2)
         END DO
      END DO
      DO jk = 3, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) /zwt(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 1, jpj                ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
         DO ji = 1, jpi
            pt_out(ji,jj,2) = zwrm(ji,jj,2)
         END DO
      END DO
      DO jk = 3, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               pt_out(ji,jj,jk) = zwrm(ji,jj,jk) - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *pt_out(ji,jj,jk-1)             
            END DO
         END DO
      END DO

      DO jj = 1, jpj                ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
         DO ji = 1, jpi
            pt_out(ji,jj,jpkm1) = pt_out(ji,jj,jpkm1) / zwt(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 2, -1
         DO jj = 1, jpj
            DO ji = 1, jpi
               pt_out(ji,jj,jk) = ( pt_out(ji,jj,jk) - zws(ji,jj,jk) * pt_out(ji,jj,jk+1) ) / zwt(ji,jj,jk)
            END DO
         END DO
      END DO
      !    
   END SUBROUTINE interp_4th_cpt_org
   

   SUBROUTINE interp_4th_cpt( pt_in, pt_out )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interp_4th_cpt  ***
      !! 
      !! **  Purpose :   Compute the interpolation of tracer at w-point
      !!
      !! **  Method  :   4th order compact interpolation
      !!----------------------------------------------------------------------
      REAL(wp),DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pt_in    ! field at t-point
      REAL(wp),DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   pt_out   ! field interpolated at w-point
      !
      INTEGER ::   ji, jj, jk   ! dummy loop integers
      INTEGER ::   ikt, ikb     ! local integers
      REAL(wp),DIMENSION(jpi,jpj,jpk) :: zwd, zwi, zws, zwrm, zwt
      !!----------------------------------------------------------------------
      !
      !                      !==  build the three diagonal matrix & the RHS  ==!
      !
      DO jk = 3, jpkm1                 ! interior (from jk=3 to jpk-1)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zwd (ji,jj,jk) = 3._wp * wmask(ji,jj,jk) + 1._wp                 !       diagonal
               zwi (ji,jj,jk) =         wmask(ji,jj,jk)                         ! lower diagonal
               zws (ji,jj,jk) =         wmask(ji,jj,jk)                         ! upper diagonal
               zwrm(ji,jj,jk) = 3._wp * wmask(ji,jj,jk)                     &   ! RHS
                  &           *       ( pt_in(ji,jj,jk) + pt_in(ji,jj,jk-1) )
            END DO
         END DO
      END DO
      !
!!gm
!      SELECT CASE( kbc )               !* boundary condition
!      CASE( np_NH   )   ! Neumann homogeneous at top & bottom
!      CASE( np_CEN2 )   ! 2nd order centered  at top & bottom
!      END SELECT
!!gm  
      !
      DO jj = 2, jpjm1                 ! 2nd order centered at top & bottom
         DO ji = fs_2, fs_jpim1
            ikt = mikt(ji,jj) + 1            ! w-point below the 1st  wet point
            ikb = mbkt(ji,jj)                !     -   above the last wet point
            !
            zwd (ji,jj,ikt) = 1._wp          ! top
            zwi (ji,jj,ikt) = 0._wp
            zws (ji,jj,ikt) = 0._wp
            zwrm(ji,jj,ikt) = 0.5_wp * ( pt_in(ji,jj,jk-1) + pt_in(ji,jj,jk) )
            !
            zwd (ji,jj,ikb) = 1._wp          ! bottom
            zwi (ji,jj,ikb) = 0._wp
            zws (ji,jj,ikb) = 0._wp
            zwrm(ji,jj,ikb) = 0.5_wp * ( pt_in(ji,jj,jk-1) + pt_in(ji,jj,jk) )            
         END DO
      END DO   
      !
      !                       !==  tridiagonal solver  ==!
      !
      DO jj = 2, jpjm1              !* 1st recurrence:   Tk = Dk - Ik Sk-1 / Tk-1
         DO ji = fs_2, fs_jpim1
            zwt(ji,jj,2) = zwd(ji,jj,2)
         END DO
      END DO
      DO jk = 3, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) /zwt(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1              !* 2nd recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
         DO ji = fs_2, fs_jpim1
            pt_out(ji,jj,2) = zwrm(ji,jj,2)
         END DO
      END DO
      DO jk = 3, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               pt_out(ji,jj,jk) = zwrm(ji,jj,jk) - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *pt_out(ji,jj,jk-1)             
            END DO
         END DO
      END DO

      DO jj = 2, jpjm1              !* 3d recurrence:    Xk = (Zk - Sk Xk+1 ) / Tk
         DO ji = fs_2, fs_jpim1
            pt_out(ji,jj,jpkm1) = pt_out(ji,jj,jpkm1) / zwt(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 2, -1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               pt_out(ji,jj,jk) = ( pt_out(ji,jj,jk) - zws(ji,jj,jk) * pt_out(ji,jj,jk+1) ) / zwt(ji,jj,jk)
            END DO
         END DO
      END DO
      !    
   END SUBROUTINE interp_4th_cpt


   SUBROUTINE tridia_solver( pD, pU, pL, pRHS, pt_out , klev )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tridia_solver  ***
      !! 
      !! **  Purpose :   solve a symmetric 3diagonal system
      !!
      !! **  Method  :   solve M.t_out = RHS(t)  where M is a tri diagonal matrix ( jpk*jpk )
      !!     
      !!             ( D_1 U_1  0   0   0  )( t_1 )   ( RHS_1 )
      !!             ( L_2 D_2 U_2  0   0  )( t_2 )   ( RHS_2 )
      !!             (  0  L_3 D_3 U_3  0  )( t_3 ) = ( RHS_3 )
      !!             (        ...          )( ... )   ( ...  )
      !!             (  0   0   0  L_k D_k )( t_k )   ( RHS_k )
      !!     
      !!        M is decomposed in the product of an upper and lower triangular matrix.
      !!        The tri-diagonals matrix is given as input 3D arrays:   pD, pU, pL 
      !!        (i.e. the Diagonal, the Upper diagonal, and the Lower diagonal).
      !!        The solution is pta.
      !!        The 3d array zwt is used as a work space array.
      !!----------------------------------------------------------------------
      REAL(wp),DIMENSION(:,:,:), INTENT(in   ) ::   pD, pU, PL    ! 3-diagonal matrix
      REAL(wp),DIMENSION(:,:,:), INTENT(in   ) ::   pRHS          ! Right-Hand-Side
      REAL(wp),DIMENSION(:,:,:), INTENT(  out) ::   pt_out        !!gm field at level=F(klev)
      INTEGER                  , INTENT(in   ) ::   klev          ! =1 pt_out at w-level 
      !                                                           ! =0 pt at t-level
      INTEGER ::   ji, jj, jk   ! dummy loop integers
      INTEGER ::   kstart       ! local indices
      REAL(wp),DIMENSION(jpi,jpj,jpk) ::   zwt   ! 3D work array
      !!----------------------------------------------------------------------
      !
      kstart =  1  + klev
      !
      DO jj = 2, jpjm1              !* 1st recurrence:   Tk = Dk - Ik Sk-1 / Tk-1
         DO ji = fs_2, fs_jpim1
            zwt(ji,jj,kstart) = pD(ji,jj,kstart)
         END DO
      END DO
      DO jk = kstart+1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zwt(ji,jj,jk) = pD(ji,jj,jk) - pL(ji,jj,jk) * pU(ji,jj,jk-1) /zwt(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1              !* 2nd recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
         DO ji = fs_2, fs_jpim1
            pt_out(ji,jj,kstart) = pRHS(ji,jj,kstart)
         END DO
      END DO
      DO jk = kstart+1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               pt_out(ji,jj,jk) = pRHS(ji,jj,jk) - pL(ji,jj,jk) / zwt(ji,jj,jk-1) *pt_out(ji,jj,jk-1)             
            END DO
         END DO
      END DO

      DO jj = 2, jpjm1              !* 3d recurrence:    Xk = (Zk - Sk Xk+1 ) / Tk
         DO ji = fs_2, fs_jpim1
            pt_out(ji,jj,jpkm1) = pt_out(ji,jj,jpkm1) / zwt(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, kstart, -1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               pt_out(ji,jj,jk) = ( pt_out(ji,jj,jk) - pU(ji,jj,jk) * pt_out(ji,jj,jk+1) ) / zwt(ji,jj,jk)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE tridia_solver

   !!======================================================================
END MODULE traadv_fct
