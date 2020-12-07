MODULE limadv_umx
   !!==============================================================================
   !!                       ***  MODULE  limadv_umx  ***
   !! LIM sea-ice model : sea-ice advection using the ULTIMATE-MACHO scheme
   !!==============================================================================
   !! History :  3.5  !  2014-11  (C. Rousset, G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   lim_adv_umx   : update the tracer trend with the 3D advection trends using a TVD scheme
   !!   ultimate      : compute a tracer value at velocity points using ULTIMATE scheme at various orders
   !!   macho         : 
   !!   nonosc_2d     : compute monotonic tracer fluxes by a non-oscillatory algorithm 
   !!----------------------------------------------------------------------
   USE phycst         ! physical constant
   USE dom_oce        ! ocean domain
   USE sbc_oce        ! ocean surface boundary condition
   USE ice            ! ice variables
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! lateral boundary conditions -- MPP exchanges
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! work arrays
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_adv_umx    ! routine called by limtrp.F90
      
   REAL(wp) ::   z1_6   = 1._wp / 6._wp   ! =1/6
   REAL(wp) ::   z1_120 = 1._wp / 120._wp ! =1/120

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: limadv_umx.F90 4499 2014-02-18 15:14:31Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_adv_umx( kt, pdt, puc, pvc, pubox, pvbox, ptc )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE lim_adv_umx  ***
      !! 
      !! **  Purpose :   Compute the now trend due to total advection of 
      !!       tracers and add it to the general trend of tracer equations
      !!
      !! **  Method  :   TVD scheme, i.e. 2nd order centered scheme with
      !!       corrected flux (monotonic correction)
      !!       note: - this advection scheme needs a leap-frog time scheme
      !!
      !! ** Action : - pt  the after advective tracer
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   )           ::   kt         ! number of iteration
      REAL(wp)                    , INTENT(in   )           ::   pdt        ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   )           ::   puc, pvc   ! 2 ice velocity components => u*e2
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   )           ::   pubox, pvbox   ! upstream velocity
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout)           ::   ptc        ! tracer content field
      !
      INTEGER  ::   ji, jj           ! dummy loop indices  
      REAL(wp) ::   ztra             ! local scalar
      REAL(wp) ::   zfp_ui, zfp_vj   !   -      -
      REAL(wp) ::   zfm_ui, zfm_vj   !   -      -
      REAL(wp), POINTER, DIMENSION(:,:) :: zt_ups, zfu_ups, zfv_ups, ztrd, zfu_ho, zfv_ho, zt_u, zt_v
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('lim_adv_umx')
      !
      CALL wrk_alloc( jpi,jpj,   zt_ups, zfu_ups, zfv_ups, ztrd, zfu_ho, zfv_ho, zt_u, zt_v )
      !
      !
      !  upstream advection with initial mass fluxes & intermediate update
      ! --------------------------------------------------------------------
      DO jj = 1, jpjm1         ! upstream tracer flux in the i and j direction
         DO ji = 1, fs_jpim1   ! vector opt.
            zfp_ui = puc(ji,jj) + ABS( puc(ji,jj) )
            zfm_ui = puc(ji,jj) - ABS( puc(ji,jj) )
            zfp_vj = pvc(ji,jj) + ABS( pvc(ji,jj) )
            zfm_vj = pvc(ji,jj) - ABS( pvc(ji,jj) )
            zfu_ups(ji,jj) = 0.5_wp * ( zfp_ui * ptc(ji,jj) + zfm_ui * ptc(ji+1,jj  ) )
            zfv_ups(ji,jj) = 0.5_wp * ( zfp_vj * ptc(ji,jj) + zfm_vj * ptc(ji  ,jj+1) )
         END DO
      END DO
      
      DO jj = 2, jpjm1            ! total intermediate advective trends
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ztra = - (   zfu_ups(ji,jj) - zfu_ups(ji-1,jj  )   &
               &       + zfv_ups(ji,jj) - zfv_ups(ji  ,jj-1)   ) * r1_e1e2t(ji,jj)
            !
            ztrd(ji,jj) =                         ztra                         ! upstream trend [ -div(uh) or -div(uhT) ]  
            zt_ups (ji,jj) = ( ptc(ji,jj) + pdt * ztra ) * tmask(ji,jj,1)      ! guess after content field with monotonic scheme
         END DO
      END DO
      CALL lbc_lnk( zt_ups, 'T', 1. )        ! Lateral boundary conditions   (unchanged sign)
      
      ! High order (_ho) fluxes 
      ! -----------------------
      SELECT CASE( nn_limadv_ord )
      CASE ( 20 )                          ! centered second order
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zfu_ho(ji,jj) = 0.5 * puc(ji,jj) * ( ptc(ji,jj) + ptc(ji+1,jj) )
               zfv_ho(ji,jj) = 0.5 * pvc(ji,jj) * ( ptc(ji,jj) + ptc(ji,jj+1) )
            END DO
         END DO
         !
      CASE ( 1:5 )                      ! 1st to 5th order ULTIMATE-MACHO scheme
         CALL macho( kt, nn_limadv_ord, pdt, ptc, puc, pvc, pubox, pvbox, zt_u, zt_v )
         !
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zfu_ho(ji,jj) = puc(ji,jj) * zt_u(ji,jj)
               zfv_ho(ji,jj) = pvc(ji,jj) * zt_v(ji,jj)
            END DO
         END DO
         !
      END SELECT
         
      ! antidiffusive flux : high order minus low order
      ! --------------------------------------------------
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zfu_ho(ji,jj) = zfu_ho(ji,jj) - zfu_ups(ji,jj)
            zfv_ho(ji,jj) = zfv_ho(ji,jj) - zfv_ups(ji,jj)
         END DO
      END DO
      CALL lbc_lnk_multi( zfu_ho, 'U', -1., zfv_ho, 'V', -1. )         ! Lateral bondary conditions
      
      ! monotonicity algorithm
      ! -------------------------
      CALL nonosc_2d( ptc, zfu_ho, zfv_ho, zt_ups, pdt )
      
      ! final trend with corrected fluxes
      ! ------------------------------------
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.  
            ztra       = ztrd(ji,jj)  - (  zfu_ho(ji,jj) - zfu_ho(ji-1,jj  )   &
               &                         + zfv_ho(ji,jj) - zfv_ho(ji  ,jj-1) ) * r1_e1e2t(ji,jj)  
            ptc(ji,jj) = ptc(ji,jj) + pdt * ztra
         END DO
      END DO
      CALL lbc_lnk( ptc(:,:) , 'T',  1. )
      !
      !
      CALL wrk_dealloc( jpi,jpj,   zt_ups, zfu_ups, zfv_ups, ztrd, zfu_ho, zfv_ho, zt_u, zt_v )
      !
      IF( nn_timing == 1 )  CALL timing_stop('lim_adv_umx')
      !
   END SUBROUTINE lim_adv_umx


   SUBROUTINE macho( kt, k_order, pdt, ptc, puc, pvc, pubox, pvbox, pt_u, pt_v )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ultimate_x  ***
      !!     
      !! **  Purpose :   compute  
      !!
      !! **  Method  :   ... ???
      !!                 TIM = transient interpolation Modeling 
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   kt         ! number of iteration
      INTEGER                     , INTENT(in   ) ::   k_order    ! order of the ULTIMATE scheme
      REAL(wp)                    , INTENT(in   ) ::   pdt        ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   ptc        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   puc, pvc   ! 2 ice velocity components
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pubox, pvbox   ! upstream velocity
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pt_u, pt_v ! tracer at u- and v-points 
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) ::   zc_box    !   -      -
      REAL(wp), POINTER, DIMENSION(:,:) :: zzt
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('macho')
      !
      CALL wrk_alloc( jpi,jpj,   zzt )
      !
      IF( MOD( (kt - 1) / nn_fsbc , 2 ) == 0 ) THEN         !==  odd ice time step:  adv_x then adv_y  ==!
         !
         !                                                           !--  ultimate interpolation of pt at u-point  --!
         CALL ultimate_x( k_order, pdt, ptc, puc, pt_u )
         !
         !                                                           !--  advective form update in zzt  --!
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zzt(ji,jj) = ptc(ji,jj) - pubox(ji,jj) * pdt * ( pt_u(ji,jj) - pt_u(ji-1,jj) ) * r1_e1t(ji,jj)  &
                  &                    - ptc  (ji,jj) * pdt * ( puc (ji,jj) - puc (ji-1,jj) ) * r1_e1e2t(ji,jj)
               zzt(ji,jj) = zzt(ji,jj) * tmask(ji,jj,1)
            END DO
         END DO
         CALL lbc_lnk( zzt, 'T', 1. )
         !
         !                                                           !--  ultimate interpolation of pt at v-point  --!
         CALL ultimate_y( k_order, pdt, zzt, pvc, pt_v )
         !
      ELSE                                                  !==  even ice time step:  adv_y then adv_x  ==!
         !
         !                                                           !--  ultimate interpolation of pt at v-point  --!
         CALL ultimate_y( k_order, pdt, ptc, pvc, pt_v )
         !
         !                                                           !--  advective form update in zzt  --!
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zzt(ji,jj) = ptc(ji,jj) - pvbox(ji,jj) * pdt * ( pt_v(ji,jj) - pt_v(ji,jj-1) ) * r1_e2t(ji,jj)  &
                  &                    - ptc  (ji,jj) * pdt * ( pvc (ji,jj) - pvc (ji,jj-1) ) * r1_e1e2t(ji,jj)
               zzt(ji,jj) = zzt(ji,jj) * tmask(ji,jj,1)
            END DO
         END DO
         CALL lbc_lnk( zzt, 'T', 1. )
         !
         !                                                           !--  ultimate interpolation of pt at u-point  --!
         CALL ultimate_x( k_order, pdt, zzt, puc, pt_u )
         !      
      ENDIF      
      !
      CALL wrk_dealloc( jpi,jpj,   zzt )
      !
      IF( nn_timing == 1 )  CALL timing_stop('macho')
      !
   END SUBROUTINE macho


   SUBROUTINE ultimate_x( k_order, pdt, pt, puc, pt_u )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ultimate_x  ***
      !!     
      !! **  Purpose :   compute  
      !!
      !! **  Method  :   ... ???
      !!                 TIM = transient interpolation Modeling 
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   k_order   ! ocean time-step index
      REAL(wp)                    , INTENT(in   ) ::   pdt       ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   puc       ! ice i-velocity component
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pt_u      ! tracer at u-point 
      !
      INTEGER  ::   ji, jj       ! dummy loop indices
      REAL(wp) ::   zcu, zdx2, zdx4    !   -      -
      REAL(wp), POINTER, DIMENSION(:,:) :: ztu1, ztu2, ztu3, ztu4
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ultimate_x')
      !
      CALL wrk_alloc( jpi,jpj,   ztu1, ztu2, ztu3, ztu4 )
      !
      !                                                     !--  Laplacian in i-direction  --!
      DO jj = 2, jpjm1         ! First derivative (gradient)
         DO ji = 1, fs_jpim1
            ztu1(ji,jj) = ( pt(ji+1,jj) - pt(ji,jj) ) * r1_e1u(ji,jj) * umask(ji,jj,1)
         END DO
         !                     ! Second derivative (Laplacian)
         DO ji = fs_2, fs_jpim1
            ztu2(ji,jj) = ( ztu1(ji,jj) - ztu1(ji-1,jj) ) * r1_e1t(ji,jj)
         END DO
      END DO
      CALL lbc_lnk( ztu2, 'T', 1. )
      !
      !                                                     !--  BiLaplacian in i-direction  --!
      DO jj = 2, jpjm1         ! Third derivative
         DO ji = 1, fs_jpim1
            ztu3(ji,jj) = ( ztu2(ji+1,jj) - ztu2(ji,jj) ) * r1_e1u(ji,jj) * umask(ji,jj,1)
         END DO
         !                     ! Fourth derivative
         DO ji = fs_2, fs_jpim1
            ztu4(ji,jj) = ( ztu3(ji,jj) - ztu3(ji-1,jj) ) * r1_e1t(ji,jj)
         END DO
      END DO
      CALL lbc_lnk( ztu4, 'T', 1. )
      !
      !
      SELECT CASE (k_order )
      !
      CASE( 1 )                                                   !==  1st order central TIM  ==! (Eq. 21)
         !        
         DO jj = 1, jpj
            DO ji = 1, fs_jpim1   ! vector opt.
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (                               pt(ji+1,jj) + pt(ji,jj)   &
                  &                                    - SIGN( 1._wp, puc(ji,jj) ) * ( pt(ji+1,jj) - pt(ji,jj) ) )
            END DO
         END DO
         !
      CASE( 2 )                                                   !==  2nd order central TIM  ==! (Eq. 23)
         !
         DO jj = 1, jpj
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = puc(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (                                   pt(ji+1,jj) + pt(ji,jj)   &
                  &                                               -              zcu   * ( pt(ji+1,jj) - pt(ji,jj) ) ) 
            END DO
         END DO
         CALL lbc_lnk( pt_u(:,:) , 'U',  1. )
         !  
      CASE( 3 )                                                   !==  3rd order central TIM  ==! (Eq. 24)
         !
         DO jj = 1, jpj
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = puc(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
               zdx2 = e1u(ji,jj) * e1u(ji,jj)
!!rachid       zdx2 = e1u(ji,jj) * e1t(ji,jj)
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (         (                         pt  (ji+1,jj) + pt  (ji,jj)        &
                  &                                               -              zcu   * ( pt  (ji+1,jj) - pt  (ji,jj) )  )   &
                  &        + z1_6 * zdx2 * ( zcu*zcu - 1._wp ) * (                         ztu2(ji+1,jj) + ztu2(ji,jj)        &
                  &                                               - SIGN( 1._wp, zcu ) * ( ztu2(ji+1,jj) - ztu2(ji,jj) )  )   )
            END DO
         END DO
         !
      CASE( 4 )                                                   !==  4th order central TIM  ==! (Eq. 27)
         !
         DO jj = 1, jpj
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = puc(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
               zdx2 = e1u(ji,jj) * e1u(ji,jj)
!!rachid       zdx2 = e1u(ji,jj) * e1t(ji,jj)
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (         (                   pt  (ji+1,jj) + pt  (ji,jj)        &
                  &                                               -          zcu * ( pt  (ji+1,jj) - pt  (ji,jj) )  )   &
                  &        + z1_6 * zdx2 * ( zcu*zcu - 1._wp ) * (                   ztu2(ji+1,jj) + ztu2(ji,jj)        &
                  &                                               - 0.5_wp * zcu * ( ztu2(ji+1,jj) - ztu2(ji,jj) )  )   )
            END DO
         END DO
         !
      CASE( 5 )                                                   !==  5th order central TIM  ==! (Eq. 29)
         !
         DO jj = 1, jpj
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = puc(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
               zdx2 = e1u(ji,jj) * e1u(ji,jj)
!!rachid       zdx2 = e1u(ji,jj) * e1t(ji,jj)
               zdx4 = zdx2 * zdx2
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (               (                   pt  (ji+1,jj) + pt  (ji,jj)       &
                  &                                                     -          zcu * ( pt  (ji+1,jj) - pt  (ji,jj) ) )   &
                  &        + z1_6   * zdx2 * ( zcu*zcu - 1._wp ) *     (                   ztu2(ji+1,jj) + ztu2(ji,jj)       &
                  &                                                     - 0.5_wp * zcu * ( ztu2(ji+1,jj) - ztu2(ji,jj) ) )   &
                  &        + z1_120 * zdx4 * ( zcu*zcu - 1._wp ) * ( zcu*zcu - 4._wp ) * ( ztu4(ji+1,jj) + ztu4(ji,jj)       &
                  &                                               - SIGN( 1._wp, zcu ) * ( ztu4(ji+1,jj) - ztu4(ji,jj) ) ) )
            END DO
         END DO
         !
      END SELECT
      !
      CALL wrk_dealloc( jpi,jpj,   ztu1, ztu2, ztu3, ztu4 )
      !
      IF( nn_timing == 1 )  CALL timing_stop('ultimate_x')
      !
   END SUBROUTINE ultimate_x
   
 
   SUBROUTINE ultimate_y( k_order, pdt, pt, pvc, pt_v )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ultimate_y  ***
      !!     
      !! **  Purpose :   compute  
      !!
      !! **  Method  :   ... ???
      !!                 TIM = transient interpolation Modeling 
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   k_order   ! ocean time-step index
      REAL(wp)                    , INTENT(in   ) ::   pdt       ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pvc       ! ice j-velocity component
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pt_v      ! tracer at v-point 
      !
      INTEGER  ::   ji, jj       ! dummy loop indices
      REAL(wp) ::   zcv, zdy2, zdy4    !   -      -
      REAL(wp), POINTER, DIMENSION(:,:) :: ztv1, ztv2, ztv3, ztv4
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ultimate_y')
      !
      CALL wrk_alloc( jpi,jpj,   ztv1, ztv2, ztv3, ztv4 )
      !
      !                                                     !--  Laplacian in j-direction  --!
      DO jj = 1, jpjm1         ! First derivative (gradient)
         DO ji = fs_2, fs_jpim1
            ztv1(ji,jj) = ( pt(ji,jj+1) - pt(ji,jj) ) * r1_e2v(ji,jj) * vmask(ji,jj,1)
         END DO
      END DO
      DO jj = 2, jpjm1         ! Second derivative (Laplacian)
         DO ji = fs_2, fs_jpim1
            ztv2(ji,jj) = ( ztv1(ji,jj) - ztv1(ji,jj-1) ) * r1_e2t(ji,jj)
         END DO
      END DO
      CALL lbc_lnk( ztv2, 'T', 1. )
      !
      !                                                     !--  BiLaplacian in j-direction  --!
      DO jj = 1, jpjm1         ! First derivative
         DO ji = fs_2, fs_jpim1
            ztv3(ji,jj) = ( ztv2(ji,jj+1) - ztv2(ji,jj) ) * r1_e2v(ji,jj) * vmask(ji,jj,1)
         END DO
      END DO
      DO jj = 2, jpjm1         ! Second derivative
         DO ji = fs_2, fs_jpim1
            ztv4(ji,jj) = ( ztv3(ji,jj) - ztv3(ji,jj-1) ) * r1_e2t(ji,jj)
         END DO
      END DO
      CALL lbc_lnk( ztv4, 'T', 1. )
      !
      !
      SELECT CASE (k_order )
         !
      CASE( 1 )                                                   !==  1st order central TIM  ==! (Eq. 21)
         !        
         DO jj = 1, jpjm1
            DO ji = 1, jpi
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (                              ( pt(ji,jj+1) + pt(ji,jj) )  &
                  &                                     - SIGN( 1._wp, pvc(ji,jj) ) * ( pt(ji,jj+1) - pt(ji,jj) ) )
            END DO
         END DO
         !
      CASE( 2 )                                                   !==  2nd order central TIM  ==! (Eq. 23)
         DO jj = 1, jpjm1
            DO ji = 1, jpi
               zcv  = pvc(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (        ( pt(ji,jj+1) + pt(ji,jj) )  &
                  &                                     - zcv * ( pt(ji,jj+1) - pt(ji,jj) ) )
            END DO
         END DO
         CALL lbc_lnk( pt_v(:,:) , 'V',  1. )
         !
      CASE( 3 )                                                   !==  3rd order central TIM  ==! (Eq. 24)
         !
         DO jj = 1, jpjm1
            DO ji = 1, jpi
               zcv  = pvc(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
               zdy2 = e2v(ji,jj) * e2v(ji,jj)
!!rachid       zdy2 = e2v(ji,jj) * e2t(ji,jj)
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (                                 ( pt  (ji,jj+1) + pt  (ji,jj)       &
                  &                                     -                        zcv   * ( pt  (ji,jj+1) - pt  (ji,jj) ) )   &
                  &        + z1_6 * zdy2 * ( zcv*zcv - 1._wp ) * (                         ztv2(ji,jj+1) + ztv2(ji,jj)       &
                  &                                               - SIGN( 1._wp, zcv ) * ( ztv2(ji,jj+1) - ztv2(ji,jj) ) ) )
            END DO
         END DO
         !
      CASE( 4 )                                                   !==  4th order central TIM  ==! (Eq. 27)
         !
         DO jj = 1, jpjm1
            DO ji = 1, jpi
               zcv  = pvc(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
               zdy2 = e2v(ji,jj) * e2v(ji,jj)
!!rachid       zdy2 = e2v(ji,jj) * e2t(ji,jj)
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (                           ( pt  (ji,jj+1) + pt  (ji,jj)       &
                  &                                               -          zcv * ( pt  (ji,jj+1) - pt  (ji,jj) ) )   &
                  &        + z1_6 * zdy2 * ( zcv*zcv - 1._wp ) * (                   ztv2(ji,jj+1) + ztv2(ji,jj)       &
                  &                                               - 0.5_wp * zcv * ( ztv2(ji,jj+1) - ztv2(ji,jj) ) ) )
            END DO
         END DO
         !
      CASE( 5 )                                                   !==  5th order central TIM  ==! (Eq. 29)
         !
         DO jj = 1, jpjm1
            DO ji = 1, jpi
               zcv  = pvc(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
               zdy2 = e2v(ji,jj) * e2v(ji,jj)
!!rachid       zdy2 = e2v(ji,jj) * e2t(ji,jj)
               zdy4 = zdy2 * zdy2
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (                                 ( pt  (ji,jj+1) + pt  (ji,jj)      &
                  &                                                     -          zcv * ( pt  (ji,jj+1) - pt  (ji,jj) ) )  &
                  &        + z1_6   * zdy2 * ( zcv*zcv - 1._wp ) *     (                   ztv2(ji,jj+1) + ztv2(ji,jj)      &
                  &                                                     - 0.5_wp * zcv * ( ztv2(ji,jj+1) - ztv2(ji,jj) ) )  &
                  &        + z1_120 * zdy4 * ( zcv*zcv - 1._wp ) * ( zcv*zcv - 4._wp ) * ( ztv4(ji,jj+1) + ztv4(ji,jj)      &
                  &                                               - SIGN( 1._wp, zcv ) * ( ztv4(ji,jj+1) - ztv4(ji,jj) ) ) )
            END DO
         END DO
         !
      END SELECT
      !
      CALL wrk_dealloc( jpi,jpj,   ztv1, ztv2, ztv3, ztv4 )
      !
      IF( nn_timing == 1 )  CALL timing_stop('ultimate_y')
      !
   END SUBROUTINE ultimate_y
   
  
   SUBROUTINE nonosc_2d( pbef, paa, pbb, paft, pdt )
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
      REAL(wp)                     , INTENT(in   ) ::   pdt          ! tracer time-step
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ) ::   pbef, paft   ! before & after field
      REAL(wp), DIMENSION (jpi,jpj), INTENT(inout) ::   paa, pbb     ! monotonic fluxes in the 2 directions
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      INTEGER  ::   ikm1      ! local integer
      REAL(wp) ::   zpos, zneg, zbt, za, zb, zc, zbig, zsml, z1_dt   ! local scalars
      REAL(wp) ::   zau, zbu, zcu, zav, zbv, zcv, zup, zdo            !   -      -
      REAL(wp), POINTER, DIMENSION(:,:) :: zbetup, zbetdo, zbup, zbdo, zmsk, zdiv
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('nonosc_2d')
      !
      CALL wrk_alloc( jpi,jpj,   zbetup, zbetdo, zbup, zbdo, zmsk, zdiv )
      !
      zbig = 1.e+40_wp
      zsml = 1.e-15_wp

      ! clem test
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.  
            zdiv(ji,jj) =  - (  paa(ji,jj) - paa(ji-1,jj  )   &
               &              + pbb(ji,jj) - pbb(ji  ,jj-1) )  
         END DO
      END DO
      CALL lbc_lnk( zdiv, 'T', 1. )        ! Lateral boundary conditions   (unchanged sign)

      ! Determine ice masks for before and after tracers 
      WHERE( pbef(:,:) == 0._wp .AND. paft(:,:) == 0._wp .AND. zdiv(:,:) == 0._wp )   ;   zmsk(:,:) = 0._wp
      ELSEWHERE                                                                       ;   zmsk(:,:) = 1._wp * tmask(:,:,1)
      END WHERE

      ! Search local extrema
      ! --------------------
      ! max/min of pbef & paft with large negative/positive value (-/+zbig) inside land
!      zbup(:,:) = MAX( pbef(:,:) * tmask(:,:,1) - zbig * ( 1.e0 - tmask(:,:,1) ),   &
!         &             paft(:,:) * tmask(:,:,1) - zbig * ( 1.e0 - tmask(:,:,1) )  )
!      zbdo(:,:) = MIN( pbef(:,:) * tmask(:,:,1) + zbig * ( 1.e0 - tmask(:,:,1) ),   &
!         &             paft(:,:) * tmask(:,:,1) + zbig * ( 1.e0 - tmask(:,:,1) )  )
      zbup(:,:) = MAX( pbef(:,:) * zmsk(:,:) - zbig * ( 1.e0 - zmsk(:,:) ),   &
         &             paft(:,:) * zmsk(:,:) - zbig * ( 1.e0 - zmsk(:,:) )  )
      zbdo(:,:) = MIN( pbef(:,:) * zmsk(:,:) + zbig * ( 1.e0 - zmsk(:,:) ),   &
         &             paft(:,:) * zmsk(:,:) + zbig * ( 1.e0 - zmsk(:,:) )  )

      z1_dt = 1._wp / pdt
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            !
            zup  = MAX(   zbup(ji,jj), zbup(ji-1,jj  ), zbup(ji+1,jj  ),   &        ! search max/min in neighbourhood
               &                       zbup(ji  ,jj-1), zbup(ji  ,jj+1)    )
            zdo  = MIN(   zbdo(ji,jj), zbdo(ji-1,jj  ), zbdo(ji+1,jj  ),   &
               &                       zbdo(ji  ,jj-1), zbdo(ji  ,jj+1)    )
               !
            zpos = MAX( 0., paa(ji-1,jj  ) ) - MIN( 0., paa(ji  ,jj  ) )   &        ! positive/negative  part of the flux
               & + MAX( 0., pbb(ji  ,jj-1) ) - MIN( 0., pbb(ji  ,jj  ) )
            zneg = MAX( 0., paa(ji  ,jj  ) ) - MIN( 0., paa(ji-1,jj  ) )   &
               & + MAX( 0., pbb(ji  ,jj  ) ) - MIN( 0., pbb(ji  ,jj-1) )
               !
            zbt = e1e2t(ji,jj) * z1_dt                                   ! up & down beta terms
            zbetup(ji,jj) = ( zup         - paft(ji,jj) ) / ( zpos + zsml ) * zbt
            zbetdo(ji,jj) = ( paft(ji,jj) - zdo         ) / ( zneg + zsml ) * zbt
         END DO
      END DO
      CALL lbc_lnk_multi( zbetup, 'T', 1., zbetdo, 'T', 1. )   ! lateral boundary cond. (unchanged sign)

      ! monotonic flux in the i & j direction (paa & pbb)
      ! -------------------------------------
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zau = MIN( 1._wp , zbetdo(ji,jj) , zbetup(ji+1,jj) )
            zbu = MIN( 1._wp , zbetup(ji,jj) , zbetdo(ji+1,jj) )
            zcu = 0.5  + SIGN( 0.5 , paa(ji,jj) )
            !
            zav = MIN( 1._wp , zbetdo(ji,jj) , zbetup(ji,jj+1) )
            zbv = MIN( 1._wp , zbetup(ji,jj) , zbetdo(ji,jj+1) )
            zcv = 0.5  + SIGN( 0.5 , pbb(ji,jj) )
            !
            paa(ji,jj) = paa(ji,jj) * ( zcu * zau + ( 1._wp - zcu) * zbu )
            pbb(ji,jj) = pbb(ji,jj) * ( zcv * zav + ( 1._wp - zcv) * zbv )
            !
         END DO
      END DO
      CALL lbc_lnk_multi( paa, 'U', -1., pbb, 'V', -1. )   ! lateral boundary condition (changed sign)
      !
      CALL wrk_dealloc( jpi,jpj,   zbetup, zbetdo, zbup, zbdo, zmsk, zdiv )
      !
      IF( nn_timing == 1 )  CALL timing_stop('nonosc_2d')
      !
   END SUBROUTINE nonosc_2d

   !!======================================================================
END MODULE limadv_umx
