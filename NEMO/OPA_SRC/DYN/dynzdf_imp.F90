MODULE dynzdf_imp
   !!======================================================================
   !!                    ***  MODULE  dynzdf_imp  ***
   !! Ocean dynamics:  vertical component(s) of the momentum mixing trend, implicit scheme
   !!======================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            8.0  !  1997-05  (G. Madec)  vertical component of isopycnal
   !!   NEMO     0.5  !  2002-08  (G. Madec)  F90: Free form and module
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!            3.4  !  2012-01  (H. Liu) Semi-implicit bottom friction
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf_imp   : compute the vertical diffusion using a implicit scheme
   !!                   together with the Leap-Frog time integration.
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE domvvl         ! variable volume
   USE sbc_oce        ! surface boundary condition: ocean
   USE dynadv   , ONLY: ln_dynadv_vec ! Momentum advection form
   USE zdf_oce        ! ocean vertical physics
   USE zdfbfr         ! Bottom friction setup
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf_imp   ! called by step.F90

   REAL(wp) ::  r_vvl     ! non-linear free surface indicator: =0 if ln_linssh=T, =1 otherwise 

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynzdf_imp.F90 7753 2017-03-03 11:46:59Z mocavero $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_zdf_imp( kt, p2dt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_imp  ***
      !!                   
      !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
      !!              together with the Leap-Frog time stepping using an 
      !!              implicit scheme.
      !!
      !! ** Method  :  - Leap-Frog time stepping on all trends but the vertical mixing
      !!         ua =         ub + 2*dt *       ua             vector form or linear free surf.
      !!         ua = ( e3u_b*ub + 2*dt * e3u_n*ua ) / e3u_a   otherwise
      !!               - update the after velocity with the implicit vertical mixing.
      !!      This requires to solver the following system: 
      !!         ua = ua + 1/e3u_a dk+1[ avmu / e3uw_a dk[ua] ]
      !!      with the following surface/top/bottom boundary condition:
      !!      surface: wind stress input (averaged over kt-1/2 & kt+1/2)
      !!      top & bottom : top stress (iceshelf-ocean) & bottom stress (cf zdfbfr.F)
      !!
      !! ** Action :   (ua,va) after velocity 
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in) ::  kt     ! ocean time-step index
      REAL(wp), INTENT(in) ::  p2dt   ! vertical profile of tracer time-step
      !
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      INTEGER  ::   ikbu, ikbv    ! local integers
      REAL(wp) ::   zzwi, ze3ua   ! local scalars
      REAL(wp) ::   zzws, ze3va   !   -      -
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwi, zwd, zws
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_zdf_imp')
      !
      CALL wrk_alloc( jpi,jpj,jpk, zwi, zwd, zws ) 
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_zdf_imp : vertical momentum diffusion implicit operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         !
         If( ln_linssh ) THEN   ;    r_vvl = 0._wp    ! non-linear free surface indicator
         ELSE                   ;    r_vvl = 1._wp
         ENDIF
      ENDIF
      !
      !              !==  Time step dynamics  ==!
      !
      IF( ln_dynadv_vec .OR. ln_linssh ) THEN      ! applied on velocity
         DO jk = 1, jpkm1
            ua(:,:,jk) = ( ub(:,:,jk) + p2dt * ua(:,:,jk) ) * umask(:,:,jk)
            va(:,:,jk) = ( vb(:,:,jk) + p2dt * va(:,:,jk) ) * vmask(:,:,jk)
         END DO
      ELSE                                         ! applied on thickness weighted velocity
         DO jk = 1, jpkm1
            ua(:,:,jk) = (         e3u_b(:,:,jk) * ub(:,:,jk)  &
               &          + p2dt * e3u_n(:,:,jk) * ua(:,:,jk)  ) / e3u_a(:,:,jk) * umask(:,:,jk)
            va(:,:,jk) = (         e3v_b(:,:,jk) * vb(:,:,jk)  &
               &          + p2dt * e3v_n(:,:,jk) * va(:,:,jk)  ) / e3v_a(:,:,jk) * vmask(:,:,jk)
         END DO
      ENDIF
      !
      !              !==  Apply semi-implicit bottom friction  ==!
      !
      ! Only needed for semi-implicit bottom friction setup. The explicit
      ! bottom friction has been included in "u(v)a" which act as the R.H.S
      ! column vector of the tri-diagonal matrix equation
      !
      IF( ln_bfrimp ) THEN
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ikbu = mbku(ji,jj)       ! ocean bottom level at u- and v-points 
               ikbv = mbkv(ji,jj)       ! (deepest ocean u- and v-points)
               avmu(ji,jj,ikbu+1) = -bfrua(ji,jj) * e3uw_n(ji,jj,ikbu+1)
               avmv(ji,jj,ikbv+1) = -bfrva(ji,jj) * e3vw_n(ji,jj,ikbv+1)
            END DO
         END DO
         IF ( ln_isfcav ) THEN
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ikbu = miku(ji,jj)       ! ocean top level at u- and v-points 
                  ikbv = mikv(ji,jj)       ! (first wet ocean u- and v-points)
                  IF( ikbu >= 2 )   avmu(ji,jj,ikbu) = -tfrua(ji,jj) * e3uw_n(ji,jj,ikbu)
                  IF( ikbv >= 2 )   avmv(ji,jj,ikbv) = -tfrva(ji,jj) * e3vw_n(ji,jj,ikbv)
               END DO
            END DO
         END IF
      ENDIF
      !
      ! With split-explicit free surface, barotropic stress is treated explicitly
      ! Update velocities at the bottom.
      ! J. Chanut: The bottom stress is computed considering after barotropic velocities, which does 
      !            not lead to the effective stress seen over the whole barotropic loop. 
      ! G. Madec : in linear free surface, e3u_a = e3u_n = e3u_0, so systematic use of e3u_a
      IF( ln_bfrimp .AND. ln_dynspg_ts ) THEN
         DO jk = 1, jpkm1        ! remove barotropic velocities
            ua(:,:,jk) = ( ua(:,:,jk) - ua_b(:,:) ) * umask(:,:,jk)
            va(:,:,jk) = ( va(:,:,jk) - va_b(:,:) ) * vmask(:,:,jk)
         END DO
         DO jj = 2, jpjm1        ! Add bottom/top stress due to barotropic component only
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ikbu = mbku(ji,jj)         ! ocean bottom level at u- and v-points 
               ikbv = mbkv(ji,jj)         ! (deepest ocean u- and v-points)
               ze3ua =  ( 1._wp - r_vvl ) * e3u_n(ji,jj,ikbu) + r_vvl * e3u_a(ji,jj,ikbu)
               ze3va =  ( 1._wp - r_vvl ) * e3v_n(ji,jj,ikbv) + r_vvl * e3v_a(ji,jj,ikbv)
               ua(ji,jj,ikbu) = ua(ji,jj,ikbu) + p2dt * bfrua(ji,jj) * ua_b(ji,jj) / ze3ua
               va(ji,jj,ikbv) = va(ji,jj,ikbv) + p2dt * bfrva(ji,jj) * va_b(ji,jj) / ze3va
            END DO
         END DO
         IF( ln_isfcav ) THEN    ! Ocean cavities (ISF)
            DO jj = 2, jpjm1        
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ikbu = miku(ji,jj)         ! top ocean level at u- and v-points 
                  ikbv = mikv(ji,jj)         ! (first wet ocean u- and v-points)
                  ze3ua =  ( 1._wp - r_vvl ) * e3u_n(ji,jj,ikbu) + r_vvl * e3u_a(ji,jj,ikbu)
                  ze3va =  ( 1._wp - r_vvl ) * e3v_n(ji,jj,ikbv) + r_vvl * e3v_a(ji,jj,ikbv)
                  ua(ji,jj,ikbu) = ua(ji,jj,ikbu) + p2dt * tfrua(ji,jj) * ua_b(ji,jj) / ze3ua
                  va(ji,jj,ikbv) = va(ji,jj,ikbv) + p2dt * tfrva(ji,jj) * va_b(ji,jj) / ze3va
               END DO
            END DO
         END IF
      ENDIF
      !
      !              !==  Vertical diffusion on u  ==!
      !
      ! Matrix and second member construction
      ! bottom boundary condition: both zwi and zws must be masked as avmu can take
      ! non zero value at the ocean bottom depending on the bottom friction used.
      !
      DO jk = 1, jpkm1        ! Matrix
         DO jj = 2, jpjm1 
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ze3ua =  ( 1._wp - r_vvl ) * e3u_n(ji,jj,jk) + r_vvl * e3u_a(ji,jj,jk)   ! after scale factor at T-point
               zzwi = - p2dt * avmu(ji,jj,jk  ) / ( ze3ua * e3uw_n(ji,jj,jk  ) )
               zzws = - p2dt * avmu(ji,jj,jk+1) / ( ze3ua * e3uw_n(ji,jj,jk+1) )
               zwi(ji,jj,jk) = zzwi * wumask(ji,jj,jk  )
               zws(ji,jj,jk) = zzws * wumask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1        ! Surface boundary conditions
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO
      END DO

      ! Matrix inversion starting from the first level
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and a lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (the after velocity) is in ua
      !-----------------------------------------------------------------------
      !
      DO jk = 2, jpkm1        !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==!
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ze3ua =  ( 1._wp - r_vvl ) * e3u_n(ji,jj,1) + r_vvl * e3u_a(ji,jj,1) 
            ua(ji,jj,1) = ua(ji,jj,1) + p2dt * 0.5_wp * ( utau_b(ji,jj) + utau(ji,jj) )   &
               &                                      / ( ze3ua * rau0 ) * umask(ji,jj,1) 
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               ua(ji,jj,jk) = ua(ji,jj,jk) - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * ua(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  thrid recurrence : SOLk = ( Lk - Uk * Ek+1 ) / Dk  ==!
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ua(ji,jj,jpkm1) = ua(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               ua(ji,jj,jk) = ( ua(ji,jj,jk) - zws(ji,jj,jk) * ua(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      !              !==  Vertical diffusion on v  ==!
      !
      ! Matrix and second member construction
      ! bottom boundary condition: both zwi and zws must be masked as avmv can take
      ! non zero value at the ocean bottom depending on the bottom friction used
      !
      DO jk = 1, jpkm1        ! Matrix
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ze3va =  ( 1._wp - r_vvl ) * e3v_n(ji,jj,jk) + r_vvl * e3v_a(ji,jj,jk)   ! after scale factor at T-point
               zzwi = - p2dt * avmv (ji,jj,jk  ) / ( ze3va * e3vw_n(ji,jj,jk  ) )
               zzws = - p2dt * avmv (ji,jj,jk+1) / ( ze3va * e3vw_n(ji,jj,jk+1) )
               zwi(ji,jj,jk) = zzwi * wvmask(ji,jj,jk  )
               zws(ji,jj,jk) = zzws * wvmask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1        ! Surface boundary conditions
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO
      END DO

      ! Matrix inversion
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (after velocity) is in 2d array va
      !-----------------------------------------------------------------------
      !
      DO jk = 2, jpkm1        !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==!
         DO ji = fs_2, fs_jpim1   ! vector opt.          
            ze3va =  ( 1._wp - r_vvl ) * e3v_n(ji,jj,1) + r_vvl * e3v_a(ji,jj,1) 
            va(ji,jj,1) = va(ji,jj,1) + p2dt * 0.5_wp * ( vtau_b(ji,jj) + vtau(ji,jj) )   &
               &                                      / ( ze3va * rau0 ) * vmask(ji,jj,1) 
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               va(ji,jj,jk) = va(ji,jj,jk) - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * va(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  third recurrence : SOLk = ( Lk - Uk * SOLk+1 ) / Dk  ==!
         DO ji = fs_2, fs_jpim1   ! vector opt.
            va(ji,jj,jpkm1) = va(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               va(ji,jj,jk) = ( va(ji,jj,jk) - zws(ji,jj,jk) * va(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO
      
      ! J. Chanut: Lines below are useless ?
      !! restore bottom layer avmu(v) 
      !!gm  I almost sure it is !!!!
      IF( ln_bfrimp ) THEN
        DO jj = 2, jpjm1
           DO ji = 2, jpim1
              ikbu = mbku(ji,jj)         ! ocean bottom level at u- and v-points 
              ikbv = mbkv(ji,jj)         ! (deepest ocean u- and v-points)
              avmu(ji,jj,ikbu+1) = 0._wp
              avmv(ji,jj,ikbv+1) = 0._wp
           END DO
        END DO
        IF (ln_isfcav) THEN
           DO jj = 2, jpjm1
              DO ji = 2, jpim1
                 ikbu = miku(ji,jj)         ! ocean top level at u- and v-points 
                 ikbv = mikv(ji,jj)         ! (first wet ocean u- and v-points)
                 IF( ikbu > 1 )   avmu(ji,jj,ikbu) = 0._wp
                 IF( ikbv > 1 )   avmv(ji,jj,ikbv) = 0._wp
              END DO
           END DO
        ENDIF
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj,jpk,   zwi, zwd, zws) 
      !
      IF( nn_timing == 1 )   CALL timing_stop('dyn_zdf_imp')
      !
   END SUBROUTINE dyn_zdf_imp

   !!==============================================================================
END MODULE dynzdf_imp
