MODULE limupdate1
   !!======================================================================
   !!                     ***  MODULE  limupdate1  ***
   !!   LIM-3 : Update of sea-ice global variables at the end of the time step
   !!======================================================================
   !! History :  3.0  !  2006-04  (M. Vancoppenolle) Original code
   !!            3.5  !  2014-06  (C. Rousset)       Complete rewriting/cleaning
   !!----------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3'                                      LIM3 sea-ice model
   !!----------------------------------------------------------------------
   !!    lim_update1   : computes update of sea-ice global variables from trend terms
   !!----------------------------------------------------------------------
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE dom_oce
   USE phycst          ! physical constants
   USE ice
   USE thd_ice         ! LIM thermodynamic sea-ice variables
   USE limitd_th
   USE limvar
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing
   USE limcons         ! conservation tests
   USE limctl          ! control prints
   USE lib_mpp         ! MPP library
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_update1

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/LIM3 4.0 , UCL - NEMO Consortium (2011)
   !! $Id: limupdate1.F90 7753 2017-03-03 11:46:59Z mocavero $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_update1( kt )
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE lim_update1  ***
      !!               
      !! ** Purpose :  Computes update of sea-ice global variables at 
      !!               the end of the dynamics.
      !!                
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! number of iteration
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) ::   zsal
      REAL(wp) ::   zvi_b, zsmv_b, zei_b, zfs_b, zfw_b, zft_b 
      !!-------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('limupdate1')

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)'' 
         WRITE(numout,*)' lim_update1 ' 
         WRITE(numout,*)' ~~~~~~~~~~~ '
      ENDIF

      ! conservation test
      IF( ln_limdiachk ) CALL lim_cons_hsm(0, 'limupdate1', zvi_b, zsmv_b, zei_b, zfw_b, zfs_b, zft_b)

      !----------------------------------------------------
      ! ice concentration should not exceed amax 
      !-----------------------------------------------------
      at_i(:,:) = 0._wp
      DO jl = 1, jpl
         at_i(:,:) = a_i(:,:,jl) + at_i(:,:)
      END DO

      DO jl  = 1, jpl
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( at_i(ji,jj) > rn_amax_2d(ji,jj) .AND. a_i(ji,jj,jl) > 0._wp ) THEN
                  a_i (ji,jj,jl) = a_i (ji,jj,jl) * ( 1._wp - ( 1._wp - rn_amax_2d(ji,jj) / at_i(ji,jj) ) )
                  oa_i(ji,jj,jl) = oa_i(ji,jj,jl) * ( 1._wp - ( 1._wp - rn_amax_2d(ji,jj) / at_i(ji,jj) ) )
               ENDIF
            END DO
         END DO
      END DO
    
      !---------------------
      ! Ice salinity bounds
      !---------------------
      IF (  nn_icesal == 2  ) THEN 
         DO jl = 1, jpl
            DO jj = 1, jpj 
               DO ji = 1, jpi
                  zsal            = smv_i(ji,jj,jl)
                  ! salinity stays in bounds
                  rswitch         = 1._wp - MAX( 0._wp, SIGN( 1._wp, - v_i(ji,jj,jl) ) )
                  smv_i(ji,jj,jl) = rswitch * MAX( MIN( rn_simax * v_i(ji,jj,jl), smv_i(ji,jj,jl) ), rn_simin * v_i(ji,jj,jl) )
                  ! associated salt flux
                  sfx_res(ji,jj) = sfx_res(ji,jj) - ( smv_i(ji,jj,jl) - zsal ) * rhoic * r1_rdtice
               END DO
            END DO
         END DO
      ENDIF

      !----------------------------------------------------
      ! Rebin categories with thickness out of bounds
      !----------------------------------------------------
      IF ( jpl > 1 ) CALL lim_itd_th_reb(1, jpl)

      !-----------------
      ! zap small values
      !-----------------
      CALL lim_var_zapsmall

      ! -------------------------------------------------
      ! Diagnostics
      ! -------------------------------------------------
      DO jl  = 1, jpl
         afx_dyn(:,:) = afx_dyn(:,:) + ( a_i(:,:,jl) - a_i_b(:,:,jl) ) * r1_rdtice
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi            
            ! heat content variation (W.m-2)
            diag_heat(ji,jj) = - ( SUM( e_i(ji,jj,1:nlay_i,:) - e_i_b(ji,jj,1:nlay_i,:) ) +  & 
               &                   SUM( e_s(ji,jj,1:nlay_s,:) - e_s_b(ji,jj,1:nlay_s,:) )    &
               &                 ) * r1_rdtice
            ! salt, volume
            diag_smvi(ji,jj) = SUM( smv_i(ji,jj,:) - smv_i_b(ji,jj,:) ) * rhoic * r1_rdtice
            diag_vice(ji,jj) = SUM( v_i  (ji,jj,:) - v_i_b  (ji,jj,:) ) * rhoic * r1_rdtice
            diag_vsnw(ji,jj) = SUM( v_s  (ji,jj,:) - v_s_b  (ji,jj,:) ) * rhosn * r1_rdtice
         END DO
      END DO

      ! conservation test
      IF( ln_limdiachk ) CALL lim_cons_hsm(1, 'limupdate1', zvi_b, zsmv_b, zei_b, zfw_b, zfs_b, zft_b)

      ! control prints
      IF( ln_ctl )       CALL lim_prt3D( 'limupdate1' )
   
      IF( nn_timing == 1 )  CALL timing_stop('limupdate1')

   END SUBROUTINE lim_update1

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty Module               No sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_update1     ! Empty routine
   END SUBROUTINE lim_update1

#endif

END MODULE limupdate1
