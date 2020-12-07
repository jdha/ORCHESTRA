MODULE dynzdf_exp
   !!==============================================================================
   !!                     ***  MODULE  dynzdf_exp  ***
   !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
   !!==============================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            8.0  !  1997-05  (G. Madec)  vertical component of isopycnal
   !!   NEMO     0.5  !  2002-08  (G. Madec)  F90: Free form and module
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!            3.7  !  2015-11  (J. Chanut) output velocities instead of trends
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf_exp   : update the momentum trend with the vertical diffusion using a split-explicit scheme
   !!                   and perform the Leap-Frog time integration.
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE zdf_oce        ! ocean vertical physics
   USE dynadv   , ONLY: ln_dynadv_vec ! Momentum advection form
   USE sbc_oce        ! surface boundary condition: ocean
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf_exp   ! called by step.F90
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2015)
   !! $Id: dynzdf_exp.F90 6140 2015-12-21 11:35:23Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_zdf_exp( kt, p2dt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_exp  ***
      !!                   
      !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
      !!              and perform the Leap-Frog time stepping.
      !!
      !! ** Method  : - Split-explicit forward time stepping.
      !!      The vertical mixing of momentum is given by:
      !!         diffu = dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ub) )
      !!      Surface boundary conditions: wind stress input (averaged over kt-1/2 & kt+1/2)
      !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F90)
      !!      Add this trend to the general trend ua :
      !!         ua = ua + dz( avmu dz(u) )
      !!              - Leap-Frog time stepping (Asselin filter will be applied in dyn_nxt) 
      !!         ua =         ub + 2*dt *       ua             vector form or linear free surf.
      !!         ua = ( e3u_b*ub + 2*dt * e3u_n*ua ) / e3u_a   otherwise
      !!
      !! ** Action : - (ua,va) after velocity
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in) ::   kt     ! ocean time-step index
      REAL(wp), INTENT(in) ::   p2dt   ! time-step 
      !
      INTEGER  ::   ji, jj, jk, jl     ! dummy loop indices
      REAL(wp) ::   zlavmr, zua, zva   ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwx, zwy, zwz, zww
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('dyn_zdf_exp')
      !
      CALL wrk_alloc( jpi,jpj,jpk,   zwx, zwy, zwz, zww ) 
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_zdf_exp : vertical momentum diffusion - explicit operator'
         WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF
      !
      !                 !==  vertical mixing trend  ==!
      !
      zlavmr = 1. / REAL( nn_zdfexp )
      !
      DO jj = 2, jpjm1           ! Surface boundary condition
         DO ji = 2, jpim1
            zwy(ji,jj,1) = ( utau_b(ji,jj) + utau(ji,jj) ) * r1_rau0
            zww(ji,jj,1) = ( vtau_b(ji,jj) + vtau(ji,jj) ) * r1_rau0
         END DO  
      END DO  
      DO jk = 1, jpk             ! Initialization of x, z and contingently trends array
         DO jj = 2, jpjm1 
            DO ji = 2, jpim1
               zwx(ji,jj,jk) = ub(ji,jj,jk)
               zwz(ji,jj,jk) = vb(ji,jj,jk)
            END DO  
         END DO  
      END DO  
      !
      DO jl = 1, nn_zdfexp       ! Time splitting loop
         !
         DO jk = 2, jpk                ! First vertical derivative
            DO jj = 2, jpjm1 
               DO ji = 2, jpim1
                  zwy(ji,jj,jk) = avmu(ji,jj,jk) * ( zwx(ji,jj,jk-1) - zwx(ji,jj,jk) ) / e3uw_n(ji,jj,jk) 
                  zww(ji,jj,jk) = avmv(ji,jj,jk) * ( zwz(ji,jj,jk-1) - zwz(ji,jj,jk) ) / e3vw_n(ji,jj,jk)
               END DO  
            END DO  
         END DO  
         DO jk = 1, jpkm1              ! Second vertical derivative and trend estimation at kt+l*rdt/nn_zdfexp
            DO jj = 2, jpjm1 
               DO ji = 2, jpim1
                  zua = zlavmr * ( zwy(ji,jj,jk) - zwy(ji,jj,jk+1) ) / e3u_n(ji,jj,jk)
                  zva = zlavmr * ( zww(ji,jj,jk) - zww(ji,jj,jk+1) ) / e3v_n(ji,jj,jk)
                  ua(ji,jj,jk) = ua(ji,jj,jk) + zua
                  va(ji,jj,jk) = va(ji,jj,jk) + zva
                  !
                  zwx(ji,jj,jk) = zwx(ji,jj,jk) + p2dt * zua * umask(ji,jj,jk)
                  zwz(ji,jj,jk) = zwz(ji,jj,jk) + p2dt * zva * vmask(ji,jj,jk)
               END DO  
            END DO  
         END DO  
      END DO                     ! End of time splitting
      !
      !
      !                 !==  Leap-Frog time integration  ==!
      !
      IF( ln_dynadv_vec .OR. ln_linssh ) THEN   ! applied on velocity
         DO jk = 1, jpkm1
            ua(:,:,jk) = ( ub(:,:,jk) + p2dt * ua(:,:,jk) ) * umask(:,:,jk)
            va(:,:,jk) = ( vb(:,:,jk) + p2dt * va(:,:,jk) ) * vmask(:,:,jk)
         END DO
      ELSE                                      ! applied on thickness weighted velocity
         DO jk = 1, jpkm1
            ua(:,:,jk) = (          e3u_b(:,:,jk) * ub(:,:,jk)    &
               &           + p2dt * e3u_n(:,:,jk) * ua(:,:,jk)  ) / e3u_a(:,:,jk) * umask(:,:,jk)
            va(:,:,jk) = (          e3v_b(:,:,jk) * vb(:,:,jk)    &
               &           + p2dt * e3v_n(:,:,jk) * va(:,:,jk)  ) / e3v_a(:,:,jk) * vmask(:,:,jk)
         END DO
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj,jpk,   zwx, zwy, zwz, zww ) 
      !
      IF( nn_timing == 1 )   CALL timing_stop('dyn_zdf_exp')
      !
   END SUBROUTINE dyn_zdf_exp

   !!==============================================================================
END MODULE dynzdf_exp
