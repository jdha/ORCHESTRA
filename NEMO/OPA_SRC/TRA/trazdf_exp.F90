MODULE trazdf_exp
   !!==============================================================================
   !!                    ***  MODULE  trazdf_exp  ***
   !! Ocean  tracers:  vertical component of the tracer mixing trend using
   !!                  a split-explicit time-stepping 
   !!==============================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            7.0  !  1991-11  (G. Madec)
   !!                 !  1992-06  (M. Imbard)  correction on tracer trend loops
   !!                 !  1996-01  (G. Madec)  statement function for e3
   !!                 !  1997-05  (G. Madec)  vertical component of isopycnal
   !!                 !  1997-07  (G. Madec)  geopotential diffusion in s-coord
   !!                 !  2000-08  (G. Madec)  double diffusive mixing
   !!   NEMO     1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2004-08  (C. Talandier) New trends organisation
   !!             -   !  2005-11  (G. Madec)  New organisation
   !!            3.0  !  2008-04  (G. Madec)  leap-frog time stepping done in trazdf
   !!            3.3  !  2010-06  (C. Ethe, G. Madec) Merge TRA-TRC
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_zdf_exp   : compute the tracer the vertical diffusion trend using a
   !!                   split-explicit time stepping and provide the after tracer
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers 
   USE dom_oce        ! ocean space and time domain 
   USE domvvl         ! variable volume levels
   USE zdf_oce        ! ocean vertical physics
   USE zdfddm         ! ocean vertical physics: double diffusion
   USE trc_oce        ! share passive tracers/Ocean variables
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_zdf_exp   ! routine called by step.F90

   !! * Substitutions
#  include "zdfddm_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: trazdf_exp.F90 6140 2015-12-21 11:35:23Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_zdf_exp( kt, kit000, cdtype, p2dt, ksts,   &
      &                                        ptb , pta , kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_exp  ***
      !!                   
      !! ** Purpose :   Compute the after tracer fields due to the vertical
      !!      tracer mixing alone, and then due to the whole tracer trend.
      !!
      !! ** Method  : - The after tracer fields due to the vertical diffusion
      !!      of tracers alone is given by:
      !!                ztb = ptb + p2dt difft
      !!      where difft = dz( avt dz(ptb) ) = 1/e3t dk+1( avt/e3w dk(ptb) )
      !!           (if lk_zdfddm=T use avs on salinity and passive tracers instead of avt)
      !!      difft is evaluated with an Euler split-explit scheme using a
      !!      no flux boundary condition at both surface and bottomi boundaries.
      !!      (N.B. bottom condition is applied through the masked field avt).
      !!              - the after tracer fields due to the whole trend is 
      !!      obtained in leap-frog environment applied on thickness weighted tracer by :
      !!          pta = [ ptb*e3tb + e3tn*( ztb - ptb + p2dt pta ) ] / e3tn
      !!
      !! ** Action : - after tracer fields pta
      !!---------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000   ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt     ! number of tracers
      INTEGER                              , INTENT(in   ) ::   ksts     ! number of sub-time step
      REAL(wp)                             , INTENT(in   ) ::   p2dt     ! vertical profile of tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb      ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta      ! in: tracer trend ; out: after tracer field 
      !
      INTEGER  ::  ji, jj, jk, jn, jl   ! dummy loop indices
      REAL(wp) ::  z1_ksts, ze3tr       ! local scalars
      REAL(wp) ::  ztra, ze3tb    !   -      -
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztb, zwf
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_zdf_exp')
      !
      CALL wrk_alloc( jpi,jpj,jpk,   ztb, zwf ) 
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_zdf_exp : explicit vertical mixing on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      ! Initializations
      ! ---------------
      z1_ksts = 1._wp / REAL( ksts, wp )
      zwf(:,:, 1 ) = 0._wp    ! no flux at the surface and at bottom level
      zwf(:,:,jpk) = 0._wp
      !
      !
      DO jn = 1, kjpt         !==  loop over tracers  ==!
         !
         ztb(:,:,:) = ptb(:,:,:,jn)    ! initial before value for tracer
         ! 
         DO jl = 1, ksts         !==  Split-explicit loop  ==!
            !              
            DO jk = 2, jpk             ! 1st vertical derivative (w-flux)
               DO jj = 2, jpjm1 
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     IF( cdtype == 'TRA' .AND. jn == jp_tem ) THEN  ! temperature : use of avt
                        zwf(ji,jj,jk) =   avt(ji,jj,jk) * ( ztb(ji,jj,jk-1) - ztb(ji,jj,jk) ) / e3w_b(ji,jj,jk)
                     ELSE                                           ! salinity or pass. tracer : use of avs
                        zwf(ji,jj,jk) = fsavs(ji,jj,jk) * ( ztb(ji,jj,jk-1) - ztb(ji,jj,jk) ) / e3w_b(ji,jj,jk)
                     END IF
                  END DO
               END DO
            END DO
            !
            DO jk = 1, jpkm1           ! 2nd vertical derivative   ==> tracer at kt+l*2*rdt/nn_zdfexp
               DO jj = 2, jpjm1 
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     ztb(ji,jj,jk) = ztb(ji,jj,jk) + p2dt * ( zwf(ji,jj,jk) - zwf(ji,jj,jk+1) ) / e3t_n(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
         END DO                  ! end sub-time stepping

         DO jk = 1, jpkm1        !==  After tracer due to all trends
            DO jj = 2, jpjm1 
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ze3tb = e3t_b(ji,jj,jk) / e3t_n(ji,jj,jk)
                  ztra  = ( ztb(ji,jj,jk) - ptb(ji,jj,jk,jn) ) + p2dt * pta(ji,jj,jk,jn)  ! total trend * 2dt 
                  pta(ji,jj,jk,jn) = ( ze3tb * ptb(ji,jj,jk,jn) + ztra ) * tmask(ji,jj,jk)    ! after tracer
               END DO
            END DO
         END DO
         !
      END DO                     ! end of tracer loop
      !
      CALL wrk_dealloc( jpi,jpj,jpk,   ztb, zwf ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_zdf_exp')
      !
   END SUBROUTINE tra_zdf_exp

   !!==============================================================================
END MODULE trazdf_exp
