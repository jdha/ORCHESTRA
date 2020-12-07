MODULE dynzdf
   !!==============================================================================
   !!                 ***  MODULE  dynzdf  ***
   !! Ocean dynamics :  vertical component of the momentum mixing trend
   !!==============================================================================
   !! History :  1.0  !  2005-11  (G. Madec)  Original code
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf       : Update the momentum trend with the vertical diffusion
   !!   dyn_zdf_init  : initializations of the vertical diffusion scheme
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables 
   USE zdf_oce        ! ocean vertical physics variables
   USE dynzdf_exp     ! vertical diffusion: explicit (dyn_zdf_exp     routine)
   USE dynzdf_imp     ! vertical diffusion: implicit (dyn_zdf_imp     routine)
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf       !  routine called by step.F90
   PUBLIC   dyn_zdf_init  !  routine called by opa.F90

   INTEGER  ::   nzdf = 0   ! type vertical diffusion algorithm used, defined from ln_zdf... namlist logicals

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynzdf.F90 7753 2017-03-03 11:46:59Z mocavero $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE dyn_zdf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf  ***
      !!
      !! ** Purpose :   compute the vertical ocean dynamics physics.
      !!---------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('dyn_zdf')
      !
      !                                          ! set time step
      IF( neuler == 0 .AND. kt == nit000     ) THEN   ;   r2dt =      rdt   ! = rdt (restart with Euler time stepping)
      ELSEIF(               kt <= nit000 + 1 ) THEN   ;   r2dt = 2. * rdt   ! = 2 rdt (leapfrog)
      ENDIF

      IF( l_trddyn )   THEN                      ! temporary save of ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdu, ztrdv ) 
         ztrdu(:,:,:) = ua(:,:,:)
         ztrdv(:,:,:) = va(:,:,:)
      ENDIF

      SELECT CASE ( nzdf )                       ! compute lateral mixing trend and add it to the general trend
      !
      CASE ( 0 )   ;   CALL dyn_zdf_exp( kt, r2dt )      ! explicit scheme
      CASE ( 1 )   ;   CALL dyn_zdf_imp( kt, r2dt )      ! implicit scheme
      !
      END SELECT

      IF( l_trddyn )   THEN                      ! save the vertical diffusive trends for further diagnostics
         ztrdu(:,:,:) = ( ua(:,:,:) - ub(:,:,:) ) / r2dt - ztrdu(:,:,:)
         ztrdv(:,:,:) = ( va(:,:,:) - vb(:,:,:) ) / r2dt - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_zdf, kt )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdu, ztrdv ) 
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask,               &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         !
      IF( nn_timing == 1 )   CALL timing_stop('dyn_zdf')
      !
   END SUBROUTINE dyn_zdf


   SUBROUTINE dyn_zdf_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_zdf_init  ***
      !!
      !! ** Purpose :   initializations of the vertical diffusion scheme
      !!
      !! ** Method  :   implicit (euler backward) scheme (default)
      !!                explicit (time-splitting) scheme if ln_zdfexp=T
      !!----------------------------------------------------------------------
      USE zdftke
      USE zdfgls
      !!----------------------------------------------------------------------
      !
      ! Choice from ln_zdfexp read in namelist in zdfini
      IF( ln_zdfexp ) THEN   ;   nzdf = 0           ! use explicit scheme
      ELSE                   ;   nzdf = 1           ! use implicit scheme
      ENDIF
      !
      ! Force implicit schemes
      IF( lk_zdftke .OR. lk_zdfgls   )   nzdf = 1   ! TKE or GLS physics
      IF( ln_dynldf_iso              )   nzdf = 1   ! iso-neutral lateral physics
      IF( ln_dynldf_hor .AND. ln_sco )   nzdf = 1   ! horizontal lateral physics in s-coordinate
      !
      IF(lwp) THEN                                  ! Print the choice
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_zdf_init : vertical dynamics physics scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         IF( nzdf ==  0 )   WRITE(numout,*) '      ===>>   Explicit time-splitting scheme'
         IF( nzdf ==  1 )   WRITE(numout,*) '      ===>>   Implicit (euler backward) scheme'
      ENDIF
      !
   END SUBROUTINE dyn_zdf_init

   !!==============================================================================
END MODULE dynzdf
