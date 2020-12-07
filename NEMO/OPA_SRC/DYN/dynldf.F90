MODULE dynldf
   !!======================================================================
   !!                       ***  MODULE  dynldf  ***
   !! Ocean physics:  lateral diffusivity trends 
   !!=====================================================================
   !! History :  2.0  ! 2005-11  (G. Madec)  Original code (new step architecture)
   !!            3.7  ! 2014-01  (F. Lemarie, G. Madec)  restructuration/simplification of ahm specification,
   !!                 !                                  add velocity dependent coefficient and optional read in file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf      : update the dynamics trend with the lateral diffusion
   !!   dyn_ldf_init : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE ldfslp         ! lateral diffusion: slopes of mixing orientation
   USE dynldf_lap_blp ! lateral mixing   (dyn_ldf_lap & dyn_ldf_blp routines)
   USE dynldf_iso     ! lateral mixing                 (dyn_ldf_iso routine )
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics   (trd_dyn      routine)
   !
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_ldf       ! called by step module 
   PUBLIC   dyn_ldf_init  ! called by opa  module 

   !                      ! Parameter to control the type of lateral viscous operator
   INTEGER, PARAMETER, PUBLIC ::   np_ERROR  =-10   ! error in setting the operator
   INTEGER, PARAMETER, PUBLIC ::   np_no_ldf = 00   ! without operator (i.e. no lateral viscous trend)
   !                          !!      laplacian     !    bilaplacian    !
   INTEGER, PARAMETER, PUBLIC ::   np_lap    = 10   ,   np_blp    = 20  ! iso-level operator
   INTEGER, PARAMETER, PUBLIC ::   np_lap_i  = 11                       ! iso-neutral or geopotential operator

   INTEGER ::   nldf   ! type of lateral diffusion used defined from ln_dynldf_... (namlist logicals)

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2015)
   !! $Id: dynldf.F90 7753 2017-03-03 11:46:59Z mocavero $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf  ***
      !! 
      !! ** Purpose :   compute the lateral ocean dynamics physics.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_ldf')
      !
      IF( l_trddyn )   THEN                      ! temporary save of momentum trends
         CALL wrk_alloc( jpi,jpj,jpk,   ztrdu, ztrdv )
         ztrdu(:,:,:) = ua(:,:,:) 
         ztrdv(:,:,:) = va(:,:,:) 
      ENDIF

      SELECT CASE ( nldf )                       ! compute lateral mixing trend and add it to the general trend
      !
      CASE ( np_lap   )    ;   CALL dyn_ldf_lap  ( kt, ub, vb, ua, va, 1 )      ! iso-level    laplacian
      CASE ( np_lap_i )    ;   CALL dyn_ldf_iso  ( kt )                         ! rotated      laplacian
      CASE ( np_blp   )    ;   CALL dyn_ldf_blp  ( kt, ub, vb, ua, va    )      ! iso-level bi-laplacian
      !
      END SELECT

      IF( l_trddyn ) THEN                        ! save the horizontal diffusive trends for further diagnostics
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_ldf, kt )
         CALL wrk_dealloc( jpi,jpj,jpk,   ztrdu, ztrdv )
      ENDIF
      !                                          ! print sum trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_ldf')
      !
   END SUBROUTINE dyn_ldf


   SUBROUTINE dyn_ldf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf_init  ***
      !! 
      !! ** Purpose :   initializations of the horizontal ocean dynamics physics
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ierr         ! temporary integers 
      !!----------------------------------------------------------------------
      !
      !                                   ! Namelist nam_dynldf: already read in ldfdyn module
      !
      IF(lwp) THEN                        ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_ldf_init : Choice of the lateral diffusive operator on dynamics'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist nam_dynldf : set lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '      laplacian operator          ln_dynldf_lap = ', ln_dynldf_lap
         WRITE(numout,*) '      bilaplacian operator        ln_dynldf_blp = ', ln_dynldf_blp
         WRITE(numout,*) '      iso-level                   ln_dynldf_lev = ', ln_dynldf_lev
         WRITE(numout,*) '      horizontal (geopotential)   ln_dynldf_hor = ', ln_dynldf_hor
         WRITE(numout,*) '      iso-neutral                 ln_dynldf_iso = ', ln_dynldf_iso
      ENDIF
      !                                   ! use of lateral operator or not
      nldf = np_ERROR
      ioptio = 0
      IF( ln_dynldf_lap )   ioptio = ioptio + 1
      IF( ln_dynldf_blp )   ioptio = ioptio + 1
      IF( ioptio >  1   )   CALL ctl_stop( 'dyn_ldf_init: use ONE or NONE of the 2 lap/bilap operator type on momentum' )
      IF( ioptio == 0   )   nldf = np_no_ldf     ! No lateral mixing operator
      !
      IF( nldf /= np_no_ldf ) THEN        ! direction ==>> type of operator  
         ioptio = 0
         IF( ln_dynldf_lev )   ioptio = ioptio + 1
         IF( ln_dynldf_hor )   ioptio = ioptio + 1
         IF( ln_dynldf_iso )   ioptio = ioptio + 1
         IF( ioptio >  1   )   CALL ctl_stop( '          use only ONE direction (level/hor/iso)' )
         IF( ioptio == 0   )   CALL ctl_stop( '          use at least ONE direction (level/hor/iso)' )
         !
         !                                   ! Set nldf, the type of lateral diffusion, from ln_dynldf_... logicals
         ierr = 0
         IF ( ln_dynldf_lap ) THEN      ! laplacian operator
            IF ( ln_zco ) THEN                ! z-coordinate
               IF ( ln_dynldf_lev )   nldf = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_dynldf_hor )   nldf = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_dynldf_iso )   nldf = np_lap_i   ! iso-neutral            (   rotation)
            ENDIF
            IF ( ln_zps ) THEN             ! z-coordinate with partial step
               IF ( ln_dynldf_lev )   nldf = np_lap     ! iso-level              (no rotation)
               IF ( ln_dynldf_hor )   nldf = np_lap     ! iso-level              (no rotation)
               IF ( ln_dynldf_iso )   nldf = np_lap_i   ! iso-neutral            (   rotation)
            ENDIF
            IF ( ln_sco ) THEN             ! s-coordinate
               IF ( ln_dynldf_lev )   nldf = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_dynldf_hor )   nldf = np_lap_i   ! horizontal             (   rotation)
               IF ( ln_dynldf_iso )   nldf = np_lap_i   ! iso-neutral            (   rotation)
            ENDIF
         ENDIF
         !
         IF( ln_dynldf_blp ) THEN          ! bilaplacian operator
            IF ( ln_zco ) THEN                ! z-coordinate
               IF ( ln_dynldf_lev )   nldf = np_blp     ! iso-level = horizontal (no rotation)
               IF ( ln_dynldf_hor )   nldf = np_blp     ! iso-level = horizontal (no rotation)
               IF ( ln_dynldf_iso )   ierr = 2          ! iso-neutral            (   rotation)
            ENDIF
            IF ( ln_zps ) THEN             ! z-coordinate with partial step
               IF ( ln_dynldf_lev )   nldf = np_blp     ! iso-level              (no rotation)
               IF ( ln_dynldf_hor )   nldf = np_blp     ! iso-level              (no rotation)
               IF ( ln_dynldf_iso )   ierr = 2          ! iso-neutral            (   rotation)
            ENDIF
            IF ( ln_sco ) THEN             ! s-coordinate
               IF ( ln_dynldf_lev )   nldf = np_blp     ! iso-level              (no rotation)
               IF ( ln_dynldf_hor )   ierr = 2          ! horizontal             (   rotation)
               IF ( ln_dynldf_iso )   ierr = 2          ! iso-neutral            (   rotation)
            ENDIF
         ENDIF
         !
         IF( ierr == 2 )   CALL ctl_stop( 'rotated bi-laplacian operator does not exist' )
         !
         IF( nldf == np_lap_i )   l_ldfslp = .TRUE.      ! rotation require the computation of the slopes
         !
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*)
         IF( nldf == np_no_ldf )   WRITE(numout,*) '      ===>>   NO lateral viscosity'
         IF( nldf == np_lap    )   WRITE(numout,*) '      ===>>   iso-level laplacian operator'
         IF( nldf == np_lap_i  )   WRITE(numout,*) '      ===>>   rotated laplacian operator with iso-level background'
         IF( nldf == np_blp    )   WRITE(numout,*) '      ===>>   iso-level bi-laplacian operator'
      ENDIF
      !
   END SUBROUTINE dyn_ldf_init

   !!======================================================================
END MODULE dynldf
