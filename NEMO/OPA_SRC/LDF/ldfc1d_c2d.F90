MODULE ldfc1d_c2d
   !!======================================================================
   !!                    ***  MODULE  ldfc1d_c2d  ***
   !! Ocean physics:  profile and horizontal shape of lateral eddy coefficients 
   !!=====================================================================
   !! History :  3.7  ! 2013-12  (G. Madec)  restructuration/simplification of aht/aeiv specification,
   !!                 !                      add velocity dependent coefficient and optional read in file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_c1d       : ah reduced by 1/4 on the vertical (tanh profile, inflection at 300m) 
   !!   ldf_c2d       : ah = F(e1,e2) (laplacian or = F(e1^3,e2^3) (bilaplacian)
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_c1d   ! called by ldftra and ldfdyn modules
   PUBLIC   ldf_c2d   ! called by ldftra and ldfdyn modules

 
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2015)
   !! $Id: $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_c1d( cd_type, prat, pahs1, pahs2, pah1, pah2 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_c1d  ***
      !!              
      !! ** Purpose :   1D eddy diffusivity/viscosity coefficients
      !!
      !! ** Method  :   1D eddy diffusivity coefficients F( depth )
      !!                Reduction by prat from surface to bottom 
      !!                hyperbolic tangent profile with inflection point 
      !!                at zh=500m and a width of zw=200m
      !!
      !!   cd_type = TRA      pah1, pah2 defined at U- and V-points
      !!             DYN      pah1, pah2 defined at T- and F-points
      !!----------------------------------------------------------------------
      CHARACTER(len=2)                , INTENT(in   ) ::   cd_type        ! DYNamique or TRAcers
      REAL(wp)                        , INTENT(in   ) ::   prat           ! ratio surface/deep values           [-]
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(in   ) ::   pahs1, pahs2   ! surface value of eddy coefficient   [m2/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pah1 , pah2    ! eddy coefficient                    [m2/s]
      !
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   zh, zc, zdep1   ! local scalars
      REAL(wp) ::   zw    , zdep2   !   -      -
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   ldf_c1d : set a given profile to eddy diffusivity/viscosity coefficients'
         WRITE(numout,*) '   ~~~~~~~'
      ENDIF

      ! initialization of the profile
      zh =  500._wp              ! depth    of the inflection point [m]
      zw =  1._wp / 200._wp      ! width^-1     -        -      -   [1/m]
      !                          ! associated coefficient           [-]
      zc = ( 1._wp - prat ) / ( 1._wp + TANH( zh * zw) )
      !
      !
      SELECT CASE( cd_type )        ! point of calculation
      !
      CASE( 'DYN' )                     ! T- and F-points
         DO jk = 1, jpk                      ! pah1 at T-point
            pah1(:,:,jk) = pahs1(:,:) * (  prat + zc * ( 1._wp + TANH( - ( gdept_n(:,:,jk) - zh ) * zw) )  ) * tmask(:,:,jk)
         END DO
         DO jk = 1, jpk                      ! pah2 at F-point (zdep2 is an approximation in zps-coord.)
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  zdep2 = (  gdept_n(ji,jj+1,jk) + gdept_n(ji+1,jj+1,jk)   &
                     &     + gdept_n(ji,jj  ,jk) + gdept_n(ji+1,jj  ,jk)  ) * 0.25_wp
                  pah2(ji,jj,jk) = pahs2(ji,jj) * (  prat + zc * ( 1._wp + TANH( - ( zdep2 - zh ) * zw) )  ) * fmask(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL lbc_lnk( pah2, 'F', 1. )   ! Lateral boundary conditions
         !
      CASE( 'TRA' )                     ! U- and V-points (zdep1 & 2 are an approximation in zps-coord.)
         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  zdep1 = (  gdept_n(ji,jj,jk) + gdept_n(ji+1,jj,jk)  ) * 0.5_wp
                  zdep2 = (  gdept_n(ji,jj,jk) + gdept_n(ji,jj+1,jk)  ) * 0.5_wp
                  pah1(ji,jj,jk) = pahs1(ji,jj) * (  prat + zc * ( 1._wp + TANH( - ( zdep1 - zh ) * zw) )  ) * umask(ji,jj,jk)
                  pah2(ji,jj,jk) = pahs2(ji,jj) * (  prat + zc * ( 1._wp + TANH( - ( zdep2 - zh ) * zw) )  ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL lbc_lnk( pah1, 'U', 1. )   ! Lateral boundary conditions
         CALL lbc_lnk( pah2, 'V', 1. )   
         !
      END SELECT
      !
   END SUBROUTINE ldf_c1d


   SUBROUTINE ldf_c2d( cd_type, cd_op, pah0, pah1, pah2 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_c2d  ***
      !!              
      !! ** Purpose :   2D eddy diffusivity/viscosity coefficients
      !!
      !! ** Method  :   2D eddy diffusivity coefficients F( e1 , e2 )
      !!       laplacian   operator :   ah proportional to the scale factor      [m2/s]
      !!       bilaplacian operator :   ah proportional to the (scale factor)^3  [m4/s]
      !!       In both cases, pah0 is the maximum value reached by the coefficient 
      !!       at the Equator in case of e1=ra*rad= ~111km, not over the whole domain.
      !!
      !!   cd_type = TRA      pah1, pah2 defined at U- and V-points
      !!             DYN      pah1, pah2 defined at T- and F-points
      !!----------------------------------------------------------------------
      CHARACTER(len=3)                , INTENT(in   ) ::   cd_type      ! DYNamique or TRAcers
      CHARACTER(len=3)                , INTENT(in   ) ::   cd_op        ! operator: LAPlacian BiLaPlacian
      REAL(wp)                        , INTENT(in   ) ::   pah0         ! eddy coefficient   [m2/s or m4/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   pah1, pah2   ! eddy coefficient   [m2/s or m4/s]
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   za00, zd_max, zemax1, zemax2   ! local scalar
      !!----------------------------------------------------------------------
      !
      zd_max = ra * rad       ! = 1 degree at the equator in meters
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   ldf_c2d :     aht = rn_aht0 *  max(e1,e2)/e_equator     (  laplacian) '
         WRITE(numout,*) '   ~~~~~~~       or  = rn_bht0 * [max(e1,e2)/e_equator]**3 (bilaplacian)'
         WRITE(numout,*)
      ENDIF
      !
      SELECT CASE( cd_type )        !==  surface values  ==!  (depending on DYN/TRA)
      !
      CASE( 'DYN' )                     ! T- and F-points
         IF( cd_op == 'LAP' ) THEN            ! laplacian operator
            IF(lwp) WRITE(numout,*) '              momentum laplacian coeffcients = rn_aht0/e_equ * max(e1,e2)'
            za00 = pah0 / zd_max
            DO jj = 1, jpj 
               DO ji = 1, jpi 
                  zemax1 = MAX( e1t(ji,jj), e2t(ji,jj) ) * tmask(ji,jj,1)
                  zemax2 = MAX( e1f(ji,jj), e2f(ji,jj) ) * fmask(ji,jj,1)
                  pah1(ji,jj,1) = za00 * zemax1
                  pah2(ji,jj,1) = za00 * zemax2
               END DO
            END DO
         ELSEIF( cd_op == 'BLP' ) THEN     ! bilaplacian operator
            IF(lwp) WRITE(numout,*) '              momentum bilaplacian coeffcients = rn_bht0/e_equ * max(e1,e2)**3'
            za00 = pah0 / ( zd_max * zd_max * zd_max )
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zemax1 = MAX( e1t(ji,jj), e2t(ji,jj) ) * tmask(ji,jj,1)
                  zemax2 = MAX( e1f(ji,jj), e2f(ji,jj) ) * fmask(ji,jj,1)
                  pah1(ji,jj,1) = za00 * zemax1 * zemax1 * zemax1 
                  pah2(ji,jj,1) = za00 * zemax2 * zemax2 * zemax2 
               END DO
            END DO
         ELSE                                ! NO diffusion/viscosity
            CALL ctl_stop( 'ldf_c2d: ', cd_op, ' case. Unknown lateral operator ' )
         ENDIF
         !                                !  deeper values  (LAP and BLP cases)
         DO jk = 2, jpk
            pah1(:,:,jk) = pah1(:,:,1) * tmask(:,:,jk) 
            pah2(:,:,jk) = pah2(:,:,1) * fmask(:,:,jk) 
         END DO
         !
      CASE( 'TRA' )                     ! U- and V-points (approximation in zps-coord.)
         IF( cd_op == 'LAP' ) THEN            ! laplacian operator
            IF(lwp) WRITE(numout,*) '              tracer laplacian coeffcients = rn_aht0/e_equ * max(e1,e2)'
            za00 = pah0 / zd_max
            DO jj = 1, jpj 
               DO ji = 1, jpi 
                  zemax1 = MAX( e1u(ji,jj), e2u(ji,jj) ) * umask(ji,jj,1)
                  zemax2 = MAX( e1v(ji,jj), e2v(ji,jj) ) * vmask(ji,jj,1)
                  pah1(ji,jj,1) = za00 * zemax1
                  pah2(ji,jj,1) = za00 * zemax2
               END DO
            END DO
         ELSEIF( cd_op == 'BLP' ) THEN      ! bilaplacian operator (NB: square root of the coeff)
            IF(lwp) WRITE(numout,*) '              tracer bilaplacian coeffcients = rn_bht0/e_equ * max(e1,e2)**3'
            za00 = pah0 / ( zd_max * zd_max * zd_max )
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zemax1 = MAX( e1u(ji,jj), e2u(ji,jj) ) * umask(ji,jj,1) 
                  zemax2 = MAX( e1v(ji,jj), e2v(ji,jj) ) * vmask(ji,jj,1) 
                  pah1(ji,jj,1) = za00 * zemax1 * zemax1 * zemax1 
                  pah2(ji,jj,1) = za00 * zemax2 * zemax2 * zemax2 
               END DO
            END DO
         ELSE                                ! NO diffusion/viscosity
            CALL ctl_stop( 'ldf_c2d: ', cd_op, ' case. Unknown lateral operator ' )
         ENDIF
         !                                !  deeper values  (LAP and BLP cases)
         DO jk = 2, jpk
            pah1(:,:,jk) = pah1(:,:,1) * umask(:,:,jk) 
            pah2(:,:,jk) = pah2(:,:,1) * vmask(:,:,jk) 
         END DO
         !
      END SELECT
      !
   END SUBROUTINE ldf_c2d

   !!======================================================================
END MODULE ldfc1d_c2d
