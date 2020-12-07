MODULE usrdef_closea
   !!======================================================================
   !!                   ***  MODULE  usrdef_closea  ***
   !!
   !!                      ===  ORCA configuration  ===
   !!                         (2, 1 and 1/4 degrees)
   !!
   !! User define : specific treatments associated with closed seas
   !!======================================================================
   !! History :   8.2  !  2000-05  (O. Marti)  Original code
   !!   NEMO      1.0  !  2002-06  (E. Durand, G. Madec)  F90
   !!             3.0  !  2006-07  (G. Madec)  add clo_rnf, clo_ups, clo_bat
   !!             3.4  !  2014-12  (P.G. Fogli) sbc_clo bug fix & mpp reproducibility
   !!             4.0  !  2016-06  (G. Madec)  move to usrdef_closea, remove clo_ups
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_clo    : modification of the ocean domain for closed seas cases
   !!   sbc_clo    : Special handling of closed seas
   !!   clo_rnf    : set close sea outflows as river mouths (see sbcrnf)
   !!   clo_bat    : set to zero a field over closed sea (see domzrg)
   !!----------------------------------------------------------------------
   USE oce             ! dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! ocean surface boundary conditions
   !
   USE in_out_manager  ! I/O manager
   USE lib_fortran,    ONLY: glob_sum, DDPDD
   USE lbclnk          ! lateral boundary condition - MPP exchanges
   USE lib_mpp         ! MPP library
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC dom_clo      ! called by domain module
   PUBLIC sbc_clo      ! called by step module
   PUBLIC clo_rnf      ! called by sbcrnf module
   PUBLIC clo_bat      ! called in domzgr module

   INTEGER, PUBLIC, PARAMETER          ::   jpncs   = 4      !: number of closed sea
   INTEGER, PUBLIC, DIMENSION(jpncs)   ::   ncstt            !: Type of closed sea
   INTEGER, PUBLIC, DIMENSION(jpncs)   ::   ncsi1, ncsj1     !: south-west closed sea limits (i,j)
   INTEGER, PUBLIC, DIMENSION(jpncs)   ::   ncsi2, ncsj2     !: north-east closed sea limits (i,j)
   INTEGER, PUBLIC, DIMENSION(jpncs)   ::   ncsnr            !: number of point where run-off pours
   INTEGER, PUBLIC, DIMENSION(jpncs,4) ::   ncsir, ncsjr     !: Location of runoff

   REAL(wp), DIMENSION (jpncs+1)       ::   surf             ! closed sea surface

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id: usrdef_closea.F90 7754 2017-03-03 11:51:06Z mocavero $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_clo( cd_cfg, kcfg )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dom_clo  ***
      !!        
      !! ** Purpose :   Closed sea domain initialization
      !!
      !! ** Method  :   if a closed sea is located only in a model grid point
      !!                just the thermodynamic processes are applied.
      !!
      !! ** Action  :   ncsi1(), ncsj1() : south-west closed sea limits (i,j)
      !!                ncsi2(), ncsj2() : north-east Closed sea limits (i,j)
      !!                ncsir(), ncsjr() : Location of runoff
      !!                ncsnr            : number of point where run-off pours
      !!                ncstt            : Type of closed sea
      !!                                   =0 spread over the world ocean
      !!                                   =2 put at location runoff
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   ) ::   cd_cfg   ! configuration name
      INTEGER         , INTENT(in   ) ::   kcfg     ! configuration identifier 
      !
      INTEGER ::   jc      ! dummy loop indices
      INTEGER ::   isrow   ! local index
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*)'dom_clo : closed seas '
      IF(lwp) WRITE(numout,*)'~~~~~~~'
      !
      ! initial values
      ncsnr(:) = 1  ;  ncsi1(:) = 1  ;  ncsi2(:) = 1  ;  ncsir(:,:) = 1
      ncstt(:) = 0  ;  ncsj1(:) = 1  ;  ncsj2(:) = 1  ;  ncsjr(:,:) = 1
      !
      ! set the closed seas (in data domain indices)
      ! -------------------
      !
      IF( cd_cfg == "orca" ) THEN      !==  ORCA configuration  ==!
         !
         SELECT CASE ( kcfg )
         !                                           ! =======================
         CASE ( 1 )                                  !  ORCA_R1 configuration
            !                                        ! =======================
            IF(lwp) WRITE(numout,*)'   ORCA_R1 closed seas :  only the Caspian Sea'
            ! This dirty section will be suppressed by simplification process:
            ! all this will come back in input files
            ! Currently these hard-wired indices relate to configuration with
            ! extend grid (jpjglo=332)
            isrow = 332 - jpjglo
            !
            ncsnr(1)   = 1    ; ncstt(1)   = 0           ! Caspian Sea  (spread over the globe)
            ncsi1(1)   = 332  ; ncsj1(1)   = 243 - isrow
            ncsi2(1)   = 344  ; ncsj2(1)   = 275 - isrow
            ncsir(1,1) = 1    ; ncsjr(1,1) = 1
            !                                        
            !                                        ! =======================
         CASE ( 2 )                                  !  ORCA_R2 configuration
            !                                        ! =======================
            IF(lwp) WRITE(numout,*)'   ORCA_R2 closed seas and lakes : '
            !                                            ! Caspian Sea
            IF(lwp) WRITE(numout,*)'      Caspian Sea  '
            ncsnr(1)   =   1  ;  ncstt(1)   =   0           ! spread over the globe
            ncsi1(1)   =  11  ;  ncsj1(1)   = 103
            ncsi2(1)   =  17  ;  ncsj2(1)   = 112
            ncsir(1,1) =   1  ;  ncsjr(1,1) =   1 
            !                                            ! Great North American Lakes
            IF(lwp) WRITE(numout,*)'      Great North American Lakes  '
            ncsnr(2)   =   1  ;  ncstt(2)   =   2           ! put at St Laurent mouth
            ncsi1(2)   =  97  ;  ncsj1(2)   = 107
            ncsi2(2)   = 103  ;  ncsj2(2)   = 111
            ncsir(2,1) = 110  ;  ncsjr(2,1) = 111           
            !                                            ! Black Sea (crossed by the cyclic boundary condition)
            IF(lwp) WRITE(numout,*)'      Black Sea  '
            ncsnr(3:4) =   4  ;  ncstt(3:4) =   2           ! put in Med Sea (north of Aegean Sea)
            ncsir(3:4,1) = 171;  ncsjr(3:4,1) = 106         !
            ncsir(3:4,2) = 170;  ncsjr(3:4,2) = 106 
            ncsir(3:4,3) = 171;  ncsjr(3:4,3) = 105 
            ncsir(3:4,4) = 170;  ncsjr(3:4,4) = 105 
            ncsi1(3)   = 174  ;  ncsj1(3)   = 107           ! 1 : west part of the Black Sea      
            ncsi2(3)   = 181  ;  ncsj2(3)   = 112           !            (ie west of the cyclic b.c.)
            ncsi1(4)   =   2  ;  ncsj1(4)   = 107           ! 2 : east part of the Black Sea 
            ncsi2(4)   =   6  ;  ncsj2(4)   = 112           !           (ie east of the cyclic b.c.)
            !
            !                                        ! =========================
         CASE ( 025 )                                !  ORCA_R025 configuration
            !                                        ! =========================
            IF(lwp) WRITE(numout,*)'   ORCA_R025 closed seas : '
            !                                            ! Caspian Sea
            IF(lwp) WRITE(numout,*)'      Caspian Sea  '
            ncsnr(1)   = 1    ; ncstt(1)   = 0               ! Caspian + Aral sea
            ncsi1(1)   = 1330 ; ncsj1(1)   = 645
            ncsi2(1)   = 1400 ; ncsj2(1)   = 795
            ncsir(1,1) = 1    ; ncsjr(1,1) = 1
            !                                        
            IF(lwp) WRITE(numout,*)'      Azov Sea  '
            ncsnr(2)   = 1    ; ncstt(2)   = 0               ! Azov Sea 
            ncsi1(2)   = 1284 ; ncsj1(2)   = 722
            ncsi2(2)   = 1304 ; ncsj2(2)   = 747
            ncsir(2,1) = 1    ; ncsjr(2,1) = 1
            !
         END SELECT
         !
      ELSE                             !==  No closed sea in the configuration  ==!
         !
         IF(lwp) WRITE(numout,*)'   No closed seas or lakes in the configuration '
         !
      ENDIF

      ! convert the position in local domain indices
      ! --------------------------------------------
      DO jc = 1, jpncs
         ncsi1(jc)   = mi0( ncsi1(jc) )
         ncsj1(jc)   = mj0( ncsj1(jc) )
         !
         ncsi2(jc)   = mi1( ncsi2(jc) )   
         ncsj2(jc)   = mj1( ncsj2(jc) )  
      END DO
      !
   END SUBROUTINE dom_clo


   SUBROUTINE sbc_clo( kt, cd_cfg, kcfg )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_clo  ***
      !!                    
      !! ** Purpose :   Special handling of closed seas
      !!
      !! ** Method  :   Water flux is forced to zero over closed sea
      !!      Excess is shared between remaining ocean, or
      !!      put as run-off in open ocean.
      !!
      !! ** Action  :   emp updated surface freshwater fluxes and associated heat content at kt
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kt       ! ocean model time step
      CHARACTER(len=*), INTENT(in   ) ::   cd_cfg   ! configuration name
      INTEGER         , INTENT(in   ) ::   kcfg     ! configuration identifier 
      !
      INTEGER             ::   ji, jj, jc, jn   ! dummy loop indices
      REAL(wp), PARAMETER ::   rsmall = 1.e-20_wp    ! Closed sea correction epsilon
      REAL(wp)            ::   zze2, ztmp, zcorr     ! 
      REAL(wp)            ::   zcoef, zcoef1         ! 
      COMPLEX(wp)         ::   ctmp 
      REAL(wp), DIMENSION(jpncs) ::   zfwf   ! 1D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_clo')
      !                                                   !------------------!
      IF( kt == nit000 ) THEN                             !  Initialisation  !
         !                                                !------------------!
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'sbc_clo : closed seas '
         IF(lwp) WRITE(numout,*)'~~~~~~~'

         surf(:) = 0._wp
         !
         surf(jpncs+1) = glob_sum( e1e2t(:,:) )   ! surface of the global ocean
         !
         !                                        ! surface of closed seas 
         DO jc = 1, jpncs
            ctmp = CMPLX( 0.e0, 0.e0, wp )
            DO jj = ncsj1(jc), ncsj2(jc)
               DO ji = ncsi1(jc), ncsi2(jc)
                  ztmp = e1e2t(ji,jj) * tmask_i(ji,jj)
                  CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
               END DO 
            END DO 
            IF( lk_mpp )   CALL mpp_sum( ctmp )
            surf(jc) = REAL(ctmp,wp)
         END DO

         IF(lwp) WRITE(numout,*)'     Closed sea surfaces'
         DO jc = 1, jpncs
            IF(lwp)WRITE(numout,FMT='(1I3,4I4,5X,F16.2)') jc, ncsi1(jc), ncsi2(jc), ncsj1(jc), ncsj2(jc), surf(jc)
         END DO

         ! jpncs+1 : surface of sea, closed seas excluded
         DO jc = 1, jpncs
            surf(jpncs+1) = surf(jpncs+1) - surf(jc)
         END DO           
         !
      ENDIF
      !                                                   !--------------------!
      !                                                   !  update emp        !
      zfwf = 0.e0_wp                                      !--------------------!
      DO jc = 1, jpncs
         ctmp = CMPLX( 0.e0, 0.e0, wp )
         DO jj = ncsj1(jc), ncsj2(jc)
            DO ji = ncsi1(jc), ncsi2(jc)
               ztmp = e1e2t(ji,jj) * ( emp(ji,jj)-rnf(ji,jj) ) * tmask_i(ji,jj)
               CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
            END DO  
         END DO 
         IF( lk_mpp )   CALL mpp_sum( ctmp )
         zfwf(jc) = REAL(ctmp,wp)
      END DO

      IF( cd_cfg == "orca" .AND. kcfg == 2 ) THEN      ! Black Sea case for ORCA_R2 configuration
         zze2    = ( zfwf(3) + zfwf(4) ) * 0.5_wp
         zfwf(3) = zze2
         zfwf(4) = zze2
      ENDIF

      zcorr = 0._wp

      DO jc = 1, jpncs
         !
         ! The following if avoids the redistribution of the round off
         IF ( ABS(zfwf(jc) / surf(jpncs+1) ) > rsmall) THEN
            !
            IF( ncstt(jc) == 0 ) THEN           ! water/evap excess is shared by all open ocean
               zcoef    = zfwf(jc) / surf(jpncs+1)
               zcoef1   = rcp * zcoef
               emp(:,:) = emp(:,:) + zcoef
               qns(:,:) = qns(:,:) - zcoef1 * sst_m(:,:)
               ! accumulate closed seas correction
               zcorr    = zcorr    + zcoef
               !
            ELSEIF( ncstt(jc) == 1 ) THEN       ! Excess water in open sea, at outflow location, excess evap shared
               IF ( zfwf(jc) <= 0.e0_wp ) THEN 
                   DO jn = 1, ncsnr(jc)
                     ji = mi0(ncsir(jc,jn))
                     jj = mj0(ncsjr(jc,jn)) ! Location of outflow in open ocean
                     IF (      ji > 1 .AND. ji < jpi   &
                         .AND. jj > 1 .AND. jj < jpj ) THEN 
                         zcoef      = zfwf(jc) / ( REAL(ncsnr(jc)) * e1e2t(ji,jj) )
                         zcoef1     = rcp * zcoef
                         emp(ji,jj) = emp(ji,jj) + zcoef
                         qns(ji,jj) = qns(ji,jj) - zcoef1 * sst_m(ji,jj)
                     ENDIF 
                   END DO 
               ELSE 
                   zcoef    = zfwf(jc) / surf(jpncs+1)
                   zcoef1   = rcp * zcoef
                   emp(:,:) = emp(:,:) + zcoef
                   qns(:,:) = qns(:,:) - zcoef1 * sst_m(:,:)
                   ! accumulate closed seas correction
                   zcorr    = zcorr    + zcoef
               ENDIF
            ELSEIF( ncstt(jc) == 2 ) THEN       ! Excess e-p-r (either sign) goes to open ocean, at outflow location
               DO jn = 1, ncsnr(jc)
                  ji = mi0(ncsir(jc,jn))
                  jj = mj0(ncsjr(jc,jn)) ! Location of outflow in open ocean
                  IF(      ji > 1 .AND. ji < jpi    &
                     .AND. jj > 1 .AND. jj < jpj ) THEN 
                     zcoef      = zfwf(jc) / ( REAL(ncsnr(jc)) *  e1e2t(ji,jj) )
                     zcoef1     = rcp * zcoef
                     emp(ji,jj) = emp(ji,jj) + zcoef
                     qns(ji,jj) = qns(ji,jj) - zcoef1 * sst_m(ji,jj)
                  ENDIF 
               END DO 
            ENDIF 
            !
            DO jj = ncsj1(jc), ncsj2(jc)
               DO ji = ncsi1(jc), ncsi2(jc)
                  zcoef      = zfwf(jc) / surf(jc)
                  zcoef1     = rcp * zcoef
                  emp(ji,jj) = emp(ji,jj) - zcoef
                  qns(ji,jj) = qns(ji,jj) + zcoef1 * sst_m(ji,jj)
               END DO  
            END DO 
            !
         END IF
      END DO 

      IF ( ABS(zcorr) > rsmall ) THEN      ! remove the global correction from the closed seas
         DO jc = 1, jpncs                  ! only if it is large enough
            DO jj = ncsj1(jc), ncsj2(jc)
               DO ji = ncsi1(jc), ncsi2(jc)
                  emp(ji,jj) = emp(ji,jj) - zcorr
                  qns(ji,jj) = qns(ji,jj) + rcp * zcorr * sst_m(ji,jj)
               END DO  
             END DO 
          END DO
      ENDIF
      !
      emp (:,:) = emp (:,:) * tmask(:,:,1)
      !
      CALL lbc_lnk( emp , 'T', 1._wp )
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_clo')
      !
   END SUBROUTINE sbc_clo


   SUBROUTINE clo_rnf( p_rnfmsk )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf  ***
      !!                    
      !! ** Purpose :   allow the treatment of closed sea outflow grid-points
      !!                to be the same as river mouth grid-points
      !!
      !! ** Method  :   set to 1 the runoff mask (mskrnf, see sbcrnf module)
      !!                at the closed sea outflow grid-point.
      !!
      !! ** Action  :   update (p_)mskrnf (set 1 at closed sea outflow)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   p_rnfmsk   ! river runoff mask (rnfmsk array)
      !
      INTEGER  ::   jc, jn, ji, jj      ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      DO jc = 1, jpncs
         IF( ncstt(jc) >= 1 ) THEN            ! runoff mask set to 1 at closed sea outflows
             DO jn = 1, 4
                DO jj =    mj0( ncsjr(jc,jn) ), mj1( ncsjr(jc,jn) )
                   DO ji = mi0( ncsir(jc,jn) ), mi1( ncsir(jc,jn) )
                      p_rnfmsk(ji,jj) = MAX( p_rnfmsk(ji,jj), 1.0_wp )
                   END DO
                END DO
            END DO 
         ENDIF 
      END DO 
      !
   END SUBROUTINE clo_rnf
   
      
   SUBROUTINE clo_bat( k_top, k_bot )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE clo_bat  ***
      !!                    
      !! ** Purpose :   suppress closed sea from the domain
      !!
      !! ** Method  :   set first and last ocean level to 0 over the closed seas.
      !!
      !! ** Action  :   set pbat=0 and kbat=0 over closed seas
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(:,:), INTENT(inout) ::   k_top, k_bot   ! ocean first and last level indices
      !
      INTEGER  ::   jc, ji, jj      ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      DO jc = 1, jpncs
         DO jj = ncsj1(jc), ncsj2(jc)
            DO ji = ncsi1(jc), ncsi2(jc)
               k_top(ji,jj) = 0   
               k_bot(ji,jj) = 0   
            END DO 
         END DO 
       END DO 
       !
   END SUBROUTINE clo_bat

   !!======================================================================
END MODULE usrdef_closea

