MODULE trcldf
   !!======================================================================
   !!                       ***  MODULE  trcldf  ***
   !! Ocean Passive tracers : lateral diffusive trends
   !!=====================================================================
   !! History :  1.0  ! 2005-11  (G. Madec)  Original code
   !!            3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA
   !!            3.7  ! 2014-03  (G. Madec)  LDF simplification
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_ldf       : update the tracer trend with the lateral diffusion
   !!   trc_ldf_ini   : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE trc            ! ocean passive tracers variables
   USE oce_trc        ! ocean dynamics and active tracers
   USE ldfslp         ! lateral diffusion: iso-neutral slope
   USE traldf_lap_blp ! lateral diffusion: lap/bilaplacian iso-level      operator  (tra_ldf_lap/_blp   routine)
   USE traldf_iso     ! lateral diffusion: laplacian iso-neutral standard operator  (tra_ldf_iso        routine)
   USE traldf_triad   ! lateral diffusion: laplacian iso-neutral triad    operator  (tra_ldf_     triad routine)
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers
   !
   USE prtctl_trc     ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ldf    
   PUBLIC   trc_ldf_ini   
   !
   LOGICAL , PUBLIC ::   ln_trcldf_lap       !:   laplacian operator
   LOGICAL , PUBLIC ::   ln_trcldf_blp       !: bilaplacian operator
   LOGICAL , PUBLIC ::   ln_trcldf_lev       !: iso-level   direction
   LOGICAL , PUBLIC ::   ln_trcldf_hor       !: horizontal  direction (rotation to geopotential)
   LOGICAL , PUBLIC ::   ln_trcldf_iso       !: iso-neutral direction (standard)
   LOGICAL , PUBLIC ::   ln_trcldf_triad     !: iso-neutral direction (triad)
   REAL(wp), PUBLIC ::   rn_ahtrc_0          !:   laplacian diffusivity coefficient for passive tracer [m2/s]
   REAL(wp), PUBLIC ::   rn_bhtrc_0          !: bilaplacian      -          --     -       -   [m4/s]
   REAL(wp), PUBLIC ::   rn_fact_lap         !: Enhanced zonal diffusivity coefficent in the equatorial domain
   !
   !                      !!: ** lateral mixing namelist (nam_trcldf) **
   REAL(wp) ::  rldf       ! ratio between active and passive tracers diffusive coefficient
   
   INTEGER  ::  nldf = 0   ! type of lateral diffusion used defined from ln_trcldf_... namlist logicals)
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.7 , NEMO Consortium (2014)
   !! $Id: trcldf.F90 7753 2017-03-03 11:46:59Z mocavero $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ldf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf  ***
      !!
      !! ** Purpose :   compute the lateral ocean tracer physics.
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !
      INTEGER            :: ji, jj, jk, jn
      REAL(wp)           :: zdep
      CHARACTER (len=22) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zahu, zahv
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   ztrtrd
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_ldf')
      !
      IF( l_trdtrc )  THEN
         CALL wrk_alloc( jpi,jpj,jpk,jptra,   ztrtrd )
         ztrtrd(:,:,:,:)  = tra(:,:,:,:)
      ENDIF
      !                                  !* set the lateral diffusivity coef. for passive tracer      
      CALL wrk_alloc( jpi,jpj,jpk,   zahu, zahv )
      zahu(:,:,:) = rldf * ahtu(:,:,:) 
      zahv(:,:,:) = rldf * ahtv(:,:,:)
      !                                  !* Enhanced zonal diffusivity coefficent in the equatorial domain
      DO jk= 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( gdept_n(ji,jj,jk) > 200. .AND. gphit(ji,jj) < 5. .AND. gphit(ji,jj) > -5. ) THEN
                  zdep = MAX( gdept_n(ji,jj,jk) - 1000., 0. ) / 1000.
                  zahu(ji,jj,jk) = zahu(ji,jj,jk) * MAX( 1., rn_fact_lap * EXP( -zdep ) )
               ENDIF
            END DO
         END DO
      END DO
      !
      SELECT CASE ( nldf )                     !* compute lateral mixing trend and add it to the general trend
      !
      CASE ( np_lap   )                               ! iso-level laplacian
         CALL tra_ldf_lap  ( kt, nittrc000,'TRC', zahu, zahv, gtru, gtrv, gtrui, gtrvi, trb,      tra, jptra,  1   )
         !
      CASE ( np_lap_i )                               ! laplacian : standard iso-neutral operator (Madec)
         CALL tra_ldf_iso  ( kt, nittrc000,'TRC', zahu, zahv, gtru, gtrv, gtrui, gtrvi, trb, trb, tra, jptra,  1   )
         !
      CASE ( np_lap_it )                              ! laplacian : triad iso-neutral operator (griffies)
         CALL tra_ldf_triad( kt, nittrc000,'TRC', zahu, zahv, gtru, gtrv, gtrui, gtrvi, trb, trb, tra, jptra,  1   )
         !
      CASE ( np_blp , np_blp_i , np_blp_it )          ! bilaplacian: all operator (iso-level, -neutral)
         CALL tra_ldf_blp  ( kt, nittrc000,'TRC', zahu, zahv, gtru, gtrv, gtrui, gtrvi, trb     , tra, jptra, nldf )
         !
      END SELECT
      !
      IF( l_trdtrc )   THEN                    ! send the trends for further diagnostics
        DO jn = 1, jptra
           ztrtrd(:,:,:,jn) = tra(:,:,:,jn) - ztrtrd(:,:,:,jn)
           CALL trd_tra( kt, 'TRC', jn, jptra_ldf, ztrtrd(:,:,:,jn) )
        END DO
        CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztrtrd )
      ENDIF
      !                
      IF( ln_ctl ) THEN                        ! print mean trends (used for debugging)
         WRITE(charout, FMT="('ldf ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj,jpk,   zahu, zahv )
      !
      IF( nn_timing == 1 )   CALL timing_stop('trc_ldf')
      !
   END SUBROUTINE trc_ldf


   SUBROUTINE trc_ldf_ini
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_ctl  ***
      !!
      !! ** Purpose :   Define the operator for the lateral diffusion
      !!
      !! ** Method  :   set nldf from the namtra_ldf logicals
      !!      nldf ==  0   laplacian operator
      !!      nldf ==  1   Rotated laplacian operator
      !!      nldf ==  2   bilaplacian operator
      !!      nldf ==  3   Rotated bilaplacian
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ierr   ! temporary integers
      INTEGER ::   ios            ! Local integer output status for namelist read
      !!
      NAMELIST/namtrc_ldf/ ln_trcldf_lap, ln_trcldf_blp,                                  &
         &                 ln_trcldf_lev, ln_trcldf_hor, ln_trcldf_iso, ln_trcldf_triad,  &
         &                 rn_ahtrc_0   , rn_bhtrc_0, rn_fact_lap  
      !!----------------------------------------------------------------------
      !
      REWIND( numnat_ref )             !  namtrc_ldf in reference namelist 
      READ  ( numnat_ref, namtrc_ldf, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrc_ldf in reference namelist', lwp )
      !
      REWIND( numnat_cfg )             !  namtrc_ldf in configuration namelist 
      READ  ( numnat_cfg, namtrc_ldf, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrc_ldf in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_ldf )
      !
      IF(lwp) THEN                     ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'trc_ldf_ini : lateral tracer diffusive operator'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtrc_ldf : set lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '      operator'
         WRITE(numout,*) '           laplacian                 ln_trcldf_lap   = ', ln_trcldf_lap
         WRITE(numout,*) '         bilaplacian                 ln_trcldf_blp   = ', ln_trcldf_blp
         WRITE(numout,*) '      direction of action'
         WRITE(numout,*) '         iso-level                   ln_trcldf_lev   = ', ln_trcldf_lev
         WRITE(numout,*) '         horizontal (geopotential)   ln_trcldf_hor   = ', ln_trcldf_hor
         WRITE(numout,*) '         iso-neutral (standard)      ln_trcldf_iso   = ', ln_trcldf_iso
         WRITE(numout,*) '         iso-neutral (triad)         ln_trcldf_triad = ', ln_trcldf_triad
         WRITE(numout,*) '      diffusivity coefficient'
         WRITE(numout,*) '           laplacian                 rn_ahtrc_0      = ', rn_ahtrc_0
         WRITE(numout,*) '         bilaplacian                 rn_bhtrc_0      = ', rn_bhtrc_0
         WRITE(numout,*) '      enhanced zonal diffusivity     rn_fact_lap     = ', rn_fact_lap

      ENDIF
      !      
      !                                ! control the namelist parameters
      ioptio = 0
      IF( ln_trcldf_lap )   ioptio = ioptio + 1
      IF( ln_trcldf_blp )   ioptio = ioptio + 1
      IF( ioptio >  1   )   CALL ctl_stop( 'trc_ldf_ctl: use ONE or NONE of the 2 lap/bilap operator type on tracer' )
      IF( ioptio == 0   )   nldf = np_no_ldf   ! No lateral diffusion
      
      IF( ln_trcldf_lap .AND. ln_trcldf_blp )   CALL ctl_stop( 'trc_ldf_ctl: bilaplacian should be used on both TRC and TRA' )
      IF( ln_trcldf_blp .AND. ln_trcldf_lap )   CALL ctl_stop( 'trc_ldf_ctl:   laplacian should be used on both TRC and TRA' )
      !
      ioptio = 0
      IF( ln_trcldf_lev )   ioptio = ioptio + 1
      IF( ln_trcldf_hor )   ioptio = ioptio + 1
      IF( ln_trcldf_iso )   ioptio = ioptio + 1
      IF( ioptio /= 1   )   CALL ctl_stop( 'trc_ldf_ctl: use only ONE direction (level/hor/iso)' )
      !
      ! defined the type of lateral diffusion from ln_trcldf_... logicals
      ! CAUTION : nldf = 1 is used in trazdf_imp, change it carefully
      ierr = 0
      IF( ln_trcldf_lap ) THEN      !==  laplacian operator  ==!
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_trcldf_lev   )   nldf = np_lap     ! iso-level = horizontal (no rotation)
            IF ( ln_trcldf_hor   )   nldf = np_lap     ! iso-level = horizontal (no rotation)
            IF ( ln_trcldf_iso   )   nldf = np_lap_i   ! iso-neutral: standard  (   rotation)
            IF ( ln_trcldf_triad )   nldf = np_lap_it  ! iso-neutral: triad     (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate with partial step
            IF ( ln_trcldf_lev   )   ierr = 1         ! iso-level not allowed 
            IF ( ln_trcldf_hor   )   nldf = np_lap     ! horizontal (no rotation)
            IF ( ln_trcldf_iso   )   nldf = np_lap_i   ! iso-neutral: standard (rotation)
            IF ( ln_trcldf_triad )   nldf = np_lap_it  ! iso-neutral: triad    (rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! s-coordinate
            IF ( ln_trcldf_lev   )   nldf = np_lap     ! iso-level  (no rotation)
            IF ( ln_trcldf_hor   )   nldf = np_lap_it  ! horizontal (   rotation)       !!gm   a checker....
            IF ( ln_trcldf_iso   )   nldf = np_lap_i   ! iso-neutral: standard (rotation)
            IF ( ln_trcldf_triad )   nldf = np_lap_it  ! iso-neutral: triad    (rotation)
         ENDIF
         !                                ! diffusivity ratio: passive / active tracers 
         IF( ABS(rn_aht_0) < 2._wp*TINY(1._wp) ) THEN
            IF( ABS(rn_ahtrc_0) < 2._wp*TINY(1._wp) ) THEN
               rldf = 1.0_wp
            ELSE
               CALL ctl_stop( 'trc_ldf_ctl : cannot define rldf, rn_aht_0==0, rn_ahtrc_0 /=0' )
            ENDIF
         ELSE
            rldf = rn_ahtrc_0 / rn_aht_0
         ENDIF
      ENDIF
      !
      IF( ln_trcldf_blp ) THEN      !==  bilaplacian operator  ==!
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_trcldf_lev   )   nldf = np_blp     ! iso-level = horizontal (no rotation)
            IF ( ln_trcldf_hor   )   nldf = np_blp     ! iso-level = horizontal (no rotation)
            IF ( ln_trcldf_iso   )   nldf = np_blp_i   ! iso-neutral: standard (rotation)
            IF ( ln_trcldf_triad )   nldf = np_blp_it  ! iso-neutral: triad    (rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate with partial step
            IF ( ln_trcldf_lev   )   ierr = 1         ! iso-level not allowed 
            IF ( ln_trcldf_hor   )   nldf = np_blp     ! horizontal (no rotation)
            IF ( ln_trcldf_iso   )   nldf = np_blp_i   ! iso-neutral: standard (rotation)
            IF ( ln_trcldf_triad )   nldf = np_blp_it  ! iso-neutral: triad    (rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! s-coordinate
            IF ( ln_trcldf_lev   )   nldf = np_blp     ! iso-level  (no rotation)
            IF ( ln_trcldf_hor   )   nldf = np_blp_it  ! horizontal (   rotation)       !!gm   a checker....
            IF ( ln_trcldf_iso   )   nldf = np_blp_i   ! iso-neutral: standard (rotation)
            IF ( ln_trcldf_triad )   nldf = np_blp_it  ! iso-neutral: triad    (rotation)
         ENDIF
         !                                ! diffusivity ratio: passive / active tracers 
         IF( ABS(rn_bht_0) < 2._wp*TINY(1._wp) ) THEN
            IF( ABS(rn_bhtrc_0) < 2._wp*TINY(1._wp) ) THEN
               rldf = 1.0_wp
            ELSE
               CALL ctl_stop( 'trc_ldf_ctl : cannot define rldf, rn_aht_0==0, rn_ahtrc_0 /=0' )
            ENDIF
         ELSE
            rldf = SQRT(  ABS( rn_bhtrc_0 / rn_bht_0 )  )
         ENDIF
      ENDIF
      !
      IF( ierr == 1 )   CALL ctl_stop( 'trc_ldf_ctl: iso-level in z-partial step, not allowed' )
      IF( ln_ldfeiv .AND. .NOT.ln_trcldf_iso )   CALL ctl_stop( 'trc_ldf_ctl: eiv requires isopycnal laplacian diffusion' )
      IF( nldf == 1 .OR. nldf == 3 )   l_ldfslp = .TRUE.    ! slope of neutral surfaces required 
      !
      IF(lwp) THEN
         WRITE(numout,*)
         SELECT CASE( nldf )
         CASE( np_no_ldf )   ;   WRITE(numout,*) '          NO lateral diffusion'
         CASE( np_lap    )   ;   WRITE(numout,*) '          laplacian iso-level operator'
         CASE( np_lap_i  )   ;   WRITE(numout,*) '          Rotated laplacian operator (standard)'
         CASE( np_lap_it )   ;   WRITE(numout,*) '          Rotated laplacian operator (triad)'
         CASE( np_blp    )   ;   WRITE(numout,*) '          bilaplacian iso-level operator'
         CASE( np_blp_i  )   ;   WRITE(numout,*) '          Rotated bilaplacian operator (standard)'
         CASE( np_blp_it )   ;   WRITE(numout,*) '          Rotated bilaplacian operator (triad)'
         END SELECT
      ENDIF
      !
   END SUBROUTINE trc_ldf_ini
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ldf( kt )
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_ldf: You should not have seen this print! error?', kt
   END SUBROUTINE trc_ldf
#endif
   !!======================================================================
END MODULE trcldf
