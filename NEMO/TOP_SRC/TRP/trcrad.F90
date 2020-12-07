MODULE trcrad
   !!======================================================================
   !!                       ***  MODULE  trcrad  ***
   !! Ocean passive tracers:  correction of negative concentrations
   !!======================================================================
   !! History :   -   !  01-01  (O. Aumont & E. Kestenare)  Original code
   !!            1.0  !  04-03  (C. Ethe)  free form F90
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_rad    : correction of negative concentrations
   !!----------------------------------------------------------------------
   USE oce_trc             ! ocean dynamics and tracers variables
   USE trc                 ! ocean passive tracers variables
   USE trd_oce
   USE trdtra
   USE prtctl_trc          ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_rad     
   PUBLIC trc_rad_ini  

   LOGICAL , PUBLIC ::   ln_trcrad           !: flag to artificially correct negative concentrations

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcrad.F90 7753 2017-03-03 11:46:59Z mocavero $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   
CONTAINS

   SUBROUTINE trc_rad( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_rad  ***
      !!
      !! ** Purpose :   "crappy" routine to correct artificial negative
      !!              concentrations due to isopycnal scheme
      !!
      !! ** Method  : - PISCES or LOBSTER: Set negative concentrations to zero
      !!                while computing the corresponding tracer content that
      !!                is added to the tracers. Then, adjust the tracer 
      !!                concentration using a multiplicative factor so that 
      !!                the total tracer concentration is preserved.
      !!              - CFC: simply set to zero the negative CFC concentration
      !!                (the total CFC content is not strictly preserved)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index      
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_rad')
      !
      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_rad : Correct artificial negative concentrations '
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF

      IF( ln_age     )   CALL trc_rad_sms( kt, trb, trn, jp_age , jp_age               )  !  AGE
      IF( ll_cfc     )   CALL trc_rad_sms( kt, trb, trn, jp_cfc0, jp_cfc1               )  !  CFC model
      IF( ln_c14     )   CALL trc_rad_sms( kt, trb, trn, jp_c14 , jp_c14               )  !  C14
      IF( ln_pisces  )   CALL trc_rad_sms( kt, trb, trn, jp_pcs0, jp_pcs1, cpreserv='Y' )  !  PISCES model
      IF( ln_my_trc  )   CALL trc_rad_sms( kt, trb, trn, jp_myt0, jp_myt1               )  !  MY_TRC model

      !
      IF(ln_ctl) THEN      ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rad')")
         CALL prt_ctl_trc_info( charout )
         CALL prt_ctl_trc( tab4d=trn, mask=tmask, clinfo=ctrcnm )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_rad')
      !
   END SUBROUTINE trc_rad

   SUBROUTINE trc_rad_ini
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trc _rad_ini ***
      !!
      !! ** Purpose : read  namelist options 
      !!----------------------------------------------------------------------
      INTEGER ::  ios                 ! Local integer output status for namelist read
      NAMELIST/namtrc_rad/ ln_trcrad
      !!----------------------------------------------------------------------

      !
      REWIND( numnat_ref )              ! namtrc_rad in reference namelist 
      READ  ( numnat_ref, namtrc_rad, IOSTAT = ios, ERR = 907)
907   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_rad in reference namelist', lwp )

      REWIND( numnat_cfg )              ! namtrc_rad in configuration namelist 
      READ  ( numnat_cfg, namtrc_rad, IOSTAT = ios, ERR = 908 )
908   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_rad in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_rad )

      IF(lwp) THEN                     !   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namtrc_rad : treatment of negative concentrations'
         WRITE(numout,*) '      correct artificially negative concen. or not ln_trcrad = ', ln_trcrad
      ENDIF
      !
   END SUBROUTINE trc_rad_ini

   SUBROUTINE trc_rad_sms( kt, ptrb, ptrn, jp_sms0, jp_sms1, cpreserv )
      !!-----------------------------------------------------------------------------
      !!                  ***  ROUTINE trc_rad_sms  ***
      !!
      !! ** Purpose :   "crappy" routine to correct artificial negative
      !!              concentrations due to isopycnal scheme
      !!
      !! ** Method  : 2 cases :
      !!                - Set negative concentrations to zero while computing
      !!                  the corresponding tracer content that is added to the
      !!                  tracers. Then, adjust the tracer concentration using
      !!                  a multiplicative factor so that the total tracer 
      !!                  concentration is preserved.
      !!                - simply set to zero the negative CFC concentration
      !!                  (the total content of concentration is not strictly preserved)
      !!--------------------------------------------------------------------------------
      !! Arguments
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      INTEGER  , INTENT( in ) ::  &
         jp_sms0, &       !: First index of the passive tracer model
         jp_sms1          !: Last  index of  the passive tracer model

      REAL(wp), DIMENSION (jpi,jpj,jpk,jptra), INTENT( inout )  :: &
         ptrb, ptrn       !: before and now traceur concentration

      CHARACTER( len = 1) , INTENT(in), OPTIONAL  :: &
         cpreserv          !: flag to preserve content or not
      
      ! Local declarations
      INTEGER  :: ji, jj, jk, jn     ! dummy loop indices
      REAL(wp) :: ztrcorb, ztrmasb   ! temporary scalars
      REAL(wp) :: zcoef, ztrcorn, ztrmasn   !    "         "
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ztrtrdb, ztrtrdn   ! workspace arrays
      REAL(wp) :: zs2rdt
      LOGICAL ::   lldebug = .FALSE.
      !!----------------------------------------------------------------------

 
      IF( l_trdtrc )  CALL wrk_alloc( jpi, jpj, jpk, ztrtrdb, ztrtrdn )
      
      IF( PRESENT( cpreserv )  ) THEN   !  total tracer concentration is preserved 
      
         DO jn = jp_sms0, jp_sms1
            !                                                        ! ===========
            ztrcorb = 0.e0   ;   ztrmasb = 0.e0
            ztrcorn = 0.e0   ;   ztrmasn = 0.e0

            IF( l_trdtrc ) THEN
               ztrtrdb(:,:,:) = ptrb(:,:,:,jn)                        ! save input trb for trend computation
               ztrtrdn(:,:,:) = ptrn(:,:,:,jn)                        ! save input trn for trend computation
            ENDIF
            !                                                         ! sum over the global domain 
            ztrcorb = glob_sum( MIN( 0., ptrb(:,:,:,jn) ) * cvol(:,:,:) )
            ztrcorn = glob_sum( MIN( 0., ptrn(:,:,:,jn) ) * cvol(:,:,:) )

            ztrmasb = glob_sum( MAX( 0., ptrb(:,:,:,jn) ) * cvol(:,:,:) )
            ztrmasn = glob_sum( MAX( 0., ptrn(:,:,:,jn) ) * cvol(:,:,:) )

            IF( ztrcorb /= 0 ) THEN
               zcoef = 1. + ztrcorb / ztrmasb
               DO jk = 1, jpkm1
                  ptrb(:,:,jk,jn) = MAX( 0., ptrb(:,:,jk,jn) )
                  ptrb(:,:,jk,jn) = ptrb(:,:,jk,jn) * zcoef * tmask(:,:,jk)
               END DO
            ENDIF

            IF( ztrcorn /= 0 ) THEN
               zcoef = 1. + ztrcorn / ztrmasn
               DO jk = 1, jpkm1
                  ptrn(:,:,jk,jn) = MAX( 0., ptrn(:,:,jk,jn) )
                  ptrn(:,:,jk,jn) = ptrn(:,:,jk,jn) * zcoef * tmask(:,:,jk)
               END DO
            ENDIF
            !
            IF( l_trdtrc ) THEN
               !
               zs2rdt = 1. / ( 2. * rdt )
               ztrtrdb(:,:,:) = ( ptrb(:,:,:,jn) - ztrtrdb(:,:,:) ) * zs2rdt
               ztrtrdn(:,:,:) = ( ptrn(:,:,:,jn) - ztrtrdn(:,:,:) ) * zs2rdt 
               CALL trd_tra( kt, 'TRC', jn, jptra_radb, ztrtrdb )       ! Asselin-like trend handling
               CALL trd_tra( kt, 'TRC', jn, jptra_radn, ztrtrdn )       ! standard     trend handling
              !
            ENDIF

         END DO
         !
         !
      ELSE  ! total CFC content is not strictly preserved

         DO jn = jp_sms0, jp_sms1  

           IF( l_trdtrc ) THEN
              ztrtrdb(:,:,:) = ptrb(:,:,:,jn)                        ! save input trb for trend computation
              ztrtrdn(:,:,:) = ptrn(:,:,:,jn)                        ! save input trn for trend computation
           ENDIF

            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     ptrn(ji,jj,jk,jn) = MAX( 0. , ptrn(ji,jj,jk,jn) )
                     ptrb(ji,jj,jk,jn) = MAX( 0. , ptrb(ji,jj,jk,jn) )
                  END DO
               END DO
            END DO
         
            IF( l_trdtrc ) THEN
               !
               zs2rdt = 1. / ( 2. * rdt * REAL( nn_dttrc, wp ) )
               ztrtrdb(:,:,:) = ( ptrb(:,:,:,jn) - ztrtrdb(:,:,:) ) * zs2rdt
               ztrtrdn(:,:,:) = ( ptrn(:,:,:,jn) - ztrtrdn(:,:,:) ) * zs2rdt 
               CALL trd_tra( kt, 'TRC', jn, jptra_radb, ztrtrdb )       ! Asselin-like trend handling
               CALL trd_tra( kt, 'TRC', jn, jptra_radn, ztrtrdn )       ! standard     trend handling
              !
            ENDIF
            !
         ENDDO

      ENDIF

      IF( l_trdtrc )  CALL wrk_dealloc( jpi, jpj, jpk, ztrtrdb, ztrtrdn )

   END SUBROUTINE trc_rad_sms
#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                                         NO TOP model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rad( kt )              ! Empty routine
      INTEGER, INTENT(in) ::   kt
      WRITE(*,*) 'trc_rad: You should not have seen this print! error?', kt
   END SUBROUTINE trc_rad
#endif
   
   !!======================================================================
END MODULE trcrad
