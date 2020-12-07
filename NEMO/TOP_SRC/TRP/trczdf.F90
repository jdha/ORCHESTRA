MODULE trczdf
   !!==============================================================================
   !!                 ***  MODULE  trczdf  ***
   !! Ocean Passive tracers : vertical diffusive trends 
   !!=====================================================================
   !! History :  9.0  ! 2005-11 (G. Madec)  Original code
   !!       NEMO 3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA 
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_zdf      : update the tracer trend with the lateral diffusion
   !!   trc_zdf_ini  : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE trc           ! ocean passive tracers variables
   USE oce_trc       ! ocean dynamics and active tracers
   USE trd_oce       ! trends: ocean variables
   USE trazdf_exp    ! vertical diffusion: explicit (tra_zdf_exp     routine)
   USE trazdf_imp    ! vertical diffusion: implicit (tra_zdf_imp     routine)
   USE trcldf        ! passive tracers: lateral diffusion
   USE trdtra        ! trends manager: tracers 
   USE prtctl_trc    ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_zdf         ! called by step.F90 
   PUBLIC   trc_zdf_ini     ! called by nemogcm.F90 
   
   !                                        !!** Vertical diffusion (nam_trczdf) **
   LOGICAL , PUBLIC ::   ln_trczdf_exp       !: explicit vertical diffusion scheme flag
   INTEGER , PUBLIC ::   nn_trczdf_exp       !: number of sub-time step (explicit time stepping)

   INTEGER ::   nzdf = 0               ! type vertical diffusion algorithm used
      !                                ! defined from ln_zdf...  namlist logicals)
   !! * Substitutions
#  include "zdfddm_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.7 , NEMO Consortium (2015)
   !! $Id: trczdf.F90 7753 2017-03-03 11:46:59Z mocavero $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_zdf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_zdf  ***
      !!
      !! ** Purpose :   compute the vertical ocean tracer physics.
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt      ! ocean time-step index
      !
      INTEGER               ::  jk, jn
      CHARACTER (len=22)    :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   ztrtrd   ! 4D workspace
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_zdf')
      !
      IF( l_trdtrc )  THEN
         CALL wrk_alloc( jpi, jpj, jpk, jptra, ztrtrd )
         ztrtrd(:,:,:,:)  = tra(:,:,:,:)
      ENDIF

      SELECT CASE ( nzdf )                       ! compute lateral mixing trend and add it to the general trend
      CASE ( 0 ) ;  CALL tra_zdf_exp( kt, nittrc000, 'TRC', r2dttrc, nn_trczdf_exp, trb, tra, jptra )    !   explicit scheme 
      CASE ( 1 ) ;  CALL tra_zdf_imp( kt, nittrc000, 'TRC', r2dttrc,                trb, tra, jptra )    !   implicit scheme          
      END SELECT

      IF( l_trdtrc )   THEN                      ! save the vertical diffusive trends for further diagnostics
         DO jn = 1, jptra
            DO jk = 1, jpkm1
               ztrtrd(:,:,jk,jn) = ( ( tra(:,:,jk,jn) - trb(:,:,jk,jn) ) / r2dttrc ) - ztrtrd(:,:,jk,jn)
            END DO
            CALL trd_tra( kt, 'TRC', jn, jptra_zdf, ztrtrd(:,:,:,jn) )
         END DO
         CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztrtrd )
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF( ln_ctl )   THEN
         WRITE(charout, FMT="('zdf ')") ;  CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      END IF
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_zdf')
      !
   END SUBROUTINE trc_zdf


   SUBROUTINE trc_zdf_ini
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE trc_zdf_ini  ***
      !!
      !! ** Purpose :   Choose the vertical mixing scheme
      !!
      !! ** Method  :   Set nzdf from ln_zdfexp
      !!      nzdf = 0   explicit (time-splitting) scheme (ln_trczdf_exp=T)
      !!           = 1   implicit (euler backward) scheme (ln_trczdf_exp=F)
      !!      NB: The implicit scheme is required when using : 
      !!             - rotated lateral mixing operator
      !!             - TKE, GLS vertical mixing scheme
      !!----------------------------------------------------------------------
      INTEGER ::  ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/namtrc_zdf/ ln_trczdf_exp  , nn_trczdf_exp
      !!----------------------------------------------------------------------
      !
      REWIND( numnat_ref )             ! namtrc_zdf in reference namelist 
      READ  ( numnat_ref, namtrc_zdf, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_zdf in reference namelist', lwp )
      !
      REWIND( numnat_cfg )             ! namtrc_zdf in configuration namelist 
      READ  ( numnat_cfg, namtrc_zdf, IOSTAT = ios, ERR = 906 )
906   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_zdf in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_zdf )
      !
      IF(lwp) THEN                     ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namtrc_zdf : set vertical diffusion  parameters'
         WRITE(numout,*) '      time splitting / backward scheme ln_trczdf_exp = ', ln_trczdf_exp
         WRITE(numout,*) '      number of time step              nn_trczdf_exp = ', nn_trczdf_exp
      ENDIF

      !                                ! Define the vertical tracer physics scheme
      IF( ln_trczdf_exp ) THEN   ;   nzdf = 0     ! explicit scheme
      ELSE                       ;   nzdf = 1     ! implicit scheme
      ENDIF

      !                                ! Force implicit schemes
      IF( ln_trcldf_iso              )   nzdf = 1      ! iso-neutral lateral physics
      IF( ln_trcldf_hor .AND. ln_sco )   nzdf = 1      ! horizontal lateral physics in s-coordinate
#if defined key_zdftke || defined key_zdfgls 
                                         nzdf = 1      ! TKE or GLS physics       
#endif
      IF( ln_trczdf_exp .AND. nzdf == 1 )  & 
         CALL ctl_stop( 'trc_zdf : If using the rotated lateral mixing operator or TKE, GLS vertical scheme ', &
            &           '          the implicit scheme is required, set ln_trczdf_exp = .false.' )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc:zdf_ctl : vertical passive tracer physics scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         IF( nzdf ==  0 )   WRITE(numout,*) '              Explicit time-splitting scheme'
         IF( nzdf ==  1 )   WRITE(numout,*) '              Implicit (euler backward) scheme'
      ENDIF
      !
   END SUBROUTINE trc_zdf_ini
   
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_zdf( kt )
      INTEGER, INTENT(in) :: kt  
      WRITE(*,*) 'trc_zdf: You should not have seen this print! error?', kt
   END SUBROUTINE trc_zdf
#endif
   !!==============================================================================
END MODULE trczdf
