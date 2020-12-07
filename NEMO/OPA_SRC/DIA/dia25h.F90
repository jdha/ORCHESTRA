MODULE dia25h 
   !!======================================================================
   !!                       ***  MODULE  diaharm  ***
   !! Harmonic analysis of tidal constituents 
   !!======================================================================
   !! History :  3.6  !  2014  (E O'Dea)  Original code
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library
   USE wrk_nemo        ! working arrays
#if defined key_zdftke 
   USE zdf_oce, ONLY: en
#endif
   USE zdf_oce, ONLY: avt, avm
#if defined key_zdfgls
   USE zdf_oce, ONLY: en
   USE zdfgls, ONLY: mxln
#endif

   IMPLICIT NONE
   PRIVATE

   LOGICAL , PUBLIC ::   ln_dia25h     !:  25h mean output
   PUBLIC   dia_25h_init               ! routine called by nemogcm.F90
   PUBLIC   dia_25h                    ! routine called by diawri.F90

  !! * variables for calculating 25-hourly means
   REAL(wp),SAVE, ALLOCATABLE,   DIMENSION(:,:,:) ::   tn_25h  , sn_25h
   REAL(wp),SAVE, ALLOCATABLE,   DIMENSION(:,:)   ::   sshn_25h 
   REAL(wp),SAVE, ALLOCATABLE,   DIMENSION(:,:,:) ::   un_25h  , vn_25h  , wn_25h
   REAL(wp),SAVE, ALLOCATABLE,   DIMENSION(:,:,:) ::   avt_25h , avm_25h
#if defined key_zdfgls || key_zdftke
   REAL(wp),SAVE, ALLOCATABLE,   DIMENSION(:,:,:) ::   en_25h
#endif
#if defined key_zdfgls 
   REAL(wp),SAVE, ALLOCATABLE,   DIMENSION(:,:,:) ::   rmxln_25h
#endif
   INTEGER, SAVE :: cnt_25h     ! Counter for 25 hour means



   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.6 , NEMO Consortium (2014)
   !! $Id:$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_25h_init 
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_25h_init  ***
      !!     
      !! ** Purpose: Initialization of 25h mean namelist 
      !!        
      !! ** Method : Read namelist
      !!   History
      !!   3.6  !  08-14  (E. O'Dea) Routine to initialize dia_25h
      !!---------------------------------------------------------------------------
      !!
      INTEGER ::   ios                 ! Local integer output status for namelist read
      INTEGER ::   ierror              ! Local integer for memory allocation
      !
      NAMELIST/nam_dia25h/ ln_dia25h
      !!----------------------------------------------------------------------
      !
      REWIND ( numnam_ref )              ! Read Namelist nam_dia25h in reference namelist : 25hour mean diagnostics
      READ   ( numnam_ref, nam_dia25h, IOSTAT=ios, ERR= 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_dia25h in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nam_dia25h in configuration namelist  25hour diagnostics
      READ  ( numnam_cfg, nam_dia25h, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_dia25h in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_dia25h )

      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_25h_init : Output 25 hour mean diagnostics'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) 'Namelist nam_dia25h : set 25h outputs '
         WRITE(numout,*) 'Switch for 25h diagnostics (T) or not (F)  ln_dia25h  = ', ln_dia25h
      ENDIF
      IF( .NOT. ln_dia25h )   RETURN
      ! ------------------- !
      ! 1 - Allocate memory !
      ! ------------------- !
      ALLOCATE( tn_25h(jpi,jpj,jpk), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate tn_25h' )   ;   RETURN
      ENDIF
      ALLOCATE( sn_25h(jpi,jpj,jpk), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate sn_25h' )   ;   RETURN
      ENDIF
      ALLOCATE( un_25h(jpi,jpj,jpk), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate un_25h' )   ;   RETURN
      ENDIF
      ALLOCATE( vn_25h(jpi,jpj,jpk), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate vn_25h' )   ;   RETURN
      ENDIF
      ALLOCATE( wn_25h(jpi,jpj,jpk), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate wn_25h' )   ;   RETURN
      ENDIF
      ALLOCATE( avt_25h(jpi,jpj,jpk), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate avt_25h' )   ;   RETURN
      ENDIF
      ALLOCATE( avm_25h(jpi,jpj,jpk), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate avm_25h' )   ;   RETURN
      ENDIF
# if defined key_zdfgls || defined key_zdftke
      ALLOCATE( en_25h(jpi,jpj,jpk), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate en_25h' )   ;   RETURN
      ENDIF
#endif
# if defined key_zdfgls 
      ALLOCATE( rmxln_25h(jpi,jpj,jpk), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate rmxln_25h' )   ;   RETURN
      ENDIF
#endif
      ALLOCATE( sshn_25h(jpi,jpj), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate sshn_25h' )   ;   RETURN
      ENDIF
      ! ------------------------- !
      ! 2 - Assign Initial Values !
      ! ------------------------- !
      cnt_25h = 1  ! sets the first value of sum at timestep 1 (note - should strictly be at timestep zero so before values used where possible) 
      tn_25h(:,:,:) = tsb(:,:,:,jp_tem)
      sn_25h(:,:,:) = tsb(:,:,:,jp_sal)
      sshn_25h(:,:) = sshb(:,:)
      un_25h(:,:,:) = ub(:,:,:)
      vn_25h(:,:,:) = vb(:,:,:)
      wn_25h(:,:,:) = wn(:,:,:)
      avt_25h(:,:,:) = avt(:,:,:)
      avm_25h(:,:,:) = avm(:,:,:)
# if defined key_zdfgls || defined key_zdftke
         en_25h(:,:,:) = en(:,:,:)
#endif
# if defined key_zdfgls
         rmxln_25h(:,:,:) = mxln(:,:,:)
#endif
#if defined key_lim3 || defined key_lim2
         CALL ctl_stop('STOP', 'dia_25h not setup yet to do tidemean ice')
#endif 

      ! -------------------------- !
      ! 3 - Return to dia_wri      !
      ! -------------------------- !


   END SUBROUTINE dia_25h_init


   SUBROUTINE dia_25h( kt )  
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_25h  ***
      !!         
      !!
      !!--------------------------------------------------------------------
      !!                   
      !! ** Purpose :   Write diagnostics with M2/S2 tide removed
      !!
      !! ** Method  :   
      !!      25hr mean outputs for shelf seas
      !!
      !! History :
      !!   ?.0  !  07-04  (A. Hines) New routine, developed from dia_wri_foam
      !!   3.4  !  02-13  (J. Siddorn) Routine taken from old dia_wri_foam
      !!   3.6  !  08-14  (E. O'Dea) adapted for VN3.6
      !!----------------------------------------------------------------------
      !! * Modules used

      IMPLICIT NONE

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index


      !! * Local declarations
      INTEGER ::   ji, jj, jk

      LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
      REAL(wp)                         ::   zsto, zout, zmax, zjulian, zmdi       ! temporary reals
      INTEGER                          ::   i_steps                               ! no of timesteps per hour
      REAL(wp), DIMENSION(jpi,jpj    ) ::   zw2d, un_dm, vn_dm                    ! temporary workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zw3d                                  ! temporary workspace
      REAL(wp), DIMENSION(jpi,jpj,3)   ::   zwtmb                                 ! temporary workspace
      INTEGER                          ::   iyear0, nimonth0,iday0                ! start year,imonth,day

      !!----------------------------------------------------------------------

      ! 0. Initialisation
      ! -----------------
      ! Define frequency of summing to create 25 h mean
      IF( MOD( 3600,INT(rdt) ) == 0 ) THEN
         i_steps = 3600/INT(rdt)
      ELSE
         CALL ctl_stop('STOP', 'dia_wri_tide: timestep must give MOD(3600,rdt) = 0 otherwise no hourly values are possible')
      ENDIF

#if defined key_lim3 || defined key_lim2
      CALL ctl_stop('STOP', 'dia_wri_tide not setup yet to do tidemean ice')
#endif

      ! local variable for debugging
      ll_print = ll_print .AND. lwp

      ! Sum of 25 hourly instantaneous values to give a 25h mean from 24hours
      ! every day
      IF( MOD( kt, i_steps ) == 0  .and. kt .ne. nn_it000 ) THEN

         IF (lwp) THEN
              WRITE(numout,*) 'dia_wri_tide : Summing instantaneous hourly diagnostics at timestep ',kt
              WRITE(numout,*) '~~~~~~~~~~~~ '
         ENDIF

         tn_25h(:,:,:)        = tn_25h(:,:,:) + tsn(:,:,:,jp_tem)
         sn_25h(:,:,:)        = sn_25h(:,:,:) + tsn(:,:,:,jp_sal)
         sshn_25h(:,:)        = sshn_25h(:,:) + sshn (:,:)
         un_25h(:,:,:)        = un_25h(:,:,:) + un(:,:,:)
         vn_25h(:,:,:)        = vn_25h(:,:,:) + vn(:,:,:)
         wn_25h(:,:,:)        = wn_25h(:,:,:) + wn(:,:,:)
         avt_25h(:,:,:)       = avt_25h(:,:,:) + avt(:,:,:)
         avm_25h(:,:,:)       = avm_25h(:,:,:) + avm(:,:,:)
# if defined key_zdfgls || defined key_zdftke
         en_25h(:,:,:)        = en_25h(:,:,:) + en(:,:,:)
#endif
# if defined key_zdfgls
         rmxln_25h(:,:,:)      = rmxln_25h(:,:,:) + mxln(:,:,:)
#endif
         cnt_25h = cnt_25h + 1

         IF (lwp) THEN
            WRITE(numout,*) 'dia_tide : Summed the following number of hourly values so far',cnt_25h
         ENDIF

      ENDIF ! MOD( kt, i_steps ) == 0

         ! Write data for 25 hour mean output streams
      IF( cnt_25h .EQ. 25 .AND.  MOD( kt, i_steps*24) == 0 .AND. kt .NE. nn_it000 ) THEN

            IF(lwp) THEN
               WRITE(numout,*) 'dia_wri_tide : Writing 25 hour mean tide diagnostics at timestep', kt
               WRITE(numout,*) '~~~~~~~~~~~~ '
            ENDIF

            tn_25h(:,:,:)        = tn_25h(:,:,:) / 25.0_wp
            sn_25h(:,:,:)        = sn_25h(:,:,:) / 25.0_wp
            sshn_25h(:,:)        = sshn_25h(:,:) / 25.0_wp
            un_25h(:,:,:)        = un_25h(:,:,:) / 25.0_wp
            vn_25h(:,:,:)        = vn_25h(:,:,:) / 25.0_wp
            wn_25h(:,:,:)        = wn_25h(:,:,:) / 25.0_wp
            avt_25h(:,:,:)       = avt_25h(:,:,:) / 25.0_wp
            avm_25h(:,:,:)       = avm_25h(:,:,:) / 25.0_wp
# if defined key_zdfgls || defined key_zdftke
            en_25h(:,:,:)        = en_25h(:,:,:) / 25.0_wp
#endif
# if defined key_zdfgls
            rmxln_25h(:,:,:)       = rmxln_25h(:,:,:) / 25.0_wp
#endif

            IF (lwp)  WRITE(numout,*) 'dia_wri_tide : Mean calculated by dividing 25 hour sums and writing output'
            zmdi=1.e+20 !missing data indicator for masking
            ! write tracers (instantaneous)
            zw3d(:,:,:) = tn_25h(:,:,:)*tmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
            CALL iom_put("temper25h", zw3d)   ! potential temperature
            zw3d(:,:,:) = sn_25h(:,:,:)*tmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
            CALL iom_put( "salin25h", zw3d  )   ! salinity
            zw2d(:,:) = sshn_25h(:,:)*tmask(:,:,1) + zmdi*(1.0-tmask(:,:,1))
            CALL iom_put( "ssh25h", zw2d )   ! sea surface 


            ! Write velocities (instantaneous)
            zw3d(:,:,:) = un_25h(:,:,:)*umask(:,:,:) + zmdi*(1.0-umask(:,:,:))
            CALL iom_put("vozocrtx25h", zw3d)    ! i-current
            zw3d(:,:,:) = vn_25h(:,:,:)*vmask(:,:,:) + zmdi*(1.0-vmask(:,:,:))
            CALL iom_put("vomecrty25h", zw3d  )   ! j-current

            zw3d(:,:,:) = wn_25h(:,:,:)*tmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
            CALL iom_put("vomecrtz25h", zw3d )   ! k-current
            zw3d(:,:,:) = avt_25h(:,:,:)*tmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
            CALL iom_put("avt25h", zw3d )   ! diffusivity
            zw3d(:,:,:) = avm_25h(:,:,:)*tmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
            CALL iom_put("avm25h", zw3d)   ! viscosity
#if defined key_zdftke || defined key_zdfgls 
            zw3d(:,:,:) = en_25h(:,:,:)*tmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
            CALL iom_put("tke25h", zw3d)   ! tke
#endif
#if defined key_zdfgls 
            zw3d(:,:,:) = rmxln_25h(:,:,:)*tmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
            CALL iom_put( "mxln25h",zw3d)
#endif

            ! After the write reset the values to cnt=1 and sum values equal current value 
            tn_25h(:,:,:) = tsn(:,:,:,jp_tem)
            sn_25h(:,:,:) = tsn(:,:,:,jp_sal)
            sshn_25h(:,:) = sshn (:,:)
            un_25h(:,:,:) = un(:,:,:)
            vn_25h(:,:,:) = vn(:,:,:)
            wn_25h(:,:,:) = wn(:,:,:)
            avt_25h(:,:,:) = avt(:,:,:)
            avm_25h(:,:,:) = avm(:,:,:)
# if defined key_zdfgls || defined key_zdftke
            en_25h(:,:,:) = en(:,:,:)
#endif
# if defined key_zdfgls
            rmxln_25h(:,:,:) = mxln(:,:,:)
#endif
            cnt_25h = 1
            IF (lwp)  WRITE(numout,*) 'dia_wri_tide :   &
        &    After 25hr mean write, reset sum to current value and cnt_25h to one for overlapping average',cnt_25h

      ENDIF !  cnt_25h .EQ. 25 .AND.  MOD( kt, i_steps * 24) == 0 .AND. kt .NE. nn_it000


   END SUBROUTINE dia_25h 

   !!======================================================================
END MODULE dia25h
