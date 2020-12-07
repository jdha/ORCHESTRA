MODULE limistate
   !!======================================================================
   !!                     ***  MODULE  limistate  ***
   !!              Initialisation of diagnostics ice variables
   !!======================================================================
   !! History :  2.0  ! 2004-01 (C. Ethe, G. Madec)  Original code
   !!            3.0  ! 2011-02 (G. Madec) dynamical allocation
   !!             -   ! 2014    (C. Rousset) add N/S initializations
   !!----------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3' :                                    LIM3 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_istate      :  Initialisation of diagnostics ice variables
   !!   lim_istate_init :  initialization of ice state and namelist read
   !!----------------------------------------------------------------------
   USE phycst           ! physical constant
   USE oce              ! dynamics and tracers variables
   USE dom_oce          ! ocean domain
   USE sbc_oce          ! Surface boundary condition: ocean fields
   USE sbc_ice          ! Surface boundary condition: ice fields
   USE eosbn2           ! equation of state
   USE ice              ! sea-ice variables
   USE par_oce          ! ocean parameters
   USE limvar           ! lim_var_salprof
   !
   USE in_out_manager   ! I/O manager
   USE lib_mpp          ! MPP library
   USE lib_fortran      ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
   USE wrk_nemo         ! work arrays
   USE fldread          ! read input fields
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_istate      ! routine called by lim_init.F90

   INTEGER , PARAMETER ::   jpfldi = 6           ! maximum number of files to read
   INTEGER , PARAMETER ::   jp_hti = 1           ! index of ice thickness (m)    at T-point
   INTEGER , PARAMETER ::   jp_hts = 2           ! index of snow thicknes (m)    at T-point
   INTEGER , PARAMETER ::   jp_ati = 3           ! index of ice fraction (%) at T-point
   INTEGER , PARAMETER ::   jp_tsu = 4           ! index of ice surface temp (K)    at T-point
   INTEGER , PARAMETER ::   jp_tmi = 5           ! index of ice temp at T-point
   INTEGER , PARAMETER ::   jp_smi = 6           ! index of ice sali at T-point
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   si  ! structure of input fields (file informations, fields read)
   !!----------------------------------------------------------------------
   !!   LIM 3.0,  UCL-LOCEAN-IPSL (2008)
   !! $Id: limistate.F90 7761 2017-03-06 17:58:35Z clem $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_istate
      !!-------------------------------------------------------------------
      !!                    ***  ROUTINE lim_istate  ***
      !!
      !! ** Purpose :   defined the sea-ice initial state
      !!
      !! ** Method  :   This routine will put some ice where ocean
      !!                is at the freezing point, then fill in ice 
      !!                state variables using prescribed initial 
      !!                values in the namelist            
      !!
      !! ** Steps   :   1) Read namelist
      !!                2) Basal temperature; ice and hemisphere masks
      !!                3) Fill in the ice thickness distribution using gaussian
      !!                4) Fill in space-dependent arrays for state variables
      !!                5) Diagnostic arrays
      !!                6) Lateral boundary conditions
      !!
      !! ** Notes   : o_i, t_su, t_s, t_i, s_i must be filled everywhere, even
      !!              where there is no ice (clem: I do not know why, is it mandatory?) 
      !!
      !! History :
      !!   2.0  !  01-04  (C. Ethe, G. Madec)  Original code
      !!   3.0  !  2007   (M. Vancoppenolle)   Rewrite for ice cats
      !!   4.0  !  09-11  (M. Vancoppenolle)   Enhanced version for ice cats
      !!--------------------------------------------------------------------
      INTEGER  :: ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) :: ztmelts, zdh
      INTEGER  :: i_hemis, i_fill, jl0  
      REAL(wp)   :: zarg, zV, zconv, zdv 
      REAL(wp), POINTER, DIMENSION(:,:)   :: zswitch    ! ice indicator
      REAL(wp), POINTER, DIMENSION(:,:)   :: zht_i_ini, zat_i_ini, zvt_i_ini            !data from namelist or nc file
      REAL(wp), POINTER, DIMENSION(:,:)   :: zts_u_ini, zht_s_ini, zsm_i_ini, ztm_i_ini !data from namelist or nc file
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zh_i_ini, za_i_ini                         !data by cattegories to fill
      INTEGER , POINTER, DIMENSION(:)     :: itest
      !--------------------------------------------------------------------

      CALL wrk_alloc( jpi, jpj, jpl, zh_i_ini,  za_i_ini )
      CALL wrk_alloc( jpi, jpj,      zht_i_ini, zat_i_ini, zvt_i_ini, zts_u_ini, zht_s_ini, zsm_i_ini, ztm_i_ini )
      CALL wrk_alloc( jpi, jpj,      zswitch )
      Call wrk_alloc( 4,             itest )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'lim_istate : sea-ice initialization '
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~ '

      !--------------------------------------------------------------------
      ! 1) Read namelist
      !--------------------------------------------------------------------
      !
      CALL lim_istate_init

      ! init surface temperature
      DO jl = 1, jpl
         t_su  (:,:,jl) = rt0 * tmask(:,:,1)
         tn_ice(:,:,jl) = rt0 * tmask(:,:,1)
      END DO

      ! init basal temperature (considered at freezing point)
      CALL eos_fzp( sss_m(:,:), t_bo(:,:) )
      t_bo(:,:) = ( t_bo(:,:) + rt0 ) * tmask(:,:,1) 


      !--------------------------------------------------------------------
      ! 2) Initialization of sea ice state variables
      !--------------------------------------------------------------------
      IF( ln_limini ) THEN
         !
         IF( ln_limini_file )THEN
         !
            zht_i_ini(:,:)  = si(jp_hti)%fnow(:,:,1)
            zht_s_ini(:,:)  = si(jp_hts)%fnow(:,:,1)
            zat_i_ini(:,:)  = si(jp_ati)%fnow(:,:,1)
            zts_u_ini(:,:)  = si(jp_tsu)%fnow(:,:,1)
            ztm_i_ini(:,:)  = si(jp_tmi)%fnow(:,:,1)
            zsm_i_ini(:,:)  = si(jp_smi)%fnow(:,:,1)
            !
            WHERE( zat_i_ini(:,:) > 0._wp ) ; zswitch(:,:) = tmask(:,:,1) 
            ELSEWHERE                       ; zswitch(:,:) = 0._wp
            END WHERE
            !
         ELSE ! ln_limini_file = F

            !--------------------------------------------------------------------
            ! 3) Basal temperature, ice mask
            !--------------------------------------------------------------------
            ! no ice if sst <= t-freez + ttest
            WHERE( ( sst_m(:,:) - (t_bo(:,:) - rt0) ) * tmask(:,:,1) >= rn_thres_sst ) ; zswitch(:,:) = 0._wp 
            ELSEWHERE                                                                  ; zswitch(:,:) = tmask(:,:,1)
            END WHERE

            !-----------------------------
            ! 3.1) Hemisphere-dependent arrays
            !-----------------------------
            ! assign initial thickness, concentration, snow depth and salinity to an hemisphere-dependent array
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( ff_t(ji,jj) >= 0._wp ) THEN
                     zht_i_ini(ji,jj) = rn_hti_ini_n * zswitch(ji,jj)
                     zht_s_ini(ji,jj) = rn_hts_ini_n * zswitch(ji,jj)
                     zat_i_ini(ji,jj) = rn_ati_ini_n * zswitch(ji,jj)
                     zts_u_ini(ji,jj) = rn_tmi_ini_n * zswitch(ji,jj)
                     zsm_i_ini(ji,jj) = rn_smi_ini_n * zswitch(ji,jj)
                     ztm_i_ini(ji,jj) = rn_tmi_ini_n * zswitch(ji,jj)
                  ELSE
                     zht_i_ini(ji,jj) = rn_hti_ini_s * zswitch(ji,jj)
                     zht_s_ini(ji,jj) = rn_hts_ini_s * zswitch(ji,jj)
                     zat_i_ini(ji,jj) = rn_ati_ini_s * zswitch(ji,jj)
                     zts_u_ini(ji,jj) = rn_tmi_ini_s * zswitch(ji,jj)
                     zsm_i_ini(ji,jj) = rn_smi_ini_s * zswitch(ji,jj)
                     ztm_i_ini(ji,jj) = rn_tmi_ini_s * zswitch(ji,jj)
                  ENDIF
               END DO
            END DO
            !
         ENDIF ! ln_limini_file
         
         zvt_i_ini(:,:) = zht_i_ini(:,:) * zat_i_ini(:,:)   ! ice volume
         !---------------------------------------------------------------------
         ! 3.2) Distribute ice concentration and thickness into the categories
         !---------------------------------------------------------------------
         ! a gaussian distribution for ice concentration is used
         ! then we check whether the distribution fullfills
         ! volume and area conservation, positivity and ice categories bounds
         zh_i_ini(:,:,:) = 0._wp 
         za_i_ini(:,:,:) = 0._wp
         !
         DO jj = 1, jpj
            DO ji = 1, jpi
               !
               IF( zat_i_ini(ji,jj) > 0._wp .AND. zht_i_ini(ji,jj) > 0._wp )THEN

                  !--- jl0: most likely index where cc will be maximum
                  jl0 = jpl
                  DO jl = 1, jpl
                     IF ( ( zht_i_ini(ji,jj) >  hi_max(jl-1) ) .AND. ( zht_i_ini(ji,jj) <= hi_max(jl) ) ) THEN
                        jl0 = jl
                        CYCLE
                     ENDIF
                  END DO
                  !
                  ! initialisation of tests
                  itest(:)  = 0
                  
                  i_fill = jpl + 1                                             !====================================
                  DO WHILE ( ( SUM( itest(:) ) /= 4 ) .AND. ( i_fill >= 2 ) )  ! iterative loop on i_fill categories
                     ! iteration                                               !====================================
                     i_fill = i_fill - 1

                     ! initialisation of ice variables for each try
                     zh_i_ini(ji,jj,:) = 0._wp 
                     za_i_ini(ji,jj,:) = 0._wp
                     itest(:) = 0
                     !
                     ! *** case very thin ice: fill only category 1
                     IF ( i_fill == 1 ) THEN
                        zh_i_ini(ji,jj,1) = zht_i_ini(ji,jj)
                        za_i_ini(ji,jj,1) = zat_i_ini(ji,jj)

                     ! *** case ice is thicker: fill categories >1
                     ELSE

                        ! Fill ice thicknesses in the (i_fill-1) cat by hmean 
                        DO jl = 1, i_fill-1
                           zh_i_ini(ji,jj,jl) = hi_mean(jl)
                        END DO
                        !
                        !--- Concentrations
                        za_i_ini(ji,jj,jl0) = zat_i_ini(ji,jj) / SQRT(REAL(jpl))
                        DO jl = 1, i_fill - 1
                           IF( jl /= jl0 )THEN
                              zarg               = ( zh_i_ini(ji,jj,jl) - zht_i_ini(ji,jj) ) / ( 0.5_wp * zht_i_ini(ji,jj) )
                              za_i_ini(ji,jj,jl) = za_i_ini(ji,jj,jl0) * EXP(-zarg**2)
                           ENDIF
                        END DO
                        !
                        ! Concentration in the last (i_fill) category
                        za_i_ini(ji,jj,i_fill) = zat_i_ini(ji,jj) - SUM( za_i_ini(ji,jj,1:i_fill-1) )

                        ! Ice thickness in the last (i_fill) category
                        zV = SUM( za_i_ini(ji,jj,1:i_fill-1) * zh_i_ini(ji,jj,1:i_fill-1) )
                        zh_i_ini(ji,jj,i_fill) = ( zvt_i_ini(ji,jj) - zV ) / MAX( za_i_ini(ji,jj,i_fill), epsi10 ) 

                        ! clem: correction if concentration of upper cat is greater than lower cat
                        !       (it should be a gaussian around jl0 but sometimes it is not)
                        IF ( jl0 /= jpl ) THEN
                           DO jl = jpl, jl0+1, -1
                              IF ( za_i_ini(ji,jj,jl) > za_i_ini(ji,jj,jl-1) ) THEN
                                 zdv = zh_i_ini(ji,jj,jl) * za_i_ini(ji,jj,jl)
                                 zh_i_ini(ji,jj,jl    ) = 0._wp
                                 za_i_ini(ji,jj,jl    ) = 0._wp
                                 za_i_ini(ji,jj,1:jl-1) = za_i_ini(ji,jj,1:jl-1)  &
                                    &                     + zdv / MAX( REAL(jl-1) * zht_i_ini(ji,jj), epsi10 )
                              END IF
                           ENDDO
                        ENDIF
                        !
                     ENDIF ! case ice is thick or thin

                     !---------------------
                     ! Compatibility tests
                     !---------------------
                     ! Test 1: area conservation
                     zconv = ABS( zat_i_ini(ji,jj) - SUM( za_i_ini(ji,jj,1:jpl) ) )
                     IF ( zconv < epsi06 ) itest(1) = 1
                     
                     ! Test 2: volume conservation
                     zconv = ABS(       zat_i_ini(ji,jj)       * zht_i_ini(ji,jj)   &
                        &        - SUM( za_i_ini (ji,jj,1:jpl) * zh_i_ini (ji,jj,1:jpl) ) )
                     IF ( zconv < epsi06 ) itest(2) = 1
                     
                     ! Test 3: thickness of the last category is in-bounds ?
                     IF ( zh_i_ini(ji,jj,i_fill) >= hi_max(i_fill-1) ) itest(3) = 1
                     
                     ! Test 4: positivity of ice concentrations
                     itest(4) = 1
                     DO jl = 1, i_fill
                        IF ( za_i_ini(ji,jj,jl) < 0._wp ) itest(4) = 0
                     END DO
                     !                                      !============================
                  END DO                                    ! end iteration on categories
                  !                                         !============================
                  !
                  IF( lwp .AND. SUM(itest) /= 4 ) THEN 
                     WRITE(numout,*)
                     WRITE(numout,*) ' !!!! ALERT itest is not equal to 4      !!! '
                     WRITE(numout,*) ' !!!! Something is wrong in the LIM3 initialization procedure '
                     WRITE(numout,*)
                     WRITE(numout,*) ' *** itest_i (i=1,4) = ', itest(:)
                     WRITE(numout,*) ' zat_i_ini : ', zat_i_ini(ji,jj)
                     WRITE(numout,*) ' zht_i_ini : ', zht_i_ini(ji,jj)
                  ENDIF
               
               ENDIF !  zat_i_ini(ji,jj) > 0._wp .AND. zht_i_ini(ji,jj) > 0._wp
               !
            END DO   
         END DO   

         !---------------------------------------------------------------------
         ! 3.3) Space-dependent arrays for ice state variables
         !---------------------------------------------------------------------

         ! Ice concentration, thickness and volume, ice salinity, ice age, surface temperature
         DO jl = 1, jpl ! loop over categories
            DO jj = 1, jpj
               DO ji = 1, jpi
                  a_i(ji,jj,jl)   = zswitch(ji,jj) * za_i_ini(ji,jj,jl)                       ! concentration
                  ht_i(ji,jj,jl)  = zswitch(ji,jj) * zh_i_ini(ji,jj,jl)                       ! ice thickness
                  sm_i(ji,jj,jl)  = zswitch(ji,jj) * zsm_i_ini(ji,jj)                         ! salinity
                  o_i(ji,jj,jl)   = zswitch(ji,jj) * 1._wp                                    ! age (1 day)
                  t_su(ji,jj,jl)  = zswitch(ji,jj) * zts_u_ini(ji,jj) + ( 1._wp - zswitch(ji,jj) ) * rt0 ! surf temp

                  IF( zht_i_ini(ji,jj) > 0._wp )THEN
                    ht_s(ji,jj,jl)= ht_i(ji,jj,jl) * ( zht_s_ini(ji,jj) / zht_i_ini(ji,jj) )  ! snow depth
                  ELSE
                    ht_s(ji,jj,jl)= 0._wp
                  ENDIF

                  ! This case below should not be used if (ht_s/ht_i) is ok in namelist
                  ! In case snow load is in excess that would lead to transformation from snow to ice
                  ! Then, transfer the snow excess into the ice (different from limthd_dh)
                  zdh = MAX( 0._wp, ( rhosn * ht_s(ji,jj,jl) + ( rhoic - rau0 ) * ht_i(ji,jj,jl) ) * r1_rau0 ) 
                  ! recompute ht_i, ht_s avoiding out of bounds values
                  ht_i(ji,jj,jl) = MIN( hi_max(jl), ht_i(ji,jj,jl) + zdh )
                  ht_s(ji,jj,jl) = MAX( 0._wp, ht_s(ji,jj,jl) - zdh * rhoic * r1_rhosn )

                  ! ice volume, salt content, age content
                  v_i(ji,jj,jl)   = ht_i(ji,jj,jl) * a_i(ji,jj,jl)              ! ice volume
                  v_s(ji,jj,jl)   = ht_s(ji,jj,jl) * a_i(ji,jj,jl)              ! snow volume
                  smv_i(ji,jj,jl) = MIN( sm_i(ji,jj,jl) , sss_m(ji,jj) ) * v_i(ji,jj,jl) ! salt content
                  oa_i(ji,jj,jl)  = o_i(ji,jj,jl) * a_i(ji,jj,jl)               ! age content
               END DO
            END DO
         END DO

         ! for constant salinity in time
         IF( nn_icesal == 1 .OR. nn_icesal == 3 )  THEN
            CALL lim_var_salprof
            smv_i = sm_i * v_i
         ENDIF
            
         ! Snow temperature and heat content
         DO jk = 1, nlay_s
            DO jl = 1, jpl ! loop over categories
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     t_s(ji,jj,jk,jl) = zswitch(ji,jj) * ztm_i_ini(ji,jj) + ( 1._wp - zswitch(ji,jj) ) * rt0
                     ! Snow energy of melting
                     e_s(ji,jj,jk,jl) = zswitch(ji,jj) * rhosn * ( cpic * ( rt0 - t_s(ji,jj,jk,jl) ) + lfus )

                     ! Mutliply by volume, and divide by number of layers to get heat content in J/m2
                     e_s(ji,jj,jk,jl) = e_s(ji,jj,jk,jl) * v_s(ji,jj,jl) * r1_nlay_s
                  END DO
               END DO
            END DO
         END DO

         ! Ice salinity, temperature and heat content
         DO jk = 1, nlay_i
            DO jl = 1, jpl ! loop over categories
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     t_i(ji,jj,jk,jl) = zswitch(ji,jj) * ztm_i_ini(ji,jj) + ( 1._wp - zswitch(ji,jj) ) * rt0 
                     s_i(ji,jj,jk,jl) = zswitch(ji,jj) * zsm_i_ini(ji,jj) + ( 1._wp - zswitch(ji,jj) ) * rn_simin
                     ztmelts          = - tmut * s_i(ji,jj,jk,jl) + rt0 !Melting temperature in K

                     ! heat content per unit volume
                     e_i(ji,jj,jk,jl) = zswitch(ji,jj) * rhoic * (   cpic    * ( ztmelts - t_i(ji,jj,jk,jl) ) &
                        +   lfus    * ( 1._wp - (ztmelts-rt0) / MIN((t_i(ji,jj,jk,jl)-rt0),-epsi20) ) &
                        -   rcp     * ( ztmelts - rt0 ) )

                     ! Mutliply by ice volume, and divide by number of layers to get heat content in J/m2
                     e_i(ji,jj,jk,jl) = e_i(ji,jj,jk,jl) * v_i(ji,jj,jl) * r1_nlay_i
                  END DO
               END DO
            END DO
         END DO

         tn_ice (:,:,:) = t_su (:,:,:)

      ELSE ! if ln_limini=false
         a_i  (:,:,:) = 0._wp
         v_i  (:,:,:) = 0._wp
         v_s  (:,:,:) = 0._wp
         smv_i(:,:,:) = 0._wp
         oa_i (:,:,:) = 0._wp
         ht_i (:,:,:) = 0._wp
         ht_s (:,:,:) = 0._wp
         sm_i (:,:,:) = 0._wp
         o_i  (:,:,:) = 0._wp

         e_i(:,:,:,:) = 0._wp
         e_s(:,:,:,:) = 0._wp

         DO jl = 1, jpl
            DO jk = 1, nlay_i
               t_i(:,:,jk,jl) = rt0 * tmask(:,:,1)
            END DO
            DO jk = 1, nlay_s
               t_s(:,:,jk,jl) = rt0 * tmask(:,:,1)
            END DO
         END DO

      ENDIF ! ln_limini
      
      at_i (:,:) = 0.0_wp
      DO jl = 1, jpl
         at_i (:,:) = at_i (:,:) + a_i (:,:,jl)
      END DO
      !
      !--------------------------------------------------------------------
      ! 4) Global ice variables for output diagnostics                    | 
      !--------------------------------------------------------------------
      u_ice (:,:)     = 0._wp
      v_ice (:,:)     = 0._wp
      stress1_i(:,:)  = 0._wp
      stress2_i(:,:)  = 0._wp
      stress12_i(:,:) = 0._wp

      !--------------------------------------------------------------------
      ! 5) Moments for advection
      !--------------------------------------------------------------------

      sxopw (:,:) = 0._wp 
      syopw (:,:) = 0._wp 
      sxxopw(:,:) = 0._wp 
      syyopw(:,:) = 0._wp 
      sxyopw(:,:) = 0._wp

      sxice (:,:,:)  = 0._wp   ;   sxsn (:,:,:)  = 0._wp   ;   sxa  (:,:,:)  = 0._wp
      syice (:,:,:)  = 0._wp   ;   sysn (:,:,:)  = 0._wp   ;   sya  (:,:,:)  = 0._wp
      sxxice(:,:,:)  = 0._wp   ;   sxxsn(:,:,:)  = 0._wp   ;   sxxa (:,:,:)  = 0._wp
      syyice(:,:,:)  = 0._wp   ;   syysn(:,:,:)  = 0._wp   ;   syya (:,:,:)  = 0._wp
      sxyice(:,:,:)  = 0._wp   ;   sxysn(:,:,:)  = 0._wp   ;   sxya (:,:,:)  = 0._wp

      sxc0  (:,:,:)  = 0._wp   ;   sxe  (:,:,:,:)= 0._wp   
      syc0  (:,:,:)  = 0._wp   ;   sye  (:,:,:,:)= 0._wp   
      sxxc0 (:,:,:)  = 0._wp   ;   sxxe (:,:,:,:)= 0._wp   
      syyc0 (:,:,:)  = 0._wp   ;   syye (:,:,:,:)= 0._wp   
      sxyc0 (:,:,:)  = 0._wp   ;   sxye (:,:,:,:)= 0._wp   

      sxsal  (:,:,:)  = 0._wp
      sysal  (:,:,:)  = 0._wp
      sxxsal (:,:,:)  = 0._wp
      syysal (:,:,:)  = 0._wp
      sxysal (:,:,:)  = 0._wp

      sxage  (:,:,:)  = 0._wp
      syage  (:,:,:)  = 0._wp
      sxxage (:,:,:)  = 0._wp
      syyage (:,:,:)  = 0._wp
      sxyage (:,:,:)  = 0._wp

      !------------------------------------
      ! 6) store fields at before time-step
      !------------------------------------
      ! it is only necessary for the 1st interpolation by Agrif
      a_i_b  (:,:,:)   = a_i  (:,:,:)
      e_i_b  (:,:,:,:) = e_i  (:,:,:,:)
      v_i_b  (:,:,:)   = v_i  (:,:,:)
      v_s_b  (:,:,:)   = v_s  (:,:,:)
      e_s_b  (:,:,:,:) = e_s  (:,:,:,:)
      smv_i_b(:,:,:)   = smv_i(:,:,:)
      oa_i_b (:,:,:)   = oa_i (:,:,:)
      u_ice_b(:,:)     = u_ice(:,:)
      v_ice_b(:,:)     = v_ice(:,:)

!!!clem
!!      ! Output the initial state and forcings
!!      CALL dia_wri_state( 'output.init', nit000 )
!!!      

      CALL wrk_dealloc( jpi, jpj, jpl, zh_i_ini,  za_i_ini )
      CALL wrk_dealloc( jpi, jpj,      zht_i_ini, zat_i_ini, zvt_i_ini, zts_u_ini, zht_s_ini, zsm_i_ini, ztm_i_ini )
      CALL wrk_dealloc( jpi, jpj,      zswitch )
      Call wrk_dealloc( 4,             itest )

   END SUBROUTINE lim_istate

   SUBROUTINE lim_istate_init
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE lim_istate_init  ***
      !!        
      !! ** Purpose : Definition of initial state of the ice 
      !!
      !! ** Method : Read the namiceini namelist and check the parameter 
      !!       values called at the first timestep (nit000)
      !!
      !! ** input : 
      !!        Namelist namiceini
      !!
      !! history :
      !!  8.5  ! 03-08 (C. Ethe) original code 
      !!  8.5  ! 07-11 (M. Vancoppenolle) rewritten initialization
      !!-----------------------------------------------------------------------------
      !
      INTEGER :: ios,ierr,inum_ice                 ! Local integer output status for namelist read
      INTEGER :: ji,jj
      INTEGER :: ifpr, ierror
      !
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ice files
      TYPE(FLD_N)                    ::   sn_hti, sn_hts, sn_ati, sn_tsu, sn_tmi, sn_smi
      TYPE(FLD_N), DIMENSION(jpfldi) ::   slf_i                 ! array of namelist informations on the fields to read
      !
      NAMELIST/namiceini/ ln_limini, ln_limini_file, rn_thres_sst, rn_hts_ini_n, rn_hts_ini_s,  &
         &                rn_hti_ini_n, rn_hti_ini_s, rn_ati_ini_n, rn_ati_ini_s, rn_smi_ini_n, &
         &                rn_smi_ini_s, rn_tmi_ini_n, rn_tmi_ini_s,                             &
         &                sn_hti, sn_hts, sn_ati, sn_tsu, sn_tmi, sn_smi, cn_dir
      !!-----------------------------------------------------------------------------
      !
      REWIND( numnam_ice_ref )              ! Namelist namiceini in reference namelist : Ice initial state
      READ  ( numnam_ice_ref, namiceini, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namiceini in reference namelist', lwp )

      REWIND( numnam_ice_cfg )              ! Namelist namiceini in configuration namelist : Ice initial state
      READ  ( numnam_ice_cfg, namiceini, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namiceini in configuration namelist', lwp )
      IF(lwm) WRITE ( numoni, namiceini )

      slf_i(jp_hti) = sn_hti  ;  slf_i(jp_hts) = sn_hts
      slf_i(jp_ati) = sn_ati  ;  slf_i(jp_tsu) = sn_tsu
      slf_i(jp_tmi) = sn_tmi  ;  slf_i(jp_smi) = sn_smi

      ! Define the initial parameters
      ! -------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lim_istate_init : ice parameters inititialisation '
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   initialization with ice (T) or not (F)       ln_limini     = ', ln_limini
         WRITE(numout,*) '   ice initialization from a netcdf file      ln_limini_file  = ', ln_limini_file
         WRITE(numout,*) '   threshold water temp. for initial sea-ice    rn_thres_sst  = ', rn_thres_sst
         WRITE(numout,*) '   initial snow thickness in the north          rn_hts_ini_n  = ', rn_hts_ini_n
         WRITE(numout,*) '   initial snow thickness in the south          rn_hts_ini_s  = ', rn_hts_ini_s 
         WRITE(numout,*) '   initial ice thickness  in the north          rn_hti_ini_n  = ', rn_hti_ini_n
         WRITE(numout,*) '   initial ice thickness  in the south          rn_hti_ini_s  = ', rn_hti_ini_s
         WRITE(numout,*) '   initial ice concentr.  in the north          rn_ati_ini_n  = ', rn_ati_ini_n
         WRITE(numout,*) '   initial ice concentr.  in the north          rn_ati_ini_s  = ', rn_ati_ini_s
         WRITE(numout,*) '   initial  ice salinity  in the north          rn_smi_ini_n  = ', rn_smi_ini_n
         WRITE(numout,*) '   initial  ice salinity  in the south          rn_smi_ini_s  = ', rn_smi_ini_s
         WRITE(numout,*) '   initial  ice/snw temp  in the north          rn_tmi_ini_n  = ', rn_tmi_ini_n
         WRITE(numout,*) '   initial  ice/snw temp  in the south          rn_tmi_ini_s  = ', rn_tmi_ini_s
      ENDIF

      IF( ln_limini_file ) THEN                      ! Ice initialization using input file
         !
         ! set si structure
         ALLOCATE( si(jpfldi), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'Ice_ini in limistate: unable to allocate si structure' )   ;   RETURN
         ENDIF

         DO ifpr = 1, jpfldi
            ALLOCATE( si(ifpr)%fnow(jpi,jpj,1) )
            ALLOCATE( si(ifpr)%fdta(jpi,jpj,1,2) )
         END DO

         ! fill si with slf_i and control print
         CALL fld_fill( si, slf_i, cn_dir, 'lim_istate', 'lim istate ini', 'numnam_ice' )

         CALL fld_read( nit000, 1, si )                ! input fields provided at the current time-step

      ENDIF

   END SUBROUTINE lim_istate_init

#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module          NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_istate          ! Empty routine
   END SUBROUTINE lim_istate
#endif

   !!======================================================================
END MODULE limistate
