MODULE dynspg_ts
   !!======================================================================
   !!                   ***  MODULE  dynspg_ts  ***
   !! Ocean dynamics:  surface pressure gradient trend, split-explicit scheme
   !!======================================================================
   !! History :   1.0  ! 2004-12  (L. Bessieres, G. Madec)  Original code
   !!              -   ! 2005-11  (V. Garnier, G. Madec)  optimization
   !!              -   ! 2006-08  (S. Masson)  distributed restart using iom
   !!             2.0  ! 2007-07  (D. Storkey) calls to BDY routines
   !!              -   ! 2008-01  (R. Benshila)  change averaging method
   !!             3.2  ! 2009-07  (R. Benshila, G. Madec) Complete revisit associated to vvl reactivation
   !!             3.3  ! 2010-09  (D. Storkey, E. O'Dea) update for BDY for Shelf configurations
   !!             3.3  ! 2011-03  (R. Benshila, R. Hordoir, P. Oddo) update calculation of ub_b
   !!             3.5  ! 2013-07  (J. Chanut) Switch to Forward-backward time stepping
   !!             3.6  ! 2013-11  (A. Coward) Update for z-tilde compatibility
   !!             3.7  ! 2015-11  (J. Chanut) free surface simplification
   !!              -   ! 2016-12  (G. Madec, E. Clementi) update for Stoke-Drift divergence
   !!---------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_spg_ts     : compute surface pressure gradient trend using a time-splitting scheme 
   !!   dyn_spg_ts_init: initialisation of the time-splitting scheme
   !!   ts_wgt         : set time-splitting weights for temporal averaging (or not)
   !!   ts_rst         : read/write time-splitting fields in restart file
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE zdf_oce         ! Bottom friction coefts
   USE sbcisf          ! ice shelf variable (fwfisf)
   USE sbcapr          ! surface boundary condition: atmospheric pressure
   USE dynadv    , ONLY: ln_dynadv_vec
   USE phycst          ! physical constants
   USE dynvor          ! vorticity term
   USE wet_dry         ! wetting/drying flux limter
   USE bdy_oce         ! open boundary
   USE bdytides        ! open boundary condition data
   USE bdydyn2d        ! open boundary conditions on barotropic variables
   USE sbctide         ! tides
   USE updtide         ! tide potential
   USE sbcwave         ! surface wave
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE iom             ! IOM library
   USE restart         ! only for lrst_oce
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing    
   USE diatmb          ! Top,middle,bottom output
#if defined key_agrif
   USE agrif_opa_interp ! agrif
#endif
#if defined key_asminc   
   USE asminc          ! Assimilation increment
#endif


   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_spg_ts        ! routine called in dynspg.F90 
   PUBLIC dyn_spg_ts_alloc  !    "      "     "    "
   PUBLIC dyn_spg_ts_init   !    "      "     "    "
   PUBLIC ts_rst            !    "      "     "    "

   INTEGER, SAVE :: icycle  ! Number of barotropic sub-steps for each internal step nn_baro <= 2.5 nn_baro
   REAL(wp),SAVE :: rdtbt   ! Barotropic time step

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::   wgtbtp1, wgtbtp2   !: 1st & 2nd weights used in time filtering of barotropic fields

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  zwz          !: ff_f/h at F points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ftnw, ftne   !: triad of coriolis parameter
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ftsw, ftse   !: (only used with een vorticity scheme)

   !! Time filtered arrays at baroclinic time step:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   un_adv , vn_adv     !: Advection vel. at "now" barocl. step

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.5 , NEMO Consortium (2013)
   !! $Id: dynspg_ts.F90 7831 2017-03-24 10:44:55Z jamesharle $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dyn_spg_ts_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_ts_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(3)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( wgtbtp1(3*nn_baro), wgtbtp2(3*nn_baro), zwz(jpi,jpj), STAT=ierr(1) )
      !
      IF( ln_dynvor_een )   ALLOCATE( ftnw(jpi,jpj) , ftne(jpi,jpj) , & 
         &                            ftsw(jpi,jpj) , ftse(jpi,jpj) , STAT=ierr(2) )
         !
      ALLOCATE( un_adv(jpi,jpj), vn_adv(jpi,jpj)                    , STAT=ierr(3) )
      !
      dyn_spg_ts_alloc = MAXVAL( ierr(:) )
      !
      IF( lk_mpp                )   CALL mpp_sum( dyn_spg_ts_alloc )
      IF( dyn_spg_ts_alloc /= 0 )   CALL ctl_warn('dyn_spg_ts_alloc: failed to allocate arrays')
      !
   END FUNCTION dyn_spg_ts_alloc


   SUBROUTINE dyn_spg_ts( kt )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose : - Compute the now trend due to the explicit time stepping
      !!              of the quasi-linear barotropic system, and add it to the
      !!              general momentum trend. 
      !!
      !! ** Method  : - split-explicit schem (time splitting) :
      !!      Barotropic variables are advanced from internal time steps
      !!      "n"   to "n+1" if ln_bt_fw=T
      !!      or from 
      !!      "n-1" to "n+1" if ln_bt_fw=F
      !!      thanks to a generalized forward-backward time stepping (see ref. below).
      !!
      !! ** Action :
      !!      -Update the filtered free surface at step "n+1"      : ssha
      !!      -Update filtered barotropic velocities at step "n+1" : ua_b, va_b
      !!      -Compute barotropic advective velocities at step "n" : un_adv, vn_adv
      !!      These are used to advect tracers and are compliant with discrete
      !!      continuity equation taken at the baroclinic time steps. This 
      !!      ensures tracers conservation.
      !!      - (ua, va) momentum trend updated with barotropic component.
      !!
      !! References : Shchepetkin and McWilliams, Ocean Modelling, 2005. 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  ::   kt   ! ocean time-step index
      !
      LOGICAL  ::   ll_fw_start			! if true, forward integration 
      LOGICAL  ::   ll_init			    ! if true, special startup of 2d equations
      LOGICAL  ::   ll_tmp1, ll_tmp2            ! local logical variables used in W/D
      INTEGER  ::   ji, jj, jk, jn   		! dummy loop indices
      INTEGER  ::   ikbu, ikbv, noffset    	! local integers
      INTEGER  ::   iktu, iktv             	! local integers
      REAL(wp) ::   zmdi
      REAL(wp) ::   zraur, z1_2dt_b, z2dt_bf  	! local scalars
      REAL(wp) ::   zx1, zy1, zx2, zy2          !   -      -
      REAL(wp) ::   z1_12, z1_8, z1_4, z1_2 	!   -      -
      REAL(wp) ::   zu_spg, zv_spg              !   -      -
      REAL(wp) ::   zhura, zhvra         	!   -      -
      REAL(wp) ::   za0, za1, za2, za3		!   -      -
      !
      REAL(wp), POINTER, DIMENSION(:,:) :: zsshp2_e
      REAL(wp), POINTER, DIMENSION(:,:) :: zu_trd, zv_trd, zu_frc, zv_frc, zssh_frc
      REAL(wp), POINTER, DIMENSION(:,:) :: zwx, zwy, zhdiv
      REAL(wp), POINTER, DIMENSION(:,:) :: zhup2_e, zhvp2_e, zhust_e, zhvst_e
      REAL(wp), POINTER, DIMENSION(:,:) :: zsshu_a, zsshv_a
      REAL(wp), POINTER, DIMENSION(:,:) :: zhf
      REAL(wp), POINTER, DIMENSION(:,:) :: zcpx, zcpy                 ! Wetting/Dying gravity filter coef.
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_ts')
      !
      !                                         !* Allocate temporary arrays
      CALL wrk_alloc( jpi,jpj,   zsshp2_e, zhdiv )
      CALL wrk_alloc( jpi,jpj,   zu_trd, zv_trd)
      CALL wrk_alloc( jpi,jpj,   zwx, zwy, zssh_frc, zu_frc, zv_frc)
      CALL wrk_alloc( jpi,jpj,   zhup2_e, zhvp2_e, zhust_e, zhvst_e)
      CALL wrk_alloc( jpi,jpj,   zsshu_a, zsshv_a                  )
      CALL wrk_alloc( jpi,jpj,   zhf )
      IF( ln_wd ) CALL wrk_alloc( jpi, jpj, zcpx, zcpy )
      !
      zmdi=1.e+20                               !  missing data indicator for masking
      !                                         !* Local constant initialization
      z1_12 = 1._wp / 12._wp 
      z1_8  = 0.125_wp                                   
      z1_4  = 0.25_wp
      z1_2  = 0.5_wp     
      zraur = 1._wp / rau0
      !  	                                       ! reciprocal of baroclinic time step 
      IF( kt == nit000 .AND. neuler == 0 ) THEN   ;   z2dt_bf =          rdt
      ELSE                                        ;   z2dt_bf = 2.0_wp * rdt
      ENDIF
      z1_2dt_b = 1.0_wp / z2dt_bf 
      !
      ll_init     = ln_bt_av                       ! if no time averaging, then no specific restart 
      ll_fw_start = .FALSE.
      !                              	            ! time offset in steps for bdy data update
      IF( .NOT.ln_bt_fw ) THEN   ;   noffset = - nn_baro
      ELSE                       ;   noffset =   0 
      ENDIF
      !
      IF( kt == nit000 ) THEN             	!* initialisation
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_ts : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~   free surface with time splitting'
         IF(lwp) WRITE(numout,*)
         !
         IF( neuler == 0 )   ll_init=.TRUE.
         !
         IF( ln_bt_fw .OR. neuler == 0 ) THEN
            ll_fw_start =.TRUE.
            noffset     = 0
         ELSE
            ll_fw_start =.FALSE.
         ENDIF
         !
         ! Set averaging weights and cycle length:
         CALL ts_wgt( ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2 )
         !
      ENDIF
      !
      ! Set arrays to remove/compute coriolis trend.
      ! Do it once at kt=nit000 if volume is fixed, else at each long time step.
      ! Note that these arrays are also used during barotropic loop. These are however frozen
      ! although they should be updated in the variable volume case. Not a big approximation.
      ! To remove this approximation, copy lines below inside barotropic loop
      ! and update depths at T-F points (ht and zhf resp.) at each barotropic time step
      !
      IF( kt == nit000 .OR. .NOT.ln_linssh ) THEN
         IF( ln_dynvor_een ) THEN               !==  EEN scheme  ==!
            SELECT CASE( nn_een_e3f )              !* ff_f/e3 at F-point
            CASE ( 0 )                                   ! original formulation  (masked averaging of e3t divided by 4)
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     zwz(ji,jj) =   ( ht_n(ji  ,jj+1) + ht_n(ji+1,jj+1) +                    &
                        &             ht_n(ji  ,jj  ) + ht_n(ji+1,jj  )   ) * 0.25_wp  
                     IF( zwz(ji,jj) /= 0._wp )   zwz(ji,jj) = ff_f(ji,jj) / zwz(ji,jj)
                  END DO
               END DO
            CASE ( 1 )                                   ! new formulation  (masked averaging of e3t divided by the sum of mask)
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     zwz(ji,jj) =   ( ht_n(ji  ,jj+1) + ht_n(ji+1,jj+1) +                     &
                        &             ht_n(ji  ,jj  ) + ht_n(ji+1,jj  )   )                   &
                        &       / ( MAX( 1._wp, tmask(ji  ,jj+1, 1) + tmask(ji+1,jj+1, 1) +    &
                        &                       tmask(ji  ,jj  , 1) + tmask(ji+1,jj  , 1) ) )
                     IF( zwz(ji,jj) /= 0._wp )   zwz(ji,jj) = ff_f(ji,jj) / zwz(ji,jj)
                  END DO
               END DO
            END SELECT
            CALL lbc_lnk( zwz, 'F', 1._wp )
            !
            ftne(1,:) = 0._wp ; ftnw(1,:) = 0._wp ; ftse(1,:) = 0._wp ; ftsw(1,:) = 0._wp
            DO jj = 2, jpj
               DO ji = 2, jpi
                  ftne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
                  ftnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
                  ftse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
                  ftsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
               END DO
            END DO
            !
         ELSE                                !== all other schemes (ENE, ENS, MIX)
            zwz(:,:) = 0._wp
            zhf(:,:) = 0._wp
            
!!gm  assume 0 in both cases (xhich is almost surely WRONG ! ) as hvatf has been removed 
!!gm    A priori a better value should be something like :
!!gm          zhf(i,j) = masked sum of  ht(i,j) , ht(i+1,j) , ht(i,j+1) , (i+1,j+1) 
!!gm                     divided by the sum of the corresponding mask 
!!gm 
!!            
              IF ( .not. ln_sco ) THEN
  
   !!gm  agree the JC comment  : this should be done in a much clear way
  
   ! JC: It not clear yet what should be the depth at f-points over land in z-coordinate case
   !     Set it to zero for the time being 
   !              IF( rn_hmin < 0._wp ) THEN    ;   jk = - INT( rn_hmin )                                      ! from a nb of level
   !              ELSE                          ;   jk = MINLOC( gdepw_0, mask = gdepw_0 > rn_hmin, dim = 1 )  ! from a depth
   !              ENDIF
   !              zhf(:,:) = gdepw_0(:,:,jk+1)
               ELSE
                 !zhf(:,:) = hbatf(:,:)
                 DO jj = 1, jpjm1
                   DO ji = 1, jpim1
                     zhf(ji,jj) = MAX( 0._wp,                                &
                                & ( ht_0(ji  ,jj  )*tmask(ji  ,jj  ,1) +     &
                                &   ht_0(ji+1,jj  )*tmask(ji+1,jj  ,1) +     &
                                &   ht_0(ji  ,jj+1)*tmask(ji  ,jj+1,1) +     &
                                &   ht_0(ji+1,jj+1)*tmask(ji+1,jj+1,1) ) /   &
                                & ( tmask(ji  ,jj  ,1) + tmask(ji+1,jj  ,1) +&
                                &   tmask(ji  ,jj+1,1) + tmask(ji+1,jj+1,1) +&
                                &   rsmall  ) )
                   END DO
                 END DO
              END IF
  
              DO jj = 1, jpjm1
                 zhf(:,jj) = zhf(:,jj) * (1._wp- umask(:,jj,1) * umask(:,jj+1,1))
              END DO
!!gm end

            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  zhf(:,jj) = zhf(:,jj) + e3f_n(:,jj,jk) * umask(:,jj,jk) * umask(:,jj+1,jk)
               END DO
            END DO
            CALL lbc_lnk( zhf, 'F', 1._wp )
            ! JC: TBC. hf should be greater than 0 
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( zhf(ji,jj) /= 0._wp )   zwz(ji,jj) = 1._wp / zhf(ji,jj) ! zhf is actually hf here but it saves an array
               END DO
            END DO
            zwz(:,:) = ff_f(:,:) * zwz(:,:)
         ENDIF
      ENDIF
      !
      ! If forward start at previous time step, and centered integration, 
      ! then update averaging weights:
      IF (.NOT.ln_bt_fw .AND.( neuler==0 .AND. kt==nit000+1 ) ) THEN
         ll_fw_start=.FALSE.
         CALL ts_wgt(ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2)
      ENDIF
                          
      ! -----------------------------------------------------------------------------
      !  Phase 1 : Coupling between general trend and barotropic estimates (1st step)
      ! -----------------------------------------------------------------------------
      !      
      !
      !                                   !* e3*d/dt(Ua) (Vertically integrated)
      !                                   ! --------------------------------------------------
      zu_frc(:,:) = 0._wp
      zv_frc(:,:) = 0._wp
      !
      DO jk = 1, jpkm1
         zu_frc(:,:) = zu_frc(:,:) + e3u_n(:,:,jk) * ua(:,:,jk) * umask(:,:,jk)
         zv_frc(:,:) = zv_frc(:,:) + e3v_n(:,:,jk) * va(:,:,jk) * vmask(:,:,jk)         
      END DO
      !
      zu_frc(:,:) = zu_frc(:,:) * r1_hu_n(:,:)
      zv_frc(:,:) = zv_frc(:,:) * r1_hv_n(:,:)
      !
      !
      !                                   !* baroclinic momentum trend (remove the vertical mean trend)
      DO jk = 1, jpkm1                    ! -----------------------------------------------------------
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - zu_frc(ji,jj) * umask(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) - zv_frc(ji,jj) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      
!!gm  Question here when removing the Vertically integrated trends, we remove the vertically integrated NL trends on momentum....
!!gm  Is it correct to do so ?   I think so...
      
      
      !                                   !* barotropic Coriolis trends (vorticity scheme dependent)
      !                                   ! --------------------------------------------------------
      zwx(:,:) = un_b(:,:) * hu_n(:,:) * e2u(:,:)        ! now fluxes 
      zwy(:,:) = vn_b(:,:) * hv_n(:,:) * e1v(:,:)
      !
      IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) * r1_e1u(ji,jj)
               zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) * r1_e1u(ji,jj)
               zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) * r1_e2v(ji,jj)
               zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) * r1_e2v(ji,jj)
               ! energy conserving formulation for planetary vorticity term
               zu_trd(ji,jj) = z1_4 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
               zv_trd(ji,jj) =-z1_4 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
            END DO
         END DO
         !
      ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zy1 =   z1_8 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1) &
                 &            + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) * r1_e1u(ji,jj)
               zx1 = - z1_8 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1) &
                 &            + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) * r1_e2v(ji,jj)
               zu_trd(ji,jj)  = zy1 * ( zwz(ji  ,jj-1) + zwz(ji,jj) )
               zv_trd(ji,jj)  = zx1 * ( zwz(ji-1,jj  ) + zwz(ji,jj) )
            END DO
         END DO
         !
      ELSEIF ( ln_dynvor_een ) THEN  ! enstrophy and energy conserving scheme
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zu_trd(ji,jj) = + z1_12 * r1_e1u(ji,jj) * (  ftne(ji,jj  ) * zwy(ji  ,jj  ) &
                &                                         + ftnw(ji+1,jj) * zwy(ji+1,jj  ) &
                &                                         + ftse(ji,jj  ) * zwy(ji  ,jj-1) &
                &                                         + ftsw(ji+1,jj) * zwy(ji+1,jj-1) )
               zv_trd(ji,jj) = - z1_12 * r1_e2v(ji,jj) * (  ftsw(ji,jj+1) * zwx(ji-1,jj+1) &
                &                                         + ftse(ji,jj+1) * zwx(ji  ,jj+1) &
                &                                         + ftnw(ji,jj  ) * zwx(ji-1,jj  ) &
                &                                         + ftne(ji,jj  ) * zwx(ji  ,jj  ) )
            END DO
         END DO
         !
      ENDIF 
      !
      !                                   !* Right-Hand-Side of the barotropic momentum equation
      !                                   ! ----------------------------------------------------
      IF( .NOT.ln_linssh ) THEN                 ! Variable volume : remove surface pressure gradient
        IF( ln_wd ) THEN                        ! Calculating and applying W/D gravity filters
           DO jj = 2, jpjm1
              DO ji = 2, jpim1 
                ll_tmp1 = MIN(   sshn(ji,jj)               ,   sshn(ji+1,jj) ) >                &
                     &    MAX( -ht_wd(ji,jj)               , -ht_wd(ji+1,jj) ) .AND.            &
                     &    MAX(   sshn(ji,jj) + ht_wd(ji,jj),   sshn(ji+1,jj) + ht_wd(ji+1,jj) ) &
                     &                                                         > rn_wdmin1 + rn_wdmin2
                ll_tmp2 = ( ABS( sshn(ji+1,jj)             -   sshn(ji  ,jj))  > 1.E-12 ).AND.( &
                     &    MAX(   sshn(ji,jj)               ,   sshn(ji+1,jj) ) >                &
                     &    MAX( -ht_wd(ji,jj)               , -ht_wd(ji+1,jj) ) + rn_wdmin1 + rn_wdmin2 )
   
                IF(ll_tmp1) THEN
                  zcpx(ji,jj) = 1.0_wp
                ELSE IF(ll_tmp2) THEN
                  ! no worries about  sshn(ji+1,jj) -  sshn(ji  ,jj) = 0, it won't happen ! here
                  zcpx(ji,jj) = ABS( (sshn(ji+1,jj) + ht_wd(ji+1,jj) - sshn(ji,jj) - ht_wd(ji,jj)) &
                              &    / (sshn(ji+1,jj) -  sshn(ji  ,jj)) )
                ELSE
                  zcpx(ji,jj) = 0._wp
                END IF
         
                ll_tmp1 = MIN(   sshn(ji,jj)               ,   sshn(ji,jj+1) ) >                &
                     &    MAX( -ht_wd(ji,jj)               , -ht_wd(ji,jj+1) ) .AND.            &
                     &    MAX(   sshn(ji,jj) + ht_wd(ji,jj),   sshn(ji,jj+1) + ht_wd(ji,jj+1) ) &
                     &                                                         > rn_wdmin1 + rn_wdmin2
                ll_tmp2 = ( ABS( sshn(ji,jj)               -   sshn(ji,jj+1))  > 1.E-12 ).AND.( &
                     &    MAX(   sshn(ji,jj)               ,   sshn(ji,jj+1) ) >                &
                     &    MAX( -ht_wd(ji,jj)               , -ht_wd(ji,jj+1) ) + rn_wdmin1 + rn_wdmin2 )
   
                IF(ll_tmp1) THEN
                  zcpy(ji,jj) = 1.0_wp
                ELSE IF(ll_tmp2) THEN
                  ! no worries about  sshn(ji,jj+1) -  sshn(ji,jj  ) = 0, it won't happen ! here
                  zcpy(ji,jj) = ABS( (sshn(ji,jj+1) + ht_wd(ji,jj+1) - sshn(ji,jj) - ht_wd(ji,jj)) &
                              &    / (sshn(ji,jj+1) -  sshn(ji,jj  )) )
                ELSE
                  zcpy(ji,jj) = 0._wp
                END IF
              END DO
           END DO
 
           DO jj = 2, jpjm1
              DO ji = 2, jpim1
                 zu_trd(ji,jj) = zu_trd(ji,jj) - grav * ( sshn(ji+1,jj  ) - sshn(ji  ,jj ) )   &
                        &                        * r1_e1u(ji,jj) * zcpx(ji,jj)
                 zv_trd(ji,jj) = zv_trd(ji,jj) - grav * ( sshn(ji  ,jj+1) - sshn(ji  ,jj ) )   &
                        &                        * r1_e2v(ji,jj) * zcpy(ji,jj)
              END DO
           END DO

         ELSE

           DO jj = 2, jpjm1
              DO ji = fs_2, fs_jpim1   ! vector opt.
                 zu_trd(ji,jj) = zu_trd(ji,jj) - grav * (  sshn(ji+1,jj  ) - sshn(ji  ,jj  )  ) * r1_e1u(ji,jj)
                 zv_trd(ji,jj) = zv_trd(ji,jj) - grav * (  sshn(ji  ,jj+1) - sshn(ji  ,jj  )  ) * r1_e2v(ji,jj) 
              END DO
           END DO
        ENDIF

      ENDIF

      DO jj = 2, jpjm1                          ! Remove coriolis term (and possibly spg) from barotropic trend
         DO ji = fs_2, fs_jpim1
             zu_frc(ji,jj) = zu_frc(ji,jj) - zu_trd(ji,jj) * ssumask(ji,jj)
             zv_frc(ji,jj) = zv_frc(ji,jj) - zv_trd(ji,jj) * ssvmask(ji,jj)
          END DO
      END DO 
      !
      !						! Add bottom stress contribution from baroclinic velocities:      
      IF (ln_bt_fw) THEN
         DO jj = 2, jpjm1                          
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ikbu = mbku(ji,jj)       
               ikbv = mbkv(ji,jj)    
               zwx(ji,jj) = un(ji,jj,ikbu) - un_b(ji,jj) ! NOW bottom baroclinic velocities
               zwy(ji,jj) = vn(ji,jj,ikbv) - vn_b(ji,jj)
            END DO
         END DO
      ELSE
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ikbu = mbku(ji,jj)       
               ikbv = mbkv(ji,jj)    
               zwx(ji,jj) = ub(ji,jj,ikbu) - ub_b(ji,jj) ! BEFORE bottom baroclinic velocities
               zwy(ji,jj) = vb(ji,jj,ikbv) - vb_b(ji,jj)
            END DO
         END DO
      ENDIF
      !
      ! Note that the "unclipped" bottom friction parameter is used even with explicit drag
      IF( ln_wd ) THEN
        zu_frc(:,:) = zu_frc(:,:) + MAX(r1_hu_n(:,:) * bfrua(:,:),-1._wp / rdtbt) * zwx(:,:)
        zv_frc(:,:) = zv_frc(:,:) + MAX(r1_hv_n(:,:) * bfrva(:,:),-1._wp / rdtbt) * zwy(:,:)
      ELSE
        zu_frc(:,:) = zu_frc(:,:) + r1_hu_n(:,:) * bfrua(:,:) * zwx(:,:)
        zv_frc(:,:) = zv_frc(:,:) + r1_hv_n(:,:) * bfrva(:,:) * zwy(:,:)
      END IF
      !
      !                                         ! Add top stress contribution from baroclinic velocities:      
      IF( ln_bt_fw ) THEN
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               iktu = miku(ji,jj)
               iktv = mikv(ji,jj)
               zwx(ji,jj) = un(ji,jj,iktu) - un_b(ji,jj) ! NOW top baroclinic velocities
               zwy(ji,jj) = vn(ji,jj,iktv) - vn_b(ji,jj)
            END DO
         END DO
      ELSE
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               iktu = miku(ji,jj)
               iktv = mikv(ji,jj)
               zwx(ji,jj) = ub(ji,jj,iktu) - ub_b(ji,jj) ! BEFORE top baroclinic velocities
               zwy(ji,jj) = vb(ji,jj,iktv) - vb_b(ji,jj)
            END DO
         END DO
      ENDIF
      !
      ! Note that the "unclipped" top friction parameter is used even with explicit drag
      zu_frc(:,:) = zu_frc(:,:) + r1_hu_n(:,:) * tfrua(:,:) * zwx(:,:)
      zv_frc(:,:) = zv_frc(:,:) + r1_hv_n(:,:) * tfrva(:,:) * zwy(:,:)
      !       
      IF (ln_bt_fw) THEN                        ! Add wind forcing
         zu_frc(:,:) =  zu_frc(:,:) + zraur * utau(:,:) * r1_hu_n(:,:)
         zv_frc(:,:) =  zv_frc(:,:) + zraur * vtau(:,:) * r1_hv_n(:,:)
      ELSE
         zu_frc(:,:) =  zu_frc(:,:) + zraur * z1_2 * ( utau_b(:,:) + utau(:,:) ) * r1_hu_n(:,:)
         zv_frc(:,:) =  zv_frc(:,:) + zraur * z1_2 * ( vtau_b(:,:) + vtau(:,:) ) * r1_hv_n(:,:)
      ENDIF  
      !
      IF ( ln_apr_dyn ) THEN                    ! Add atm pressure forcing
         IF (ln_bt_fw) THEN
            DO jj = 2, jpjm1              
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_spg =  grav * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj) ) * r1_e1u(ji,jj)
                  zv_spg =  grav * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj) ) * r1_e2v(ji,jj)
                  zu_frc(ji,jj) = zu_frc(ji,jj) + zu_spg
                  zv_frc(ji,jj) = zv_frc(ji,jj) + zv_spg
               END DO
            END DO
         ELSE
            DO jj = 2, jpjm1              
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_spg =  grav * z1_2 * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj)    &
                      &                    + ssh_ibb(ji+1,jj  ) - ssh_ibb(ji,jj)  ) * r1_e1u(ji,jj)
                  zv_spg =  grav * z1_2 * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj)    &
                      &                    + ssh_ibb(ji  ,jj+1) - ssh_ibb(ji,jj)  ) * r1_e2v(ji,jj)
                  zu_frc(ji,jj) = zu_frc(ji,jj) + zu_spg
                  zv_frc(ji,jj) = zv_frc(ji,jj) + zv_spg
               END DO
            END DO
         ENDIF 
      ENDIF
      !                                   !* Right-Hand-Side of the barotropic ssh equation
      !                                   ! -----------------------------------------------
      !                                         ! Surface net water flux and rivers
      IF (ln_bt_fw) THEN
         zssh_frc(:,:) = zraur * ( emp(:,:) - rnf(:,:) + fwfisf(:,:) )
      ELSE
         zssh_frc(:,:) = zraur * z1_2 * (  emp(:,:) + emp_b(:,:) - rnf(:,:) - rnf_b(:,:)   &
                &                        + fwfisf(:,:) + fwfisf_b(:,:)                     )
      ENDIF
      !
      IF( ln_sdw ) THEN                         ! Stokes drift divergence added if necessary
         zssh_frc(:,:) = zssh_frc(:,:) + div_sd(:,:)
      ENDIF
      !
#if defined key_asminc
      !                                         ! Include the IAU weighted SSH increment
      IF( lk_asminc .AND. ln_sshinc .AND. ln_asmiau ) THEN
         zssh_frc(:,:) = zssh_frc(:,:) - ssh_iau(:,:)
      ENDIF
#endif
      !                                   !* Fill boundary data arrays for AGRIF
      !                                   ! ------------------------------------
#if defined key_agrif
         IF( .NOT.Agrif_Root() ) CALL agrif_dta_ts( kt )
#endif
      !
      ! -----------------------------------------------------------------------
      !  Phase 2 : Integration of the barotropic equations 
      ! -----------------------------------------------------------------------
      !
      !                                             ! ==================== !
      !                                             !    Initialisations   !
      !                                             ! ==================== !  
      ! Initialize barotropic variables:      
      IF( ll_init )THEN
         sshbb_e(:,:) = 0._wp
         ubb_e  (:,:) = 0._wp
         vbb_e  (:,:) = 0._wp
         sshb_e (:,:) = 0._wp
         ub_e   (:,:) = 0._wp
         vb_e   (:,:) = 0._wp
      ENDIF

      !
      IF (ln_bt_fw) THEN                  ! FORWARD integration: start from NOW fields                    
         sshn_e(:,:) =    sshn(:,:)            
         un_e  (:,:) =    un_b(:,:)            
         vn_e  (:,:) =    vn_b(:,:)
         !
         hu_e  (:,:) =    hu_n(:,:)       
         hv_e  (:,:) =    hv_n(:,:) 
         hur_e (:,:) = r1_hu_n(:,:)    
         hvr_e (:,:) = r1_hv_n(:,:)
      ELSE                                ! CENTRED integration: start from BEFORE fields
         sshn_e(:,:) =    sshb(:,:)
         un_e  (:,:) =    ub_b(:,:)         
         vn_e  (:,:) =    vb_b(:,:)
         !
         hu_e  (:,:) =    hu_b(:,:)       
         hv_e  (:,:) =    hv_b(:,:) 
         hur_e (:,:) = r1_hu_b(:,:)    
         hvr_e (:,:) = r1_hv_b(:,:)
      ENDIF
      !
      !
      !
      ! Initialize sums:
      ua_b  (:,:) = 0._wp       ! After barotropic velocities (or transport if flux form)          
      va_b  (:,:) = 0._wp
      ssha  (:,:) = 0._wp       ! Sum for after averaged sea level
      un_adv(:,:) = 0._wp       ! Sum for now transport issued from ts loop
      vn_adv(:,:) = 0._wp
      !                                             ! ==================== !
      DO jn = 1, icycle                             !  sub-time-step loop  !
         !                                          ! ==================== !
         !                                                !* Update the forcing (BDY and tides)
         !                                                !  ------------------
         ! Update only tidal forcing at open boundaries
         IF( ln_bdy      .AND. ln_tide )   CALL bdy_dta_tides( kt, kit=jn, time_offset= noffset+1 )
         IF( ln_tide_pot .AND. ln_tide )   CALL upd_tide     ( kt, kit=jn, time_offset= noffset   )
         !
         ! Set extrapolation coefficients for predictor step:
         IF ((jn<3).AND.ll_init) THEN      ! Forward           
           za1 = 1._wp                                          
           za2 = 0._wp                        
           za3 = 0._wp                        
         ELSE                              ! AB3-AM4 Coefficients: bet=0.281105 
           za1 =  1.781105_wp              ! za1 =   3/2 +   bet
           za2 = -1.06221_wp               ! za2 = -(1/2 + 2*bet)
           za3 =  0.281105_wp              ! za3 = bet
         ENDIF

         ! Extrapolate barotropic velocities at step jit+0.5:
         ua_e(:,:) = za1 * un_e(:,:) + za2 * ub_e(:,:) + za3 * ubb_e(:,:)
         va_e(:,:) = za1 * vn_e(:,:) + za2 * vb_e(:,:) + za3 * vbb_e(:,:)

         IF( .NOT.ln_linssh ) THEN                        !* Update ocean depth (variable volume case only)
            !                                             !  ------------------
            ! Extrapolate Sea Level at step jit+0.5:
            zsshp2_e(:,:) = za1 * sshn_e(:,:)  + za2 * sshb_e(:,:) + za3 * sshbb_e(:,:)
            !
            DO jj = 2, jpjm1                                    ! Sea Surface Height at u- & v-points
               DO ji = 2, fs_jpim1   ! Vector opt.
                  zwx(ji,jj) = z1_2 * ssumask(ji,jj)  * r1_e1e2u(ji,jj)     &
                     &              * ( e1e2t(ji  ,jj) * zsshp2_e(ji  ,jj)  &
                     &              +   e1e2t(ji+1,jj) * zsshp2_e(ji+1,jj) )
                  zwy(ji,jj) = z1_2 * ssvmask(ji,jj)  * r1_e1e2v(ji,jj)     &
                     &              * ( e1e2t(ji,jj  ) * zsshp2_e(ji,jj  )  &
                     &              +   e1e2t(ji,jj+1) * zsshp2_e(ji,jj+1) )
               END DO
            END DO
            CALL lbc_lnk_multi( zwx, 'U', 1._wp, zwy, 'V', 1._wp )
            !
            zhup2_e (:,:) = hu_0(:,:) + zwx(:,:)                ! Ocean depth at U- and V-points
            zhvp2_e (:,:) = hv_0(:,:) + zwy(:,:)
         ELSE
            zhup2_e (:,:) = hu_n(:,:)
            zhvp2_e (:,:) = hv_n(:,:)
         ENDIF
         !                                                !* after ssh
         !                                                !  -----------
         ! One should enforce volume conservation at open boundaries here
         ! considering fluxes below:
         !
         zwx(:,:) = e2u(:,:) * ua_e(:,:) * zhup2_e(:,:)         ! fluxes at jn+0.5
         zwy(:,:) = e1v(:,:) * va_e(:,:) * zhvp2_e(:,:)
         !
#if defined key_agrif
         ! Set fluxes during predictor step to ensure volume conservation
         IF( .NOT.Agrif_Root() .AND. ln_bt_fw ) THEN
            IF((nbondi == -1).OR.(nbondi == 2)) THEN
               DO jj=1,jpj
                  zwx(2,jj) = ubdy_w(jj) * e2u(2,jj)
               END DO
            ENDIF
            IF((nbondi ==  1).OR.(nbondi == 2)) THEN
               DO jj=1,jpj
                  zwx(nlci-2,jj) = ubdy_e(jj) * e2u(nlci-2,jj)
               END DO
            ENDIF
            IF((nbondj == -1).OR.(nbondj == 2)) THEN
               DO ji=1,jpi
                  zwy(ji,2) = vbdy_s(ji) * e1v(ji,2)
               END DO
            ENDIF
            IF((nbondj ==  1).OR.(nbondj == 2)) THEN
               DO ji=1,jpi
                  zwy(ji,nlcj-2) = vbdy_n(ji) * e1v(ji,nlcj-2)
               END DO
            ENDIF
         ENDIF
#endif
         IF( ln_wd ) CALL wad_lmt_bt(zwx, zwy, sshn_e, zssh_frc, rdtbt)
         !
         ! Sum over sub-time-steps to compute advective velocities
         za2 = wgtbtp2(jn)
         un_adv(:,:) = un_adv(:,:) + za2 * zwx(:,:) * r1_e2u(:,:)
         vn_adv(:,:) = vn_adv(:,:) + za2 * zwy(:,:) * r1_e1v(:,:)
         !
         ! Set next sea level:
         DO jj = 2, jpjm1                                 
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zhdiv(ji,jj) = (   zwx(ji,jj) - zwx(ji-1,jj)   &
                  &             + zwy(ji,jj) - zwy(ji,jj-1)   ) * r1_e1e2t(ji,jj)
            END DO
         END DO
         ssha_e(:,:) = (  sshn_e(:,:) - rdtbt * ( zssh_frc(:,:) + zhdiv(:,:) )  ) * ssmask(:,:)
         
         CALL lbc_lnk( ssha_e, 'T',  1._wp )

         ! Duplicate sea level across open boundaries (this is only cosmetic if linssh=T)
         IF( ln_bdy )   CALL bdy_ssh( ssha_e )
#if defined key_agrif
         IF( .NOT.Agrif_Root() )   CALL agrif_ssh_ts( jn )
#endif
         !  
         ! Sea Surface Height at u-,v-points (vvl case only)
         IF( .NOT.ln_linssh ) THEN                                
            DO jj = 2, jpjm1
               DO ji = 2, jpim1      ! NO Vector Opt.
                  zsshu_a(ji,jj) = z1_2 * ssumask(ji,jj) * r1_e1e2u(ji,jj)    &
                     &              * ( e1e2t(ji  ,jj  )  * ssha_e(ji  ,jj  ) &
                     &              +   e1e2t(ji+1,jj  )  * ssha_e(ji+1,jj  ) )
                  zsshv_a(ji,jj) = z1_2 * ssvmask(ji,jj) * r1_e1e2v(ji,jj)    &
                     &              * ( e1e2t(ji  ,jj  )  * ssha_e(ji  ,jj  ) &
                     &              +   e1e2t(ji  ,jj+1)  * ssha_e(ji  ,jj+1) )
               END DO
            END DO
            CALL lbc_lnk_multi( zsshu_a, 'U', 1._wp, zsshv_a, 'V', 1._wp )
         ENDIF   
         !                                 
         ! Half-step back interpolation of SSH for surface pressure computation:
         !----------------------------------------------------------------------
         IF ((jn==1).AND.ll_init) THEN
           za0=1._wp                        ! Forward-backward
           za1=0._wp                           
           za2=0._wp
           za3=0._wp
         ELSEIF ((jn==2).AND.ll_init) THEN  ! AB2-AM3 Coefficients; bet=0 ; gam=-1/6 ; eps=1/12
           za0= 1.0833333333333_wp          ! za0 = 1-gam-eps
           za1=-0.1666666666666_wp          ! za1 = gam
           za2= 0.0833333333333_wp          ! za2 = eps
           za3= 0._wp              
         ELSE                               ! AB3-AM4 Coefficients; bet=0.281105 ; eps=0.013 ; gam=0.0880 
           za0=0.614_wp                     ! za0 = 1/2 +   gam + 2*eps    
           za1=0.285_wp                     ! za1 = 1/2 - 2*gam - 3*eps 
           za2=0.088_wp                     ! za2 = gam
           za3=0.013_wp                     ! za3 = eps
         ENDIF
         !
         zsshp2_e(:,:) = za0 *  ssha_e(:,:) + za1 *  sshn_e (:,:) &
          &            + za2 *  sshb_e(:,:) + za3 *  sshbb_e(:,:)
         IF( ln_wd ) THEN                   ! Calculating and applying W/D gravity filters
           DO jj = 2, jpjm1
              DO ji = 2, jpim1 
                ll_tmp1 = MIN( zsshp2_e(ji,jj)               , zsshp2_e(ji+1,jj) ) >                &
                     &    MAX(   -ht_wd(ji,jj)               ,   -ht_wd(ji+1,jj) ) .AND.            &
                     &    MAX( zsshp2_e(ji,jj) + ht_wd(ji,jj), zsshp2_e(ji+1,jj) + ht_wd(ji+1,jj) ) &
                     &                                                             > rn_wdmin1 + rn_wdmin2
                ll_tmp2 = (ABS(zsshp2_e(ji,jj)               - zsshp2_e(ji+1,jj))  > 1.E-12 ).AND.( &
                     &    MAX( zsshp2_e(ji,jj)               , zsshp2_e(ji+1,jj) ) >                &
                     &    MAX(   -ht_wd(ji,jj)               ,   -ht_wd(ji+1,jj) ) + rn_wdmin1 + rn_wdmin2 )
   
                IF(ll_tmp1) THEN
                  zcpx(ji,jj) = 1.0_wp
                ELSE IF(ll_tmp2) THEN
                  ! no worries about  zsshp2_e(ji+1,jj) - zsshp2_e(ji  ,jj) = 0, it won't happen ! here
                  zcpx(ji,jj) = ABS( (zsshp2_e(ji+1,jj) +    ht_wd(ji+1,jj) - zsshp2_e(ji,jj) - ht_wd(ji,jj)) &
                              &    / (zsshp2_e(ji+1,jj) - zsshp2_e(ji  ,jj)) )
                ELSE
                  zcpx(ji,jj) = 0._wp
                END IF
         
                ll_tmp1 = MIN( zsshp2_e(ji,jj)               , zsshp2_e(ji,jj+1) ) >                &
                     &    MAX(   -ht_wd(ji,jj)               ,   -ht_wd(ji,jj+1) ) .AND.            &
                     &    MAX( zsshp2_e(ji,jj) + ht_wd(ji,jj), zsshp2_e(ji,jj+1) + ht_wd(ji,jj+1) ) &
                     &                                                             > rn_wdmin1 + rn_wdmin2
                ll_tmp2 = (ABS(zsshp2_e(ji,jj)               - zsshp2_e(ji,jj+1))  > 1.E-12 ).AND.( &
                     &    MAX( zsshp2_e(ji,jj)               , zsshp2_e(ji,jj+1) ) >                &
                     &    MAX(   -ht_wd(ji,jj)               ,   -ht_wd(ji,jj+1) ) + rn_wdmin1 + rn_wdmin2 )
   
                IF(ll_tmp1) THEN
                  zcpy(ji,jj) = 1.0_wp
                ELSE IF(ll_tmp2) THEN
                  ! no worries about  zsshp2_e(ji,jj+1) - zsshp2_e(ji,jj  ) = 0, it won't happen ! here
                  zcpy(ji,jj) = ABS( (zsshp2_e(ji,jj+1) +    ht_wd(ji,jj+1) - zsshp2_e(ji,jj) - ht_wd(ji,jj)) &
                              &    / (zsshp2_e(ji,jj+1) - zsshp2_e(ji,jj  )) )
                ELSE
                  zcpy(ji,jj) = 0._wp
                END IF
              END DO
           END DO
         END IF
         !
         ! Compute associated depths at U and V points:
         IF( .NOT.ln_linssh  .AND. .NOT.ln_dynadv_vec ) THEN 	!* Vector form
            !                                        
            DO jj = 2, jpjm1                            
               DO ji = 2, jpim1
                  zx1 = z1_2 * ssumask(ji  ,jj) *  r1_e1e2u(ji  ,jj)    &
                     &      * ( e1e2t(ji  ,jj  ) * zsshp2_e(ji  ,jj)    &
                     &      +   e1e2t(ji+1,jj  ) * zsshp2_e(ji+1,jj  ) )
                  zy1 = z1_2 * ssvmask(ji  ,jj) *  r1_e1e2v(ji  ,jj  )  &
                     &       * ( e1e2t(ji ,jj  ) * zsshp2_e(ji  ,jj  )  &
                     &       +   e1e2t(ji ,jj+1) * zsshp2_e(ji  ,jj+1) )
                  zhust_e(ji,jj) = hu_0(ji,jj) + zx1 
                  zhvst_e(ji,jj) = hv_0(ji,jj) + zy1
               END DO
            END DO

         ENDIF
         !
         ! Add Coriolis trend:
         ! zwz array below or triads normally depend on sea level with ln_linssh=F and should be updated
         ! at each time step. We however keep them constant here for optimization.
         ! Recall that zwx and zwy arrays hold fluxes at this stage:
         ! zwx(:,:) = e2u(:,:) * ua_e(:,:) * zhup2_e(:,:)   ! fluxes at jn+0.5
         ! zwy(:,:) = e1v(:,:) * va_e(:,:) * zhvp2_e(:,:)
         !
         IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN     !==  energy conserving or mixed scheme  ==!
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zy1 = ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1) ) * r1_e1u(ji,jj)
                  zy2 = ( zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) * r1_e1u(ji,jj)
                  zx1 = ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1) ) * r1_e2v(ji,jj)
                  zx2 = ( zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) * r1_e2v(ji,jj)
                  zu_trd(ji,jj) = z1_4 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
                  zv_trd(ji,jj) =-z1_4 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
               END DO
            END DO
            !
         ELSEIF ( ln_dynvor_ens ) THEN                   !==  enstrophy conserving scheme  ==!
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zy1 =   z1_8 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1) &
                   &             + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) * r1_e1u(ji,jj)
                  zx1 = - z1_8 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1) &
                   &             + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) * r1_e2v(ji,jj)
                  zu_trd(ji,jj)  = zy1 * ( zwz(ji  ,jj-1) + zwz(ji,jj) )
                  zv_trd(ji,jj)  = zx1 * ( zwz(ji-1,jj  ) + zwz(ji,jj) )
               END DO
            END DO
            !
         ELSEIF ( ln_dynvor_een ) THEN                   !==  energy and enstrophy conserving scheme  ==!
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_trd(ji,jj) = + z1_12 * r1_e1u(ji,jj) * (  ftne(ji,jj  ) * zwy(ji  ,jj  ) &
                     &                                       + ftnw(ji+1,jj) * zwy(ji+1,jj  ) &
                     &                                       + ftse(ji,jj  ) * zwy(ji  ,jj-1) & 
                     &                                       + ftsw(ji+1,jj) * zwy(ji+1,jj-1) )
                  zv_trd(ji,jj) = - z1_12 * r1_e2v(ji,jj) * (  ftsw(ji,jj+1) * zwx(ji-1,jj+1) & 
                     &                                       + ftse(ji,jj+1) * zwx(ji  ,jj+1) &
                     &                                       + ftnw(ji,jj  ) * zwx(ji-1,jj  ) & 
                     &                                       + ftne(ji,jj  ) * zwx(ji  ,jj  ) )
               END DO
            END DO
            ! 
         ENDIF
         !
         ! Add tidal astronomical forcing if defined
         IF ( ln_tide .AND. ln_tide_pot ) THEN
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_spg = grav * ( pot_astro(ji+1,jj) - pot_astro(ji,jj) ) * r1_e1u(ji,jj)
                  zv_spg = grav * ( pot_astro(ji,jj+1) - pot_astro(ji,jj) ) * r1_e2v(ji,jj)
                  zu_trd(ji,jj) = zu_trd(ji,jj) + zu_spg
                  zv_trd(ji,jj) = zv_trd(ji,jj) + zv_spg
               END DO
            END DO
         ENDIF
         !
         ! Add bottom stresses:
         zu_trd(:,:) = zu_trd(:,:) + bfrua(:,:) * un_e(:,:) * hur_e(:,:)
         zv_trd(:,:) = zv_trd(:,:) + bfrva(:,:) * vn_e(:,:) * hvr_e(:,:)
         !
         ! Add top stresses:
         zu_trd(:,:) = zu_trd(:,:) + tfrua(:,:) * un_e(:,:) * hur_e(:,:)
         zv_trd(:,:) = zv_trd(:,:) + tfrva(:,:) * vn_e(:,:) * hvr_e(:,:)
         !
         ! Surface pressure trend:

         IF( ln_wd ) THEN
           DO jj = 2, jpjm1
              DO ji = 2, jpim1 
                 ! Add surface pressure gradient
                 zu_spg = - grav * ( zsshp2_e(ji+1,jj) - zsshp2_e(ji,jj) ) * r1_e1u(ji,jj)
                 zv_spg = - grav * ( zsshp2_e(ji,jj+1) - zsshp2_e(ji,jj) ) * r1_e2v(ji,jj)
                 zwx(ji,jj) = zu_spg * zcpx(ji,jj) 
                 zwy(ji,jj) = zv_spg * zcpy(ji,jj)
              END DO
           END DO
         ELSE
           DO jj = 2, jpjm1
              DO ji = fs_2, fs_jpim1   ! vector opt.
                 ! Add surface pressure gradient
                 zu_spg = - grav * ( zsshp2_e(ji+1,jj) - zsshp2_e(ji,jj) ) * r1_e1u(ji,jj)
                 zv_spg = - grav * ( zsshp2_e(ji,jj+1) - zsshp2_e(ji,jj) ) * r1_e2v(ji,jj)
                 zwx(ji,jj) = zu_spg
                 zwy(ji,jj) = zv_spg
              END DO
           END DO
         END IF

         !
         ! Set next velocities:
         IF( ln_dynadv_vec .OR. ln_linssh ) THEN 	!* Vector form
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ua_e(ji,jj) = (                                 un_e(ji,jj)   & 
                            &     + rdtbt * (                      zwx(ji,jj)   &
                            &                                 + zu_trd(ji,jj)   &
                            &                                 + zu_frc(ji,jj) ) & 
                            &   ) * ssumask(ji,jj)

                  va_e(ji,jj) = (                                 vn_e(ji,jj)   &
                            &     + rdtbt * (                      zwy(ji,jj)   &
                            &                                 + zv_trd(ji,jj)   &
                            &                                 + zv_frc(ji,jj) ) &
                            &   ) * ssvmask(ji,jj)
               END DO
            END DO
            !
         ELSE					 	                     !* Flux form
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.

                  IF( ln_wd ) THEN
                    zhura = MAX(hu_0(ji,jj) + zsshu_a(ji,jj), rn_wdmin1)
                    zhvra = MAX(hv_0(ji,jj) + zsshv_a(ji,jj), rn_wdmin1)
                  ELSE
                    zhura = hu_0(ji,jj) + zsshu_a(ji,jj)
                    zhvra = hv_0(ji,jj) + zsshv_a(ji,jj)
                  END IF
                  zhura = ssumask(ji,jj)/(zhura + 1._wp - ssumask(ji,jj))
                  zhvra = ssvmask(ji,jj)/(zhvra + 1._wp - ssvmask(ji,jj))

                  ua_e(ji,jj) = (                hu_e(ji,jj)  *   un_e(ji,jj)   & 
                            &     + rdtbt * ( zhust_e(ji,jj)  *    zwx(ji,jj)   & 
                            &               + zhup2_e(ji,jj)  * zu_trd(ji,jj)   &
                            &               +    hu_n(ji,jj)  * zu_frc(ji,jj) ) &
                            &   ) * zhura

                  va_e(ji,jj) = (                hv_e(ji,jj)  *   vn_e(ji,jj)   &
                            &     + rdtbt * ( zhvst_e(ji,jj)  *    zwy(ji,jj)   &
                            &               + zhvp2_e(ji,jj)  * zv_trd(ji,jj)   &
                            &               +    hv_n(ji,jj)  * zv_frc(ji,jj) ) &
                            &   ) * zhvra
               END DO
            END DO
         ENDIF
         !
         IF( .NOT.ln_linssh ) THEN                     !* Update ocean depth (variable volume case only)
            IF( ln_wd ) THEN
              hu_e (:,:) = MAX(hu_0(:,:) + zsshu_a(:,:), rn_wdmin1)
              hv_e (:,:) = MAX(hv_0(:,:) + zsshv_a(:,:), rn_wdmin1)
            ELSE
              hu_e (:,:) = hu_0(:,:) + zsshu_a(:,:)
              hv_e (:,:) = hv_0(:,:) + zsshv_a(:,:)
            END IF
            hur_e(:,:) = ssumask(:,:) / ( hu_e(:,:) + 1._wp - ssumask(:,:) )
            hvr_e(:,:) = ssvmask(:,:) / ( hv_e(:,:) + 1._wp - ssvmask(:,:) )
            !
         ENDIF
         !                                             !* domain lateral boundary
         CALL lbc_lnk_multi( ua_e, 'U', -1._wp, va_e , 'V', -1._wp )
         !
         !                                                 ! open boundaries
         IF( ln_bdy )   CALL bdy_dyn2d( jn, ua_e, va_e, un_e, vn_e, hur_e, hvr_e, ssha_e )
#if defined key_agrif                                                           
         IF( .NOT.Agrif_Root() )  CALL agrif_dyn_ts( jn )  ! Agrif
#endif
         !                                             !* Swap
         !                                             !  ----
         ubb_e  (:,:) = ub_e  (:,:)
         ub_e   (:,:) = un_e  (:,:)
         un_e   (:,:) = ua_e  (:,:)
         !
         vbb_e  (:,:) = vb_e  (:,:)
         vb_e   (:,:) = vn_e  (:,:)
         vn_e   (:,:) = va_e  (:,:)
         !
         sshbb_e(:,:) = sshb_e(:,:)
         sshb_e (:,:) = sshn_e(:,:)
         sshn_e (:,:) = ssha_e(:,:)

         !                                             !* Sum over whole bt loop
         !                                             !  ----------------------
         za1 = wgtbtp1(jn)                                    
         IF( ln_dynadv_vec .OR. ln_linssh ) THEN    ! Sum velocities
            ua_b  (:,:) = ua_b  (:,:) + za1 * ua_e  (:,:) 
            va_b  (:,:) = va_b  (:,:) + za1 * va_e  (:,:) 
         ELSE                                              ! Sum transports
            ua_b  (:,:) = ua_b  (:,:) + za1 * ua_e  (:,:) * hu_e (:,:)
            va_b  (:,:) = va_b  (:,:) + za1 * va_e  (:,:) * hv_e (:,:)
         ENDIF
         !                       			   ! Sum sea level
         ssha(:,:) = ssha(:,:) + za1 * ssha_e(:,:)
         !                                                 ! ==================== !
      END DO                                               !        end loop      !
      !                                                    ! ==================== !
      ! -----------------------------------------------------------------------------
      ! Phase 3. update the general trend with the barotropic trend
      ! -----------------------------------------------------------------------------
      !
      ! Set advection velocity correction:
      zwx(:,:) = un_adv(:,:)
      zwy(:,:) = vn_adv(:,:)
      IF( ( kt == nit000 .AND. neuler==0 ) .OR. .NOT.ln_bt_fw ) THEN     
         un_adv(:,:) = zwx(:,:) * r1_hu_n(:,:)
         vn_adv(:,:) = zwy(:,:) * r1_hv_n(:,:)
      ELSE
         un_adv(:,:) = z1_2 * ( ub2_b(:,:) + zwx(:,:) ) * r1_hu_n(:,:)
         vn_adv(:,:) = z1_2 * ( vb2_b(:,:) + zwy(:,:) ) * r1_hv_n(:,:)
      END IF

      IF( ln_bt_fw ) THEN ! Save integrated transport for next computation
         ub2_b(:,:) = zwx(:,:)
         vb2_b(:,:) = zwy(:,:)
      ENDIF
      !
      ! Update barotropic trend:
      IF( ln_dynadv_vec .OR. ln_linssh ) THEN
         DO jk=1,jpkm1
            ua(:,:,jk) = ua(:,:,jk) + ( ua_b(:,:) - ub_b(:,:) ) * z1_2dt_b
            va(:,:,jk) = va(:,:,jk) + ( va_b(:,:) - vb_b(:,:) ) * z1_2dt_b
         END DO
      ELSE
         ! At this stage, ssha has been corrected: compute new depths at velocity points
         DO jj = 1, jpjm1
            DO ji = 1, jpim1      ! NO Vector Opt.
               zsshu_a(ji,jj) = z1_2 * umask(ji,jj,1)  * r1_e1e2u(ji,jj) &
                  &              * ( e1e2t(ji  ,jj) * ssha(ji  ,jj)    &
                  &              +   e1e2t(ji+1,jj) * ssha(ji+1,jj) )
               zsshv_a(ji,jj) = z1_2 * vmask(ji,jj,1)  * r1_e1e2v(ji,jj) &
                  &              * ( e1e2t(ji,jj  ) * ssha(ji,jj  )    &
                  &              +   e1e2t(ji,jj+1) * ssha(ji,jj+1) )
            END DO
         END DO
         CALL lbc_lnk_multi( zsshu_a, 'U', 1._wp, zsshv_a, 'V', 1._wp ) ! Boundary conditions
         !
         DO jk=1,jpkm1
            ua(:,:,jk) = ua(:,:,jk) + r1_hu_n(:,:) * ( ua_b(:,:) - ub_b(:,:) * hu_b(:,:) ) * z1_2dt_b
            va(:,:,jk) = va(:,:,jk) + r1_hv_n(:,:) * ( va_b(:,:) - vb_b(:,:) * hv_b(:,:) ) * z1_2dt_b
         END DO
         ! Save barotropic velocities not transport:
         ua_b(:,:) =  ua_b(:,:) / ( hu_0(:,:) + zsshu_a(:,:) + 1._wp - ssumask(:,:) )
         va_b(:,:) =  va_b(:,:) / ( hv_0(:,:) + zsshv_a(:,:) + 1._wp - ssvmask(:,:) )
      ENDIF
      !
      DO jk = 1, jpkm1
         ! Correct velocities:
         un(:,:,jk) = ( un(:,:,jk) + un_adv(:,:) - un_b(:,:) ) * umask(:,:,jk)
         vn(:,:,jk) = ( vn(:,:,jk) + vn_adv(:,:) - vn_b(:,:) ) * vmask(:,:,jk)
         !
      END DO
      !
      CALL iom_put(  "ubar", un_adv(:,:)      )    ! barotropic i-current
      CALL iom_put(  "vbar", vn_adv(:,:)      )    ! barotropic i-current
      !
#if defined key_agrif
      ! Save time integrated fluxes during child grid integration
      ! (used to update coarse grid transports at next time step)
      !
      IF( .NOT.Agrif_Root() .AND. ln_bt_fw ) THEN
         IF( Agrif_NbStepint() == 0 ) THEN
            ub2_i_b(:,:) = 0._wp
            vb2_i_b(:,:) = 0._wp
         END IF
         !
         za1 = 1._wp / REAL(Agrif_rhot(), wp)
         ub2_i_b(:,:) = ub2_i_b(:,:) + za1 * ub2_b(:,:)
         vb2_i_b(:,:) = vb2_i_b(:,:) + za1 * vb2_b(:,:)
      ENDIF
#endif      
      !                                   !* write time-spliting arrays in the restart
      IF( lrst_oce .AND.ln_bt_fw )   CALL ts_rst( kt, 'WRITE' )
      !
      CALL wrk_dealloc( jpi,jpj,   zsshp2_e, zhdiv )
      CALL wrk_dealloc( jpi,jpj,   zu_trd, zv_trd )
      CALL wrk_dealloc( jpi,jpj,   zwx, zwy, zssh_frc, zu_frc, zv_frc )
      CALL wrk_dealloc( jpi,jpj,   zhup2_e, zhvp2_e, zhust_e, zhvst_e )
      CALL wrk_dealloc( jpi,jpj,   zsshu_a, zsshv_a                                   )
      CALL wrk_dealloc( jpi,jpj,   zhf )
      IF( ln_wd ) CALL wrk_dealloc( jpi, jpj, zcpx, zcpy )
      !
      IF ( ln_diatmb ) THEN
         CALL iom_put( "baro_u" , un_b*umask(:,:,1)+zmdi*(1-umask(:,:,1 ) ) )  ! Barotropic  U Velocity
         CALL iom_put( "baro_v" , vn_b*vmask(:,:,1)+zmdi*(1-vmask(:,:,1 ) ) )  ! Barotropic  V Velocity
      ENDIF
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg_ts')
      !
   END SUBROUTINE dyn_spg_ts


   SUBROUTINE ts_wgt( ll_av, ll_fw, jpit, zwgt1, zwgt2)
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_wgt  ***
      !!
      !! ** Purpose : Set time-splitting weights for temporal averaging (or not)
      !!----------------------------------------------------------------------
      LOGICAL, INTENT(in) ::   ll_av      ! temporal averaging=.true.
      LOGICAL, INTENT(in) ::   ll_fw      ! forward time splitting =.true.
      INTEGER, INTENT(inout) :: jpit      ! cycle length    
      REAL(wp), DIMENSION(3*nn_baro), INTENT(inout) ::   zwgt1, & ! Primary weights
                                                         zwgt2    ! Secondary weights
      
      INTEGER ::  jic, jn, ji                      ! temporary integers
      REAL(wp) :: za1, za2
      !!----------------------------------------------------------------------

      zwgt1(:) = 0._wp
      zwgt2(:) = 0._wp

      ! Set time index when averaged value is requested
      IF (ll_fw) THEN 
         jic = nn_baro
      ELSE
         jic = 2 * nn_baro
      ENDIF

      ! Set primary weights:
      IF (ll_av) THEN
           ! Define simple boxcar window for primary weights 
           ! (width = nn_baro, centered around jic)     
         SELECT CASE ( nn_bt_flt )
              CASE( 0 )  ! No averaging
                 zwgt1(jic) = 1._wp
                 jpit = jic

              CASE( 1 )  ! Boxcar, width = nn_baro
                 DO jn = 1, 3*nn_baro
                    za1 = ABS(float(jn-jic))/float(nn_baro) 
                    IF (za1 < 0.5_wp) THEN
                      zwgt1(jn) = 1._wp
                      jpit = jn
                    ENDIF
                 ENDDO

              CASE( 2 )  ! Boxcar, width = 2 * nn_baro
                 DO jn = 1, 3*nn_baro
                    za1 = ABS(float(jn-jic))/float(nn_baro) 
                    IF (za1 < 1._wp) THEN
                      zwgt1(jn) = 1._wp
                      jpit = jn
                    ENDIF
                 ENDDO
              CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt' )
         END SELECT

      ELSE ! No time averaging
         zwgt1(jic) = 1._wp
         jpit = jic
      ENDIF
    
      ! Set secondary weights
      DO jn = 1, jpit
        DO ji = jn, jpit
             zwgt2(jn) = zwgt2(jn) + zwgt1(ji)
        END DO
      END DO

      ! Normalize weigths:
      za1 = 1._wp / SUM(zwgt1(1:jpit))
      za2 = 1._wp / SUM(zwgt2(1:jpit))
      DO jn = 1, jpit
        zwgt1(jn) = zwgt1(jn) * za1
        zwgt2(jn) = zwgt2(jn) * za2
      END DO
      !
   END SUBROUTINE ts_wgt


   SUBROUTINE ts_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_rst  ***
      !!
      !! ** Purpose : Read or write time-splitting arrays in restart file
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      !
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN
         CALL iom_get( numror, jpdom_autoglo, 'ub2_b'  , ub2_b  (:,:) )   
         CALL iom_get( numror, jpdom_autoglo, 'vb2_b'  , vb2_b  (:,:) ) 
         IF( .NOT.ln_bt_av ) THEN
            CALL iom_get( numror, jpdom_autoglo, 'sshbb_e'  , sshbb_e(:,:) )   
            CALL iom_get( numror, jpdom_autoglo, 'ubb_e'    ,   ubb_e(:,:) )   
            CALL iom_get( numror, jpdom_autoglo, 'vbb_e'    ,   vbb_e(:,:) )
            CALL iom_get( numror, jpdom_autoglo, 'sshb_e'   ,  sshb_e(:,:) ) 
            CALL iom_get( numror, jpdom_autoglo, 'ub_e'     ,    ub_e(:,:) )   
            CALL iom_get( numror, jpdom_autoglo, 'vb_e'     ,    vb_e(:,:) )
         ENDIF
#if defined key_agrif
         ! Read time integrated fluxes
         IF ( .NOT.Agrif_Root() ) THEN
            CALL iom_get( numror, jpdom_autoglo, 'ub2_i_b'  , ub2_i_b(:,:) )   
            CALL iom_get( numror, jpdom_autoglo, 'vb2_i_b'  , vb2_i_b(:,:) )
         ENDIF
#endif
      !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN
         CALL iom_rstput( kt, nitrst, numrow, 'ub2_b'   , ub2_b  (:,:) )
         CALL iom_rstput( kt, nitrst, numrow, 'vb2_b'   , vb2_b  (:,:) )
         !
         IF (.NOT.ln_bt_av) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'sshbb_e'  , sshbb_e(:,:) ) 
            CALL iom_rstput( kt, nitrst, numrow, 'ubb_e'    ,   ubb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vbb_e'    ,   vbb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'sshb_e'   ,  sshb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'ub_e'     ,    ub_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vb_e'     ,    vb_e(:,:) )
         ENDIF
#if defined key_agrif
         ! Save time integrated fluxes
         IF ( .NOT.Agrif_Root() ) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'ub2_i_b'  , ub2_i_b(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vb2_i_b'  , vb2_i_b(:,:) )
         ENDIF
#endif
      ENDIF
      !
   END SUBROUTINE ts_rst


   SUBROUTINE dyn_spg_ts_init
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_spg_ts_init  ***
      !!
      !! ** Purpose : Set time splitting options
      !!----------------------------------------------------------------------
      INTEGER  ::   ji ,jj              ! dummy loop indices
      REAL(wp) ::   zxr2, zyr2, zcmax   ! local scalar
      REAL(wp), POINTER, DIMENSION(:,:) ::   zcu
      !!----------------------------------------------------------------------
      !
      ! Max courant number for ext. grav. waves
      !
      CALL wrk_alloc( jpi,jpj,   zcu )
      !
      DO jj = 1, jpj
         DO ji =1, jpi
            zxr2 = r1_e1t(ji,jj) * r1_e1t(ji,jj)
            zyr2 = r1_e2t(ji,jj) * r1_e2t(ji,jj)
            zcu(ji,jj) = SQRT( grav * MAX(ht_0(ji,jj),0._wp) * (zxr2 + zyr2) )
         END DO
      END DO
      !
      zcmax = MAXVAL( zcu(:,:) )
      IF( lk_mpp )   CALL mpp_max( zcmax )

      ! Estimate number of iterations to satisfy a max courant number= rn_bt_cmax
      IF( ln_bt_auto )   nn_baro = CEILING( rdt / rn_bt_cmax * zcmax)
      
      rdtbt = rdt / REAL( nn_baro , wp )
      zcmax = zcmax * rdtbt
							! Print results
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dyn_spg_ts : split-explicit free surface'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      IF( ln_bt_auto ) THEN
         IF(lwp) WRITE(numout,*) '     ln_ts_auto=.true. Automatically set nn_baro '
         IF(lwp) WRITE(numout,*) '     Max. courant number allowed: ', rn_bt_cmax
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_ts_auto=.false.: Use nn_baro in namelist '
      ENDIF

      IF(ln_bt_av) THEN
         IF(lwp) WRITE(numout,*) '     ln_bt_av=.true.  => Time averaging over nn_baro time steps is on '
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_bt_av=.false. => No time averaging of barotropic variables '
      ENDIF
      !
      !
      IF(ln_bt_fw) THEN
         IF(lwp) WRITE(numout,*) '     ln_bt_fw=.true.  => Forward integration of barotropic variables '
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_bt_fw =.false.=> Centred integration of barotropic variables '
      ENDIF
      !
#if defined key_agrif
      ! Restrict the use of Agrif to the forward case only
      IF( .NOT.ln_bt_fw .AND. .NOT.Agrif_Root() )   CALL ctl_stop( 'AGRIF not implemented if ln_bt_fw=.FALSE.' )
#endif
      !
      IF(lwp) WRITE(numout,*)    '     Time filter choice, nn_bt_flt: ', nn_bt_flt
      SELECT CASE ( nn_bt_flt )
         CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '           Dirac'
         CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '           Boxcar: width = nn_baro'
         CASE( 2 )      ;   IF(lwp) WRITE(numout,*) '           Boxcar: width = 2*nn_baro' 
         CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt: should 0,1,2' )
      END SELECT
      !
      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) '     nn_baro = ', nn_baro
      IF(lwp) WRITE(numout,*) '     Barotropic time step [s] is :', rdtbt
      IF(lwp) WRITE(numout,*) '     Maximum Courant number is   :', zcmax
      !
      IF( .NOT.ln_bt_av .AND. .NOT.ln_bt_fw ) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: No time averaging => only forward integration is possible' )
      ENDIF
      IF( zcmax>0.9_wp ) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: Maximum Courant number is greater than 0.9: Inc. nn_baro !' )          
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj,   zcu )
      !
   END SUBROUTINE dyn_spg_ts_init

   !!======================================================================
END MODULE dynspg_ts
