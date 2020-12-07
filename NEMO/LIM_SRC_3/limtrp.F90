MODULE limtrp
   !!======================================================================
   !!                       ***  MODULE limtrp   ***
   !! LIM transport ice model : sea-ice advection/diffusion
   !!======================================================================
   !! History : LIM-2 ! 2000-01 (M.A. Morales Maqueda, H. Goosse, and T. Fichefet)  Original code
   !!            3.0  ! 2005-11 (M. Vancoppenolle)   Multi-layer sea ice, salinity variations
   !!            4.0  ! 2011-02 (G. Madec) dynamical allocation
   !!----------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3'                                      LIM3 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_trp      : advection/diffusion process of sea ice
   !!----------------------------------------------------------------------
   USE phycst         ! physical constant
   USE dom_oce        ! ocean domain
   USE sbc_oce        ! ocean surface boundary condition
   USE ice            ! ice variables
   USE limhdf         ! ice horizontal diffusion
   USE limvar         ! 
   USE limadv_prather ! advection scheme (Prather)
   USE limadv_umx     ! advection scheme (ultimate-macho)
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! lateral boundary conditions -- MPP exchanges
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! work arrays
   USE prtctl         ! Print control
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
   USE timing         ! Timing
   USE limcons        ! conservation tests
   USE limctl         ! control prints

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_trp    ! called by sbcice_lim

   INTEGER  ::   ncfl                 ! number of ice time step with CFL>1/2  

   !! * Substitution
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/LIM3 4.0 , UCL - NEMO Consortium (2011)
   !! $Id: limtrp.F90 7753 2017-03-03 11:46:59Z mocavero $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_trp( kt ) 
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE lim_trp ***
      !!                    
      !! ** purpose : advection/diffusion process of sea ice
      !!
      !! ** method  : variables included in the process are scalar,   
      !!     other values are considered as second order. 
      !!     For advection, one can choose between
      !!     a) an Ultimate-Macho scheme (whose order is defined by nn_limadv_ord) => nn_limadv=0
      !!     b) and a second order Prather scheme => nn_limadv=-1
      !!
      !! ** action :
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! number of iteration
      !
      INTEGER  ::   ji, jj, jk, jm, jl, jt  ! dummy loop indices
      INTEGER  ::   initad                  ! number of sub-timestep for the advection
      REAL(wp) ::   zcfl , zusnit           !   -      -
      CHARACTER(len=80) :: cltmp
      !
      REAL(wp) ::    zvi_b, zsmv_b, zei_b, zfs_b, zfw_b, zft_b
      REAL(wp) ::    zdv, zda
      REAL(wp), POINTER, DIMENSION(:,:)      ::   zatold, zeiold, zesold, zsmvold 
      REAL(wp), POINTER, DIMENSION(:,:,:)    ::   zhimax, zviold, zvsold
      ! --- diffusion --- !
      REAL(wp), POINTER, DIMENSION(:,:,:)    ::   zhdfptab
      INTEGER , PARAMETER                    ::   ihdf_vars  = 6 ! Number of variables in which we apply horizontal diffusion
                                                                 !  inside limtrp for each ice category , not counting the 
                                                                 !  variables corresponding to ice_layers 
      ! --- ultimate macho only --- !
      REAL(wp)                               ::   zdt
      REAL(wp), POINTER, DIMENSION(:,:)      ::   zudy, zvdx, zcu_box, zcv_box
      ! --- prather only --- !
      REAL(wp), POINTER, DIMENSION(:,:)      ::   zarea
      REAL(wp), POINTER, DIMENSION(:,:,:)    ::   z0opw
      REAL(wp), POINTER, DIMENSION(:,:,:)    ::   z0ice, z0snw, z0ai, z0es , z0smi , z0oi
      REAL(wp), POINTER, DIMENSION(:,:,:,:)  ::   z0ei
      !!
      !!---------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('limtrp')

      CALL wrk_alloc( jpi,jpj,                            zatold, zeiold, zesold, zsmvold )
      CALL wrk_alloc( jpi,jpj,jpl,                        zhimax, zviold, zvsold )
      CALL wrk_alloc( jpi,jpj,jpl*(ihdf_vars + nlay_i)+1, zhdfptab)
 
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)''
         WRITE(numout,*)'limtrp'
         WRITE(numout,*)'~~~~~~'
         ncfl = 0                ! nb of time step with CFL > 1/2
      ENDIF
      
      CALL lim_var_agg( 1 ) ! integrated values + ato_i

      !-------------------------------------!
      !   Advection of sea ice properties   !
      !-------------------------------------!

      ! conservation test
      IF( ln_limdiachk )   CALL lim_cons_hsm(0, 'limtrp', zvi_b, zsmv_b, zei_b, zfw_b, zfs_b, zft_b)
      
      ! store old values for diag
      zviold = v_i
      zvsold = v_s
      zsmvold(:,:) = SUM( smv_i(:,:,:), dim=3 )
      zeiold (:,:) = et_i
      zesold (:,:) = et_s 

      !--- Thickness correction init. --- !
      zatold(:,:) = at_i
      DO jl = 1, jpl
         DO jj = 1, jpj
            DO ji = 1, jpi
               rswitch          = MAX( 0._wp , SIGN( 1._wp, a_i(ji,jj,jl) - epsi20 ) )
               ht_i  (ji,jj,jl) = v_i (ji,jj,jl) / MAX( a_i(ji,jj,jl) , epsi20 ) * rswitch
               ht_s  (ji,jj,jl) = v_s (ji,jj,jl) / MAX( a_i(ji,jj,jl) , epsi20 ) * rswitch
            END DO
         END DO
      END DO
      ! --- Record max of the surrounding ice thicknesses for correction in case advection creates ice too thick --- !
      zhimax(:,:,:) = ht_i(:,:,:) + ht_s(:,:,:)
      DO jl = 1, jpl
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zhimax(ji,jj,jl) = MAXVAL( ht_i(ji-1:ji+1,jj-1:jj+1,jl) + ht_s(ji-1:ji+1,jj-1:jj+1,jl) )
            END DO
         END DO
         CALL lbc_lnk(zhimax(:,:,jl),'T',1.)
      END DO
         
      ! --- If ice drift field is too fast, use an appropriate time step for advection --- !        
      zcfl  =            MAXVAL( ABS( u_ice(:,:) ) * rdt_ice * r1_e1u(:,:) )    ! CFL test for stability
      zcfl  = MAX( zcfl, MAXVAL( ABS( v_ice(:,:) ) * rdt_ice * r1_e2v(:,:) ) )
      IF( lk_mpp )   CALL mpp_max( zcfl )
      
      IF( zcfl > 0.5 ) THEN   ;   initad = 2   ;   zusnit = 0.5_wp
      ELSE                    ;   initad = 1   ;   zusnit = 1.0_wp
      ENDIF
      
!!      IF( zcfl > 0.5_wp .AND. lwp ) THEN
!!         ncfl = ncfl + 1
!!         IF( ncfl > 0 ) THEN   
!!            WRITE(cltmp,'(i6.1)') ncfl
!!            CALL ctl_warn( 'lim_trp: ncfl= ', TRIM(cltmp), 'advective ice time-step using a split in sub-time-step ')
!!         ENDIF
!!      ENDIF

      SELECT CASE ( nn_limadv )
         
                       !=============================!
      CASE ( 0 )       !==  Ultimate-MACHO scheme  ==!                   
                       !=============================!
      
         CALL wrk_alloc( jpi,jpj, zudy, zvdx, zcu_box, zcv_box )
      
         IF( kt == nit000 .AND. lwp ) THEN
            WRITE(numout,*)''
            WRITE(numout,*)'lim_adv_umx : Ultimate-MACHO advection scheme'
            WRITE(numout,*)'~~~~~~~~~~~'
         ENDIF
         !
         zdt = rdt_ice / REAL(initad)
         
         ! transport
         zudy(:,:) = u_ice(:,:) * e2u(:,:)
         zvdx(:,:) = v_ice(:,:) * e1v(:,:)
         
         ! define velocity for advection: u*grad(H)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               IF    ( u_ice(ji,jj) * u_ice(ji-1,jj) <= 0._wp ) THEN   ;   zcu_box(ji,jj) = 0._wp
               ELSEIF( u_ice(ji,jj)                  >  0._wp ) THEN   ;   zcu_box(ji,jj) = u_ice(ji-1,jj)
               ELSE                                                    ;   zcu_box(ji,jj) = u_ice(ji  ,jj)
               ENDIF
               
               IF    ( v_ice(ji,jj) * v_ice(ji,jj-1) <= 0._wp ) THEN   ;   zcv_box(ji,jj) = 0._wp
               ELSEIF( v_ice(ji,jj)                  >  0._wp ) THEN   ;   zcv_box(ji,jj) = v_ice(ji,jj-1)
               ELSE                                                    ;   zcv_box(ji,jj) = v_ice(ji,jj  )
               ENDIF
            END DO
         END DO
         
         ! advection
         DO jt = 1, initad
            CALL lim_adv_umx( kt, zdt, zudy, zvdx, zcu_box, zcv_box, ato_i(:,:) )       ! Open water area 
            DO jl = 1, jpl
               CALL lim_adv_umx( kt, zdt, zudy, zvdx, zcu_box, zcv_box, a_i(:,:,jl) )      ! Ice area
               CALL lim_adv_umx( kt, zdt, zudy, zvdx, zcu_box, zcv_box, v_i(:,:,jl) )      ! Ice  volume
               CALL lim_adv_umx( kt, zdt, zudy, zvdx, zcu_box, zcv_box, smv_i(:,:,jl) )    ! Salt content
               CALL lim_adv_umx( kt, zdt, zudy, zvdx, zcu_box, zcv_box, oa_i (:,:,jl) )    ! Age content
               DO jk = 1, nlay_i
                  CALL lim_adv_umx( kt, zdt, zudy, zvdx, zcu_box, zcv_box, e_i(:,:,jk,jl) )   ! Ice  heat content
               END DO
               CALL lim_adv_umx( kt, zdt, zudy, zvdx, zcu_box, zcv_box, v_s(:,:,jl) )      ! Snow volume
               CALL lim_adv_umx( kt, zdt, zudy, zvdx, zcu_box, zcv_box, e_s(:,:,1,jl) )    ! Snow heat content
            END DO
         END DO
         !
         at_i(:,:) = a_i(:,:,1)      ! total ice fraction
         DO jl = 2, jpl
            at_i(:,:) = at_i(:,:) + a_i(:,:,jl)
         END DO
         !
         CALL wrk_dealloc( jpi,jpj, zudy, zvdx, zcu_box, zcv_box )
         
                       !=============================!
      CASE ( -1 )      !==     Prather scheme      ==!                   
                       !=============================!

         CALL wrk_alloc( jpi,jpj,            zarea )
         CALL wrk_alloc( jpi,jpj,1,          z0opw )
         CALL wrk_alloc( jpi,jpj,jpl,        z0ice, z0snw, z0ai, z0es , z0smi , z0oi )
         CALL wrk_alloc( jpi,jpj,nlay_i,jpl, z0ei )
         
         IF( kt == nit000 .AND. lwp ) THEN
            WRITE(numout,*)''
            WRITE(numout,*)'lim_adv_xy : Prather advection scheme'
            WRITE(numout,*)'~~~~~~~~~~~'
         ENDIF
         
         zarea(:,:) = e1e2t(:,:)
         
         !-------------------------
         ! transported fields                                        
         !-------------------------
         z0opw(:,:,1) = ato_i(:,:) * e1e2t(:,:)             ! Open water area 
         DO jl = 1, jpl
            z0snw (:,:,jl)  = v_s  (:,:,  jl) * e1e2t(:,:)  ! Snow volume
            z0ice(:,:,jl)   = v_i  (:,:,  jl) * e1e2t(:,:)  ! Ice  volume
            z0ai  (:,:,jl)  = a_i  (:,:,  jl) * e1e2t(:,:)  ! Ice area
            z0smi (:,:,jl)  = smv_i(:,:,  jl) * e1e2t(:,:)  ! Salt content
            z0oi (:,:,jl)   = oa_i (:,:,  jl) * e1e2t(:,:)  ! Age content
            z0es (:,:,jl)   = e_s  (:,:,1,jl) * e1e2t(:,:)  ! Snow heat content
            DO jk = 1, nlay_i
               z0ei  (:,:,jk,jl) = e_i  (:,:,jk,jl) * e1e2t(:,:) ! Ice  heat content
            END DO
         END DO


         IF( MOD( ( kt - 1) / nn_fsbc , 2 ) == 0 ) THEN       !==  odd ice time step:  adv_x then adv_y  ==!
            DO jt = 1, initad
               CALL lim_adv_x( zusnit, u_ice, 1._wp, zarea, z0opw (:,:,1), sxopw(:,:),   &             !--- ice open water area
                  &                                         sxxopw(:,:)  , syopw(:,:), syyopw(:,:), sxyopw(:,:)  )
               CALL lim_adv_y( zusnit, v_ice, 0._wp, zarea, z0opw (:,:,1), sxopw(:,:),   &
                  &                                         sxxopw(:,:)  , syopw(:,:), syyopw(:,:), sxyopw(:,:)  )
               DO jl = 1, jpl
                  CALL lim_adv_x( zusnit, u_ice, 1._wp, zarea, z0ice (:,:,jl), sxice(:,:,jl),   &    !--- ice volume  ---
                     &                                         sxxice(:,:,jl), syice(:,:,jl), syyice(:,:,jl), sxyice(:,:,jl)  )
                  CALL lim_adv_y( zusnit, v_ice, 0._wp, zarea, z0ice (:,:,jl), sxice(:,:,jl),   &
                     &                                         sxxice(:,:,jl), syice(:,:,jl), syyice(:,:,jl), sxyice(:,:,jl)  )
                  CALL lim_adv_x( zusnit, u_ice, 1._wp, zarea, z0snw (:,:,jl), sxsn (:,:,jl),   &    !--- snow volume  ---
                     &                                         sxxsn (:,:,jl), sysn (:,:,jl), syysn (:,:,jl), sxysn (:,:,jl)  )
                  CALL lim_adv_y( zusnit, v_ice, 0._wp, zarea, z0snw (:,:,jl), sxsn (:,:,jl),   &
                     &                                         sxxsn (:,:,jl), sysn (:,:,jl), syysn (:,:,jl), sxysn (:,:,jl)  )
                  CALL lim_adv_x( zusnit, u_ice, 1._wp, zarea, z0smi (:,:,jl), sxsal(:,:,jl),   &    !--- ice salinity ---
                     &                                         sxxsal(:,:,jl), sysal(:,:,jl), syysal(:,:,jl), sxysal(:,:,jl)  )
                  CALL lim_adv_y( zusnit, v_ice, 0._wp, zarea, z0smi (:,:,jl), sxsal(:,:,jl),   &
                     &                                         sxxsal(:,:,jl), sysal(:,:,jl), syysal(:,:,jl), sxysal(:,:,jl)  )
                  CALL lim_adv_x( zusnit, u_ice, 1._wp, zarea, z0oi  (:,:,jl), sxage(:,:,jl),   &    !--- ice age      ---     
                     &                                         sxxage(:,:,jl), syage(:,:,jl), syyage(:,:,jl), sxyage(:,:,jl)  )
                  CALL lim_adv_y( zusnit, v_ice, 0._wp, zarea, z0oi  (:,:,jl), sxage(:,:,jl),   &
                     &                                         sxxage(:,:,jl), syage(:,:,jl), syyage(:,:,jl), sxyage(:,:,jl)  )
                  CALL lim_adv_x( zusnit, u_ice, 1._wp, zarea, z0ai  (:,:,jl), sxa  (:,:,jl),   &    !--- ice concentrations ---
                     &                                         sxxa  (:,:,jl), sya  (:,:,jl), syya  (:,:,jl), sxya  (:,:,jl)  )
                  CALL lim_adv_y( zusnit, v_ice, 0._wp, zarea, z0ai  (:,:,jl), sxa  (:,:,jl),   & 
                     &                                         sxxa  (:,:,jl), sya  (:,:,jl), syya  (:,:,jl), sxya  (:,:,jl)  )
                  CALL lim_adv_x( zusnit, u_ice, 1._wp, zarea, z0es  (:,:,jl), sxc0 (:,:,jl),   &    !--- snow heat contents ---
                     &                                         sxxc0 (:,:,jl), syc0 (:,:,jl), syyc0 (:,:,jl), sxyc0 (:,:,jl)  )
                  CALL lim_adv_y( zusnit, v_ice, 0._wp, zarea, z0es  (:,:,jl), sxc0 (:,:,jl),   &
                     &                                         sxxc0 (:,:,jl), syc0 (:,:,jl), syyc0 (:,:,jl), sxyc0 (:,:,jl)  )
                  DO jk = 1, nlay_i                                                                !--- ice heat contents ---
                     CALL lim_adv_x( zusnit, u_ice, 1._wp, zarea, z0ei(:,:,jk,jl), sxe (:,:,jk,jl),   & 
                        &                                         sxxe(:,:,jk,jl), sye (:,:,jk,jl),   &
                        &                                         syye(:,:,jk,jl), sxye(:,:,jk,jl) )
                     CALL lim_adv_y( zusnit, v_ice, 0._wp, zarea, z0ei(:,:,jk,jl), sxe (:,:,jk,jl),   & 
                        &                                         sxxe(:,:,jk,jl), sye (:,:,jk,jl),   &
                        &                                         syye(:,:,jk,jl), sxye(:,:,jk,jl) )
                  END DO
               END DO
            END DO
         ELSE
            DO jt = 1, initad
               CALL lim_adv_y( zusnit, v_ice, 1._wp, zarea, z0opw (:,:,1), sxopw(:,:),   &             !--- ice open water area
                  &                                         sxxopw(:,:)  , syopw(:,:), syyopw(:,:), sxyopw(:,:)  )
               CALL lim_adv_x( zusnit, u_ice, 0._wp, zarea, z0opw (:,:,1), sxopw(:,:),   &
                  &                                         sxxopw(:,:)  , syopw(:,:), syyopw(:,:), sxyopw(:,:)  )
               DO jl = 1, jpl
                  CALL lim_adv_y( zusnit, v_ice, 1._wp, zarea, z0ice (:,:,jl), sxice(:,:,jl),   &    !--- ice volume  ---
                     &                                         sxxice(:,:,jl), syice(:,:,jl), syyice(:,:,jl), sxyice(:,:,jl)  )
                  CALL lim_adv_x( zusnit, u_ice, 0._wp, zarea, z0ice (:,:,jl), sxice(:,:,jl),   &
                     &                                         sxxice(:,:,jl), syice(:,:,jl), syyice(:,:,jl), sxyice(:,:,jl)  )
                  CALL lim_adv_y( zusnit, v_ice, 1._wp, zarea, z0snw (:,:,jl), sxsn (:,:,jl),   &    !--- snow volume  ---
                     &                                         sxxsn (:,:,jl), sysn (:,:,jl), syysn (:,:,jl), sxysn (:,:,jl)  )
                  CALL lim_adv_x( zusnit, u_ice, 0._wp, zarea, z0snw (:,:,jl), sxsn (:,:,jl),   &
                     &                                         sxxsn (:,:,jl), sysn (:,:,jl), syysn (:,:,jl), sxysn (:,:,jl)  )
                  CALL lim_adv_y( zusnit, v_ice, 1._wp, zarea, z0smi (:,:,jl), sxsal(:,:,jl),   &    !--- ice salinity ---
                     &                                         sxxsal(:,:,jl), sysal(:,:,jl), syysal(:,:,jl), sxysal(:,:,jl)  )
                  CALL lim_adv_x( zusnit, u_ice, 0._wp, zarea, z0smi (:,:,jl), sxsal(:,:,jl),   &
                     &                                         sxxsal(:,:,jl), sysal(:,:,jl), syysal(:,:,jl), sxysal(:,:,jl)  )
                  CALL lim_adv_y( zusnit, v_ice, 1._wp, zarea, z0oi  (:,:,jl), sxage(:,:,jl),   &   !--- ice age      ---
                     &                                         sxxage(:,:,jl), syage(:,:,jl), syyage(:,:,jl), sxyage(:,:,jl)  )
                  CALL lim_adv_x( zusnit, u_ice, 0._wp, zarea, z0oi  (:,:,jl), sxage(:,:,jl),   &
                     &                                         sxxage(:,:,jl), syage(:,:,jl), syyage(:,:,jl), sxyage(:,:,jl)  )
                  CALL lim_adv_y( zusnit, v_ice, 1._wp, zarea, z0ai  (:,:,jl), sxa  (:,:,jl),   &   !--- ice concentrations ---
                     &                                         sxxa  (:,:,jl), sya  (:,:,jl), syya  (:,:,jl), sxya  (:,:,jl)  )
                  CALL lim_adv_x( zusnit, u_ice, 0._wp, zarea, z0ai  (:,:,jl), sxa  (:,:,jl),   &
                     &                                         sxxa  (:,:,jl), sya  (:,:,jl), syya  (:,:,jl), sxya  (:,:,jl)  )
                  CALL lim_adv_y( zusnit, v_ice, 1._wp, zarea, z0es  (:,:,jl), sxc0 (:,:,jl),   &  !--- snow heat contents ---
                     &                                         sxxc0 (:,:,jl), syc0 (:,:,jl), syyc0 (:,:,jl), sxyc0 (:,:,jl)  )
                  CALL lim_adv_x( zusnit, u_ice, 0._wp, zarea, z0es  (:,:,jl), sxc0 (:,:,jl),   &
                     &                                         sxxc0 (:,:,jl), syc0 (:,:,jl), syyc0 (:,:,jl), sxyc0 (:,:,jl)  )
                  DO jk = 1, nlay_i                                                           !--- ice heat contents ---
                     CALL lim_adv_y( zusnit, v_ice, 1._wp, zarea, z0ei(:,:,jk,jl), sxe (:,:,jk,jl),   & 
                        &                                         sxxe(:,:,jk,jl), sye (:,:,jk,jl),   &
                        &                                         syye(:,:,jk,jl), sxye(:,:,jk,jl) )
                     CALL lim_adv_x( zusnit, u_ice, 0._wp, zarea, z0ei(:,:,jk,jl), sxe (:,:,jk,jl),   & 
                        &                                         sxxe(:,:,jk,jl), sye (:,:,jk,jl),   &
                        &                                         syye(:,:,jk,jl), sxye(:,:,jk,jl) )
                  END DO
               END DO
            END DO
         ENDIF

         !-------------------------------------------
         ! Recover the properties from their contents
         !-------------------------------------------
         ato_i(:,:) = z0opw(:,:,1) * r1_e1e2t(:,:)
         DO jl = 1, jpl
            v_i  (:,:,  jl) = z0ice(:,:,jl) * r1_e1e2t(:,:)
            v_s  (:,:,  jl) = z0snw(:,:,jl) * r1_e1e2t(:,:)
            smv_i(:,:,  jl) = z0smi(:,:,jl) * r1_e1e2t(:,:)
            oa_i (:,:,  jl) = z0oi (:,:,jl) * r1_e1e2t(:,:)
            a_i  (:,:,  jl) = z0ai (:,:,jl) * r1_e1e2t(:,:)
            e_s  (:,:,1,jl) = z0es (:,:,jl) * r1_e1e2t(:,:)
            DO jk = 1, nlay_i
               e_i(:,:,jk,jl) = z0ei(:,:,jk,jl) * r1_e1e2t(:,:)
            END DO
         END DO

         at_i(:,:) = a_i(:,:,1)      ! total ice fraction
         DO jl = 2, jpl
            at_i(:,:) = at_i(:,:) + a_i(:,:,jl)
         END DO
         
         CALL wrk_dealloc( jpi,jpj,            zarea )
         CALL wrk_dealloc( jpi,jpj,1,          z0opw )
         CALL wrk_dealloc( jpi,jpj,jpl,        z0ice, z0snw, z0ai, z0es , z0smi , z0oi )
         CALL wrk_dealloc( jpi,jpj,nlay_i,jpl, z0ei )

      END SELECT
      
      !------------------------------!
      ! Diffusion of Ice fields                  
      !------------------------------!
      IF( nn_ahi0 /= -1 .AND. nn_limdyn == 2 ) THEN
         !
         ! --- Prepare diffusion for variables with categories --- !
         !     mask eddy diffusivity coefficient at ocean U- and V-points
         jm=1
         DO jl = 1, jpl
            DO jj = 1, jpjm1                 ! NB: has not to be defined on jpj line and jpi row
               DO ji = 1 , fs_jpim1
                  pahu3D(ji,jj,jl) = ( 1._wp - MAX( 0._wp, SIGN( 1._wp, -a_i(ji  ,jj,  jl ) ) ) )   &
                  &                * ( 1._wp - MAX( 0._wp, SIGN( 1._wp, -a_i(ji+1,jj,  jl ) ) ) ) * ahiu(ji,jj)
                  pahv3D(ji,jj,jl) = ( 1._wp - MAX( 0._wp, SIGN( 1._wp, -a_i(ji,  jj,  jl ) ) ) )   &
                  &                * ( 1._wp - MAX( 0._wp, SIGN( 1._wp,- a_i(ji,  jj+1,jl ) ) ) ) * ahiv(ji,jj)
               END DO
            END DO

            zhdfptab(:,:,jm)= a_i  (:,:,  jl); jm = jm + 1
            zhdfptab(:,:,jm)= v_i  (:,:,  jl); jm = jm + 1
            zhdfptab(:,:,jm)= v_s  (:,:,  jl); jm = jm + 1
            zhdfptab(:,:,jm)= smv_i(:,:,  jl); jm = jm + 1
            zhdfptab(:,:,jm)= oa_i (:,:,  jl); jm = jm + 1
            zhdfptab(:,:,jm)= e_s  (:,:,1,jl); jm = jm + 1
            ! Sample of adding more variables to apply lim_hdf (ihdf_vars must be increased)
            !   zhdfptab(:,:,jm) = variable_1 (:,:,1,jl); jm = jm + 1  
            !   zhdfptab(:,:,jm) = variable_2 (:,:,1,jl); jm = jm + 1 
            DO jk = 1, nlay_i
              zhdfptab(:,:,jm)=e_i(:,:,jk,jl); jm= jm+1
            END DO
         END DO

         ! --- Prepare diffusion for open water area --- !
         !     mask eddy diffusivity coefficient at ocean U- and V-points
         DO jj = 1, jpjm1                    ! NB: has not to be defined on jpj line and jpi row
            DO ji = 1 , fs_jpim1
               pahu3D(ji,jj,jpl+1) = ( 1._wp - MAX( 0._wp, SIGN( 1._wp, -at_i(ji  ,jj) ) ) )   &
                  &                * ( 1._wp - MAX( 0._wp, SIGN( 1._wp, -at_i(ji+1,jj) ) ) ) * ahiu(ji,jj)
               pahv3D(ji,jj,jpl+1) = ( 1._wp - MAX( 0._wp, SIGN( 1._wp, -at_i(ji,jj  ) ) ) )   &
                  &                * ( 1._wp - MAX( 0._wp, SIGN( 1._wp,- at_i(ji,jj+1) ) ) ) * ahiv(ji,jj)
            END DO
         END DO
         !
         zhdfptab(:,:,jm)= ato_i  (:,:);

         ! --- Apply diffusion --- !
         CALL lim_hdf( zhdfptab, ihdf_vars )

         ! --- Recover properties --- !
         jm=1
         DO jl = 1, jpl
            a_i  (:,:,  jl) = zhdfptab(:,:,jm); jm = jm + 1
            v_i  (:,:,  jl) = zhdfptab(:,:,jm); jm = jm + 1
            v_s  (:,:,  jl) = zhdfptab(:,:,jm); jm = jm + 1
            smv_i(:,:,  jl) = zhdfptab(:,:,jm); jm = jm + 1
            oa_i (:,:,  jl) = zhdfptab(:,:,jm); jm = jm + 1
            e_s  (:,:,1,jl) = zhdfptab(:,:,jm); jm = jm + 1
            ! Sample of adding more variables to apply lim_hdf
            !   variable_1  (:,:,1,jl) = zhdfptab(:,:, jm  ) ; jm + 1 
            !   variable_2  (:,:,1,jl) = zhdfptab(:,:, jm  ) ; jm + 1
            DO jk = 1, nlay_i
               e_i(:,:,jk,jl) = zhdfptab(:,:,jm);jm= jm + 1
            END DO
         END DO
         ato_i  (:,:) = zhdfptab(:,:,jm)
              
      ENDIF

      ! --- diags ---
      DO jj = 1, jpj
         DO ji = 1, jpi
            diag_trp_ei (ji,jj) = ( SUM( e_i  (ji,jj,1:nlay_i,:) ) -  zeiold(ji,jj) ) * r1_rdtice
            diag_trp_es (ji,jj) = ( SUM( e_s  (ji,jj,1:nlay_s,:) ) -  zesold(ji,jj) ) * r1_rdtice
            diag_trp_smv(ji,jj) = ( SUM( smv_i(ji,jj,:)          ) - zsmvold(ji,jj) ) * r1_rdtice
            diag_trp_vi (ji,jj) =   SUM(   v_i(ji,jj,:)            -  zviold(ji,jj,:) ) * r1_rdtice
            diag_trp_vs (ji,jj) =   SUM(   v_s(ji,jj,:)            -  zvsold(ji,jj,:) ) * r1_rdtice
         END DO
      END DO
      
      IF( nn_limdyn == 2) THEN

         ! zap small areas
         CALL lim_var_zapsmall
           
         !--- Thickness correction in case too high --- !
         DO jl = 1, jpl
            DO jj = 1, jpj
               DO ji = 1, jpi
                  
                  IF ( v_i(ji,jj,jl) > 0._wp ) THEN
                     
                     rswitch          = MAX( 0._wp , SIGN( 1._wp, a_i(ji,jj,jl) - epsi20 ) )
                     ht_i  (ji,jj,jl) = v_i (ji,jj,jl) / MAX( a_i(ji,jj,jl) , epsi20 ) * rswitch
                     ht_s  (ji,jj,jl) = v_s (ji,jj,jl) / MAX( a_i(ji,jj,jl) , epsi20 ) * rswitch
                     
                     zdv  = v_i(ji,jj,jl) + v_s(ji,jj,jl) - zviold(ji,jj,jl) - zvsold(ji,jj,jl)  
                     
                     IF ( ( zdv >  0.0 .AND. (ht_i(ji,jj,jl)+ht_s(ji,jj,jl)) > zhimax(ji,jj,jl) .AND. zatold(ji,jj) < 0.80 ) .OR. &
                        & ( zdv <= 0.0 .AND. (ht_i(ji,jj,jl)+ht_s(ji,jj,jl)) > zhimax(ji,jj,jl) ) ) THEN
                        
                        rswitch        = MAX( 0._wp, SIGN( 1._wp, zhimax(ji,jj,jl) - epsi20 ) )
                        a_i(ji,jj,jl)  = rswitch * ( v_i(ji,jj,jl) + v_s(ji,jj,jl) ) / MAX( zhimax(ji,jj,jl), epsi20 )
                        
                        ! small correction due to *rswitch for a_i
                        v_i  (ji,jj,jl)        = rswitch * v_i  (ji,jj,jl)
                        v_s  (ji,jj,jl)        = rswitch * v_s  (ji,jj,jl)
                        smv_i(ji,jj,jl)        = rswitch * smv_i(ji,jj,jl)
                        e_s(ji,jj,1,jl)        = rswitch * e_s(ji,jj,1,jl)
                        e_i(ji,jj,1:nlay_i,jl) = rswitch * e_i(ji,jj,1:nlay_i,jl)
                                                
                     ENDIF
                     
                  ENDIF
                
               END DO
            END DO
         END DO
         ! -------------------------------------------------
         
         ! Force the upper limit of ht_i to always be < hi_max (99 m).
         DO jj = 1, jpj
            DO ji = 1, jpi
               rswitch         = MAX( 0._wp , SIGN( 1._wp, ht_i(ji,jj,jpl) - epsi20 ) )
               ht_i(ji,jj,jpl) = MIN( ht_i(ji,jj,jpl) , hi_max(jpl) )
               a_i (ji,jj,jpl) = v_i(ji,jj,jpl) / MAX( ht_i(ji,jj,jpl) , epsi20 ) * rswitch
            END DO
         END DO

      ENDIF
         
      !------------------------------------------------------------
      ! Impose a_i < amax if no ridging/rafting or in mono-category
      !------------------------------------------------------------
      !
      at_i(:,:) = SUM( a_i(:,:,:), dim=3 )
      IF ( nn_limdyn == 1 .OR. ( ( nn_monocat == 2 ) .AND. ( jpl == 1 ) ) ) THEN ! simple conservative piling, comparable with LIM2
         DO jl = 1, jpl
            DO jj = 1, jpj
               DO ji = 1, jpi
                  rswitch       = MAX( 0._wp, SIGN( 1._wp, at_i(ji,jj) - epsi20 ) )
                  zda           = rswitch * MIN( rn_amax_2d(ji,jj) - at_i(ji,jj), 0._wp )  &
                     &                    * a_i(ji,jj,jl) / MAX( at_i(ji,jj), epsi20 )
                  a_i(ji,jj,jl) = a_i(ji,jj,jl) + zda
               END DO
            END DO
         END DO
      ENDIF
      
      ! --- agglomerate variables -----------------
      vt_i(:,:) = SUM( v_i(:,:,:), dim=3 )
      vt_s(:,:) = SUM( v_s(:,:,:), dim=3 )
      at_i(:,:) = SUM( a_i(:,:,:), dim=3 )
      
      ! --- open water = 1 if at_i=0 --------------------------------
      WHERE( at_i == 0._wp ) ato_i = 1._wp 
      
      ! conservation test
      IF( ln_limdiachk )   CALL lim_cons_hsm(1, 'limtrp', zvi_b, zsmv_b, zei_b, zfw_b, zfs_b, zft_b)
        
      ! -------------------------------------------------
      ! control prints
      ! -------------------------------------------------
      IF( ln_limctl )   CALL lim_prt( kt, iiceprt, jiceprt,-1, ' - ice dyn & trp - ' )
      !
      CALL wrk_dealloc( jpi,jpj,                            zatold, zeiold, zesold, zsmvold )
      CALL wrk_dealloc( jpi,jpj,jpl,                        zhimax, zviold, zvsold )
      CALL wrk_dealloc( jpi,jpj,jpl*(ihdf_vars + nlay_i)+1, zhdfptab)
      !
      IF( nn_timing == 1 )  CALL timing_stop('limtrp')
      !
   END SUBROUTINE lim_trp

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty Module                No sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_trp        ! Empty routine
   END SUBROUTINE lim_trp
#endif
   !!======================================================================
END MODULE limtrp

