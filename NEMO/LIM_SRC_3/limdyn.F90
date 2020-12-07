MODULE limdyn
   !!======================================================================
   !!                     ***  MODULE  limdyn  ***
   !!   Sea-Ice dynamics :  
   !!======================================================================
   !! history :  1.0  ! 2002-08  (C. Ethe, G. Madec)  original VP code 
   !!            3.0  ! 2007-03  (MA Morales Maqueda, S. Bouillon, M. Vancoppenolle)  LIM3: EVP-Cgrid
   !!            3.5  ! 2011-02  (G. Madec) dynamical allocation
   !!----------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3' :                                 LIM3 sea-ice model
   !!----------------------------------------------------------------------
   !!    lim_dyn      : computes ice velocities
   !!    lim_dyn_init : initialization and namelist read
   !!----------------------------------------------------------------------
   USE phycst           ! physical constants
   USE dom_oce          ! ocean space and time domain
   USE sbc_ice          ! Surface boundary condition: ice   fields
   USE ice              ! LIM-3 variables
   USE limrhg           ! LIM-3 rheology
   USE lbclnk           ! lateral boundary conditions - MPP exchanges
   USE lib_mpp          ! MPP library
   USE wrk_nemo         ! work arrays
   USE in_out_manager   ! I/O manager
   USE lib_fortran      ! glob_sum
   USE timing           ! Timing
   USE limcons          ! conservation tests
   USE limctl           ! control prints
   USE limvar

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_dyn        ! routine called by sbcice_lim.F90
   PUBLIC   lim_dyn_init   ! routine called by sbcice_lim.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/LIM3 4.0 , UCL - NEMO Consortium (2011)
   !! $Id: limdyn.F90 7753 2017-03-03 11:46:59Z mocavero $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_dyn( kt )
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE lim_dyn  ***
      !!               
      !! ** Purpose :   compute ice velocity
      !!                
      !! ** Method  : 
      !!
      !! ** Action  : - Initialisation
      !!              - Call of the dynamic routine for each hemisphere
      !!------------------------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! number of iteration
      !!
      INTEGER  :: jl, jk ! dummy loop indices
      REAL(wp) :: zvi_b, zsmv_b, zei_b, zfs_b, zfw_b, zft_b 
     !!---------------------------------------------------------------------

      IF( nn_timing == 1 )  CALL timing_start('limdyn')

      CALL lim_var_agg(1)                      ! aggregate ice categories
      !
      ! conservation test
      IF( ln_limdiachk ) CALL lim_cons_hsm(0, 'limdyn', zvi_b, zsmv_b, zei_b, zfw_b, zfs_b, zft_b)
      
      ! ice velocities before rheology
      u_ice_b(:,:) = u_ice(:,:) * umask(:,:,1)
      v_ice_b(:,:) = v_ice(:,:) * vmask(:,:,1)
      
      ! Landfast ice parameterization: define max bottom friction
      tau_icebfr(:,:) = 0._wp
      IF( ln_landfast ) THEN
         DO jl = 1, jpl
            WHERE( ht_i(:,:,jl) > ht_n(:,:) * rn_gamma )  tau_icebfr(:,:) = tau_icebfr(:,:) + a_i(:,:,jl) * rn_icebfr
         END DO
      ENDIF
      
      ! Rheology (ice dynamics)
      ! ========     
      CALL lim_rhg
      !
      ! conservation test
      IF( ln_limdiachk ) CALL lim_cons_hsm(1, 'limdyn', zvi_b, zsmv_b, zei_b, zfw_b, zfs_b, zft_b)

      ! Control prints
      IF( ln_ctl )       CALL lim_prt3D( 'limdyn' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('limdyn')

   END SUBROUTINE lim_dyn


   SUBROUTINE lim_dyn_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE lim_dyn_init  ***
      !!
      !! ** Purpose : Physical constants and parameters linked to the ice
      !!      dynamics
      !!
      !! ** Method  :  Read the namicedyn namelist and check the ice-dynamic
      !!       parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namicedyn
      !!-------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      NAMELIST/namicedyn/ nn_limadv, nn_limadv_ord,  &
         &                nn_icestr, ln_icestr_bvf, rn_pe_rdg, rn_pstar, rn_crhg, rn_cio, rn_creepl, rn_ecc, &
         &                nn_nevp, rn_relast, ln_landfast, rn_gamma, rn_icebfr, rn_lfrelax
      !!-------------------------------------------------------------------

      REWIND( numnam_ice_ref )              ! Namelist namicedyn in reference namelist : Ice dynamics
      READ  ( numnam_ice_ref, namicedyn, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicedyn in reference namelist', lwp )

      REWIND( numnam_ice_cfg )              ! Namelist namicedyn in configuration namelist : Ice dynamics
      READ  ( numnam_ice_cfg, namicedyn, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicedyn in configuration namelist', lwp )
      IF(lwm) WRITE ( numoni, namicedyn )
      
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'lim_dyn_init : ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~'
         ! limtrp
         WRITE(numout,*)'    choose the advection scheme (-1=Prather, 0=Ulimate-Macho)   nn_limadv     = ', nn_limadv 
         WRITE(numout,*)'    choose the order of the scheme (if ultimate)                nn_limadv_ord = ', nn_limadv_ord  
         ! limrhg
         WRITE(numout,*)'    ice strength parameterization (0=Hibler 1=Rothrock)         nn_icestr     = ', nn_icestr 
         WRITE(numout,*)'    Including brine volume in ice strength comp.                ln_icestr_bvf = ', ln_icestr_bvf
         WRITE(numout,*)'    Ratio of ridging work to PotEner change in ridging          rn_pe_rdg     = ', rn_pe_rdg 
         WRITE(numout,*) '   drag coefficient for oceanic stress                         rn_cio        = ', rn_cio
         WRITE(numout,*) '   first bulk-rheology parameter                               rn_pstar      = ', rn_pstar
         WRITE(numout,*) '   second bulk-rhelogy parameter                               rn_crhg       = ', rn_crhg
         WRITE(numout,*) '   creep limit                                                 rn_creepl     = ', rn_creepl
         WRITE(numout,*) '   eccentricity of the elliptical yield curve                  rn_ecc        = ', rn_ecc
         WRITE(numout,*) '   number of iterations for subcycling                         nn_nevp       = ', nn_nevp
         WRITE(numout,*) '   ratio of elastic timescale over ice time step               rn_relast     = ', rn_relast
         WRITE(numout,*) '   Landfast: param (T or F)                                    ln_landfast   = ', ln_landfast
         WRITE(numout,*) '   Landfast: fraction of ocean depth that ice must reach       rn_gamma      = ', rn_gamma
         WRITE(numout,*) '   Landfast: maximum bottom stress per unit area of contact    rn_icebfr     = ', rn_icebfr
         WRITE(numout,*) '   Landfast: relax time scale (s-1) to reach static friction   rn_lfrelax    = ', rn_lfrelax
      ENDIF
      !
   END SUBROUTINE lim_dyn_init

#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module           NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_dyn         ! Empty routine
   END SUBROUTINE lim_dyn
#endif 

   !!======================================================================
END MODULE limdyn
