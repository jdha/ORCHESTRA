MODULE zdfqiao
   !!======================================================================
   !!                       ***  MODULE  zdfqiao  ***
   !! Qiao module      : vertical mixing enhancement due to surface waves 
   !!======================================================================
   !! History :  3.6  !  2014-10  (E. Clementi)  Original code
   !!----------------------------------------------------------------------
   !!   zdf_qiao        : compute Qiao parameters
   !!----------------------------------------------------------------------

   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE zdf_oce
   USE sbcwave         ! wave module
   USE dom_oce
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)  
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC zdf_qiao    ! routine called in step

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: qbv, qbvu, qbvv

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011) 
   !! $Id: $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE zdf_qiao( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE zdf_qiao ***
      !!
      !! ** Purpose :Compute the Qiao term (qbv) to be added to
      !!             vertical viscosity and diffusivity coeffs.  
      !!
      !! ** Method  :qbv = alpha * A * Us(0) * exp (3 * k * z)
      !!             
      !! ** action  :Compute the Qiao wave dependent term 
      !!             only if ln_wave=.true.
      !!               
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in  ) ::  kt   ! ocean time step
      !
      INTEGER :: jj, ji, jk   ! dummy loop indices
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN                   ! First call kt=nit000 !
         IF( .NOT. ( ln_wave .AND. ln_sdw ) )   &
            &   CALL ctl_stop ( 'Ask for wave Qiao enhanced turbulence but ln_wave   &
            &                    and ln_sdw have to be activated')
         IF( zdf_qiao_alloc() /= 0 )   &
            &   CALL ctl_stop( 'STOP', 'zdf_qiao : unable to allocate arrays' )
      ENDIF

      !
      ! Compute the Qiao term Bv (qbv) to be added to
      ! vertical viscosity and diffusivity
      ! qbv = alpha * A * Us(0) * exp (3 * k * z)
      ! alpha here is set to 1
      !---------------------------------------------------------------------------------
      !
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               qbv(ji,jj,jk) = 1.0 * 0.353553 * hsw(ji,jj) * tsd2d(ji,jj) *             &
            &                  EXP(3.0 * wnum(ji,jj) *                                  &                     
            &                  (-MIN( gdepw_n(ji  ,jj  ,jk), gdepw_n(ji+1,jj  ,jk),     &
            &                         gdepw_n(ji  ,jj+1,jk), gdepw_n(ji+1,jj+1,jk))))   &
            &                          * wmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      CALL lbc_lnk( qbv, 'W', 1. )   ! Lateral boundary conditions
         
      !
      ! Interpolate Qiao parameter qbv into the grid_U and grid_V
      !----------------------------------------------------------
      !
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               qbvu(ji,jj,jk) = 0.5 * wumask(ji,jj,jk)  *              &  
            &                  ( qbv(ji,jj,jk) + qbv(ji+1,jj  ,jk) )
               qbvv(ji,jj,jk) = 0.5 * wvmask(ji,jj,jk)  *              &
            &                  ( qbv(ji,jj,jk) + qbv(ji  ,jj+1,jk) )
            END DO
         END DO
      END DO
      ! 
      CALL lbc_lnk( qbvu, 'U', 1. ) ; CALL lbc_lnk( qbvv, 'V', 1. )   ! Lateral boundary conditions

      ! Enhance vertical mixing coeff.         
      !-------------------------------
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               avmu(ji,jj,jk) = ( avmu(ji,jj,jk) + qbvu(ji,jj,jk) ) * umask(ji,jj,jk)
               avmv(ji,jj,jk) = ( avmv(ji,jj,jk) + qbvv(ji,jj,jk) ) * vmask(ji,jj,jk)
               avt (ji,jj,jk) = ( avt (ji,jj,jk) + qbv (ji,jj,jk) ) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE zdf_qiao

   INTEGER FUNCTION zdf_qiao_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION zdf_qiao_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( qbv(jpi,jpj,jpk), qbvu(jpi,jpj,jpk), qbvv(jpi,jpj,jpk),   &
         &      STAT = zdf_qiao_alloc )
      !
      IF( lk_mpp             )  CALL mpp_sum ( zdf_qiao_alloc )
      IF( zdf_qiao_alloc > 0 )  CALL ctl_warn('zdf_qiao_alloc: allocation of arrays failed.')
      !
   END FUNCTION zdf_qiao_alloc
      
   !!======================================================================
END MODULE zdfqiao
