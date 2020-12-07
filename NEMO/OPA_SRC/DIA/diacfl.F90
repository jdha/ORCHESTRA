MODULE diacfl
   !!==============================================================================
   !!                       ***  MODULE  diacfl  ***
   !! Output CFL diagnostics to ascii file
   !!==============================================================================
   !! History :  1.0  !  2010-03  (E. Blockley)  Original code
   !!                 !  2014-06  (T Graham) Removed CPP key & Updated to vn3.6
   !! 
   !!----------------------------------------------------------------------
   !!   dia_cfl        : Compute and output Courant numbers at each timestep
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE lib_mpp         ! distribued memory computing
   USE lbclnk          ! ocean lateral boundary condition (or mpp link)
   USE in_out_manager  ! I/O manager
   USE domvvl     
   USE timing          ! Performance output

   IMPLICIT NONE
   PRIVATE

   REAL(wp) :: cu_max, cv_max, cw_max                      ! Run max U Courant number 
   INTEGER, DIMENSION(3) :: cu_loc, cv_loc, cw_loc         ! Run max locations
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: zcu_cfl           ! Courant number arrays
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: zcv_cfl           ! Courant number arrays
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: zcw_cfl           ! Courant number arrays

   INTEGER  :: numcfl                                       ! outfile unit
   CHARACTER(LEN=50) :: clname="cfl_diagnostics.ascii"      ! ascii filename

   PUBLIC   dia_cfl       ! routine called by step.F90
   PUBLIC   dia_cfl_init  ! routine called by nemogcm

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LOCEAN-IPSL (2009) 
   !! $Id: diacfl.F90 8329 2017-07-13 14:14:54Z rblod $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------


CONTAINS


   SUBROUTINE dia_cfl ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_cfl  ***
      !!
      !! ** Purpose :  Compute the Courant numbers Cu=u*dt/dx and Cv=v*dt/dy
      !!               and output to ascii file 'cfl_diagnostics.ascii'
      !!----------------------------------------------------------------------

      INTEGER, INTENT(in) ::  kt                            ! ocean time-step index

      REAL(wp) :: zcu_max, zcv_max, zcw_max                 ! max Courant numbers per timestep
      INTEGER, DIMENSION(3) :: zcu_loc, zcv_loc, zcw_loc    ! max Courant number locations

      REAL(wp) :: dt                                        ! temporary scalars
      INTEGER, DIMENSION(3) :: zlocu, zlocv, zlocw          ! temporary arrays 
      INTEGER  :: ji, jj, jk                                ! dummy loop indices

      
      IF( nn_diacfl == 1) THEN
         IF( nn_timing == 1 )   CALL timing_start('dia_cfl')
         ! setup timestep multiplier to account for initial Eulerian timestep
         IF( neuler == 0 .AND. kt == nit000 ) THEN   ;    dt = rdt
         ELSE                                        ;    dt = rdt * 2.0
         ENDIF

             ! calculate Courant numbers
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, fs_jpim1   ! vector opt.

                  ! Courant number for x-direction (zonal current)
                  zcu_cfl(ji,jj,jk) = ABS(un(ji,jj,jk))*dt/e1u(ji,jj)

                  ! Courant number for y-direction (meridional current)
                  zcv_cfl(ji,jj,jk) = ABS(vn(ji,jj,jk))*dt/e2v(ji,jj)

                  ! Courant number for z-direction (vertical current)
                  zcw_cfl(ji,jj,jk) = ABS(wn(ji,jj,jk))*dt/e3w_n(ji,jj,jk)
               END DO
            END DO         
         END DO

         ! calculate maximum values and locations
         IF( lk_mpp ) THEN
            CALL mpp_maxloc(zcu_cfl,umask,zcu_max, zcu_loc(1), zcu_loc(2), zcu_loc(3))
            CALL mpp_maxloc(zcv_cfl,vmask,zcv_max, zcv_loc(1), zcv_loc(2), zcv_loc(3))
            CALL mpp_maxloc(zcw_cfl,tmask,zcw_max, zcw_loc(1), zcw_loc(2), zcw_loc(3))
         ELSE
            zlocu = MAXLOC( ABS( zcu_cfl(:,:,:) ) )
            zcu_loc(1) = zlocu(1) + nimpp - 1
            zcu_loc(2) = zlocu(2) + njmpp - 1
            zcu_loc(3) = zlocu(3)
            zcu_max = zcu_cfl(zcu_loc(1),zcu_loc(2),zcu_loc(3))

            zlocv = MAXLOC( ABS( zcv_cfl(:,:,:) ) )
            zcv_loc(1) = zlocv(1) + nimpp - 1
            zcv_loc(2) = zlocv(2) + njmpp - 1
            zcv_loc(3) = zlocv(3)
            zcv_max = zcv_cfl(zcv_loc(1),zcv_loc(2),zcv_loc(3))

            zlocw = MAXLOC( ABS( zcw_cfl(:,:,:) ) )
            zcw_loc(1) = zlocw(1) + nimpp - 1
            zcw_loc(2) = zlocw(2) + njmpp - 1
            zcw_loc(3) = zlocw(3)
            zcw_max = zcw_cfl(zcw_loc(1),zcw_loc(2),zcw_loc(3))
         ENDIF
      
         ! write out to file
         IF( lwp ) THEN
            WRITE(numcfl,FMT='(2x,i4,5x,a6,5x,f6.4,1x,i4,1x,i4,1x,i4)') kt, 'Max Cu', zcu_max, zcu_loc(1), zcu_loc(2), zcu_loc(3)
            WRITE(numcfl,FMT='(11x,a6,5x,f6.4,1x,i4,1x,i4,1x,i4)') 'Max Cv', zcv_max, zcv_loc(1), zcv_loc(2), zcv_loc(3)
            WRITE(numcfl,FMT='(11x,a6,5x,f6.4,1x,i4,1x,i4,1x,i4)') 'Max Cw', zcw_max, zcw_loc(1), zcw_loc(2), zcw_loc(3)
         ENDIF

         ! update maximum Courant numbers from whole run if applicable
         IF( zcu_max > cu_max ) THEN
            cu_max = zcu_max
            cu_loc = zcu_loc
         ENDIF
         IF( zcv_max > cv_max ) THEN
            cv_max = zcv_max
            cv_loc = zcv_loc
         ENDIF
         IF( zcw_max > cw_max ) THEN
            cw_max = zcw_max
            cw_loc = zcw_loc
         ENDIF

         ! at end of run output max Cu and Cv and close ascii file
         IF( kt == nitend .AND. lwp ) THEN
            ! to ascii file
            WRITE(numcfl,*) '******************************************'
            WRITE(numcfl,FMT='(3x,a12,7x,f6.4,1x,i4,1x,i4,1x,i4)') 'Run Max Cu', cu_max, cu_loc(1), cu_loc(2), cu_loc(3)
            WRITE(numcfl,FMT='(3x,a8,11x,f7.1)') ' => dt/C', dt*(1.0/cu_max)
            WRITE(numcfl,*) '******************************************'
            WRITE(numcfl,FMT='(3x,a12,7x,f6.4,1x,i4,1x,i4,1x,i4)') 'Run Max Cv', cv_max, cv_loc(1), cv_loc(2), cv_loc(3)
            WRITE(numcfl,FMT='(3x,a8,11x,f7.1)') ' => dt/C', dt*(1.0/cv_max)
            WRITE(numcfl,*) '******************************************'
            WRITE(numcfl,FMT='(3x,a12,7x,f6.4,1x,i4,1x,i4,1x,i4)') 'Run Max Cw', cw_max, cw_loc(1), cw_loc(2), cw_loc(3)
            WRITE(numcfl,FMT='(3x,a8,11x,f7.1)') ' => dt/C', dt*(1.0/cw_max)
            CLOSE( numcfl ) 

            ! to ocean output
            WRITE(numout,*)
            WRITE(numout,*) 'dia_cfl     : Maximum Courant number information for the run:'
            WRITE(numout,*) '~~~~~~~~~~~~'
            WRITE(numout,FMT='(12x,a12,7x,f6.4,5x,a16,i4,1x,i4,1x,i4,a1)') 'Run Max Cu', cu_max, 'at (i, j, k) =   &
                          &   (', cu_loc(1), cu_loc(2), cu_loc(3), ')'
            WRITE(numout,FMT='(12x,a8,11x,f7.1)') ' => dt/C', dt*(1.0/cu_max)
            WRITE(numout,FMT='(12x,a12,7x,f6.4,5x,a16,i4,1x,i4,1x,i4,a1)') 'Run Max Cv', cv_max, 'at (i, j, k) =   &
                          &   (', cv_loc(1), cv_loc(2), cv_loc(3), ')'
            WRITE(numout,FMT='(12x,a8,11x,f7.1)') ' => dt/C', dt*(1.0/cv_max)
            WRITE(numout,FMT='(12x,a12,7x,f6.4,5x,a16,i4,1x,i4,1x,i4,a1)') 'Run Max Cw', cw_max, 'at (i, j, k) =   &
                          &   (', cw_loc(1), cw_loc(2), cw_loc(3), ')'
            WRITE(numout,FMT='(12x,a8,11x,f7.1)') ' => dt/C', dt*(1.0/cw_max)

         ENDIF

         IF( nn_timing == 1 )   CALL timing_stop('dia_cfl')
      ENDIF

   END SUBROUTINE dia_cfl

   SUBROUTINE dia_cfl_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_cfl_init  ***
      !!                   
      !! ** Purpose :   create output file, initialise arrays
      !!----------------------------------------------------------------------


      IF( nn_diacfl == 1 ) THEN
         IF( nn_timing == 1 )   CALL timing_start('dia_cfl_init')

         cu_max=0.0
         cv_max=0.0
         cw_max=0.0

         ALLOCATE( zcu_cfl(jpi, jpj, jpk), zcv_cfl(jpi, jpj, jpk), zcw_cfl(jpi, jpj, jpk) )

         zcu_cfl(:,:,:)=0.0
         zcv_cfl(:,:,:)=0.0
         zcw_cfl(:,:,:)=0.0

         IF( lwp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'dia_cfl     : Outputting CFL diagnostics to '//TRIM(clname)
            WRITE(numout,*) '~~~~~~~~~~~~'
            WRITE(numout,*)

            ! create output ascii file
            CALL ctl_opn( numcfl, clname, 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, 1 )
            WRITE(numcfl,*) 'Timestep  Direction  Max C     i    j    k'
            WRITE(numcfl,*) '******************************************'
         ENDIF

         IF( nn_timing == 1 )   CALL timing_stop('dia_cfl_init')

      ENDIF

   END SUBROUTINE dia_cfl_init

END MODULE diacfl
