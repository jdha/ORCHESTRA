MODULE p4zligand
   !!======================================================================
   !!                         ***  MODULE p4zligand  ***
   !! TOP :   PISCES Compute remineralization/dissolution of organic ligands
   !!=========================================================================
   !! History :   3.6  !  2016-03  (O. Aumont, A. Tagliabue) Quota model and reorganization
   !!----------------------------------------------------------------------
   !!   p4z_ligand       :  Compute remineralization/dissolution of organic ligands
   !!   p4z_ligand_init  :  Initialisation of parameters for remineralisation
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_ligand         ! called in p4zbio.F90
   PUBLIC   p4z_ligand_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  rlgw     !: lifetime (years) of weak ligands
   REAL(wp), PUBLIC ::  rlgs     !: lifetime (years) of strong ligands
   REAL(wp), PUBLIC ::  rlig     !: Remin ligand production
   REAL(wp), PUBLIC ::  prlgw    !: Photochemical of weak ligand
   REAL(wp), PUBLIC ::  rfep     !: Dissolution rate of FeP

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zligand.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_ligand( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_ligand  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic ligands
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zlgwp, zlgwpr, zlgwr, zlablgw, zrfepa, zfepr
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_ligand')
      !
      ! ------------------------------------------------------------------
      ! Remineralization of iron ligands
      ! ------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! production from remineralisation of organic matter
               zlgwp  = orem(ji,jj,jk) * rlig
               ! decay of weak ligand
               ! This is based on the idea that as LGW is lower
               ! there is a larger fraction of refractory OM
               zlgwr = max( rlgs , rlgw * exp( -2 * (trb(ji,jj,jk,jplgw)*1e9) ) ) ! years
               zlgwr = 1. / zlgwr * tgfunc(ji,jj,jk) * ( xstep / nyear_len(1) ) * trb(ji,jj,jk,jplgw)
               ! photochem loss of weak ligand
               zlgwpr = prlgw * xstep * etot(ji,jj,jk) * trb(ji,jj,jk,jplgw) * (1. - fr_i(ji,jj))
               tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + zlgwp - zlgwr - zlgwpr
            END DO
         END DO
      END DO

      ! ----------------------------------------------------------
      ! Dissolution of nanoparticle Fe
      ! ----------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! dissolution rate is maximal in the presence of light and 
               ! lower in the aphotici zone
               ! ! 25 Wm-2 constant
               zrfepa = rfep * ( 1. - EXP( -1. * etot(ji,jj,jk) / 25. ) ) * (1.- fr_i(ji,jj))
               zrfepa = MAX( (zrfepa / 10.0), zrfepa ) ! min of 10 days lifetime
               zfepr = rfep * xstep * trb(ji,jj,jk,jpfep)
               tra(ji,jj,jk,jpfep) = tra(ji,jj,jk,jpfep) - zfepr
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zfepr
            END DO
         END DO
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('ligand1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_ligand')
      !
   END SUBROUTINE p4z_ligand


   SUBROUTINE p4z_ligand_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_ligand_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampislig namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampislig
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampislig/ rlgw, prlgw, rlgs, rfep, rlig
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampislig in reference namelist : Pisces remineralization
      READ  ( numnatp_ref, nampislig, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampislig in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampislig in configuration namelist : Pisces remineralization
      READ  ( numnatp_cfg, nampislig, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampislig in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampislig )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for ligands, nampislig'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Dissolution rate of FeP                        rfep =', rfep
         WRITE(numout,*) '    Lifetime (years) of weak ligands               rlgw =', rlgw
         WRITE(numout,*) '    Remin ligand production per unit C             rlig =', rlig
         WRITE(numout,*) '    Photolysis of weak ligand                     prlgw =', prlgw
         WRITE(numout,*) '    Lifetime (years) of strong ligands             rlgs =', rlgs
      ENDIF
      !
   END SUBROUTINE p4z_ligand_init

   !!======================================================================
END MODULE p4zligand
