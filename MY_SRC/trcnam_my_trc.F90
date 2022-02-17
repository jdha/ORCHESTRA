MODULE trcnam_my_trc
   !!======================================================================
   !!                      ***  MODULE trcnam_my_trc  ***
   !! TOP :   initialisation of some run parameters for MY_TRC bio-model
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec) Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_nam_my_trc      : MY_TRC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_my_trc   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2016)
   !! $Id: trcnam_my_trc.F90 7646 2017-02-06 09:25:03Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_my_trc
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_my_trc  ***  
      !!
      !! ** Purpose :   read MY_TRC namelist
      !!
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_my_trc : set NAge variables, no namelist'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~'
      !
      ! Variable setting
      ctrcnm    (jp_myt0:jp_myt1) = 'NAge'
      ctrcln    (jp_myt0:jp_myt1) = 'Sea water age since northern boundary'
      ctrcun    (jp_myt0:jp_myt1) = 'year'
      ln_trc_ini(jp_myt0:jp_myt1) = .false.
      ln_trc_sbc(jp_myt0:jp_myt1) = .false.
      ln_trc_cbc(jp_myt0:jp_myt1) = .false.
      ln_trc_obc(jp_myt0:jp_myt1) = .false.
      !
   END SUBROUTINE trc_nam_my_trc
   
   !!======================================================================
END MODULE trcnam_my_trc
