MODULE asmpar
   !!======================================================================
   !!                       ***  MODULE asmpar  ***
   !! Assimilation increment : Parameters for assimilation interface
   !!======================================================================

   IMPLICIT NONE
   PRIVATE

   CHARACTER(LEN=40), PUBLIC, PARAMETER ::   c_asmbkg = 'assim_background_state_Jb'   !: Filename for storing the background state
   !                                                                                  !  for use in the Jb term
   CHARACTER(LEN=40), PUBLIC, PARAMETER ::   c_asmdin = 'assim_background_state_DI'   !: Filename for storing the background state
   !                                                                                  !  for direct initialization
   CHARACTER(LEN=40), PUBLIC, PARAMETER ::   c_asmtrj = 'assim_trj'                   !: Filename for storing the reference trajectory
   CHARACTER(LEN=40), PUBLIC, PARAMETER ::   c_asminc = 'assim_background_increments' !: Filename for storing the increments 
   !                                                                                  !  to the background state

   INTEGER, PUBLIC ::   nitbkg_r      !: Background time step referenced to nit000
   INTEGER, PUBLIC ::   nitdin_r      !: Direct Initialization time step referenced to nit000
   INTEGER, PUBLIC ::   nitiaustr_r   !: IAU starting time step referenced to nit000
   INTEGER, PUBLIC ::   nitiaufin_r   !: IAU final time step referenced to nit000
   INTEGER, PUBLIC ::   nittrjfrq     !: Frequency of trajectory output for 4D-VAR

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: asmpar.F90 5836 2015-10-26 14:49:40Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE asmpar
