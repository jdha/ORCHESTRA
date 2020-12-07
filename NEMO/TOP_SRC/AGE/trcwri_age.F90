MODULE trcwri_age
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    age :   Output of age tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top &&  defined key_iomput
   !!----------------------------------------------------------------------
   !! trc_wri_age   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE par_age     
   USE trc         
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_age 

CONTAINS

   SUBROUTINE trc_wri_age
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      CHARACTER (len=20)   :: cltra
      INTEGER              :: jn
      !!---------------------------------------------------------------------

      ! write the tracer concentrations in the file

      cltra = TRIM( ctrcnm(jp_age) )                  ! short title for tracer
      CALL iom_put( cltra, trn(:,:,:,jp_age) )

      !
   END SUBROUTINE trc_wri_age

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_age
CONTAINS
   SUBROUTINE trc_wri_age                     ! Empty routine  
   END SUBROUTINE trc_wri_age
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcwri_age.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_age
