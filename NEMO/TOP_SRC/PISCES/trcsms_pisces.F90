MODULE trcsms_pisces
   !!======================================================================
   !!                         ***  MODULE trcsms_pisces  ***
   !! TOP :   PISCES Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
   !!   trcsms_pisces        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE par_pisces
   USE sms_pisces
   USE p4zsms
   USE p2zsms

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_pisces    ! called in trcsms.F90
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_pisces.F90 7646 2017-02-06 09:25:03Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ini_pisces ***
      !!
      !! ** Purpose :   Initialisation of the PISCES biochemical model
      !!----------------------------------------------------------------------


   SUBROUTINE trc_sms_pisces( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms_pisces  ***
      !!
      !! ** Purpose :   Managment of the call to Biological sources and sinks 
      !!                routines of PISCES or LOBSTER bio-model
      !!
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!---------------------------------------------------------------------
      !
      IF( ln_p4z .OR. ln_p5z ) THEN  ;   CALL p4z_sms( kt )   !  PISCES
      ELSE                           ;   CALL p2z_sms( kt )   !  LOBSTER
      ENDIF

      !
   END SUBROUTINE trc_sms_pisces

   !!======================================================================
END MODULE trcsms_pisces 
