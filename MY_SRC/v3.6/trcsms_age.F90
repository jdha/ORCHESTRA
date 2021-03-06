MODULE trcsms_age
   !!======================================================================
   !!                         ***  MODULE trcsms_age  ***
   !! TOP :   Main module of the AGE tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
   !! trc_sms_age       : AGE model main routine
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_age       ! called by trcsms.F90 module

   INTEGER , PUBLIC :: nl_age             ! T level surrounding age_depth
   INTEGER , PUBLIC :: nla_age            ! T level wholly above age_depth
   INTEGER , PUBLIC :: nlb_age            ! T level wholly below age_depth

   REAL(wp), PUBLIC :: rn_age_depth       ! = 10       depth over which age tracer reset to zero
   REAL(wp), PUBLIC :: rn_age_kill_rate   ! = -1./7200  recip of relaxation timescale (s) for  age tracer shallower than age_depth
   
   REAL(wp), PUBLIC :: rryear          !: recip number of seconds in one year
   REAL(wp), PUBLIC :: frac_kill_age   !: fraction of level nl_age above age_depth where it is relaxed towards zero
   REAL(wp), PUBLIC :: frac_add_age    !: fraction of level nl_age below age_depth where it is incremented


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_age( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_age  ***
      !!
      !! ** Purpose :   main routine of AGE model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER ::   jn, jk   ! dummy loop index
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_age')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_age:  AGE model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'


      DO jk = 1, nla_age
         tra(:,:,jk,jp_age) = rn_age_kill_rate * trb(:,:,jk,jp_age)
      ENDDO
      !
      tra(:,:,nl_age,jp_age) = frac_kill_age * rn_age_kill_rate * trb(:,:,nl_age,jp_age)  &
          &                   + frac_add_age  * rryear * tmask(:,:,nl_age)
      !
      DO jk = nlb_age, jpk
         tra(:,:,jk,jp_age) = tmask(:,:,jk) * rryear
      ENDDO
      !
      ! DRM - include aging of northern boundary tracer at every depth level.
      DO jk = 1, jpk
         tra(:,:,jk,jp_age+1) = tmask(:,:,jk) * rryear
      ENDDO
      !
      ! DRM - save northern boundary tracer trends.
      !
      IF( l_trdtrc )THEN
         CALL trd_trc( tra(:,:,:,jp_age), jp_age, jptra_sms, kt )   ! save trends
         CALL trd_trc( tra(:,:,:,jp_age+1), jp_age+1, jptra_sms, kt )   ! save trends
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_age')
      !
   END SUBROUTINE trc_sms_age

   !!======================================================================
END MODULE trcsms_age
