MODULE trctrp
   !!======================================================================
   !!                       ***  MODULE trctrp  ***
   !! Ocean Physics    : manage the passive tracer transport
   !!======================================================================
   !! History :   1.0  !  2004-03 (C. Ethe) Original code
   !!             3.3  !  2010-07 (C. Ethe) Merge TRA-TRC
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_trp        : passive tracer transport
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and active tracers variables
   USE trc             ! ocean passive tracers variables 
   USE trabbl          ! bottom boundary layer               (trc_bbl routine)
   USE trcbbl          ! bottom boundary layer               (trc_bbl routine)
   USE trcdmp          ! internal damping                    (trc_dmp routine)
   USE trcldf          ! lateral mixing                      (trc_ldf routine)
   USE trcadv          ! advection                           (trc_adv routine)
   USE trczdf          ! vertical diffusion                  (trc_zdf routine)
   USE trcnxt          ! time-stepping                       (trc_nxt routine)
   USE trcrad          ! positivity                          (trc_rad routine)
   USE trcsbc          ! surface boundary condition          (trc_sbc routine)
   USE zpshde          ! partial step: hor. derivative       (zps_hde routine)
   USE bdy_oce   , ONLY: ln_bdy
   USE trcbdy          ! BDY open boundaries

#if defined key_agrif
   USE agrif_top_sponge ! tracers sponges
   USE agrif_top_update ! tracers updates
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_trp    ! called by trc_stp

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trctrp.F90 7646 2017-02-06 09:25:03Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_trp( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_trp  ***
      !!                      
      !! ** Purpose :   Management of passive tracers transport
      !! 
      !! ** Method  : - Compute the passive tracers trends 
      !!              - Update the passive tracers
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
      !! ---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_trp')
      !
      IF( .NOT. lk_c1d ) THEN
         !
                                CALL trc_sbc    ( kt )      ! surface boundary condition
         IF( lk_trabbl )        CALL trc_bbl    ( kt )      ! advective (and/or diffusive) bottom boundary layer scheme
         IF( ln_trcdmp )        CALL trc_dmp    ( kt )      ! internal damping trends
         IF( ln_bdy )           CALL trc_bdy_dmp( kt )      ! BDY damping trends
                                CALL trc_adv    ( kt )      ! horizontal & vertical advection 
         !                                                         ! Partial top/bottom cell: GRADh( trb )  
         IF( ln_zps ) THEN
           IF( ln_isfcav ) THEN ; CALL zps_hde_isf( kt, jptra, trb, pgtu=gtru, pgtv=gtrv, pgtui=gtrui, pgtvi=gtrvi )  ! both top & bottom
           ELSE                 ; CALL zps_hde    ( kt, jptra, trb, gtru, gtrv )                                      !  only bottom
           ENDIF
         ENDIF
         !                                                      
                                CALL trc_ldf    ( kt )      ! lateral mixing
#if defined key_agrif
         IF(.NOT. Agrif_Root()) CALL Agrif_Sponge_trc       ! tracers sponge
#endif
                                CALL trc_zdf    ( kt )      ! vertical mixing and after tracer fields
                                CALL trc_nxt    ( kt )      ! tracer fields at next time step     
         IF( ln_trcrad )        CALL trc_rad    ( kt )      ! Correct artificial negative concentrations
         IF( ln_trcdmp_clo )    CALL trc_dmp_clo( kt )      ! internal damping trends on closed seas only

#if defined key_agrif
         IF( .NOT.Agrif_Root()) CALL Agrif_Update_Trc( kt ) ! Update tracer at AGRIF zoom boundaries : children only
#endif
         !
      ELSE                                               ! 1D vertical configuration
                                CALL trc_sbc( kt )            ! surface boundary condition
         IF( ln_trcdmp )        CALL trc_dmp( kt )            ! internal damping trends
                                CALL trc_zdf( kt )            ! vertical mixing and after tracer fields
                                CALL trc_nxt( kt )            ! tracer fields at next time step     
          IF( ln_trcrad )       CALL trc_rad( kt )            ! Correct artificial negative concentrations
         !
      END IF
      !
      IF( nn_timing == 1 )   CALL timing_stop('trc_trp')
      !
   END SUBROUTINE trc_trp

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                                        No TOP models
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_trp( kt )              ! Empty routine
      INTEGER, INTENT(in) ::   kt
      WRITE(*,*) 'trc_trp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_trp
#endif
   
   !!======================================================================
END MODULE trctrp
