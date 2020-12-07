MODULE limthd_da
   !!======================================================================
   !!                       ***  MODULE limthd_da ***
   !! LIM-3 sea-ice :  computation of lateral melting in the ice
   !!======================================================================
   !! History :   4.0   ! 2016-03 (C. Rousset) original code
   !!---------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3'                                      LIM-3 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_thd_da   : sea ice lateral melting
   !!----------------------------------------------------------------------
   USE par_oce                ! ocean parameters
   USE phycst                 ! physical constants (ocean directory)
   USE sbc_oce, ONLY: sst_m   ! Surface boundary condition: ocean fields
   USE ice                    ! LIM variables
   USE lib_mpp                ! MPP library
   USE wrk_nemo               ! work arrays
   USE lib_fortran            ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_thd_da        ! called by limthd module

   !!----------------------------------------------------------------------
   !! NEMO/LIM3 4.0 , UCL - NEMO Consortium (2011)
   !! $Id: limthd_da.F90 5123 2015-03-04 16:06:03Z clem $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_thd_da
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE lim_thd_da  ***    
      !!   
      !! ** Purpose :   computes sea ice lateral melting
      !!
      !! ** Method  :   dA/dt = - P * W   [s-1]
      !!                   W = melting velocity [m.s-1]
      !!                   P = perimeter of ice-ocean lateral interface normalized by grid cell area [m.m-2]
      !!
      !!                   W = m1 * (Tw -Tf)**m2                    --- originally from Josberger 1979 ---
      !!                      (Tw - Tf) = elevation of water temp above freezing
      !!                      m1 and m2 = (1.6e-6 , 1.36) best fit from field experiment near the coast of Prince Patrick Island (Perovich 1983) => static ice
      !!                      m1 and m2 = (3.0e-6 , 1.36) best fit from MIZEX 84 experiment (Maykut and Perovich 1987) => moving ice
      !!
      !!                   P = N * pi * D                           --- from Rothrock and Thorndike 1984 ---
      !!                      D = mean floe caliper diameter
      !!                      N = number of floes = ice area / floe area(average) = A / (Cs * D**2)
      !!                         A = ice concentration
      !!                         Cs = deviation from a square (square:Cs=1 ; circle:Cs=pi/4 ; floe:Cs=0.66)
      !!
      !!                   D = Dmin * ( Astar / (Astar-A) )**beta   --- from Lupkes et al., 2012 (eq. 26-27) ---
      !!                                                             
      !!                      Astar = 1 / ( 1 - (Dmin/Dmax)**(1/beta) )
      !!                      Dmin = minimum floe diameter (recommended to be 8m +- 20%)
      !!                      Dmax = maximum floe diameter (recommended to be 300m, but it does not impact melting much except for Dmax<100m)
      !!                      beta = 1.0 +-20% (recommended value)
      !!                           = 0.3 best fit for western Fram Strait and Antarctica
      !!                           = 1.4 best fit for eastern Fram Strait
      !!
      !! ** Tunable parameters  :   We propose to tune the lateral melting via 2 parameters
      !!                               Dmin [6-10m]   => 6  vs 8m = +40% melting at the peak (A~0.5)
      !!                                                 10 vs 8m = -20% melting
      !!                               beta [0.8-1.2] => decrease = more melt and melt peaks toward higher concentration
      !!                                                                  (A~0.5 for beta=1 ; A~0.8 for beta=0.2)
      !!                                                 0.3 = best fit for western Fram Strait and Antarctica
      !!                                                 1.4 = best fit for eastern Fram Strait
      !!
      !! ** Note   :   Former and more simple formulations for floe diameters can be found in Mai (1995), 
      !!               Birnbaum and Lupkes (2002), Lupkes and Birnbaum (2005). They are reviewed in Lupkes et al 2012
      !!               A simpler implementation for CICE can be found in Bitz et al (2001) and Tsamados et al (2015)
      !!
      !! ** References
      !!    Bitz, C. M., Holland, M. M., Weaver, A. J., & Eby, M. (2001).
      !!              Simulating the ice‐thickness distribution in a coupled climate model.
      !!              Journal of Geophysical Research: Oceans, 106(C2), 2441-2463.
      !!    Josberger, E. G. (1979).
      !!              Laminar and turbulent boundary layers adjacent to melting vertical ice walls in salt water
      !!              (No. SCIENTIFIC-16). WASHINGTON UNIV SEATTLE DEPT OF ATMOSPHERIC SCIENCES.
      !!    Lüpkes, C., Gryanik, V. M., Hartmann, J., & Andreas, E. L. (2012).
      !!              A parametrization, based on sea ice morphology, of the neutral atmospheric drag coefficients
      !!              for weather prediction and climate models.
      !!              Journal of Geophysical Research: Atmospheres, 117(D13).
      !!    Maykut, G. A., & Perovich, D. K. (1987).
      !!              The role of shortwave radiation in the summer decay of a sea ice cover.
      !!              Journal of Geophysical Research: Oceans, 92(C7), 7032-7044.
      !!    Perovich, D. K. (1983).
      !!              On the summer decay of a sea ice cover. (Doctoral dissertation, University of Washington).
      !!    Rothrock, D. A., & Thorndike, A. S. (1984).
      !!              Measuring the sea ice floe size distribution.
      !!              Journal of Geophysical Research: Oceans, 89(C4), 6477-6486.
      !!    Tsamados, M., Feltham, D., Petty, A., Schroeder, D., & Flocco, D. (2015).
      !!              Processes controlling surface, bottom and lateral melt of Arctic sea ice in a state of the art sea ice model.
      !!              Phil. Trans. R. Soc. A, 373(2052), 20140167.
      !!---------------------------------------------------------------------
      INTEGER             ::   ji, jj, jl      ! dummy loop indices
      REAL(wp)            ::   zastar, zdfloe, zperi, zwlat, zda
      REAL(wp), PARAMETER ::   zdmax = 300._wp
      REAL(wp), PARAMETER ::   zcs   = 0.66_wp
      REAL(wp), PARAMETER ::   zm1   = 3.e-6_wp
      REAL(wp), PARAMETER ::   zm2   = 1.36_wp
      !
      REAL(wp), POINTER, DIMENSION(:,:) ::   zda_tot
      !!---------------------------------------------------------------------
      CALL wrk_alloc( jpi,jpj, zda_tot )

      !------------------------------------------------------------!
      ! --- Calculate reduction of total sea ice concentration --- !
      !------------------------------------------------------------!
      zastar = 1._wp / ( 1._wp - (rn_dmin / zdmax)**(1._wp/rn_beta) )
      
      DO jj = 1, jpj
         DO ji = 1, jpi
            
            ! Mean floe caliper diameter [m]
            zdfloe = rn_dmin * ( zastar / ( zastar - at_i(ji,jj) ) )**rn_beta
            
            ! Mean perimeter of the floe = N*pi*D = (A/cs*D^2)*pi*D [m.m-2]
            zperi = at_i(ji,jj) * rpi / ( zcs * zdfloe )
            
            ! Melt speed rate [m/s]
            zwlat = zm1 * ( MAX( 0._wp, sst_m(ji,jj) - ( t_bo(ji,jj) - rt0 ) ) )**zm2
            
            ! sea ice concentration decrease
            zda_tot(ji,jj) = - MIN( zwlat * zperi * rdt_ice, at_i(ji,jj) )
            
         END DO
      END DO
      
      !---------------------------------------------------------------------------------------------!
      ! --- Distribute reduction among ice categories and calculate associated ice-ocean fluxes --- !
      !---------------------------------------------------------------------------------------------!
      DO jl = jpl, 1, -1
         DO jj = 1, jpj
            DO ji = 1, jpi
               
               ! decrease of concentration for the category jl
               !    1st option: each category contributes to melting in proportion to its concentration
               rswitch = MAX( 0._wp , SIGN( 1._wp, at_i(ji,jj) - epsi10 ) )
               zda     = rswitch * zda_tot(ji,jj) * a_i(ji,jj,jl) / MAX( at_i(ji,jj), epsi10 )
               !    2d option: melting of the upper cat first
               !!zda = MAX( zda_tot(ji,jj), - a_i(ji,jj,jl) )
               !!zda_tot(ji,jj) = zda_tot(ji,jj) + zda
               
               ! Contribution to salt flux
               sfx_lam(ji,jj) = sfx_lam(ji,jj) - rhoic *  ht_i(ji,jj,jl) * zda * sm_i(ji,jj,jl) * r1_rdtice
               
               ! Contribution to heat flux into the ocean [W.m-2], <0  
               hfx_thd(ji,jj) = hfx_thd(ji,jj) + zda * r1_rdtice * ( ht_i(ji,jj,jl) * SUM( e_i(ji,jj,:,jl) ) * r1_nlay_i  &
                  &                                                + ht_s(ji,jj,jl) *      e_s(ji,jj,1,jl)   * r1_nlay_s )
               
               ! Contribution to mass flux
               wfx_lam(ji,jj) =  wfx_lam(ji,jj) - zda * r1_rdtice * ( rhoic * ht_i(ji,jj,jl) + rhosn * ht_s(ji,jj,jl) )
               
               ! new concentration
               a_i(ji,jj,jl) = a_i(ji,jj,jl) + zda
            END DO
         END DO
      END DO
      
      ! total concentration
      at_i(:,:) = SUM( a_i(:,:,:), dim=3 )
      
      ! --- ensure that ht_i = 0 where a_i = 0 ---
      WHERE( a_i == 0._wp ) ht_i = 0._wp
      !
      CALL wrk_dealloc( jpi,jpj, zda_tot )
      !
   END SUBROUTINE lim_thd_da
   
#else
   !!----------------------------------------------------------------------
   !!   Default option         Dummy Module          No LIM-3 sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_thd_da          ! Empty routine
   END SUBROUTINE lim_thd_da
#endif
   !!======================================================================
END MODULE limthd_da
