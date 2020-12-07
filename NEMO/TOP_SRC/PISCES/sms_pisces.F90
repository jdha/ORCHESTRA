MODULE sms_pisces   
   !!----------------------------------------------------------------------
   !!                     ***  sms_pisces.F90  ***  
   !! TOP :   PISCES Source Minus Sink variables
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc

   IMPLICIT NONE
   PUBLIC

   INTEGER ::   numnatp_ref = -1           !! Logical units for namelist pisces
   INTEGER ::   numnatp_cfg = -1           !! Logical units for namelist pisces
   INTEGER ::   numonp      = -1           !! Logical unit for namelist pisces output

   !                                                       !:  PISCES  : silicon dependant half saturation

   !!* Model used
   LOGICAL  ::  ln_p2z            !: Flag to use LOBSTER model
   LOGICAL  ::  ln_p4z            !: Flag to use PISCES  model
   LOGICAL  ::  ln_p5z            !: Flag to use PISCES  quota model
   LOGICAL  ::  ln_ligand         !: Flag to enable organic ligands

   !!*  Time variables
   INTEGER  ::   nrdttrc           !: ???
   INTEGER  ::   ndayflxtr         !: ???
   REAL(wp) ::   rfact , rfactr    !: ???
   REAL(wp) ::   rfact2, rfact2r   !: ???
   REAL(wp) ::   xstep             !: Time step duration for biology
   REAL(wp) ::   ryyss             !: number of seconds per year 
   REAL(wp) ::   r1_ryyss          !: inverse number of seconds per year 


   !!*  Biological parameters 
   INTEGER  ::   niter1max, niter2max   !: Maximum number of iterations for sinking
   REAL(wp) ::   rno3              !: ???
   REAL(wp) ::   o2ut              !: ???
   REAL(wp) ::   po4r              !: ???
   REAL(wp) ::   rdenit            !: ???
   REAL(wp) ::   rdenita           !: ???
   REAL(wp) ::   o2nit             !: ???
   REAL(wp) ::   wsbio, wsbio2     !: ???
   REAL(wp) ::   wsbio2max         !: ???
   REAL(wp) ::   wsbio2scale       !: ???
   REAL(wp) ::   xkmort            !: ???
   REAL(wp) ::   ferat3            !: ???
   REAL(wp) ::   wfep              !: ???
   REAL(wp) ::   ldocp             !: ???
   REAL(wp) ::   ldocz             !: ???
   REAL(wp) ::   lthet             !: ???
   REAL(wp) ::   no3rat3           !: ???
   REAL(wp) ::   po4rat3           !: ???


   !!*  diagnostic parameters 
   REAL(wp) ::  tpp                !: total primary production
   REAL(wp) ::  t_oce_co2_exp      !: total carbon export
   REAL(wp) ::  t_oce_co2_flx      !: Total ocean carbon flux
   REAL(wp) ::  t_oce_co2_flx_cum  !: Cumulative Total ocean carbon flux
   REAL(wp) ::  t_atm_co2_flx      !: global mean of atmospheric pco2

   !!* restoring
   LOGICAL  ::  ln_pisdmp          !: restoring or not of nutrients to a mean value
   INTEGER  ::  nn_pisdmp          !: frequency of relaxation or not of nutrients to a mean value

   !!* Mass conservation
   LOGICAL  ::  ln_check_mass      !: Flag to check mass conservation

   !!*  Biological fluxes for light : variables shared by pisces & lobster
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  neln  !: number of T-levels + 1 in the euphotic layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup  !: euphotic layer depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etot  !: par (photosynthetic available radiation)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etot_ndcy      !: PAR over 24h in case of diurnal cycle
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  enano, ediat   !: PAR for phyto, nano and diat 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  epico          !: PAR for pico
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  emoy           !: averaged PAR in the mixed layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup_01 !: Absolute euphotic layer depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  xksi  !:  LOBSTER : zooplakton closure

   !!*  Biological fluxes for primary production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)    ::   xksimax    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   biron      !: bioavailable fraction of iron
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   plig       !: proportion of iron organically complexed

   !!*  Sinking speed
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio4   !: GOC sinking speed
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wscal    !: Calcite and BSi sinking speeds
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsfep



   !!*  SMS for the organic matter
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xfracal    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   nitrfac    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   orem       !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xdiss      !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodcal    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodpoc    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   conspoc    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodgoc    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   consgoc    !: Calcite production

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizen      !: size of diatoms 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizep      !: size of diatoms 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sized      !: size of diatoms 


   !!* Variable for chemistry of the CO2 cycle
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak13       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak23       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksp       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hi         !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   excess     !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aphscale   !: 


   !!* Temperature dependancy of SMS terms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc    !: Temp. dependancy of various biological rates
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc2   !: Temp. dependancy of mesozooplankton rates

#if defined key_sed
   LOGICAL, PUBLIC, PARAMETER ::   lk_sed = .TRUE.     !: sediment flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_sed = .FALSE.     !: sediment flag
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: sms_pisces.F90 7646 2017-02-06 09:25:03Z timgraham $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS


   INTEGER FUNCTION sms_pisces_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_pisces_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_warn
      INTEGER ::   ierr(10)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !*  Biological fluxes for light : shared variables for pisces & lobster
      ALLOCATE( etot(jpi,jpj,jpk), neln(jpi,jpj), heup(jpi,jpj),    &
        &       heup_01(jpi,jpj) , xksi(jpi,jpj)               ,  STAT=ierr(1) )
      !
  
      IF( ln_p4z .OR. ln_p5z ) THEN
         !*  Biological fluxes for light 
         ALLOCATE(  enano(jpi,jpj,jpk)    , ediat(jpi,jpj,jpk) ,   &
           &        etot_ndcy(jpi,jpj,jpk), emoy(jpi,jpj,jpk)  ,  STAT=ierr(2) ) 

         !*  Biological fluxes for primary production
         ALLOCATE( xksimax(jpi,jpj)  , biron(jpi,jpj,jpk)      ,  STAT=ierr(3) )
         !
         !*  SMS for the organic matter
         ALLOCATE( xfracal (jpi,jpj,jpk), nitrfac(jpi,jpj,jpk) ,    &
            &      orem    (jpi,jpj,jpk),                           &
            &      prodcal(jpi,jpj,jpk),  xdiss   (jpi,jpj,jpk),    &
            &      prodpoc(jpi,jpj,jpk) , conspoc(jpi,jpj,jpk) ,    &
            &      prodgoc(jpi,jpj,jpk) , consgoc(jpi,jpj,jpk) ,  STAT=ierr(4) )

         !* Variable for chemistry of the CO2 cycle
         ALLOCATE( ak13  (jpi,jpj,jpk) ,                            &
            &      ak23(jpi,jpj,jpk)    , aksp  (jpi,jpj,jpk) ,     &
            &      hi  (jpi,jpj,jpk)    , excess(jpi,jpj,jpk) ,     &
            &      aphscale(jpi,jpj,jpk),                         STAT=ierr(5) )
         !
         !* Temperature dependancy of SMS terms
         ALLOCATE( tgfunc(jpi,jpj,jpk)  , tgfunc2(jpi,jpj,jpk),   STAT=ierr(6) )
         !
         !* Sinkong speed
         ALLOCATE( wsbio3 (jpi,jpj,jpk) , wsbio4 (jpi,jpj,jpk),     &
            &      wscal(jpi,jpj,jpk)                         ,   STAT=ierr(7) )   
         ! 
         IF( ln_ligand ) THEN
           ALLOCATE( plig(jpi,jpj,jpk)  , wsfep(jpi,jpj,jpk)  ,   STAT=ierr(8) )
         ENDIF
         !
      ENDIF
      !
      IF( ln_p5z ) THEN
         !       
         ALLOCATE( epico(jpi,jpj,jpk)                         ,   STAT=ierr(9) ) 

         !*  Size of phytoplankton cells
         ALLOCATE( sizen(jpi,jpj,jpk), sizep(jpi,jpj,jpk),         &
           &       sized(jpi,jpj,jpk),                            STAT=ierr(10) )
      ENDIF
      !
      sms_pisces_alloc = MAXVAL( ierr )
      !
      IF( sms_pisces_alloc /= 0 )   CALL ctl_warn('sms_pisces_alloc: failed to allocate arrays') 
      !
   END FUNCTION sms_pisces_alloc

   !!======================================================================   
END MODULE sms_pisces    
