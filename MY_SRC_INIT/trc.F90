MODULE trc
   !!======================================================================
   !!                      ***  MODULE  trc  ***
   !! Passive tracers   :  module for tracers defined
   !!======================================================================
   !! History :   OPA  !  1996-01  (M. Levy)  Original code
   !!              -   !  2000-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
   !!   NEMO      1.0  !  2004-03  (C. Ethe)  Free form and module
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc
   USE bdy_oce, only: ln_bdy, nb_bdy, OBC_DATA
   
   IMPLICIT NONE
   PUBLIC

   PUBLIC   trc_alloc   ! called by nemogcm.F90

   !! parameters for the control of passive tracers
   !! ---------------------------------------------   
   INTEGER, PUBLIC                                                 ::   numnat_ref = -1   !: logical unit for the reference passive tracer namelist_top_ref
   INTEGER, PUBLIC                                                 ::   numnat_cfg = -1   !: logical unit for the reference passive tracer namelist_top_cfg
   INTEGER, PUBLIC                                                 ::   numont     = -1   !: logical unit for the reference passive tracer namelist output output.namelist.top
   INTEGER, PUBLIC                                                 ::   numtrc_ref = -1   !: logical unit for the reference passive tracer namelist_top_ref
   INTEGER, PUBLIC                                                 ::   numtrc_cfg = -1   !: logical unit for the reference passive tracer namelist_top_cfg
   INTEGER, PUBLIC                                                 ::   numonr     = -1   !: logical unit for the reference passive tracer namelist output output.namelist.top
   INTEGER, PUBLIC                                                 ::   numstr        !: logical unit for tracer statistics
   INTEGER, PUBLIC                                                 ::   numrtr        !: logical unit for trc restart (read )
   INTEGER, PUBLIC                                                 ::   numrtw        !: logical unit for trc restart ( write )

   !! passive tracers fields (before,now,after)
   !! --------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)               ::  trai           !: initial total tracer
   REAL(wp), PUBLIC                                                ::  areatot        !: total volume 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  )         ::  cvol           !: volume correction -degrad option- 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)         ::  trn            !: tracer concentration for now time step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)         ::  tra            !: tracer concentration for next time step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)         ::  trb            !: tracer concentration for before time step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  )         ::  sbc_trc_b      !: Before sbc fluxes for tracers
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  )         ::  sbc_trc        !: Now sbc fluxes for tracers

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  )         ::  trc_i          !: prescribed tracer concentration in sea ice for SBC
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  )         ::  trc_o          !: prescribed tracer concentration in ocean for SBC
   INTEGER             , PUBLIC                                    ::  nn_ice_tr      !: handling of sea ice tracers

   !! interpolated gradient
   !!--------------------------------------------------  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)           ::  gtru           !: hor. gradient at u-points at bottom ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)           ::  gtrv           !: hor. gradient at v-points at bottom ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)           ::  gtrui          !: hor. gradient at u-points at top    ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)           ::  gtrvi          !: hor. gradient at v-points at top    ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)             ::  qsr_mean        !: daily mean qsr
   
   !! passive tracers  (input and output)
   !! ------------------------------------------  
   LOGICAL             , PUBLIC                                    ::  ln_rsttr       !: boolean term for restart i/o for passive tracers (namelist)
   LOGICAL             , PUBLIC                                    ::  lrst_trc       !: logical to control the trc restart write
   INTEGER             , PUBLIC                                    ::  nn_writetrc    !: time step frequency for concentration outputs (namelist)
   INTEGER             , PUBLIC                                    ::  nutwrs         !: output FILE for passive tracers restart
   INTEGER             , PUBLIC                                    ::  nutrst         !: logical unit for restart FILE for passive tracers
   INTEGER             , PUBLIC                                    ::  nn_rsttr       !: control of the time step ( 0 or 1 ) for pass. tr.
   CHARACTER(len = 80) , PUBLIC                                    ::  cn_trcrst_in   !: suffix of pass. tracer restart name (input)
   CHARACTER(len = 256), PUBLIC                                    ::  cn_trcrst_indir  !: restart input directory
   CHARACTER(len = 80) , PUBLIC                                    ::  cn_trcrst_out  !: suffix of pass. tracer restart name (output)
   CHARACTER(len = 256), PUBLIC                                    ::  cn_trcrst_outdir  !: restart output directory
   REAL(wp)            , PUBLIC                                    ::  rdttrc         !: passive tracer time step
   REAL(wp)            , PUBLIC                                    ::  r2dttrc        !: = 2*rdttrc except at nit000 (=rdttrc) if neuler=0
   LOGICAL             , PUBLIC                                    ::  ln_top_euler   !: boolean term for euler integration 
   LOGICAL             , PUBLIC                                    ::  ln_trcdta      !: Read inputs data from files
   LOGICAL             , PUBLIC                                    ::  ln_trcdmp      !: internal damping flag
   LOGICAL             , PUBLIC                                    ::  ln_trcdmp_clo  !: internal damping flag on closed seas
   INTEGER             , PUBLIC                                    ::  nittrc000      !: first time step of passive tracers model
   LOGICAL             , PUBLIC                                    ::  l_trcdm2dc     !: Diurnal cycle for TOP

   !! Information for the ice module for tracers
   !! ------------------------------------------
   TYPE TRC_I_NML                    !--- Ice tracer namelist structure
         REAL(wp)         :: trc_ratio  ! ice-ocean trc ratio
         REAL(wp)         :: trc_prescr ! prescribed ice trc cc
         CHARACTER(len=2) :: ctrc_o     ! choice of ocean trc cc
   END TYPE

   REAL(wp)        , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)  :: trc_ice_ratio      ! ice-ocean tracer ratio
   REAL(wp)        , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)  :: trc_ice_prescr     ! prescribed ice trc cc
   CHARACTER(len=2), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)  :: cn_trc_o ! choice of ocean tracer cc


   !! information for outputs
   !! --------------------------------------------------
   TYPE, PUBLIC :: PTRACER                                                            !: Passive tracer type
       CHARACTER(len = 20)  :: clsname  !: short name
       CHARACTER(len = 80)  :: cllname  !: long name
       CHARACTER(len = 20)  :: clunit   !: unit
       LOGICAL              :: llinit   !: read in a file or not
       LOGICAL              :: llsbc   !: read in a file or not
       LOGICAL              :: llcbc   !: read in a file or not
       LOGICAL              :: llobc   !: read in a file or not
   END TYPE PTRACER

   CHARACTER(len = 20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ctrcnm         !: tracer name 
   CHARACTER(len = 80), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ctrcln         !: trccer field long name
   CHARACTER(len = 20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ctrcun         !: tracer unit

   TYPE, PUBLIC :: DIAG                                                               !: passive trcacer ddditional diagnostic type
      CHARACTER(len = 20)  :: sname    !: short name
      CHARACTER(len = 80)  :: lname    !: long name
      CHARACTER(len = 20)  :: units    !: unit
   END TYPE DIAG

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)         ::  trc3d          !: 3D diagnostics for tracers
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)         ::  trc2d          !: 2D diagnostics for tracers

   !! information for inputs
   !! --------------------------------------------------
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ln_trc_ini     !: Initialisation from data input file
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ln_trc_obc     !: Use open boundary condition data
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ln_trc_sbc     !: Use surface boundary condition data
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ln_trc_cbc     !: Use coastal boundary condition data
   LOGICAL , PUBLIC                                     ::  ln_rnf_ctl     !: remove runoff dilution on tracers
   REAL(wp), PUBLIC                                     ::  rn_bc_time     !: Time scaling factor for SBC and CBC data (seconds in a day)


   !! variables to average over physics over passive tracer sub-steps.
   !! ----------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  un_tm       !: i-horizontal velocity average     [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  vn_tm       !: j-horizontal velocity average     [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::  tsn_tm      !: t/s average     [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  avt_tm      !: vertical diffusivity coeff. at  w-point   [m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  rhop_tm     !: 
# if defined key_zdfddm
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  avs_tm      !: vertical double diffusivity coeff. at w-point   [m/s]
# endif
#if defined key_trabbl
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  ahu_bbl_tm  !: u-, w-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  ahv_bbl_tm  !: j-direction slope at u-, w-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  utr_bbl_tm  !: j-direction slope at u-, w-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  vtr_bbl_tm  !: j-direction slope at u-, w-points
#endif
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  sshn_tm     !: average ssh for the now step [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  sshb_hold   !:hold sshb from the beginning of each sub-stepping[m]  

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  rnf_tm     !: river runoff
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  h_rnf_tm   !: depth in metres to the bottom of the relevant grid box
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  hmld_tm    !: mixed layer depth average [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  fr_i_tm    !: average ice fraction     [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  emp_tm     !: freshwater budget: volume flux [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  fmmflx_tm  !: freshwater budget: freezing/melting [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  emp_b_hold !: hold emp from the beginning of each sub-stepping[m]  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  qsr_tm     !: solar radiation average [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  wndm_tm    !: 10m wind average [m]
   !

   ! Temporary physical arrays for sub_stepping
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::  tsn_temp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  un_temp,vn_temp,wn_temp     !: hold current values of avt, un, vn, wn
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  avt_temp, rhop_temp     !: hold current values of avt, un, vn, wn
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  sshn_temp, sshb_temp, ssha_temp, rnf_temp,h_rnf_temp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  hdivn_temp, rotn_temp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  hdivb_temp, rotb_temp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  hmld_temp, qsr_temp, fr_i_temp,wndm_temp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  emp_temp, fmmflx_temp, emp_b_temp
   !
#if defined key_trabbl
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  ahu_bbl_temp, ahv_bbl_temp, utr_bbl_temp, vtr_bbl_temp !: hold current values 
#endif
   !
# if defined key_zdfddm
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  avs_temp      !: salinity vertical diffusivity coeff. at w-point   [m/s]
# endif
   !
   ! DRM, 04/07/17 - Fix the dimensions of the namelist variables to prevent segfault due to ALLOCATABLE and assumed-size-arrays.
   !
   CHARACTER(len=20), PUBLIC, DIMENSION(1)   ::  cn_trc_dflt          ! Default OBC condition for all tracers
   CHARACTER(len=20), PUBLIC, DIMENSION(1)   ::  cn_trc               ! Choice of boundary condition for tracers
   INTEGER,           PUBLIC, DIMENSION(1)   ::  nn_trcdmp_bdy        !: =T Tracer damping
   ! External data structure of BDY for TOP. Available elements: cn_obc, ll_trc, trcnow, dmp
   TYPE(OBC_DATA),    PUBLIC, DIMENSION(1,1), TARGET ::  trcdta_bdy           !: bdy external data (local process)
   !

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3.1 , NEMO Consortium (2010)
   !! $Id: trc.F90 7881 2017-04-07 08:33:13Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trc_alloc()
      !!-------------------------------------------------------------------
      !!                    *** ROUTINE trc_alloc ***
      !!-------------------------------------------------------------------
      USE lib_mpp, ONLY: ctl_warn
      !!-------------------------------------------------------------------
      INTEGER :: ierr(4)
      !!-------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( trn(jpi,jpj,jpk,jptra), trb(jpi,jpj,jpk,jptra), tra(jpi,jpj,jpk,jptra),       &  
         &      trc_i(jpi,jpj,jptra)  , trc_o(jpi,jpj,jptra)                          ,       &
         &      gtru (jpi,jpj,jptra)  , gtrv (jpi,jpj,jptra)                          ,       &
         &      gtrui(jpi,jpj,jptra)  , gtrvi(jpi,jpj,jptra)                          ,       &
         &      trc_ice_ratio(jptra)  , trc_ice_prescr(jptra) , cn_trc_o(jptra)       ,       &
         &      sbc_trc_b(jpi,jpj,jptra), sbc_trc(jpi,jpj,jptra)                      ,       &  
         &      cvol(jpi,jpj,jpk)     , trai(jptra)           , qsr_mean(jpi,jpj)     ,       &
         &      ctrcnm(jptra)         , ctrcln(jptra)         , ctrcun(jptra)         ,       &
         &      ln_trc_ini(jptra)     ,                                                       &
         &      ln_trc_sbc(jptra)     , ln_trc_cbc(jptra)     , ln_trc_obc(jptra)     ,       &
         &      STAT = ierr(1)  )
      !
      ! DRM, 04/07/17 - disable allocate statemet, size of namelist variables now fixed.
      !IF ( ln_bdy ) THEN
      !   ALLOCATE( cn_trc_dflt(nb_bdy)   , cn_trc(nb_bdy)     , nn_trcdmp_bdy(nb_bdy) ,       &
      !   &      trcdta_bdy(jptra,nb_bdy)                                              ,       &
      !   &      STAT = ierr(2)  )
      !ENDIF
      !
      IF (jp_dia3d > 0 ) ALLOCATE( trc3d(jpi,jpj,jpk,jp_dia3d), STAT = ierr(3) )
      !
      IF (jp_dia2d > 0 ) ALLOCATE( trc2d(jpi,jpj,jpk,jp_dia2d), STAT = ierr(4) )
      ! 
      trc_alloc = MAXVAL( ierr )
      IF( trc_alloc /= 0 )   CALL ctl_warn('trc_alloc: failed to allocate arrays')
      !
   END FUNCTION trc_alloc

   !!======================================================================
END MODULE trc
