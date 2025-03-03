SetSFPR = function(){
  
  SFPR = list(
    t1        = c(min = 0      , max =   100   ),# Thresholds to Sugawara Tank Model
    t2        = c(min = 0      , max =   100   ),
    a1        = c(min = 0      , max =   1     ),
    a2        = c(min = 0      , max =   1     ),
    i1        = c(min = 0      , max =   1     ),
    t3        = c(min = 0      , max =   100   ),
    t4        = c(min = 0      , max =   100   ),
    D_I        = c(min = 0     , max =    1    ), # Inflow to IR
    alpha_Qq_FR= c(min = 0.1   , max =    10   ), # alpha in Q=K*S^alpha
    m_E_FR     = c(min = 0.01  , max =     2   ), # smoothing factor for FR evaporation
    K_Qq_FR    = c(min = 1e-4  , max =    10   ), # K in Q=K*S^alpha
    K_Qb_FR    = c(min = 0     , max =     4   ), # controls exchange between FR and SR
    K_Qb_UR    = c(min = 1e-06 , max =     2   ), # controls percolation from UR to SR
    Ce         = c(min = 0.1   , max =     3   ), # evaporation multiplication factor
    Beta_Qq_UR = c(min = 0.001 , max =    10   ), # UR runoff coefficient parameter: power, logistic, linear, quadroid
    Beta_E_UR  = c(min = 0.01  , max =    10   ), # factor for UR evaporation : mkinetic
    Smax_UR    = c(min = 0.1   , max =  1e04   ), # maximum  reservoir capacity for unsaturated reservoir
    Smax_IR    = c(min = 0.1  , max =    20   ), # maximum reservoir capacity for interception reservoir
    m_QE_IR    = c(min = 0.001 , max =     1   ), # smoothing factor in interception reservoir. Power smoothing range is 0.3 to 0.5
    K_Qq_RR    = c(min = 0.05  , max =     4   ), # K in Q=K*, default=1
    Smax_CR    = c(min = 0.1   , max =  1e04   ),
    Umax_uCR   = c(min = 0.8   , max =  0.9    ),
    Smin_uCR   = c(min = 0.001 , max =     1   ),
    Beta_Qq_uCR= c(min = 0.001 , max =    10   ),
    mu_Qq_uCR  = c(min = 0.1   , max =  1   ),
    K_Qb_uCR   = c(min = 1e-06 , max =     1   ),
    Sevmax_CR  = c(min = 0.01   , max =  1e04  ),
    Beta_E_CR  = c(min = 0.01  , max =    10   ),
    Beta_Qq_sCR= c(min = 0.01  , max =    10   ),
    K_Qb_sCR   = c(min = 1e-06 , max =     2   ),
    K_Qd_sCR   = c(min = 0     , max =     4   ),
    Cp_WR      = c(min = 0.5   , max =     5   ),
    m_Q_WR     = c(min = 0.01  , max =     2   ),
    Kq_WR      = c(min = 0.01  , max =    10   ),
    Tp_WR      = c(min = 0     , max =    10   ),
    Tm_WR      = c(min = 0     , max =     4   ),
    K_Qq_SR    = c(min = 1e-07 , max =     1   ), # K in Q=K*S^alpha
    K_Qb_SR    = c(min = 1e-07 , max =     1   ), # K in Q=K*S : loss from groundwater reservoir
    Tlag       = c(min = 1     , max =    10   ), # Base of the half triangle lag function
    D_R        = c(min = 0     , max =     1   ), # Inflow to RR
    D_S        = c(min = 0     , max =     1   ), # Diverts input from FR to UR and SR (slow reservoirs)
    D_F        = c(min = 0     , max =     1   ), # Inflow to FR from UR
    D_C        = c(min = 0     , max =     1   ),
    SiniFr_UR  = c(min = 0     , max =     1   ), # state intial factor for unsaturated reservoir. Initial storage= SiniFr_UR*Smax_UR
    alpha_Qq_SR= c(min = 0.01  , max =    10  ), # alpha in Q=K*S^alpha
    m_E_SR     = c(min = 0.001 , max =    1    ), # evaporation smoothing factor for slow reservoir
    P_ED_max   = c(min = 0.1   , max =  1e07   ), # threshold for infiltration excess overland flow
    m_P_ED     = c(min = 1e-03 , max =    10   ), # Smoothing factor for infiltration excess overland flow
    mu_Qq_UR   = c(min = 0.1   , max =     1   ), # paramters to be used along with Beta_Qq_UR in modified logistic curve function
    Lag_RR     = c(min = 0     , max =     1   ),
    Lag_FR     = c(min = 0     , max =     1   ),
    Lag_SR     = c(min = 0     , max =     1   ),
    option_i   = c(min = 1     , max =     2   ), #choice of functions for calculating IR discharge
    option_u   = c(min = 0.5 + .Machine[["double.xmin"]]  , max = 4.5 - .Machine[["double.xmin"]]), #choice of functions for calculating UR discharge, TRUE range is 1:4, but this settings is needed due to scaling     
    option_f   = c(min = 1     , max =     2   ), #choice of functions for calculating FR discharge 
    option_c   = c(min = 1     , max =     2   ),
    option_s   = c(min = 1     , max =     2   ), #choice of functions for calculating SR discharge 
    u_res      = c(min = 0     , max =     1   ), #choice of unsaturated reservoir 1: on or 0:off
    s_res      = c(min = 0     , max =     1   ), #choice of Slow reservoir on or off
    i_res      = c(min = 0     , max =     1   ), #choice of interception reservoir on or off
    r_res      = c(min = 0     , max =     1   ), #choice of riparian reservoir on or off
    f_res      = c(min = 0     , max =     1   ), #choice of fast reservoir on or off
    c_res      = c(min = 0     , max =     1   ),
    w_res      = c(min = 0     , max =     1   ),
    
    ###### FUSE-Parameters(24)
    timedelay  = c(min = 0.01  , max =     5   ),
    rferr_add  = c(min = 0     , max =     0   ),
    rferr_mlt  = c(min = 1     , max =     1   ),
    maxwatr_1  = c(min = 25    , max =   500   ),
    maxwatr_2  = c(min = 50    , max =  5000   ),
    fracten  = c(min = 0.05    , max =  0.95   ),
    frchzne  = c(min = 0.05    , max =  0.95   ),
    fprimqb  = c(min = 0.05    , max =  0.95   ),
    rtfrac1  = c(min = 0.05    , max =  0.95   ),
    percrte  = c(min = 0.01    , max =  1000   ),
    percexp  = c(min = 1      , max =     20   ),
    sacpmlt  = c(min = 1     , max =     250   ),
    sacpexp  = c(min = 1       , max =     5   ),
    percfrac  = c(min = 0.05  , max =   0.95   ),
    iflwrte  = c(min = 0.01     , max = 1000   ),
    baserte  = c(min = 0.001    , max = 1000   ),
    qb_powr  = c(min = 1      , max =     10   ),
    qb_prms  = c(min = 0.001   , max =  0.25   ),
    qbrate_2a  = c(min = 0.001  , max = 0.25   ),
    qbrate_2b  = c(min = 0.001 , max =  0.25   ),
    sareamax  = c(min = 0.05   , max =  0.95   ),
    axv_bexp  = c(min = 0.001     , max =  3   ),
    loglamb  = c(min = 5     , max =  10   ),
    tishape  = c(min = 2     , max =  5  ),
    
    ###### FUSE-Options(8)
    rferr  = c(min = 11     , max =  12  ),
    arch1  = c(min = 21     , max =  23  ),
    arch2  = c(min = 31     , max =  34  ),
    qsurf  = c(min = 41     , max =  43  ),
    qperc  = c(min = 51     , max =  53  ),
    esoil  = c(min = 61     , max =  62  ),
    qintf  = c(min = 71     , max =  72  ),
    q_tdh  = c(min = 81     , max =  82  ),
    
    ##### DCT
    R_t1   = c(min = 0.01     , max =  5  ),
    R_t2   = c(min = 0.01     , max =  5  ),
    R_t3   = c(min = 0.01     , max =  5  ),
    R_1   = c(min = 624     , max =  625  ),
    R_2   = c(min = 624     , max =  625  ),
    R_3   = c(min = 624     , max =  625  ),
    
    #### MM_29 (M29_smax,M29_b,M29_a,M29_kf,M29_ks,M29_s1,M29_s2,M29_s3,M29_s4,M29_s5)
    M29_smax = c(min =1, max =2000),
    M29_b = c(min =0, max =10),
    M29_a = c(min =0, max =1),
    M29_kf = c(min =0, max =1),
    M29_ks= c(min =0, max =1),
    M29_s1 = c(min =0, max =1),
    M29_s2 = c(min =1, max =300),
    M29_s3 = c(min =1, max =300),
    M29_s4 = c(min =1, max =300),
    M29_s5 = c(min =1, max =300),
    
    #M01
    M01_smax = c(min =1, max =2000),
    M01_s1 = c(min =0, max =1),
    
    #M02
    M02_dw = c(min =0, max =5),
    M02_betaw = c(min =0, max =10),
    M02_swmax = c(min =1, max =2000),
    M02_kw = c(min =0, max =1),
    M02_s1 = c(min =0, max =1),
    
    #M03
    M03_S1max = c(min =1, max =2000),
    M03_Sfc = c(min =0.05, max =0.95),
    M03_a = c(min =0, max =1),
    M03_M = c(min =0.05, max =0.95),
    M03_s1 = c(min =0, max =1),
    
    #M04
    M04_S1max = c(min =1, max =2000),
    M04_Sfc = c(min =0.05, max =0.95),
    M04_m = c(min =0.05, max =0.95),
    M04_a = c(min =0, max =1),
    M04_b = c(min =1, max =5),
    M04_tcbf = c(min =0, max =1),
    M04_s1 = c(min =0, max =1),
    
    #M05 M05_lp,M05_d,M05_p,M05_alpha,M05_tau_q,M05_tau_s,M05_tau_d
    M05_lp = c(min =1, max =2000),
    M05_d= c(min =1, max =2000),
    M05_p= c(min =0, max =10),
    M05_alpha= c(min =0, max =1),
    M05_tau_q= c(min =1, max =700),
    M05_tau_s= c(min =1, max =700),
    M05_tau_d= c(min =0, max =119),
    M05_s1 = c(min =1, max =2000),
    
    #M06 "M06_tt","M06_ddf","M06_smax","M06_tc","M06_s1","M06_s2"
    M06_tt= c(min =(-3), max =5),
    M06_ddf= c(min =0, max =20),
    M06_smax= c(min =1, max =2000),
    M06_tc= c(min =0, max =1),
    M06_s1= c(min =1, max =2000),
    M06_s2= c(min =0, max =1),
    
    #M07 "M07_x1","M07_x2","M07_x3","M07_x4","M07_s1","M07_s2"
    M07_x1= c(min =1, max =2000),
    M07_x2= c(min =(-10), max =15),
    M07_x3= c(min =1, max =300),
    M07_x4= c(min =1, max =15),
    M07_s1= c(min =0, max =1),
    M07_s2= c(min =0, max =1),
    
    #M08 "M08_alpha_ei","M08_m","M08_smax","M08_fc","M08_alpha_ss","M08_s1","M08_s2"
    M08_alpha_ei= c(min =0, max =1),
    M08_m= c(min =0.05, max =0.95),
    M08_smax= c(min =1, max =2000),
    M08_fc= c(min =0.05, max =0.95),
    M08_alpha_ss= c(min =0, max =1),
    M08_s1= c(min =0, max =1),
    M08_s2= c(min =0, max =1),
    
    #M09 "M09_sb","M09_sfc","M09_m","M09_a","M09_b","M09_r","M09_s1","M09_s2"
    M09_sb= c(min =1, max =2000),
    M09_sfc= c(min =0.05, max =0.95),
    M09_m= c(min =0.05, max =0.95),
    M09_a= c(min =1, max =50),
    M09_b= c(min =0.2, max =1),
    M09_r= c(min =0, max =1),
    M09_s1= c(min =0, max =1),
    M09_s2= c(min =0, max =1),
    
    #M10 "M10_sb","M10_phi","M10_fc","M10_r","M10_c","M10_d","M10_s1","M10_s2"
    M10_sb= c(min =1, max =2000),
    M10_phi= c(min =0.05, max =0.95),
    M10_fc= c(min =0.05, max =0.95),
    M10_r= c(min =0, max =1),
    M10_c= c(min =0, max =1),
    M10_d= c(min =1, max =5),
    M10_s1= c(min =0, max =1),
    M10_s2= c(min =0, max =1),
    
    #M11 "M11_s1max","M11_sfc","M11_a","M11_M","M11_b","M11_lambda","M11_s1","M11_s2"
    M11_s1max= c(min =1, max =2000),
    M11_sfc= c(min =0.05, max =0.95),
    M11_a= c(min =0, max =1),
    M11_M= c(min =0.05, max =0.95),
    M11_b= c(min =1, max =5),
    M11_lambda= c(min =0, max =1),
    M11_s1= c(min =0, max =1),
    M11_s2= c(min =1, max =2000),
    
    #M12 "M12_tt","M12_ddf","M12_smax","M12_Cfc","M12_tcin","M12_tcbf","M12_s1","M12_s2"
    M12_tt= c(min =(-3), max =5),
    M12_ddf= c(min =0, max =20),
    M12_smax= c(min =1, max =2000),
    M12_Cfc= c(min =0.05, max =0.95),
    M12_tcin= c(min =0, max =1),
    M12_tcbf= c(min =0, max =1),
    M12_s1= c(min =1, max =2000),
    M12_s2= c(min =0, max =1),
    
    #M13 "M13_dw","M13_betaw","M13_swmax","M13_a","M13_th","M13_c","M13_kh","M13_s1","M13_s2"
    M13_dw= c(min =0, max =5),
    M13_betaw= c(min =0, max =10),
    M13_swmax= c(min =1, max =2000),
    M13_a= c(min =0, max =1),
    M13_th= c(min =1, max =120),
    M13_c= c(min =0, max =4),
    M13_kh= c(min =0, max =1),
    M13_s1= c(min =0, max =1),
    M13_s2= c(min =1, max =2000),
    
    #M14 "M14_suzmax","M14_st","M14_kd","M14_q0","M14_f","M14_chi","M14_phi","M14_s1","M14_s2"
    M14_suzmax= c(min =1, max =2000),
    M14_st= c(min =0.05, max =0.95),
    M14_kd= c(min =0, max =1),
    M14_q0= c(min =0.1, max =200),
    M14_f= c(min =0, max =1),
    M14_chi= c(min =1, max =7.5),
    M14_phi= c(min =0.1, max =5),
    M14_s1= c(min =0, max =1),
    M14_s2= c(min =1, max =2000),
    
    #M15 "M15_fmax","M15_dp","M15_sumax","M15_lp","M15_p","M15_tp","M15_c","M15_kp","M15_s1","M15_s2"
    M15_fmax= c(min =0, max =200),
    M15_dp= c(min =0, max =5),
    M15_sumax= c(min =1, max =2000),
    M15_lp= c(min =0.05, max =0.95),
    M15_p= c(min =0, max =1),
    M15_tp= c(min =1, max =120),
    M15_c= c(min =0, max =4),
    M15_kp= c(min =0, max =1),
    M15_s1= c(min =0, max =1),
    M15_s2= c(min =1, max =2000),
    
    #M16 "M16_s1max","M16_s2max","M16_sfc","M16_m","M16_a","M16_b","M16_tcbf","M16_d","M16_s1","M16_s2"
    M16_s1max= c(min =0, max =5),
    M16_s2max= c(min =1, max =2000),
    M16_sfc= c(min =0.05, max =0.95),
    M16_m= c(min =0.05, max =0.95),
    M16_a= c(min =0, max =1),
    M16_b= c(min =1, max =5),
    M16_tcbf= c(min =0, max =1),
    M16_d= c(min =1, max =120),
    M16_s1= c(min =0, max =1),
    M16_s2= c(min =0, max =1),
    
    #M17 "M17_smax","M17_phi","M17_gam","M17_kl","M17_s1","M17_s2","M17_s3"
    M17_smax= c(min =1, max =2000),
    M17_phi= c(min =0, max =1),
    M17_gam= c(min =0, max =1),
    M17_kl= c(min =0, max =1),
    M17_s1= c(min =1, max =2000),
    M17_s2= c(min =1, max =2000),
    M17_s3= c(min =1, max =300),
    
    #M18 "M18_insc","M18_coeff","M18_sq","M18_smsc","M18_sub","M18_crak","M18_k","M18_s1","M18_s2","M18_s3"
    M18_insc= c(min =0, max =5),
    M18_coeff= c(min =0, max =600),
    M18_sq= c(min =0, max =15),
    M18_smsc= c(min =1, max =2000),
    M18_sub= c(min =0, max =1),
    M18_crak= c(min =0, max =1),
    M18_k= c(min =0, max =1),
    M18_s1= c(min =0, max =1),
    M18_s2= c(min =0, max =1),
    M18_s3= c(min =1, max =2000),
    
    #M19 M19_sb","M19_phi","M19_fc","M19_alpha_ss","M19_beta_ss","M19_k_deep","M19_alpha_bf","M19_beta_bf","M19_s1","M19_s2","M19_s3"
    M19_sb= c(min =1, max =2000),
    M19_phi= c(min =0.05, max =0.95),
    M19_fc= c(min =0.01, max =1),
    M19_alpha_ss= c(min =0, max =1),
    M19_beta_ss= c(min =1, max =5),
    M19_k_deep= c(min =0, max =1),
    M19_alpha_bf= c(min =0, max =1),
    M19_beta_bf= c(min =1, max =5),
    M19_s1= c(min =0, max =1),
    M19_s2= c(min =0, max =1),
    M19_s3= c(min =1, max =2000),
    
    #M20 "M20_c","M20_ndc","M20_smax","M20_emax","M20_frate","M20_b","M20_dpf","M20_sdrmax","M20_s1","M20_s2","M20_s3"
    M20_c= c(min =0, max =1),
    M20_ndc= c(min =0.05, max =0.95),
    M20_smax= c(min =1, max =2000),
    M20_emax= c(min =0, max =20),
    M20_frate= c(min =0, max =200),
    M20_b= c(min =0, max =1),
    M20_dpf= c(min =0, max =1),
    M20_sdrmax= c(min =1, max =300),
    M20_s1= c(min =1, max =2000),
    M20_s2= c(min =0, max =1),
    M20_s3= c(min =0, max =1),
    
    #M21 "M21_s1max","M21_beta","M21_d","M21_percmax","M21_lp","M21_nlagf","M21_nalgs","M21_kf","M21_ks","M21_s1","M21_s2","M21_s3"
    M21_s1max= c(min =1, max =2000),
    M21_beta= c(min =0, max =10),
    M21_d= c(min =0, max =1),
    M21_percmax= c(min =0, max =20),
    M21_lp= c(min =0.05, max =0.95),
    M21_nlagf= c(min =1, max =5),
    M21_nalgs= c(min =1, max =15),
    M21_kf= c(min =0, max =1),
    M21_ks= c(min =0, max =1),
    M21_s1= c(min =0, max =1),
    M21_s2= c(min =1, max =300),
    M21_s3= c(min =1, max =300),
    
    #M22 "M22_ibar","M22_idelta","M22_ishift","M22_stot","M22_fsm","M22_b","M22_k1","M22_c1","M22_k2","M22_c2","M22_s1","M22_s2","M22_s3"
    M22_ibar= c(min =0.1, max =5),
    M22_idelta= c(min =0, max =1),
    M22_ishift= c(min =1, max =365),
    M22_stot= c(min =1, max =2000),
    M22_fsm= c(min =0.01, max =0.99),
    M22_b= c(min =0, max =10),
    M22_k1= c(min =0, max =1),
    M22_c1= c(min =0, max =10),
    M22_k2= c(min =0, max =1),
    M22_c2= c(min =1, max =5),
    M22_s1= c(min =0, max =1),
    M22_s2= c(min =0, max =1),
    M22_s3= c(min =0, max =1),
    
    #M23 "","","","","","","","","","","","","","","","","","","","","","","","","","",""
    M23_af= c(min =0, max =200),
    M23_bf= c(min =0, max =5),
    M23_stot= c(min =1, max =2000),
    M23_xa= c(min =0.01, max =0.99),
    M23_xf= c(min =0.01, max =0.99),
    M23_na= c(min =0.01, max =0.99),
    M23_ac= c(min =0, max =5),
    M23_bc= c(min =0, max =10),
    M23_ass= c(min =0, max =5),
    M23_bss= c(min =0, max =10),
    M23_c= c(min =0, max =200),
    M23_ag= c(min =0, max =5),
    M23_bg= c(min =0, max =1),
    M23_gf= c(min =0, max =1),
    M23_df= c(min =0, max =10),
    M23_td= c(min =0, max =1),
    M23_ab= c(min =0, max =1),
    M23_bb= c(min =0.01, max =200),
    M23_ga= c(min =0, max =1),
    M23_da= c(min =0, max =10),
    M23_aa= c(min =0.01, max =200),
    M23_ba= c(min =1, max =5),
    M23_gb= c(min =0, max =1),
    M23_db= c(min =0, max =10),
    M23_s1= c(min =0, max =5),
    M23_s2= c(min =0, max =1),
    M23_s3= c(min =0, max =1),
    
    #M24 "M24_s1max","M24_tw","M24_tu","M24_se","M24_tc","M24_s1","M24_s2","M24_s3","M24_s4"
    M24_s1max= c(min =1, max =2000),
    M24_tw= c(min =0, max =1),
    M24_tu= c(min =0, max =1),
    M24_se= c(min =1, max =2000),
    M24_tc= c(min =0, max =1),
    M24_s1= c(min =0, max =1),
    M24_s2= c(min =1, max =2000),
    M24_s3= c(min =1, max =300),
    M24_s4= c(min =1, max =300),
    
    #M25 "M25_phi","M25_rc","M25_gam","M25_k1","M25_fa","M25_k2","M25_s1","M25_s2","M25_s3","M25_s4"
    M25_phi= c(min =0, max =1),
    M25_rc= c(min =1, max =2000),
    M25_gam= c(min =0, max =1),
    M25_k1= c(min =0, max =1),
    M25_fa= c(min =0, max =1),
    M25_k2= c(min =0, max =1),
    M25_s1= c(min =1, max =2000),
    M25_s2= c(min =1, max =2000),
    M25_s3= c(min =1, max =2000),
    M25_s4= c(min =1, max =2000),
    
    #M26 "M26_smax","M26_beta","M26_d","M26_percmax","M26_lp","M26_nlagf","M26_nlags","M26_kf","M26_ks","M26_imax","M26_s1","M26_s2","M26_s3","M26_s4"
    M26_smax= c(min =1, max =2000),
    M26_beta= c(min =0, max =10),
    M26_d= c(min =0, max =1),
    M26_percmax= c(min =0, max =20),
    M26_lp= c(min =0.05, max =0.95),
    M26_nlagf= c(min =1, max =5),
    M26_nlags= c(min =1, max =15),
    M26_kf= c(min =0, max =1),
    M26_ks= c(min =0, max =1),
    M26_imax= c(min =0, max =5),
    M26_s1= c(min =0, max =1),
    M26_s2= c(min =0, max =1),
    M26_s3= c(min =1, max =300),
    M26_s4= c(min =1, max =300),
    
    #M27 "M27_a0","M27_b0","M27_c0","M27_a1","M27_fa","M27_fb","M27_fc","M27_fd","M27_st","M27_f2","M27_f1","M27_f3","M27_s1","M27_s2","M27_s3","M27_s4"
    M27_a0= c(min =0, max =1),
    M27_b0= c(min =0, max =1),
    M27_c0= c(min =0, max =1),
    M27_a1= c(min =0, max =1),
    M27_fa= c(min =0, max =1),
    M27_fb= c(min =0, max =1),
    M27_fc= c(min =0, max =1),
    M27_fd= c(min =0, max =1),
    M27_st= c(min =1, max =2000),
    M27_f2= c(min =0.01, max =0.99),
    M27_f1= c(min =0.01, max =0.99),
    M27_f3= c(min =0.01, max =0.99),
    M27_s1= c(min =0, max =1),
    M27_s2= c(min =0, max =1),
    M27_s3= c(min =0, max =1),
    M27_s4= c(min =0, max =1),
    
    #M28 "M28_aim","M28_a","M28_b","M28_stot","M28_fwmx","M28_flm","M28_c","M28_ex","M28_ki","M28_kg","M28_ci","M28_cg","M28_s1","M28_s2","M28_s3","M28_s4"
    M28_aim= c(min =0, max =1),
    M28_a= c(min =(-0.49), max =0.49),
    M28_b= c(min =0, max =10),
    M28_stot= c(min =1, max =2000),
    M28_fwmx= c(min =0.01, max =0.99),
    M28_flm= c(min =0.01, max =0.99),
    M28_c= c(min =0.01, max =0.99),
    M28_ex= c(min =0, max =10),
    M28_ki= c(min =0, max =1),
    M28_kg= c(min =0, max =1),
    M28_ci= c(min =0, max =1),
    M28_cg= c(min =0, max =1),
    M28_s1= c(min =0, max =1),
    M28_s2= c(min =0, max =1),
    M28_s3= c(min =0, max =300),
    M28_s4= c(min =0, max =300),
    
    #M30 "M30_tcrit","M30_ddf","M30_s2max","M30_tw","M30_tu","M30_se","M30_tc","M30_s1","M30_s2","M30_s3","M30_s4","M30_s5"
    M30_tcrit= c(min =(-3), max =3),
    M30_ddf= c(min =0, max =20),
    M30_s2max= c(min =1, max =2000),
    M30_tw= c(min =0, max =1),
    M30_tu= c(min =0, max =1),
    M30_se= c(min =1, max =2000),
    M30_tc= c(min =0, max =1),
    M30_s1= c(min =1, max =2000),
    M30_s2= c(min =0, max =1),
    M30_s3= c(min =1, max =2000),
    M30_s4= c(min =1, max =300),
    M30_s5= c(min =1, max =300),
    
    #M31 M31_tcrit","M31_ddf","M31_s2max","M31_tw","M31_tu","M31_se","M31_s3max","M31_tc","M31_s1","M31_s2","M31_s3","M31_s4","M31_s5
    M31_tcrit= c(min =(-3), max =3),
    M31_ddf= c(min =0, max =20),
    M31_s2max= c(min =1, max =2000),
    M31_tw= c(min =0, max =1),
    M31_tu= c(min =0, max =1),
    M31_se= c(min =0.05, max =0.95),
    M31_s3max= c(min =1, max =2000),
    M31_tc= c(min =0, max =1),
    M31_s1= c(min =1, max =2000),
    M31_s2= c(min =0, max =1),
    M31_s3= c(min =0, max =1),
    M31_s4= c(min =1, max =300),
    M31_s5= c(min =1, max =300),
    
    #M32 M32_tcrit","M32_ddf","M32_s2max","M32_tw","M32_i_alpha","M32_i_s","M32_tu","M32_se","M32_s3max","M32_tc","M32_s1","M32_s2","M32_s3","M32_s4","M32_s5
    M32_tcrit= c(min =(-3), max =3),
    M32_ddf= c(min =0, max =20),
    M32_s2max= c(min =1, max =2000),
    M32_tw= c(min =0, max =1),
    M32_i_alpha= c(min =0, max =1),
    M32_i_s= c(min =1, max =365),
    M32_tu= c(min =0, max =1),
    M32_se= c(min =0.05, max =0.95),
    M32_s3max= c(min =1, max =2000),
    M32_tc= c(min =0, max =1),
    M32_s1= c(min =1, max =2000),
    M32_s2= c(min =0, max =1),
    M32_s3= c(min =0, max =1),
    M32_s4= c(min =1, max =300),
    M32_s5= c(min =1, max =300),
    
    #M33 "M33_pctim","M33_smax","M33_f1","M33_f2","M33_kuz","M33_rexp","M33_f3","M33_f4","M33_pfree","M33_klzp","M33_klzs","M33_s1","M33_s2","M33_s3","M33_s4","M33_s5"
    M33_pctim= c(min =0, max =1),
    M33_smax= c(min =1, max =2000),
    M33_f1= c(min =0.005, max =0.995),
    M33_f2= c(min =0.005, max =0.995),
    M33_kuz= c(min =0, max =1),
    M33_rexp= c(min =0, max =7),
    M33_f3= c(min =0.005, max =0.995),
    M33_f4= c(min =0.005, max =0.995),
    M33_pfree= c(min =0, max =1),
    M33_klzp= c(min =0, max =1),
    M33_klzs= c(min =0, max =1),
    M33_s1= c(min =0, max =1),
    M33_s2= c(min =0, max =1),
    M33_s3= c(min =0, max =1),
    M33_s4= c(min =0, max =1),
    M33_s5= c(min =0, max =1),
    
    #M34 "M34_smax","M34_beta","M34_d","M34_percmax","M34_lp","M34_nlagf","M34_nlags","M34_kf","M34_ks","M34_imax","M34_tt","M34_ddf","M34_s1","M34_s2","M34_s3","M34_s4","M34_s5"
    M34_smax= c(min =1, max =2000),
    M34_beta= c(min =0, max =10),
    M34_d= c(min =0, max =1),
    M34_percmax= c(min =0, max =20),
    M34_lp= c(min =0.05, max =0.95),
    M34_nlagf= c(min =1, max =5),
    M34_nlags= c(min =1, max =15),
    M34_kf= c(min =0, max =1),
    M34_ks= c(min =0, max =1),
    M34_imax= c(min =0, max =5),
    M34_tt= c(min =(-3), max =5),
    M34_ddf= c(min =0, max =20),
    M34_s1= c(min =1, max =2000),
    M34_s2= c(min =0, max =1),
    M34_s3= c(min =0, max =1),
    M34_s4= c(min =1, max =300),
    M34_s5= c(min =1, max =300),
    
    #M35 "M35_tcrit","M35_ddf","M35_s2max","M35_tw","M35_i_alpha","M35_i_s","M35_tmin","M35_trange","M35_tu","M35_se","M35_s3max","M35_tc","M35_s1","M35_s2","M35_s3","M35_s4","M35_s5"
    M35_tcrit= c(min =(-3), max =3),
    M35_ddf= c(min =0, max =20),
    M35_s2max= c(min =1, max =2000),
    M35_tw= c(min =0, max =1),
    M35_i_alpha= c(min =0, max =1),
    M35_i_s= c(min =1, max =365),
    M35_tmin= c(min =(-10), max =0),
    M35_trange= c(min =1, max =20),
    M35_tu= c(min =0, max =1),
    M35_se= c(min =0.05, max =0.95),
    M35_s3max= c(min =1, max =2000),
    M35_tc= c(min =0, max =1),
    M35_s1= c(min =0, max =2000),
    M35_s2= c(min =0, max =1),
    M35_s3= c(min =0, max =1),
    M35_s4= c(min =1, max =300),
    M35_s5= c(min =1, max =300),
    
    #M36 "M36_INSC","M36_COEFF","M36_SQ","M36_SMSC","M36_SUB","M36_CRAK","M36_EM","M36_DSC","M36_ADS","M36_MD","M36_VCOND","M36_DLEV","M36_K1","M36_K2","M36_K3","M36_s1","M36_s2","M36_s3","M36_s4","M36_s5"
    M36_INSC= c(min =0, max =5),
    M36_COEFF= c(min =0, max =600),
    M36_SQ= c(min =0, max =15),
    M36_SMSC= c(min =1, max =2000),
    M36_SUB= c(min =0, max =1),
    M36_CRAK= c(min =0, max =1),
    M36_EM= c(min =0, max =20),
    M36_DSC= c(min =0, max =50),
    M36_ADS= c(min =0, max =1),
    M36_MD= c(min =0.99, max =1),
    M36_VCOND= c(min =0, max =0.5),
    M36_DLEV= c(min =(-10), max =10),
    M36_K1= c(min =0, max =1),
    M36_K2= c(min =0, max =1),
    M36_K3= c(min =0, max =100),
    M36_s1= c(min =0, max =1),
    M36_s2= c(min =0, max =1),
    M36_s3= c(min =0, max =1),
    M36_s4= c(min =1, max =2000),
    M36_s5= c(min =1, max =300),
    
    #M37 "M37_tt","M37_tti","M37_ttm","M37_cfr","M37_cfmax","M37_whc","M37_cflux","M37_fc","M37_lp","M37_beta","M37_k0","M37_alpha","M37_perc","M37_K1","M37_maxbas","M37_s1","M37_s2","M37_s3","M37_s4","M37_s5"
    M37_tt= c(min =(-3), max =5),
    M37_tti= c(min =0, max =17),
    M37_ttm= c(min =(-3), max =3),
    M37_cfr= c(min =0, max =1),
    M37_cfmax= c(min =0, max =20),
    M37_whc= c(min =0, max =1),
    M37_cflux= c(min =0, max =4),
    M37_fc= c(min =1, max =2000),
    M37_lp= c(min =0.05, max =0.95),
    M37_beta= c(min =0, max =10),
    M37_k0= c(min =0, max =1),
    M37_alpha= c(min =0, max =4),
    M37_perc= c(min =0, max =20),
    M37_K1= c(min =0, max =1),
    M37_maxbas= c(min =1, max =120),
    M37_s1= c(min =0, max =2000),
    M37_s2= c(min =0, max =1),
    M37_s3= c(min =0, max =1),
    M37_s4= c(min =1, max =2000),
    M37_s5= c(min =1, max =2000),
    
    #M38 "M38_a0","M38_b0","M38_c0","M38_a1","M38_fa","M38_fb","M38_fc","M38_fd","M38_st","M38_f2","M38_f1","M38_f3","M38_k1","M38_K2","M38_z1","M38_z2","M38_s1","M38_s2","M38_s3","M38_s4","M38_s5"
    M38_a0= c(min =0, max =1),
    M38_b0= c(min =0, max =1),
    M38_c0= c(min =0, max =1),
    M38_a1= c(min =0, max =1),
    M38_fa= c(min =0, max =1),
    M38_fb= c(min =0, max =1),
    M38_fc= c(min =0, max =1),
    M38_fd= c(min =0, max =1),
    M38_st= c(min =1, max =2000),
    M38_f2= c(min =0.01, max =0.99),
    M38_f1= c(min =0.01, max =0.99),
    M38_f3= c(min =0.01, max =0.99),
    M38_k1= c(min =0, max =4),
    M38_k2= c(min =0, max =4),
    M38_z1= c(min =0.01, max =0.99),
    M38_z2= c(min =0.01, max =0.99),
    M38_s1= c(min =0, max =1),
    M38_s2= c(min =0, max =1),
    M38_s3= c(min =0, max =1),
    M38_s4= c(min =0, max =1),
    M38_s5= c(min =0, max =1),
    
    #M39 "M39_smax","M39_cmax","M39_ct","M39_c1","M39_ce","M39_dsurp","M39_kd","M39_gamd","M39_qpmax","M39_kg","M39_tau","M39_sbf","M39_kcr","M39_gamcr","M39_kor","M39_gamor","M39_s1","M39_s2","M39_s3","M39_s4","M39_s5"
    M39_smax= c(min =0, max =5),
    M39_cmax= c(min =0.01, max =0.99),
    M39_ct= c(min =0.01, max =0.99),
    M39_c1= c(min =0, max =2),
    M39_ce= c(min =0, max =1),
    M39_dsurp= c(min =1, max =2000),
    M39_kd= c(min =0, max =1),
    M39_gamd= c(min =1, max =5),
    M39_qpmax= c(min =0, max =20),
    M39_kg= c(min =0, max =1),
    M39_tau= c(min =1, max =120),
    M39_sbf= c(min =1, max =300),
    M39_kcr= c(min =0, max =1),
    M39_gamcr= c(min =1, max =5),
    M39_kor= c(min =0, max =1),
    M39_gamor= c(min =1, max =5),
    M39_s1= c(min =0, max =1),
    M39_s2= c(min =0, max =1),
    M39_s3= c(min =1, max =2000),
    M39_s4= c(min =0, max =1),
    M39_s5= c(min =0, max =1),
    
    #M40 "M40_h","M40_y","M40_smax","M40_c","M40_g","M40_kg","M40_n","M40_nk","M40_s1","M40_s2","M40_s3","M40_s4","M40_s5","M40_s6"
    M40_h= c(min =0, max =1),
    M40_y= c(min =0, max =200),
    M40_smax= c(min =1, max =2000),
    M40_c= c(min =0, max =1),
    M40_g= c(min =0, max =1),
    M40_kg= c(min =0, max =1),
    M40_n= c(min =1, max =10),
    M40_nk= c(min =1, max =120),
    M40_s1= c(min =0, max =1),
    M40_s2= c(min =0, max =1),
    M40_s3= c(min =0, max =1),
    M40_s4= c(min =0, max =1),
    M40_s5= c(min =0, max =1),
    M40_s6= c(min =1, max =2000),
    
    #41 "M41_cs","M41_cif","M41_stot","M41_cl1","M41_f1","M41_cof","M41_cl2","M41_k0","M41_k1","M41_kb","M41_s1","M41_s2","M41_s3","M41_s4","M41_s5","M41_s6"
    M41_cs= c(min =0, max =20),
    M41_cif= c(min =0, max =1),
    M41_stot= c(min =1, max =2000),
    M41_cl1= c(min =0, max =0.99),
    M41_f1= c(min =0.01, max =0.99),
    M41_cof= c(min =0, max =1),
    M41_cl2= c(min =0, max =0.99),
    M41_k0= c(min =0, max =1),
    M41_k1= c(min =0, max =1),
    M41_kb= c(min =0, max =1),
    M41_s1= c(min =1, max =2000),
    M41_s2= c(min =0, max =1),
    M41_s3= c(min =0, max =1),
    M41_s4= c(min =1, max =300),
    M41_s5= c(min =1, max =300),
    M41_s6= c(min =1, max =2000),
    
    #42 M42_c","M42_imax","M42_a","M42_fi2","M42_kin","M42_d50","M42_fd16","M42_sbc","M42_kb","M42_pb","M42_kh","M42_kc","M42_s1","M42_s2","M42_s3","M42_s4","M42_s5","M42_s6"
    M42_c= c(min =0, max =1),
    M42_imax= c(min =0, max =5),
    M42_a= c(min =0, max =1),
    M42_fi2= c(min =0.01, max =0.99),
    M42_kin= c(min =0, max =1),
    M42_d50= c(min =1, max =2000),
    M42_fd16= c(min =0.01, max =0.99),
    M42_sbc= c(min =1, max =2000),
    M42_kb= c(min =0, max =1),
    M42_pb= c(min =1, max =5),
    M42_kh= c(min =0, max =1),
    M42_kc= c(min =0, max =1),
    M42_s1= c(min =0, max =1),
    M42_s2= c(min =0, max =1),
    M42_s3= c(min =1, max =2000),
    M42_s4= c(min =1, max =300),
    M42_s5= c(min =1, max =300),
    M42_s6= c(min =1, max =300),
    
    #43 "M43_fice","M43_t0","M43_asnow","M43_tm","M43_ks","M43_aice","M43_ki","M43_a","M43_x","M43_y","M43_ksl","M43_beta","M43_s1","M43_s2","M43_s3","M43_s4","M43_s5","M43_s6"
    M43_fice= c(min =0, max =1),
    M43_t0= c(min =(-3), max =5),
    M43_asnow= c(min =0, max =20),
    M43_tm= c(min =(-3), max =3),
    M43_ks= c(min =0, max =1),
    M43_aice= c(min =0, max =20),
    M43_ki= c(min =0, max =1),
    M43_a= c(min =1, max =2000),
    M43_x= c(min =0, max =10),
    M43_y= c(min =0, max =5),
    M43_ksl= c(min =0, max =1),
    M43_beta= c(min =0, max =1),
    M43_s1= c(min =1, max =2000),
    M43_s2= c(min =1, max =2000),
    M43_s3= c(min =1, max =2000),
    M43_s4= c(min =1, max =2000),
    M43_s5= c(min =0, max =1),
    M43_s6= c(min =1, max =300),
    
    #M44 "M44_rho","M44_ts","M44_tm","M44_as","M44_af","M44_gmax","M44_the","M44_phi","M44_smax","M44_fsm","M44_fsw","M44_ksat","M44_c","M44_lmax","M44_kf","M44_ks","M44_s1","M44_s2","M44_s3","M44_s4","M44_s5","M44_s6"
    M44_rho= c(min =0, max =5),
    M44_ts= c(min =(-3), max =5),
    M44_tm= c(min =(-3), max =3),
    M44_as= c(min =0, max =20),
    M44_af= c(min =0, max =1),
    M44_gmax= c(min =0, max =2),
    M44_the= c(min =0, max =1),
    M44_phi= c(min =0, max =200),
    M44_smax= c(min =1, max =2000),
    M44_fsm= c(min =0.05, max =0.95),
    M44_fsw= c(min =0.05, max =0.95),
    M44_ksat= c(min =0, max =1),
    M44_c= c(min =0, max =5),
    M44_lmax= c(min =0, max =20),
    M44_kf= c(min =0, max =1),
    M44_ks= c(min =0, max =1),
    M44_s1= c(min =0, max =1),
    M44_s2= c(min =1, max =2000),
    M44_s3= c(min =1, max =2000),
    M44_s4= c(min =0, max =1),
    M44_s5= c(min =1, max =300),
    M44_s6= c(min =1, max =300),
    
    #M45 "M45_tt","M45_ddf","M45_alpha","M45_beta","M45_stor","M45_retip","M45_fscn","M45_scx","M45_flz","M45_stot","M45_cgw","M45_resmax","M45_k1","M45_k2","M45_k3","M45_k4","M45_k5","M45_k6","M45_s1","M45_s2","M45_s3","M45_s4","M45_s5","M45_s6","M45_s7"
    M45_tt= c(min =(-3), max =5),
    M45_ddf= c(min =0, max =20),
    M45_alpha= c(min =0, max =1),
    M45_beta= c(min =0, max =1),
    M45_stor= c(min =0, max =5),
    M45_retip= c(min =0, max =50),
    M45_fscn= c(min =0, max =1),
    M45_scx= c(min =0, max =1),
    M45_flz= c(min =0.005, max =0.995),
    M45_stot= c(min =1, max =2000),
    M45_cgw= c(min =0, max =20),
    M45_resmax= c(min =1, max =300),
    M45_k1= c(min =0, max =1),
    M45_k2= c(min =1, max =5),
    M45_k3= c(min =0, max =1),
    M45_k4= c(min =0, max =1),
    M45_k5= c(min =0, max =1),
    M45_k6= c(min =0, max =1),
    M45_s1= c(min =0, max =2000),
    M45_s2= c(min =0, max =1),
    M45_s3= c(min =0, max =1),
    M45_s4= c(min =0, max =1),
    M45_s5= c(min =0, max =1),
    M45_s6= c(min =0, max =1),
    M45_s7= c(min =1, max =2000),
    
    #M46 "M46_fap","M46_fdp","M46_dp","M46_cq","M46_d1","M46_tf","M46_fds","M46_ds","M46_d2","M46_cxq","M46_cxs","M46_cu","M46_s1","M46_s2","M46_s3","M46_s4","M46_s5","M46_s6","M46_s7","M46_s8"
    M46_fap= c(min =0, max =1),
    M46_fdp= c(min =0.01, max =0.99),
    M46_dp= c(min =1, max =2000),
    M46_cq= c(min =0, max =1),
    M46_d1= c(min =0, max =1),
    M46_tf= c(min =0, max =1),
    M46_fds= c(min =0.01, max =0.99),
    M46_ds= c(min =1, max =2000),
    M46_d2= c(min =0, max =1),
    M46_cxq= c(min =0, max =1),
    M46_cxs= c(min =0, max =1),
    M46_cu= c(min =0, max =1),
    M46_s1= c(min =0, max =1),
    M46_s2= c(min =0, max =1),
    M46_s3= c(min =1, max =300),
    M46_s4= c(min =0, max =1),
    M46_s5= c(min =0, max =1),
    M46_s6= c(min =1, max =300),
    M46_s7= c(min =1, max =300),
    M46_s8= c(min =1, max =300)
    
    
  )
  
  SFPR
  
}

CEDprecomputation= function(x){
  # x are measured data
  
  optimum_bins = function(data){
    
    minM = 1
    maxM = 1000
    N = length(data)
    logp = rep(1, maxM)
    
    for (M in minM:maxM){
      
      w = hist(data, seq(min(data), max(data),length.out = M + 1), plot=FALSE) 
      n = w[['counts']]# Bin the data (equal width bins here)
      p = 0
      
      for (k in 1:M){
        
        p = p + lgamma(n[k] + 0.5)
        
      }
      
      logp[M] = N * log(M) + lgamma(M / 2) - M * lgamma(1 / 2) - lgamma(N + M / 2) + p
      
    }
    
    optBins = which.max(logp)
    
    optBins
    
  }
  
  Qobs_d= sort(x, decreasing = TRUE)
  
  range_high         = c( 0,   2)
  range_medium       = c( 2,  20)
  range_intermediate = c(20,  70)
  range_low          = c(70,  100)
  lenQo = length(Qobs_d)
  
  h1 = ceiling(range_high[1] / 100 * lenQo) + 1
  h2 = ceiling(range_high[2] / 100 * lenQo)
  # 
  m1 = ceiling(range_medium[1] / 100 * lenQo) + 1
  m2 = ceiling(range_medium[2] / 100 * lenQo)
  # 
  i1 = ceiling(range_intermediate[1] / 100 * lenQo) + 1
  i2 = ceiling(range_intermediate[2] / 100 * lenQo)
  # 
  l1 = ceiling(range_low[1] / 100 * lenQo) + 1 
  l2 = ceiling(range_low[2] / 100 * lenQo)
  # 
  
  Qobs_r_high = Qobs_d[h1:h2]
  Qobs_r_medium = Qobs_d[m1:m2]
  Qobs_r_intermediate = Qobs_d[i1:i2]
  Qobs_r_low = Qobs_d[l1:l2] 
  # 
  optBins_high         = optimum_bins(Qobs_r_high)
  optBins_medium       = optimum_bins(Qobs_r_medium)
  optBins_intermediate = optimum_bins(Qobs_r_intermediate)
  optBins_low          = optimum_bins(Qobs_r_low)
  
  optBins = optimum_bins(Qobs_d)
  # 
  Qobs2 =quantile(Qobs_d,0.98) # 2 percentile of observed flows
  Qobs20=quantile(Qobs_d,0.80)  # 20 percentile of observed flows
  Qobs70=quantile(Qobs_d,0.30)  # 70 percentile of observed flows
  # 
  Qobs_min=min(Qobs_d)
  Qobs_max=max(Qobs_d)
  # 
  list(OptBins = optBins, OptBins_high = optBins_high, OptBins_medium = optBins_medium, OptBins_intermediate = optBins_intermediate, OptBins_low = optBins_low,
       lenQo=lenQo, Qobs2=Qobs2, Qobs20=Qobs20, Qobs70=Qobs70, Qobs_min=Qobs_min, Qobs_max=Qobs_max)
  
}

MeaningfulCTcombs = function(){
  
  c('Zero Tank',        #0,0,0,0,0,0,0,
    'One Tank - CR',      #0,0,0,0,0,0,1,
    'One Tank - SR',      #0,0,0,0,0,1,0,
    'Two Tanks - SR,CR',  #0,0,0,0,0,1,1,
    'One Tank - FR',      #0,0,0,0,1,0,0,
    'Two Tanks - FR,CR',  #0,0,0,0,1,0,1,
    'Two Tanks - FR,SR',  #0,0,0,0,1,1,0,
    'Three Tanks - FR,SR,CR',  #0,0,0,0,1,1,1,
    'One Tank - UR',           #0,0,0,1,0,0,0,
    'Two Tanks - UR,CR', #0,0,0,1,0,0,1,
    'Two Tanks - UR,SR', #0,0,0,1,0,1,0,
    'Three Tanks - UR,SR,CR', #0,0,0,1,0,1,1,
    'Two Tanks - UR,FR', #0,0,0,1,1,0,0,
    'Three Tanks - UR,FR,CR', #0,0,0,1,1,0,1,
    'Three Tanks - UR,FR,SR', #0,0,0,1,1,1,0,
    'Four Tanks - UR,FR,SR,CR', #0,0,0,1,1,1,1,
    'One Tank - RR',            #0,0,1,0,0,0,0,
    'Two Tanks - RR,CR', #0,0,1,0,0,0,1,
    'Two Tanks - RR,SR', #0,0,1,0,0,1,0,
    'Three Tanks - RR,SR,CR', #0,0,1,0,0,1,1,
    'Two Tanks - RR,FR', #0,0,1,0,1,0,0,
    'Three Tanks - RR,FR,CR', #0,0,1,0,1,0,1,
    'Three Tanks - RR,FR,SR', #0,0,1,0,1,1,0,
    'Four Tanks - RR,FR,SR,CR', #0,0,1,0,1,1,1,
    'Two Tanks - RR,UR', #0,0,1,1,0,0,0,
    'Three Tanks - RR,UR,CR', #0,0,1,1,0,0,1,
    'Three Tanks - RR,UR,SR', #0,0,1,1,0,1,0,
    'Four Tanks - RR,UR,SR,CR', #0,0,1,1,0,1,1,
    'Three Tanks - RR,UR,FR', #0,0,1,1,1,0,0,
    'Four Tanks - RR,UR,FR,CR', #0,0,1,1,1,0,1,
    'Four Tanks - RR,UR,FR,SR', #0,0,1,1,1,1,0,
    'Five Tanks - RR,UR,FR,SR,CR', #0,0,1,1,1,1,1,
    'One Tank - IR',               #0,1,0,0,0,0,0,
    'Two Tanks - IR,CR', #0,1,0,0,0,0,1,
    'Two Tanks - IR,SR', #0,1,0,0,0,1,0,
    'Three Tanks - IR,SR,CR', #0,1,0,0,0,1,1,
    'Two Tanks - IR,FR', #0,1,0,0,1,0,0,
    'Three Tanks - IR,FR,CR', #0,1,0,0,1,0,1,
    'Three Tanks - IR,FR,SR', #0,1,0,0,1,1,0,
    'Four Tanks - IR,FR,SR,CR', #0,1,0,0,1,1,1,
    'Two Tanks - IR,UR', #0,1,0,1,0,0,0,
    'Three Tanks - IR,UR,CR', #0,1,0,1,0,0,1,
    'Three Tanks - IR,UR,SR', #0,1,0,1,0,1,0,
    'Four Tanks - IR,UR,SR,CR', #0,1,0,1,0,1,1,
    'Three Tanks - IR,UR,FR', #0,1,0,1,1,0,0,
    'Four Tanks - IR,UR,FR,CR', #0,1,0,1,1,0,1,
    'Four Tanks - IR,UR,FR,SR', #0,1,0,1,1,1,0,
    'Five Tanks - IR,UR,FR,SR,CR', #0,1,0,1,1,1,1,
    'Two Tanks - IR,RR', #0,1,1,0,0,0,0,
    'Three Tanks - IR,RR,CR', #0,1,1,0,0,0,1,
    'Three Tanks - IR,RR,SR', #0,1,1,0,0,1,0,
    'Four Tanks - IR,RR,SR,CR', #0,1,1,0,0,1,1,
    'Three Tanks - IR,RR,FR', #0,1,1,0,1,0,0,
    'Four Tanks - IR,RR,FR,CR', #0,1,1,0,1,0,1,
    'Four Tanks - IR,RR,FR,SR', #0,1,1,0,1,1,0,
    'Five Tanks - IR,RR,FR,SR,CR', #0,1,1,0,1,1,1,
    'Three Tanks - IR,RR,UR', #0,1,1,1,0,0,0,
    'Four Tanks - IR,RR,UR,CR', #0,1,1,1,0,0,1,
    'Four Tanks - IR,RR,UR,SR', #0,1,1,1,0,1,0,
    'Five Tanks - IR,RR,UR,SR,CR', #0,1,1,1,0,1,1,
    'Four Tanks - IR,RR,UR,FR', #0,1,1,1,1,0,0,
    'Five Tanks - IR,RR,UR,FR,CR', #0,1,1,1,1,0,1,
    'Five Tanks - IR,RR,UR,FR,SR', #0,1,1,1,1,1,0,
    'Six Tanks - IR,RR,UR,FR,SR,CR', #0,1,1,1,1,1,1,
    'One Tank - WR',                 #1,0,0,0,0,0,0,
    'Two Tanks - WR,CR', #1,0,0,0,0,0,1,
    'Two Tanks - WR,SR', #1,0,0,0,0,1,0,
    'Three Tanks - WR,SR,CR', #1,0,0,0,0,1,1,
    'Two Tanks - WR,FR', #1,0,0,0,1,0,0,
    'Three Tanks - WR,FR,CR', #1,0,0,0,1,0,1,
    'Three Tanks - WR,FR,SR', #1,0,0,0,1,1,0,
    'Four Tanks - WR,FR,SR,CR', #1,0,0,0,1,1,1,
    'Two Tanks - WR,UR', #1,0,0,1,0,0,0,
    'Three Tanks - WR,UR,CR', #1,0,0,1,0,0,1,
    'Three Tanks - WR,UR,SR', #1,0,0,1,0,1,0,
    'Four Tanks - IR,UR,SR,CR', #1,0,0,1,0,1,1,
    'Three Tanks - WR,UR,FR', #1,0,0,1,1,0,0,
    'Four Tanks - WR,UR,FR,CR', #1,0,0,1,1,0,1,
    'Four Tanks - WR,UR,FR,SR', #1,0,0,1,1,1,0,
    'Five Tanks - WR,UR,FR,SR,CR', #1,0,0,1,1,1,1,
    'Two Tanks - WR,RR', #1,0,1,0,0,0,0,
    'Three Tanks - WR,RR,CR', #1,0,1,0,0,0,1,
    'Three Tanks - WR,RR,SR', #1,0,1,0,0,1,0,
    'Four Tanks - WR,RR,SR,CR', #1,0,1,0,0,1,1,
    'Three Tanks - WR,RR,FR', #1,0,1,0,1,0,0,
    'Four Tanks - WR,RR,FR,CR', #1,0,1,0,1,0,1,
    'Four Tanks - WR,RR,FR,SR', #1,0,1,0,1,1,0,
    'Five Tanks - WR,RR,FR,SR,CR', #1,0,1,0,1,1,1,
    'Three Tanks - WR,RR,UR', #1,0,1,1,0,0,0,
    'Four Tanks - WR,RR,UR,CR', #1,0,1,1,0,0,1,
    'Four Tanks - WR,RR,UR,SR', #1,0,1,1,0,1,0,
    'Five Tanks - WR,RR,UR,SR,CR', #1,0,1,1,0,1,1,
    'Four Tanks - WR,RR,UR,FR', #1,0,1,1,1,0,0,
    'Five Tanks - WR,RR,UR,FR,CR', #1,0,1,1,1,0,1,
    'Five Tanks - WR,RR,UR,FR,SR', #1,0,1,1,1,1,0,
    'Six Tanks - WR,RR,UR,FR,SR,CR', #1,0,1,1,1,1,1,
    'Two Tanks - WR,IR', #1,1,0,0,0,0,0,
    'Three Tanks - WR,IR,CR', #1,1,0,0,0,0,1,
    'Three Tanks - WR,IR,SR', #1,1,0,0,0,1,0,
    'Four Tanks - WR,IR,SR,CR', #1,1,0,0,0,1,1,
    'Three Tanks - WR,IR,FR', #1,1,0,0,1,0,0,
    'Four Tanks - WR,IR,FR,CR', #1,1,0,0,1,0,1,
    'Four Tanks - WR,IR,FR,SR', #1,1,0,0,1,1,0,
    'Five Tanks - WR,IR,FR,SR,CR', #1,1,0,0,1,1,1,
    'Three Tanks - WR,IR,UR', #1,1,0,1,0,0,0,
    'Four Tanks - WR,IR,UR,CR', #1,1,0,1,0,0,1,
    'Four Tanks - WR,IR,UR,SR', #1,1,0,1,0,1,0,
    'Five Tanks - WR,IR,UR,SR,CR', #1,1,0,1,0,1,1,
    'Four Tanks - WR,IR,UR,FR', #1,1,0,1,1,0,0,
    'Five Tanks - WR,IR,UR,FR,CR', #1,1,0,1,1,0,1,
    'Five Tanks - WR,IR,UR,FR,SR', #1,1,0,1,1,1,0,
    'Six Tanks - WR,IR,UR,FR,SR,CR', #1,1,0,1,1,1,1,
    'Three Tanks - WR,IR,RR', #1,1,1,0,0,0,0,
    'Four Tanks - WR,IR,RR,CR', #1,1,1,0,0,0,1,
    'Four Tanks - WR,IR,RR,SR', #1,1,1,0,0,1,0,
    'Five Tanks - WR,IR,RR,SR,CR', #1,1,1,0,0,1,1,
    'Four Tanks - WR,IR,RR,FR', #1,1,1,0,1,0,0,
    'Five Tanks - WR,IR,RR,FR,CR', #1,1,1,0,1,0,1,
    'Five Tanks - WR,IR,RR,FR,SR', #1,1,1,0,1,1,0,
    'Six Tanks - WR,IR,RR,FR,SR,CR', #1,1,1,0,1,1,1,
    'Four Tanks - WR,IR,RR,UR', #1,1,1,1,0,0,0,
    'Five Tanks - WR,IR,RR,UR,CR', #1,1,1,1,0,0,1,
    'Five Tanks - WR,IR,RR,UR,SR', #1,1,1,1,0,1,0,
    'Six Tanks - WR,IR,RR,UR,SR,CR', #1,1,1,1,0,1,1,
    'Five Tanks - WR,IR,RR,UR,FR', #1,1,1,1,1,0,0,
    'Six Tanks - WR,IR,RR,UR,FR,CR', #1,1,1,1,1,0,1,
    'Six Tanks - WR,IR,RR,UR,FR,SR', #1,1,1,1,1,1,0,
    'Seven Tanks - WR,IR,RR,UR,FR,SR,CR' #1,1,1,1,1,1,1 
  )
  
}
MeaningfulCTcombsFun = function(){
  
  mCTmatrix = matrix(c(          0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,1,
                                 0,0,0,0,0,1,0,
                                 0,0,0,0,0,1,1,
                                 0,0,0,0,1,0,0,
                                 0,0,0,0,1,0,1,
                                 0,0,0,0,1,1,0,
                                 0,0,0,0,1,1,1,
                                 0,0,0,1,0,0,0,
                                 0,0,0,1,0,0,1,
                                 0,0,0,1,0,1,0,
                                 0,0,0,1,0,1,1,
                                 0,0,0,1,1,0,0,
                                 0,0,0,1,1,0,1,
                                 0,0,0,1,1,1,0,
                                 0,0,0,1,1,1,1,
                                 0,0,1,0,0,0,0,
                                 0,0,1,0,0,0,1,
                                 0,0,1,0,0,1,0,
                                 0,0,1,0,0,1,1,
                                 0,0,1,0,1,0,0,
                                 0,0,1,0,1,0,1,
                                 0,0,1,0,1,1,0,
                                 0,0,1,0,1,1,1,
                                 0,0,1,1,0,0,0,
                                 0,0,1,1,0,0,1,
                                 0,0,1,1,0,1,0,
                                 0,0,1,1,0,1,1,
                                 0,0,1,1,1,0,0,
                                 0,0,1,1,1,0,1,
                                 0,0,1,1,1,1,0,
                                 0,0,1,1,1,1,1,
                                 0,1,0,0,0,0,0,
                                 0,1,0,0,0,0,1,
                                 0,1,0,0,0,1,0,
                                 0,1,0,0,0,1,1,
                                 0,1,0,0,1,0,0,
                                 0,1,0,0,1,0,1,
                                 0,1,0,0,1,1,0,
                                 0,1,0,0,1,1,1,
                                 0,1,0,1,0,0,0,
                                 0,1,0,1,0,0,1,
                                 0,1,0,1,0,1,0,
                                 0,1,0,1,0,1,1,
                                 0,1,0,1,1,0,0,
                                 0,1,0,1,1,0,1,
                                 0,1,0,1,1,1,0,
                                 0,1,0,1,1,1,1,
                                 0,1,1,0,0,0,0,
                                 0,1,1,0,0,0,1,
                                 0,1,1,0,0,1,0,
                                 0,1,1,0,0,1,1,
                                 0,1,1,0,1,0,0,
                                 0,1,1,0,1,0,1,
                                 0,1,1,0,1,1,0,
                                 0,1,1,0,1,1,1,
                                 0,1,1,1,0,0,0,
                                 0,1,1,1,0,0,1,
                                 0,1,1,1,0,1,0,
                                 0,1,1,1,0,1,1,
                                 0,1,1,1,1,0,0,
                                 0,1,1,1,1,0,1,
                                 0,1,1,1,1,1,0,
                                 0,1,1,1,1,1,1,
                                 1,0,0,0,0,0,0,
                                 1,0,0,0,0,0,1,
                                 1,0,0,0,0,1,0,
                                 1,0,0,0,0,1,1,
                                 1,0,0,0,1,0,0,
                                 1,0,0,0,1,0,1,
                                 1,0,0,0,1,1,0,
                                 1,0,0,0,1,1,1,
                                 1,0,0,1,0,0,0,
                                 1,0,0,1,0,0,1,
                                 1,0,0,1,0,1,0,
                                 1,0,0,1,0,1,1,
                                 1,0,0,1,1,0,0,
                                 1,0,0,1,1,0,1,
                                 1,0,0,1,1,1,0,
                                 1,0,0,1,1,1,1,
                                 1,0,1,0,0,0,0,
                                 1,0,1,0,0,0,1,
                                 1,0,1,0,0,1,0,
                                 1,0,1,0,0,1,1,
                                 1,0,1,0,1,0,0,
                                 1,0,1,0,1,0,1,
                                 1,0,1,0,1,1,0,
                                 1,0,1,0,1,1,1,
                                 1,0,1,1,0,0,0,
                                 1,0,1,1,0,0,1,
                                 1,0,1,1,0,1,0,
                                 1,0,1,1,0,1,1,
                                 1,0,1,1,1,0,0,
                                 1,0,1,1,1,0,1,
                                 1,0,1,1,1,1,0,
                                 1,0,1,1,1,1,1,
                                 1,1,0,0,0,0,0,
                                 1,1,0,0,0,0,1,
                                 1,1,0,0,0,1,0,
                                 1,1,0,0,0,1,1,
                                 1,1,0,0,1,0,0,
                                 1,1,0,0,1,0,1,
                                 1,1,0,0,1,1,0,
                                 1,1,0,0,1,1,1,
                                 1,1,0,1,0,0,0,
                                 1,1,0,1,0,0,1,
                                 1,1,0,1,0,1,0,
                                 1,1,0,1,0,1,1,
                                 1,1,0,1,1,0,0,
                                 1,1,0,1,1,0,1,
                                 1,1,0,1,1,1,0,
                                 1,1,0,1,1,1,1,
                                 1,1,1,0,0,0,0,
                                 1,1,1,0,0,0,1,
                                 1,1,1,0,0,1,0,
                                 1,1,1,0,0,1,1,
                                 1,1,1,0,1,0,0,
                                 1,1,1,0,1,0,1,
                                 1,1,1,0,1,1,0,
                                 1,1,1,0,1,1,1,
                                 1,1,1,1,0,0,0,
                                 1,1,1,1,0,0,1,
                                 1,1,1,1,0,1,0,
                                 1,1,1,1,0,1,1,
                                 1,1,1,1,1,0,0,
                                 1,1,1,1,1,0,1,
                                 1,1,1,1,1,1,0,
                                 1,1,1,1,1,1,1),  
                     ncol = 7, byrow = TRUE)
  
  mCTdf = as.data.frame(mCTmatrix)
  names(mCTdf) = c('w', 'i', 'r', 'u', 'f', 's','c')
  mCTlst = split(mCTdf, 1:nrow(mCTdf))
  lapply(mCTlst, unlist)
  
}

DepArgsOfCTfun = function(){
  
  CTargsNames = names(formals(CT))
  list(WR = grep('_WR', CTargsNames, value = TRUE),
       IR = grep('_IR', CTargsNames, value = TRUE),
       RR = grep('_RR', CTargsNames, value = TRUE),
       UR = grep('_UR', CTargsNames, value = TRUE),
       FR = grep('_FR', CTargsNames, value = TRUE),
       SR = grep('_SR', CTargsNames, value = TRUE),
       CR = grep('_uCR', CTargsNames, value = TRUE),
       CR = grep('_sCR', CTargsNames, value = TRUE),
       CR = grep('_CR', CTargsNames, value = TRUE))
  
}

FitnessGlobal = function(referenceValues, assignIt = FALSE){
  
  lowValPositions = referenceValues <= quantile(referenceValues, probs = 0.2, na.rm = TRUE)
  highValPositions = referenceValues >= quantile(referenceValues, probs = 0.7, na.rm = TRUE) # this can be done outside of this function for faster run
  mediumValPositions = !(lowValPositions|highValPositions)
  logOriginalValues = Qlog(referenceValues + quantile(referenceValues[referenceValues > 0], 0.1))
  MeanOriginalValues = mean(referenceValues)
  indexm1=which.min(abs(referenceValues-(as.numeric(quantile(referenceValues, probs = 0.2, na.rm = TRUE)))))
  originalm1=referenceValues[indexm1]
  indexm2=which.min(abs(referenceValues-(as.numeric(quantile(referenceValues, probs = 0.7, na.rm = TRUE)))))
  originalm2=referenceValues[indexm2]
  originall2=as.numeric(quantile(referenceValues, probs = 0, na.rm = TRUE))
  orlength=length(referenceValues)
  lag_num=numeric(orlength-1)
  
  for (i in 1:(orlength-1))
  {
    lag_num[i]=(referenceValues[i]-MeanOriginalValues)*(referenceValues[i+1]-MeanOriginalValues)
  }
  fdcmOriginalValues=Qlog10(originalm2)-Qlog10(originalm1)  
  OriginalValuesLow= referenceValues[lowValPositions] 
  OriginalValuesHigh=referenceValues[highValPositions] 
  fdchOriginalValues=sum(OriginalValuesHigh)
  fdclOriginalValues=-1*(sum(Qlog10(OriginalValuesLow)-Qlog10(originall2)))  
  Lag_ac_Original=sum(lag_num)/sum((referenceValues-MeanOriginalValues)^2)  
  VarOriginalValues=sqrt(sum((referenceValues-MeanOriginalValues)^2)/(orlength-1)) 
  MeanlogOriginalValues=mean(logOriginalValues) 
  VarlogOriginalValues=sqrt(sum((logOriginalValues-MeanlogOriginalValues)^2)/(orlength-1)) 
  MedianOriginalValues=median(referenceValues) 
  PeakOriginalValues= max(referenceValues)
  
  if(assignIt){
    
    assign('MultiObjectiveFitnesses', c('Multi_Madsen', 'Dawson', 'Vis_2')    , envir = parent.env(environment()))
    assign('LowValPositions'     , lowValPositions                            , envir = parent.env(environment()))
    assign('HighValPositions'    , highValPositions                           , envir = parent.env(environment()))
    assign('MediumValPositions'  , mediumValPositions                         , envir = parent.env(environment()))
    assign('OriginalValuesLow'   , OriginalValuesLow                          , envir = parent.env(environment()))
    assign('OriginalValuesHigh'  , OriginalValuesHigh                         , envir = parent.env(environment()))
    assign('OriginalValuesMedium', referenceValues[mediumValPositions]        , envir = parent.env(environment()))
    assign('NS0denominator'      , referenceValues - mean(referenceValues)    , envir = parent.env(environment()))
    assign('logOriginalValues'   , logOriginalValues                          , envir = parent.env(environment()))
    assign('logNS0denominator'   , logOriginalValues - mean(logOriginalValues), envir = parent.env(environment()))
    assign('MeanOriginalValues'  , mean(referenceValues)                      , envir = parent.env(environment()))
    assign('originalm1'          , originalm1                                 , envir = parent.env(environment()))
    assign('originalm2'          , originalm2                                 , envir = parent.env(environment()))
    assign('originall2'          , originall2                                 , envir = parent.env(environment()))
    assign('orlength'            , orlength                                   , envir = parent.env(environment()))
    assign('lag_num'             , lag_num                                    , envir = parent.env(environment()))
    assign('fdcmOriginalValues'  , fdcmOriginalValues                         , envir = parent.env(environment()))
    assign('fdchOriginalValues'  , fdchOriginalValues                         , envir = parent.env(environment()))
    assign('fdclOriginalValues'  , fdclOriginalValues                         , envir = parent.env(environment()))
    assign('MedianOriginalValues' , MedianOriginalValues                      , envir = parent.env(environment()))
    assign('PeakOriginalValues'   , PeakOriginalValues                        , envir = parent.env(environment()))
    assign('Lag_ac_Original'      , Lag_ac_Original                           , envir = parent.env(environment()))
    assign('VarOriginalValues'    , VarOriginalValues                         , envir = parent.env(environment()))
    assign('MeanlogOriginalValues', MeanlogOriginalValues                     , envir = parent.env(environment())) 
    assign('VarlogOriginalValues' , VarlogOriginalValues                      , envir = parent.env(environment())) 
    
    # CT globals
    DepArgsOfCT = DepArgsOfCTfun()
    assign('DepArgsOfCT', DepArgsOfCT, envir = parent.env(environment()), inherits = TRUE)
    
    MeaningfullCTcombinations = MeaningfulCTcombsFun()
    assign('MeaningfullCTcombinations', MeaningfullCTcombinations, envir = parent.env(environment()), inherits = TRUE)   
    
    # CEDs globals
    if(FitnessFunction == 'SUSE'|FitnessFunction == 'CED'|FitnessFunction == 'CED_new'){
      
      ceds = CEDprecomputation(referenceValues)
      sapply(names(ceds), function(nm) assign(nm, ceds[[nm]], envir = parent.env(environment()), inherits = TRUE))
      
    }
    
  }
  
  list(lowValPositions = lowValPositions, highValPositions = highValPositions, 
       mediumValPositions = mediumValPositions, logOriginalValues = logOriginalValues, MeanOriginalValues = MeanOriginalValues,
       originalm1 = originalm1, originalm2 = originalm2, originall2 = originall2, orlength = orlength, lag_num = lag_num,
       fdcmOriginalValues=fdcmOriginalValues, fdchOriginalValues=fdchOriginalValues, fdclOriginalValues=fdclOriginalValues,
       MedianOriginalValues=MedianOriginalValues, PeakOriginalValues=PeakOriginalValues,Lag_ac_Original=Lag_ac_Original,
       VarOriginalValues=VarOriginalValues, MeanlogOriginalValues=MeanlogOriginalValues,VarlogOriginalValues=VarlogOriginalValues)
  
}


SetGlobalVariables = function(functionSet, indepVar, constRng, StrictFunsArgsRng, roundFactor, dataSet, depVar, fitFun, weights, punish1nodeInd, maxDepthRun){
  
  nArgs = function(...){
    
    sapply(list(...) , function(x) length(formals(x)))
    
  }
  
  # StrictFunctions is vector of argument-strict functions
  StrictFunctions = c( "MI", "MII", "MIII", "MIV", "MV", "MVI", "MVII", "MVIII", "MIX", "MX","MXI","MXII",
                       "CT","FUSE","DCT","CT_new","WR","IR","RR","UR","FR","SR","CR","MARRMot") # EDIT BY ADDING NEW FUNCTION!!!
  NonStrictArguments = c( MI = 0, MII = 0, MIII = 0, MIV = 0, MV = 0, MVI = 0, 
                         MVII = 0, MVIII = 0, MIX = 0, MX = 0, MXI=0, MXII = 0, CT=0, FUSE =0, DCT=0, CT_new=0,WR = 0, IR = 0, RR = 0, UR = 0, FR = 0, SR = 0, CR = 0, MARRMot = 0)
  SpecialFunctions = c("SMA", "RES", "DLY", "SSUM", StrictFunctions)
  
  functionsDefinition = data.frame(
    
    Functions = c("sqrt", "log", "exp", "log10", "sin", "cos", "tan", "tanh", "Hstep", "floor", "ceiling",
                  "+", "-", "*", "/", "^", "max", "min", "Gr", "GrEq", "Eq", "Pdiv", SpecialFunctions),
    
    Arity = c(rep(1, 11), rep(2, 15), 
              nArgs(MI, MII, MIII, MIV, MV, MVI, MVII, MVIII, MIX, MX, MXI, MXII, CT, FUSE, DCT, CT_new,WR,IR,RR,UR,FR,SR,CR,MARRMot)), # EDIT BY ADDING NEW FUNCTION!!!
    
    PrefixFun = c(rep(FALSE, 16), rep(TRUE, 6 + length(SpecialFunctions))),
    
    stringsAsFactors = FALSE
    
  )
  rownames(functionsDefinition) = functionsDefinition[['Functions']]
  
  MaxNoColumns = max(functionsDefinition[functionSet, 'Arity']) + 1
  
  variableArityRanges = data.frame(TANK = c(5, functionsDefinition[["Arity"]][functionsDefinition[["Functions"]] == "TANK"]))
  
  SFPR = SetSFPR()
  
  DepArgsOfCT = DepArgsOfCTfun()
  assign('DepArgsOfCT', DepArgsOfCT, envir = parent.env(environment()))  
  
  MeaningfullCTcombinations = MeaningfulCTcombsFun()
  assign('MeaningfullCTcombinations', MeaningfullCTcombinations, envir = parent.env(environment()))   
  
  assign('StrictFunctions', StrictFunctions, envir = parent.env(environment()))
  assign('SFPR', SFPR, envir = parent.env(environment()))
  assign('VariableArityRanges', variableArityRanges, envir = parent.env(environment()))  
  
  assign('NonStrictArguments', NonStrictArguments, envir = parent.env(environment()))
  assign('SpecialFunctions', SpecialFunctions, envir = parent.env(environment()))
  assign('FunctionsDefinition', functionsDefinition, envir = parent.env(environment()))
  assign('PrefixFunctions', functionsDefinition[['Functions']][functionsDefinition[['PrefixFun']]], envir = parent.env(environment()))
  assign('MaxNoColumns', MaxNoColumns, envir = parent.env(environment()))
  
  assign('FunctionSet', functionSet, envir = parent.env(environment()))  
  assign('StandardFunctions', functionSet[!functionSet %in% StrictFunctions], envir = parent.env(environment()))  
  
  assign('DependentVariable', depVar, envir = parent.env(environment())) 
  assign('IndependentVariables', indepVar, envir = parent.env(environment())) 
  assign('ConstantRange', constRng, envir = parent.env(environment())) 
  assign('StrictFunctionsArgumentsRange', StrictFunsArgsRng, envir = parent.env(environment())) 
  
  assign('RoundingFactor', roundFactor, envir = parent.env(environment()))
  
  assign('DataSet', dataSet, envir = parent.env(environment()))   
  assign('FitnessFunction', fitFun, envir = parent.env(environment())) 
  
  assign('Weights', weights, envir = parent.env(environment())) 
  assign('PunishOneNodeIndividuals', punish1nodeInd, envir = parent.env(environment())) 
  assign('MaxDepthRun', maxDepthRun, envir = parent.env(environment())) 
  
  
  FitGlobs = FitnessGlobal(DataSet[[depVar]], assignIt = FALSE)
  
  lowValPositions    = FitGlobs[['lowValPositions']]
  highValPositions   = FitGlobs[['highValPositions']]
  mediumValPositions = FitGlobs[['mediumValPositions']]
  logOriginalValues  = FitGlobs[['logOriginalValues']]
  MeanOriginalValues = FitGlobs[['MeanOriginalValues']]
  originalm1         = FitGlobs[['originalm1']]
  originalm2         = FitGlobs[['originalm2']]
  originall2         = FitGlobs[['originall2']]
  lag_num            = FitGlobs[['lag_num']]
  orlength           = FitGlobs[['orlength']] 
  fdcmOriginalValues = FitGlobs[['fdcmOriginalValues']] 
  fdchOriginalValues = FitGlobs[['fdchOriginalValues']] 
  fdclOriginalValues = FitGlobs[['fdclOriginalValues']] 
  MedianOriginalValues = FitGlobs[['MedianOriginalValues']] 
  PeakOriginalValues = FitGlobs[['PeakOriginalValues']] 
  Lag_ac_Original = FitGlobs[['Lag_ac_Original']] 
  VarOriginalValues= FitGlobs[['VarOriginalValues']]
  MeanlogOriginalValues=FitGlobs[['MeanlogOriginalValues']]
  VarlogOriginalValues=FitGlobs[['VarlogOriginalValues']]
  
  assign('MultiObjectiveFitnesses', c('Multi_Madsen', 'Dawson', 'Vis_2')                                                , envir = parent.env(environment()))  
  assign('LowValPositions'      , lowValPositions                                                                       , envir = parent.env(environment()))
  assign('HighValPositions'     , highValPositions                                                                      , envir = parent.env(environment()))
  assign('MediumValPositions'   , mediumValPositions                                                                    , envir = parent.env(environment()))
  assign('OriginalValuesLow'    , DataSet[[depVar]][lowValPositions]                                                    , envir = parent.env(environment()))
  assign('OriginalValuesHigh'   , DataSet[[depVar]][highValPositions]                                                   , envir = parent.env(environment()))
  assign('OriginalValuesMedium' , DataSet[[depVar]][mediumValPositions]                                                 , envir = parent.env(environment())) 
  assign('NS0denominator'       , DataSet[[depVar]] - MeanOriginalValues                                                , envir = parent.env(environment()))
  assign('logOriginalValues'    , logOriginalValues                                                                     , envir = parent.env(environment()))
  assign('MeanlogOriginalValues', MeanlogOriginalValues                                                                 , envir = parent.env(environment()))
  assign('VarlogOriginalValues' , VarlogOriginalValues                                                                  , envir = parent.env(environment()))
  assign('originalm1'           , originalm1                                                                            , envir = parent.env(environment()))
  assign('originalm2'           , originalm2                                                                            , envir = parent.env(environment()))
  assign('originall2'           , originall2                                                                            , envir = parent.env(environment()))
  assign('orlength'             , orlength                                                                              , envir = parent.env(environment()))
  assign('lag_num'              , lag_num                                                                               , envir = parent.env(environment()))
  assign('logNS0denominator'    , logOriginalValues - MeanlogOriginalValues                                             , envir = parent.env(environment()))
  assign('MeanOriginalValues'   , mean(DataSet[[depVar]])                                                               , envir = parent.env(environment()))
  assign('VarOriginalValues'    , VarOriginalValues                                                                     , envir = parent.env(environment()))   
  assign('SumP'                 , sum(DataSet[['P']])                                                                   , envir = parent.env(environment()))
  assign('orrOriginalValues'    , sum(DataSet[[depVar]])/SumP                                                           , envir = parent.env(environment()))
  assign('orrlogOriginalValues' , sum(logOriginalValues)/SumP                                                           , envir = parent.env(environment()))
  assign('Month'                , months(as.Date(DataSet[['Date']]))                                                    , envir = parent.env(environment()))
  assign('MaxMonthlyOriginal'   , max(aggregate( DataSet[[depVar]] ~ Month , DataSet , mean )[2])                       , envir = parent.env(environment()))
  assign('fdcmOriginalValues'   , fdcmOriginalValues                                                                    , envir = parent.env(environment()))
  assign('fdchOriginalValues'   , fdchOriginalValues                                                                    , envir = parent.env(environment()))
  assign('fdclOriginalValues'   , fdclOriginalValues                                                                    , envir = parent.env(environment()))
  assign('MedianOriginalValues' , MedianOriginalValues                                                                  , envir = parent.env(environment()))
  assign('PeakOriginalValues'   , PeakOriginalValues                                                                    , envir = parent.env(environment()))
  assign('Lag_ac_Original'      , Lag_ac_Original                                                                       , envir = parent.env(environment()))
  
  if(FitnessFunction == 'CED' | FitnessFunction == 'SUSE'|FitnessFunction == 'CED_new') ceds = CEDprecomputation(DataSet[[depVar]])
  #   if(FitnessFunction == 'SUSE') ceds = CEDprecomputation_all(DataSet[[depVar]])
  if(FitnessFunction == 'CED' | FitnessFunction == 'SUSE'|FitnessFunction == 'CED_new') sapply(names(ceds), function(nm) assign(nm, ceds[[nm]], envir = parent.env(environment()), inherits = TRUE))
  
}

Qlog = function(x){
  
  x[x <= 0] = 1e-6
  log(x)
  
}

Qlog10 = function(x){
  
  x[x <= 0] = 1e-6
  log10(x)
  
}

Sample = function(x, ...){
  
  x[sample.int(length(x), ...)]
  
}

NoNodesInEquation = function(equation){
  
  parseData = getParseData(parse(text = equation, keep.source = TRUE))
  uncleanNodesPositions = (parseData[['text']]) != '' & parseData[['text']] != '(' & parseData[['text']] != ')'
  uncleanNodes = parseData[uncleanNodesPositions, 'text']
  noUncleanNodes = length(uncleanNodes)
  noNegativeConstants = ifelse(uncleanNodes[1] == '-', 1, 0)  
  
  for(i in 1:(noUncleanNodes - 1)) {
    
    if(any(uncleanNodes[i] == FunctionSet) & any(uncleanNodes[i + 1] == '-')) noNegativeConstants  = noNegativeConstants + 1 
    
  }
  
  noUncleanNodes - noNegativeConstants
  
}

SaveRunResults = function(population, front_list,  runResult, generation, NumberOfGenerations){
  
  #generationLengths = sapply(population, function(x) x[['IndLength']])
  
  #generationDepths = mapply(function(ind, indLength) NodeDepth(ind, indLength), population, generationLengths) 
  
  #bestEquation = population[[which.min(generationFitnesses)]][['Equation']]
  
  #g1 = generation + 1
  
  if (generation == NumberOfGenerations) {
    runResult <- population
  } else {
    runResult <- runResult
  }
  
  # if(yacasSimplification){
  #   
  #   if(generation > 0 && runResult[['BestModel']][generation] == bestEquation) {
  #     
  #     runResult[['BestModelSimplified']][g1] = runResult[['BestModelSimplified']][generation]
  #     
  #   } else {
  #     
  #     simplifiedBestEquation = YacasAnalysis(bestEquation)
  #     
  #     runResult[['BestModelSimplified']][g1] = simplifiedBestEquation 
  #     
  #   }
    
    # Simple complexity Computation is meaningful only for simplified equations
    #runResult[['SimpleComplexity']][g1] = SimpleComplexity(runResult[['BestModelSimplified']][g1]) 
    
  #}
  
  # runResult[['MinFit']][g1] = min(generationFitnesses, na.rm = TRUE)
  # runResult[['MedianFit']][g1] = median(generationFitnesses, na.rm = TRUE)
  # runResult[['BestModelDepth']][g1] = generationDepths[[which.min(generationFitnesses)]]  
  # runResult[['MedianDepth']][g1] = median(generationDepths, na.rm = TRUE)
  # runResult[['HealthyIndividuals']][g1] = sum(!is.na(generationFitnesses))
  # runResult[['BestModel']][g1] = bestEquation
  #runResult[['Ind_front_1']][g1] = list(population[front_list[[1]]])
  
  #   if(!yacasSimplification) runResult[['ComplexityVann']][g1] = runResult[['SimpleComplexity']][g1]
  
  runResult
  
}

YacasAnalysis = function(equation){
  
  roundingConstants = function(regExp, eq){
    
    if(grepl(regExp, eq)){
      
      parsedEquationData = getParseData(parse(text = eq, keep.source = TRUE))
      parsedEquationConstants = as.numeric(parsedEquationData[parsedEquationData[["token"]] == "NUM_CONST", 'text'])
      parsedEquationData[parsedEquationData[["token"]] == "NUM_CONST", 'text'] = as.character(round(parsedEquationConstants, RoundingFactor))
      eqNew = paste0(parsedEquationData[, 'text'], collapse = '')
      
    } else {
      
      eqNew = eq
      
    }
    
    eqNew
    
  }
  
  if(any(grepl('[[:punct:]]', IndependentVariables))) stop('Independent variables names must be without punctuation for YacasAnalysis')
  
  rightSideOfEquation = substring(equation, 5)
  
  # Changing -- for - - , yacas has problem with --
  rightSideOfEquationClear = gsub('([-]){2}', '- -', rightSideOfEquation)
  
  # Simplification
  systemYacasCommand = paste0("yacas -pc --execute '[Echo(N(", rightSideOfEquationClear, "));Exit();]'")
  eqSimp1 = system(systemYacasCommand, intern = TRUE) 
  
  # Rounding constants for yacas simplification
  regexpString = paste0('(([0-9])([[:punct:]])([0-9]){', RoundingFactor + 1, ',})')
  
  eqSimp1rounded = roundingConstants(regexpString, eqSimp1)
  
  systemYacasCommand = paste0("yacas -pc --execute '[Echo(Simplify(", eqSimp1rounded, "));Exit();]'")  
  eqSimp2 = system(systemYacasCommand, intern = TRUE) 
  
  if(grepl('Abs', eqSimp2) | grepl('Gcd', eqSimp2)) {
    
    eqSimp = eqSimp1rounded 
    
  }else{
    
    eqSimp = substr(eqSimp2, 1, nchar(eqSimp2)-1) # the second yacas simplification adds a white space at the end of string
    
  }
  
  # Rounding constants from yacas simplification
  regexpString = paste0('(([0-9])([[:punct:]])([0-9]){', RoundingFactor + 1, ',})')
  
  eqSimpRounded = roundingConstants(regexpString, eqSimp)
  
  
  eqAllSimp = RegexpSimplification(eqSimpRounded)
  
  simplyfiedEquation = paste0("y = ", eqAllSimp)
  
  #   print(equation)
  #   print(simplyfiedEquation)
  #   cat ("Press [enter] to continue or b to start browser")
  #   line <- readline()
  #   if (line == 'b') browser()
  
  simplyfiedEquation
  
  
}

RegexpSimplification = function(equation){
  # Works only for range 0-1 of strict function arguments
  # Term function means here strict function
  
  NewReplacement = function(targetLst, regExpString){
    
    rangeFun = function(x){
      
      xn = as.numeric(substring(x,2))
      
      if(xn > StrictFunctionsArgumentsRange[['upper']]) xn = StrictFunctionsArgumentsRange[['upper']]
      if(xn < StrictFunctionsArgumentsRange[['lower']]) xn = StrictFunctionsArgumentsRange[['lower']]
      if(substring(x,1,1) == ',') out = paste0(',', xn)
      if(substring(x,1,1) == '(') out = paste0('(', xn)
      out
    }
    
    positions = lapply(targetLst, 
                       function(x){
                         pos = vector(mode = 'list', length(x))
                         for(i in 1:length(x)) pos[[i]] = gregexpr(regExpString, x[i])[[1]]
                         pos
                       }
    ) 
    
    outLst = vector(mode = 'list', length(targetLst))
    for(i in seq_along(targetLst)) {
      
      lstElement = targetLst[[i]]
      newLstElement = lstElement
      
      for(j in seq_along(positions[[i]][[1]])) {
        
        pos = positions[[i]][[1]][j]
        if(pos != -1) {
          const = substring(lstElement, pos, pos + attributes(positions[[i]][[1]])[['match.length']][j]-1)[[1]]
          if(substring(const,1,1) == '(') const = sub(',','',const)
          checkConst = rangeFun(const)
          if(substring(const,1,1) == '(') {
            const = sub('\\(','',const)
            const = paste0('\\(',const)
          }
          stS = strsplit(newLstElement, const)
          newLstElement = paste0(stS[[1]], collapse = checkConst)
          
        } else{
          
          newLstElement = lstElement
          
        }
        
      }
      
      outLst[[i]] = newLstElement
      
    }
    
    outLst
    
  }
  
  # Possible function strings
  regStrFun = sapply(StrictFunctions, function(x) paste0(x, '\\(.[^\\)]*\\),*[-0-9\\.,\\)]*'))
  
  # Positions of functions (while function is wrapped function, it will be solved again, firstly as parent fun and secondly as fun itself, but it doesnt matter.???)
  posFun = lapply(regStrFun, function(x) gregexpr(x, equation)[[1]])
  
  strFunIncluded = sapply(posFun, function(x) x[1]!=-1)
  
  # No strict function case
  if(!any(strFunIncluded == TRUE)) return(equation)
  
  # Number of strict function in equation and sequence of this length
  okList = posFun[strFunIncluded]
  
  # Separation of identified functions
  separFuns = lapply(okList, function(x){sepFun = character(length(x)); 
  for(i in 1:length(x)) sepFun[i] = substr(equation, x[i], x[i]+attributes(x)[["match.length"]][i]-1);
  sepFun})
  
  # Replacement of negative values     
  # Regexp for negative values
  separFuns = as.list(separFuns[[1]])
  
  regNegative = ',-[0-9]+\\.[0-9]+|,-[0-9]+|\\(-[0-9]+\\.[0-9]+,|\\(-[0-9]+,'
  eqsNoNeg = NewReplacement(separFuns, regNegative)
  
  # Replacement of values, bigger than StrictFunctionsArgumentsRange[['upper']]. Minimal value of upper bound is 1!!!
  regBigPositive = ',[1-9]+\\.[0-9]+|,[1-9]+|\\([1-9]+\\.[0-9]+,|\\([1-9]+,'
  
  eqsCleanFinal = NewReplacement(eqsNoNeg, regBigPositive)
  
  #   print(separFuns)
  #   print(eqsCleanFinal)
  #   cat ("Press [enter] to continue or b to start browser")
  #   line <- readline()
  #   if (line == 'b') browser()
  
  eqFinal = equation
  
  for(i in seq_along(separFuns)){
    
    eqFinal = sub(separFuns[i], eqsCleanFinal[i], eqFinal, fixed = TRUE)
    
  }
  eqFinal
  
}

Simulation = function(data, model){
  
  computedValues = try(eval(parse(text = model), envir = data), silent = TRUE)
  
  if(class(computedValues) == 'try-error' | any(!is.finite(computedValues))){ 
    
    computedValues = NA 
    
  }
  
  computedValues
  
}
SimulationFitness = function(data, model, referenceVariable, rounding = 3){
  
  round(FitnessComputation(Simulation(data, model), referenceVariable), rounding)
  
}
SimpleComplexity = function(equation){
  
  # Function for number of features of given identifier (regular expression form) in equation
  NrFeatures = function(idenitifier, equation){
    
    positions = gregexpr(idenitifier, equation)[[1]]
    np = length(positions)
    nrf = ifelse(positions[1] > -1, np, 0)
    
    lengths = attributes(positions)[["match.length"]] - 1
    
    SFmodels = paste0('M', as.roman(c(1:12)))
    
    if(any(SFmodels %in% substring(equation, positions, positions + lengths))){
      
      sepFun = character(np)
      
      # identification of superflex functions in equation
      for(i in 1:np) sepFun[i] = substr(equation, positions[i], positions[i] + lengths[i])
      
      # number of superflex arguments
      nrsa = sum(sapply(sepFun, function(x) length(formals(x))))
      
      # complexities of superflex functions
      SFcomplexities = numeric(length(SFmodels))
      names(SFcomplexities) = SFmodels
      
      # one tank models
      ms1 = paste0('M', as.roman(1:2))
      SFcomplexities[ms1] = rep(10, length(ms1)) + sapply(ms1, function(x) length(formals(x))) * 2
      
      # two tank models      
      ms2 = paste0('M', as.roman(c(3:5,8)))
      SFcomplexities[ms2] = rep(20, length(ms2)) + sapply(ms2, function(x) length(formals(x))) * 2
      
      # three tank models
      ms3 = paste0('M', as.roman(c(6,7,9,10,11)))
      SFcomplexities[ms3] = rep(30, length(ms3)) + sapply(ms3, function(x) length(formals(x))) * 2
      
      # four tank models
      ms4 = 'MXII'
      SFcomplexities[ms4] = rep(40, length(ms4)) + sapply(ms4, function(x) length(formals(x))) * 2
      
      cpxSF = sum(SFcomplexities[sepFun])     
      
    } else{
      
      nrsa = 0
      cpxSF = 0
      
    }
    
    list(nrFeat = nrf, nrStricArgs = nrsa, CpxSFfun = cpxSF)
    
  }
  
  # Superflex functions definition for search
  SFf = paste0('M', as.roman(c(c(1:12))), collapse = '|')
  
  # Number of constants (minus the number of constants in all R2T and R4T functions)
  SFcheck = NrFeatures(SFf, equation)
  nrConst = NrFeatures('[0-9]+\\.[0-9]+|[0-9]+', equation)[['nrFeat']] - SFcheck[['nrStricArgs']] 
  
  # Number of variables
  nrVars = NrFeatures(paste0(IndependentVariables, collapse = '|'), equation)[['nrFeat']]
  
  # Number of functions class 1 (eg. - and +)
  nrFunClass1 = NrFeatures('\\+|-', equation)[['nrFeat']] # problem with minus sign of negative values
  
  # Number of functions class 2 (eg. * and /)
  nrFunClass2 = NrFeatures('\\*|\\/', equation)[['nrFeat']]
  
  # Number of functions class 3  (eg. ^ and sqrt)
  nrFunClass3 = NrFeatures('\\^|sqrt', equation)[['nrFeat']]
  
  cpxSeries = c(1,2,3,5,7)
  cpxOfCommons = sum(c(nrConst, nrVars, nrFunClass1, nrFunClass2, nrFunClass3) * cpxSeries) 
  
  cpxOfSFfunctions = SFcheck[['CpxSFfun']]
  
  cpxOfCommons + cpxOfSFfunctions
  
}

SimplifyCT = function(eq){
  ## Works with CT equations only. eg. max tree depth 1. Problem with equations like below.
  #   y = 0.051/(CT(0,0.277,0,0.827,0,0.971,0.105,0.169,0.555,1,0.836,0,1,0.321,0.062,0.16,1,0.712,0,0.043,0.124,0.17,0.955,0.145,0.228,0.864,0,0,1,0.036,1,0.088,0.931,0.85)*CT(0.173,0.323,0.014,0.036,1,0.408,0.603,0.144,0.139,0.925,0.313,1,0.471,0.784,0.489,0.683,0.784,0.132,0.921,0,0.75,0.376,0.185,0.57,0.293,0.076,0,0.681,0.069,0.172,0.719,0.98,0.474,0.592))'  if(!grepl('CT', eq)) return(list(SimplifiedEquation = eq, ResultString = NA))
  if(!grepl('CT', eq)) return(list(SimplifiedEquation = eq, ResultString = NA))
  positions = function(target, regexString, firstCoeff){
    
    pos = gregexpr(regexString, target)[[1]]
    posDF = data.frame(position = as.vector(pos), length = attributes(pos)[["match.length"]] - firstCoeff) # -1 is here because of last character which is comma and must not be in selection
    if(posDF[["position"]][1] == -1) posDF = NULL
    posDF
    
  }
  
  ## This is related to MeaningfulCTcombsFun!!!
  combs = MeaningfulCTcombs()
  
  MeaningfullCTcombinations = MeaningfulCTcombsFun()
  #   ## Checking of negative arguments as well. Unfortunately this can happen, need to be fixed somewhere else.!!!
  
  # Position of first arguments 
  regex4DEargsFirsts  = '\\([0-9]+\\.[0-9]+,|\\([0-9]+,|\\(-[0-9]+\\.[0-9]+,|\\(-[0-9]+,'
  posArgsDfFirsts = positions(eq, regex4DEargsFirsts, 1)
  
  # Position of other arguments 
  regex4DEargs = ',[0-9]+\\.[0-9]+|,[0-9]+|,-[0-9]+\\.[0-9]+|,-[0-9]+'
  posArgsDf = positions(eq, regex4DEargs, 0)
  
  # Position of non argument constants
  regex4DEconst = '[-*/+ (^][0-9]+\\.[0-9]+|[-*/+ (^][0-9]+'
  posConstDf = positions(eq, regex4DEconst, 0)
  
  allConst = unique(rbind(posArgsDfFirsts, posArgsDf))
  # allConst = allConst[order(allConst[["position"]]),] # For the one CT function equation is not needed.
  allConst[['type']] = allConst[["position"]] %in% c(as.vector(posArgsDf[['position']]), as.vector(posArgsDfFirsts[['position']]))
  allConst[['type']] = ifelse(allConst[['type']], 'argument', 'constant')  
  
  # Separation of numeric values
  nrAllConst = nrow(allConst)
  constRaw = character(nrAllConst)
  
  for(i in seq_len(nrAllConst)){
    
    from = allConst[['position']][i]
    to = from + allConst[['length']][i] - 1
    if(i == 1 & allConst[['type']][i] == 'argument') constRaw[i] = substr(eq, from, to + 1) else constRaw[i] = substr(eq, from, to)
    
  } 
  
  numVals = sub('[*/ (,^]', '', constRaw)  
  numVals = sub('\\(|,', '', numVals) # must be done twice because of strict SF functions first argument 
  numVals = as.numeric(numVals)
  names(numVals) = names(formals(CT))
  
  # Simplification of res and option parameters  
  # numVals are changed only in res parameters, thus the final equation is still computable by CHMGPpredict
  # scaledNumVals is prepared for view of final equation with real values of parameters eg. scaled to range of SFPR. Not usable with CHMGPpredict!
  scaledNumVals = round(ScaleParams(numVals), 3)
  
  optResNam = tail(names(numVals), 15)
  scaledNumVals[optResNam] = round(scaledNumVals[optResNam])
  
  resNam = tail(names(numVals), 7)
  lagNam = head(numVals,3)
  
  numVals[lagNam][numVals[lagNam] > 1] = 1
  
  
  numVals[resNam][numVals[resNam] > 1] = 1
  numVals[resNam][numVals[resNam] < 0] = 0
  numVals[resNam] = round(numVals[resNam])
  scaledNumVals[resNam] = numVals[resNam]
  
  if(numVals[['D_R']] == 0) numVals['r_res'] = 0  
  if(numVals[['D_S']] == 1) numVals['u_res'] = 0
  #if(numVals[['option_u']] != 1 ) numVals['mu_Qq_UR'] = 0 # i.e. option_u =4, mu_Qq_UR is to be removed, how??
  #if(numVals[['alpha_Qq_FR']] == 1) numVals['option_f'] = 1
  #if(numVals[['alpha_Qq_SR']] == 1) numVals['option_s'] = 1
  scaledNumVals['r_res'] = numVals['r_res']  
  scaledNumVals['u_res'] = numVals['u_res'] 
  #scaledNumVals['option_f'] = numVals['option_f']
  #scaledNumVals['option_s'] = numVals['option_s']
  #scaledNumVals['mu_Qq_UR'] = numVals['mu_Qq_UR']
  
  # Preparing computable equation
  cnstRawNew = constRaw
  names(cnstRawNew) = names(numVals)
  cnstRawNew[resNam] = paste0(',',format(round(as.numeric(numVals[resNam]),3), nsmall = 3))
  cnstRawNew[1] = gsub(',', '', cnstRawNew[1])
  newEq = eq
  
  newEq = sub(substr(newEq, allConst[['position']][1], tail(allConst[['position']],1) + tail(allConst[['length']],1) - 1),
              paste0(cnstRawNew, collapse = ""), 
              newEq, 
              fixed = TRUE)
  
  # Preparing real parameters equation
  cRN4Show = constRaw
  names(cRN4Show) = names(scaledNumVals)
  cRN4Show = paste0(',',round(as.numeric(scaledNumVals),3))
  cRN4Show[1] = gsub(',', '(', cRN4Show[1])
  newEqForShow = eq
  
  newEqForShow = sub(substr(newEqForShow, allConst[['position']][1], tail(allConst[['position']],1) + tail(allConst[['length']],1) - 1),
                     paste0(cRN4Show, collapse = ''), 
                     newEqForShow, 
                     fixed = TRUE)  
  
  ress = numVals[tail(names(numVals), 7)]
  names(ress) = names(MeaningfullCTcombinations[[1]])
  
  r_option=2
  names(r_option)='option_r'
  #options= scaledNumVals[26:29]
  options=c(scaledNumVals[40],r_option,scaledNumVals[41:44])
  names(options)=c("w","i","r","u","f","s","c")
  print(options)
  #paste(as.character(scaledNumVals[26:29]), collapse=", ")
  
  indices=which(ress>0)
  ress1<-ress[indices]
  options_logic=options[match(names(ress1),names(options))]
  #options_logic=options_logic[!is.na(options_logic)]
  temp_i=as.character(options_logic['i'])
  newvals_i<-c(`1`='power',`2`='R.Hyperbolic')
  if(is.na(temp_i)){
    optionsstring_i=NA
  }else{
    optionsstring_i=newvals_i[temp_i]
  }
  newvals <- c(`1`='power',`2`='linear',`3`='Monod',`4`='M.Logistic')
  #optionsstring=newvals[as.character(options_logic)]
  optionsstring_others=newvals[as.character( options_logic[!names(options_logic) %in% 'i'])]
  temp=paste(optionsstring_others, collapse=", ")
  if(is.na(optionsstring_i)){
    optionsstring=temp
  } else {
    optionsstring=paste(optionsstring_i,',',temp)
  }
  
  
  testMeaningfulness = sapply(MeaningfullCTcombinations, function(x) identical(unlist(x),ress))
  
  resultString = combs[testMeaningfulness]
  
  
  if(length(resultString) < 1) resultString = 'NaMC' # NaMC - Not a Meaningful Combination
  
  #   if(length(resultString) > 1) browser()
  
  #   browser()
  
  list(SimplifiedEquation = newEq, RealParamEq = newEqForShow, ResultString = resultString,OptionString=optionsstring)
  
}

TANKidentify = function(eq){
  if(!grepl('CT', eq)) return(NA)
  positions = function(target, regexString, firstCoeff){
    
    pos = gregexpr(regexString, target)[[1]]
    posDF = data.frame(position = as.vector(pos), length = attributes(pos)[["match.length"]] - firstCoeff) # -1 is here because of last character which is comma and must not be in selection
    if(posDF[["position"]][1] == -1) posDF = NULL
    posDF
    
  }
  
  # Position of other arguments 
  regex4DEargs = ',[0-9]+\\.[0-9]+|,[0-9]+'
  posArgsDf = positions(eq, regex4DEargs, 0)
  ress = nrow(posArgsDf)
  
  
  c(Arguments = ress, Reservoirs = ress/6)
  
  
}  
CTtypesPlot = function(CTs){
  
  shortenCTtype = function(fullCTstring){
    
    combs = gsub('Tank - |Tanks - ', 'T-', fullCTstring)
    combs = gsub('One '  , '1', combs)
    combs = gsub('Two '  , '2', combs)
    combs = gsub('Three ', '3', combs)
    combs = gsub('Four ' , '4', combs)
    combs = gsub('Five ' , '5', combs)
    combs = gsub('Six ' , '6', combs)
    
    combs
    
  }
  
  #   CTs = GPresult[['CTtypes']] #smaz
  
  combsOrig = MeaningfulCTcombs()
  
  combs = shortenCTtype(combsOrig)
  
  allCTdf = data.frame(matrix(0, nrow = length(CTs), ncol = length(combs) + 1))
  colnames(allCTdf) = c(combs, 'NaMC')
  namsGen = names(CTs)
  rownames(allCTdf) = namsGen
  
  countCTperGen = sapply(CTs, function(x){ if(is.atomic(x)) y = table(x) else y = NA; y})
  
  #   names(countCTperGen[['Generation_10']])
  
  for(i in 1:nrow(allCTdf)){
    
    colNam = shortenCTtype(names(countCTperGen[[i]]))
    allCTdf[i, colNam] = countCTperGen[[i]]
    
  }
  
  library('ggplot2')
  library('reshape2')
  library('ggrepel')
  
  
  allCTdf[['Generation']] = as.factor(0:(length(namsGen)-1))
  allCTdfMelted = melt(allCTdf ,  id.vars = 'Generation', variable.name = 'CTtypes')
  
  names(allCTdfMelted) = c(head(names(allCTdfMelted), -1), 'Occurence')
  
  allCTdfMelted[['CTtypesShort']] = substr(allCTdfMelted[['CTtypes']],1,2)
  
  seqBreaks = round(seq(0,nrow(allCTdf)-1, length.out = 4))
  
  #   browser()
  
  pAll = 
    ggplot(allCTdfMelted, aes(x = Generation, y = Occurence)) + 
    geom_point(na.rm = TRUE) + 
    scale_x_discrete('Generation', breaks = seqBreaks, labels = seqBreaks ) + 
    scale_y_continuous(limits = c(0, max(allCTdf[, !colnames(allCTdf) %in% c('NaMC', 'Generation')]))) +
    facet_wrap(~CTtypes, ncol = 7, scale = 'free_x') + 
    theme(axis.text.x=element_text(angle=0))
  ggsave('AllCTcombinations.pdf', plot = pAll, height = 16, width = 25, units = 'cm')
  
  for(i in unique(allCTdfMelted[['CTtypesShort']])) { 
    
    sbst = subset(allCTdfMelted, CTtypesShort == i)
    labs_repel = sbst[['CTtypes']]
    labs_repel[sbst[['Occurence']] == 0] = NA
    ggplot(sbst, aes(x = Generation, y = Occurence, color = CTtypes, label = CTtypes)) + 
      geom_point() + 
      geom_label_repel(aes(label=labs_repel), size = 3, na.rm = TRUE, show.legend = FALSE) +
      scale_y_continuous(breaks = 0:max(sbst[['Occurence']])) + 
      guides(color=guide_legend(override.aes=list(shape=rep(16, length(unique(sbst[['CTtypes']]))), size=3, linetype=0)))
    ggsave(paste0(i, '_CT.pdf')) 
    
  }    
  
}
