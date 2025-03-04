# Clear Global Environment
rm(list = ls())

# Setting working directory
setwd('C:/Naila/Singapore/Research/Vidura codes/Final_Codes/Lumped/MARRMOT')
path <- 'C:/Naila/Singapore/Research/Vidura codes/Final_Codes/Lumped/MARRMOT'
# Reading input data
luxHue = read.table('./Syn_data/Blackwater.txt', sep=',',header = TRUE)
luxHue_train =luxHue[8493:10684,]

# set CHMGP folder 
setwd('./CHMGP')
# source('./CHMGP.R', local = TRUE)

# source('CHMGP_GlobalVar.R')

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




# source('CHMGP_FirstGeneration.R')

FirstGeneration = function(maxDepthIni, populationSize){
  method = "grow"
  
  population = vector('list', populationSize) 
  
  for (i in 1:populationSize){
    
    method = ifelse(method == "grow", "full", "grow")
    actualDepth = ifelse(method == "grow", sample(0:maxDepthIni, 1), maxDepthIni)
    population[[i]] = CreateIndividual(actualDepth, method, restricted = FALSE)
    
  }
  
  population
  
}




# source('CHMGP_Individual.R')

NodeType = function(method, functionSet, parentFun = NA, terminal = FALSE, restrictedNode = FALSE){
  
  if(terminal && restrictedNode) nodeType = "Constant"
  
  if(terminal && !restrictedNode){
    
    if(is.na(parentFun)) nodeType = Sample(c("Constant", "Variable"), 1)
    nodeType = ifelse(parentFun %in% StrictFunctions, "Variable", Sample(c("Constant", "Variable"), 1))
    
  } 
  
  if(!terminal && restrictedNode){
    # No special functions, because these functions returns a vector!
    # Random choice of function or constant in grow method
    possibleNodes = functionSet[!functionSet %in% SpecialFunctions]
    if(length(possibleNodes) == 0) possibleNodes = c("Constant", "Variable")
    
    if(method == "grow"){
      
      nodeType = Sample(c("Function", "Constant"), 1)
      if(nodeType == "Function") nodeType = Sample(possibleNodes, 1)
      
    }
    
    # Random choice of function or constant in full method
    if(method == "full"){
      
      nodeType = Sample(possibleNodes, 1)
      
    }
    
  }  
  
  if(!terminal && !restrictedNode){
    
    # Boosting probability of Strict function occurence
    stricFuncBoost = rbinom(1,1,0.5)
    if(stricFuncBoost & any(functionSet %in% StrictFunctions)) fn = Sample(functionSet[functionSet %in% StrictFunctions], 1) else fn = Sample(functionSet, 1)
    
    # Random choice of function or terminal in grow method
    if(method == "grow"){
      
      nodeType = Sample(c("Function", "Terminal"), 1)
      
      if(parentFun %in% StrictFunctions) {
        
        nodeType = ifelse(nodeType == "Function", fn, "Variable")      
        
      } else{
        
        nodeType = ifelse(nodeType == "Function", fn, Sample(c("Constant", "Variable"), 1))
        
      }
      
    }
    
    # Random choice of function or constant in full method
    if(method == "full"){
      
      nodeType = fn
      
    }
    
  }
  
  nodeType
  
}

NodeCompletion = function(new_row, nodeType, usedArities, functionSet){
  
  tuple = rep(NA, MaxNoColumns)
  tuple[1] = nodeType
  
  if(any(functionSet == nodeType)) {
    
    if(any(names(VariableArityRanges) == tuple[1])){
      
      funcsArityRanges = VariableArityRanges[[tuple[1]]]
      ## Random number of arities, sequence is created by specific range for TANK, see variableArityRanges in GlobalVar
      n_arit = sample(seq(funcsArityRanges[1], funcsArityRanges[2], funcsArityRanges[1]), 1)
      
    } else {
      
      n_arit = usedArities[which(functionSet == tuple[1])]
      
    }
    
    tuple[2:(n_arit + 1)] = (new_row + 1):(new_row + n_arit)
    
    new_row = new_row + n_arit
    
  }
  
  return(list(tuple = tuple, new_row = new_row))
  
}

InsertTerminals = function(indArray){
  
  constantsPositions = which(indArray[ ,1] == 'Constant')
  lcp = length(constantsPositions)
  
  if(lcp > 0) {
    
    indArray[constantsPositions, 1] = runif(lcp, ConstantRange[1], ConstantRange[2])
    
  }
  
  variablesPositions = which(indArray[ ,1] == 'Variable')
  lvp = length(variablesPositions)
  
  if(lvp > 0) {
    
    indArray[variablesPositions, 1] = Sample(IndependentVariables, lvp, replace = TRUE)
    
  }
  
  indArray
  
}

CreateIndividual = function(maxDepth, method, restricted){
  
  # Arities of used functions
  usedArities = FunctionsDefinition[FunctionSet, 'Arity']
  
  ##########
  ReservoirSet <-c("WR","IR","RR","UR","FR","SR","CR")
  AddusedArities = FunctionsDefinition[ReservoirSet, 'Arity']
  
  #########
  # Test of one node individual
  singleNodeIndividual = ifelse(maxDepth == 0 | is.null(usedArities), TRUE, FALSE)
  
  if(restricted){  
    
    # vector of nonrestricted functions
    funSet = FunctionSet[!FunctionSet %in% StrictFunctions]
    ##### this might be the reason for subscript out of problem, this will added to all codes ###
    usedArities = usedArities[!FunctionSet %in% StrictFunctions]
    
  } else{
    
    funSet = FunctionSet
    usedArities = usedArities
    
  }
  
  if(length(funSet) == 0) singleNodeIndividual = TRUE
  
  # Max. number of rows in individual array
  noRow = ifelse(singleNodeIndividual, 1, sum(max(usedArities)^(0:maxDepth)))
  
  # Individual array preparation
  
  indArray = array(c(NA), dim = c(noRow, MaxNoColumns))
  
  # First node
  nodeType = NodeType(method, funSet, parentFun = NA, terminal = singleNodeIndividual)
  nodeComplete = NodeCompletion(1, nodeType, usedArities, funSet)
  
  indArray[1, ] = nodeComplete[['tuple']]
  
  # Finishing one node individual
  if(singleNodeIndividual){
    
    indArray = InsertTerminals(indArray) 
    
    return(list(IndArray = indArray, Equation = NA, Fitness = c(rep(NA,length(FitnessFunction))), DriveVec = 'Terminal',  IndLength = 1, Changed = TRUE,Front = NA, Crowding_distance = NA,Sim =c())) 
    
  }
  
  # Initial settings of actual depth, row, maximal row number in given depth and determination of pointers of actual node
  actualDepth = 1
  actualRow = 1
  maxRowLevel = 1
  newRow = nodeComplete[['new_row']]
  
  # Vector of array rows which will have restriction due to occurence of strict function (positions of restricted arguments of this strict function)
  transferredRestrictions = numeric()
  
  ########################################
  
  if (indArray[1,1]=="CT_new"){
    while(actualDepth <= maxDepth) {
      if (actualDepth==1){
        actualColumn = 2
        while(actualColumn <= MaxNoColumns && !is.na(indArray[actualRow, actualColumn])){
          actFun = indArray[actualRow, 1]
          actNSA = NonStrictArguments[actFun]
          restrictedNode = (actFun %in% StrictFunctions && actualColumn > (1 + actNSA)) | actualRow %in% transferredRestrictions
          if (actualColumn==2){
            nodeType <- "WR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==3) {
            nodeType <- "IR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==4) {
            nodeType <- "RR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==5) {
            nodeType <- "UR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==6) {
            nodeType <- "FR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==7) {
            nodeType <- "SR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==8) {
            nodeType <- "CR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else {
            #actFun = indArray[actualRow, 1]
            #actNSA = NonStrictArguments[actFun]
            #restrictedNode = (actFun %in% StrictFunctions && actualColumn > (1 + actNSA)) | actualRow %in% transferredRestrictions
            
            if(actualDepth < maxDepth && !restrictedNode){
              
              nodeType = NodeType(method, funSet, parentFun = actFun, terminal = FALSE, restrictedNode = FALSE)
              
            } else if(actualDepth < maxDepth && restrictedNode) {
              
              nodeType = NodeType(method, funSet, parentFun = actFun, terminal = FALSE, restrictedNode = TRUE)
              
            } else if(actualDepth == maxDepth && !restrictedNode){
              
              nodeType = NodeType(method, funSet, parentFun = actFun, terminal = TRUE, restrictedNode = FALSE)
              
            } else if(actualDepth == maxDepth && restrictedNode){
              
              nodeType = NodeType(method, funSet, parentFun = actFun, terminal = TRUE, restrictedNode = TRUE)
              
            }
            nodeComplete = NodeCompletion(newRow, nodeType, usedArities, funSet)
          }
          #nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
          
          targetRow = as.numeric(indArray[actualRow, actualColumn])
          
          indArray[targetRow, ] = nodeComplete[['tuple']]
          
          if(restrictedNode) transferredRestrictions = c(transferredRestrictions, targetRow)
          
          newRow = nodeComplete[['new_row']]
          
          actualColumn = actualColumn + 1
        }
      } else {
        actualColumn = 2
        while(actualColumn <= MaxNoColumns && !is.na(indArray[actualRow, actualColumn])){
          actFun = indArray[actualRow, 1]
          actNSA = NonStrictArguments[actFun]
          restrictedNode = (actFun %in% StrictFunctions && actualColumn > (1 + actNSA)) | actualRow %in% transferredRestrictions
          
          if(actualDepth < maxDepth && !restrictedNode){
            
            nodeType = NodeType(method, funSet, parentFun = actFun, terminal = FALSE, restrictedNode = FALSE)
            
          } else if(actualDepth < maxDepth && restrictedNode) {
            
            nodeType = NodeType(method, funSet, parentFun = actFun, terminal = FALSE, restrictedNode = TRUE)
            
          } else if(actualDepth == maxDepth && !restrictedNode){
            
            nodeType = NodeType(method, funSet, parentFun = actFun, terminal = TRUE, restrictedNode = FALSE)
            
          } else if(actualDepth == maxDepth && restrictedNode){
            
            nodeType = NodeType(method, funSet, parentFun = actFun, terminal = TRUE, restrictedNode = TRUE)
            
          }
          
          
          nodeComplete = NodeCompletion(newRow, nodeType, usedArities, funSet)
          
          targetRow = as.numeric(indArray[actualRow, actualColumn])
          
          indArray[targetRow, ] = nodeComplete[['tuple']]
          
          if(restrictedNode) transferredRestrictions = c(transferredRestrictions, targetRow)
          
          newRow = nodeComplete[['new_row']]
          
          actualColumn = actualColumn + 1
          
        }
      }
      if(actualRow == maxRowLevel){
        
        if(any(!is.na(as.numeric(indArray[1:actualRow, 2:MaxNoColumns])))){ ## as.numeric pryc
          
          maxRowLevel = max(as.numeric(indArray[1:actualRow, 2:MaxNoColumns]), na.rm = TRUE) ## as.numeric pryc
          
        }
        
        actualDepth = actualDepth + 1
        
      }   
      
      actualRow = actualRow + 1
      
      # Termination of main cycle
      if(is.na(indArray[actualRow, 1])) {
        
        actualDepth = maxDepth + 1
        
      }
      
    }
  } else {
    
    # Construction of the individual array
    while(actualDepth <= maxDepth){
      
      actualColumn = 2
      
      # Construction of target node, eg. on which the actual node points out
      while(actualColumn <= MaxNoColumns && !is.na(indArray[actualRow, actualColumn])){
        
        # some simplification of variables from global variables
        actFun = indArray[actualRow, 1]
        # actNSA = NonStrictArguments[which(actFun == StrictFunctions)]
        actNSA = NonStrictArguments[actFun]
        
        restrictedNode = (actFun %in% StrictFunctions && actualColumn > (1 + actNSA)) | actualRow %in% transferredRestrictions
        
        if(actualDepth < maxDepth && !restrictedNode){
          
          nodeType = NodeType(method, funSet[-1], parentFun = actFun, terminal = FALSE, restrictedNode = FALSE)
          
        } else if(actualDepth < maxDepth && restrictedNode) {
          
          nodeType = NodeType(method, funSet, parentFun = actFun, terminal = FALSE, restrictedNode = TRUE)
          
        } else if(actualDepth == maxDepth && !restrictedNode){
          
          nodeType = NodeType(method, funSet[-1], parentFun = actFun, terminal = TRUE, restrictedNode = FALSE)
          
        } else if(actualDepth == maxDepth && restrictedNode){
          
          nodeType = NodeType(method, funSet, parentFun = actFun, terminal = TRUE, restrictedNode = TRUE)
          
        }
        
        
        nodeComplete = NodeCompletion(newRow, nodeType, usedArities, funSet)
        
        targetRow = as.numeric(indArray[actualRow, actualColumn])
        
        indArray[targetRow, ] = nodeComplete[['tuple']]
        
        if(restrictedNode) transferredRestrictions = c(transferredRestrictions, targetRow)
        
        newRow = nodeComplete[['new_row']]
        
        actualColumn = actualColumn + 1
        
      }
      
      # Change of actual row
      if(actualRow == maxRowLevel){
        
        if(any(!is.na(as.numeric(indArray[1:actualRow, 2:MaxNoColumns])))){ ## as.numeric pryc
          
          maxRowLevel = max(as.numeric(indArray[1:actualRow, 2:MaxNoColumns]), na.rm = TRUE) ## as.numeric pryc
          
        }
        
        actualDepth = actualDepth + 1
        
      }   
      
      actualRow = actualRow + 1
      
      # Termination of main cycle
      if(is.na(indArray[actualRow, 1])) {
        
        actualDepth = maxDepth + 1
        
      }
      
    }
  }
  
  # Shortening of individual
  IndRows = max(which(!is.na(indArray[, 1])))
  indArray = head(indArray, IndRows)
  
  # Driving vector for next use
  driveVec = indArray[, 1]
  driveVec[driveVec == 'Variable' | driveVec == 'Constant'] = 'Terminal'
  driveVec[driveVec != 'Terminal'] = 'Function'
  
  if(restricted){
    
    # Changing all terminals to constants
    indArray[driveVec == 'Terminal', 1] = 'Constant'
    
  }
  
  # Terminals adding
  indArray = InsertTerminals(indArray) 
  
  if(length(transferredRestrictions) == 0) transferredRestrictions = NA
  
  list(IndArray = indArray, Equation = NA, Fitness = c(rep(NA,length(FitnessFunction))), DriveVec = driveVec, RestrictedRows = transferredRestrictions, IndLength = IndRows, Changed = TRUE, Front = NA, Crowding_distance = NA,Sim =c())
  
}

# source('CHMGP_ComputationFunctions.R')

## Heaviside function - unit step

Hstep = function(x){
  
  as.numeric(x >= 0)
  
}


## Gr: greater than

Gr = function(x1, x2){
  
  as.numeric(x1 > x2)
  
}


## GrEq: greater or equal

GrEq = function(x1, x2){
  
  as.numeric(x1 >= x2)
  
}


## Eq: equal

Eq = function(x1, x2){
  
  as.numeric(x1 == x2)
  
}

Pdiv = function(x1, x2){
  
  if(any(x2 == 0)){
    
    x2[which(x2 == 0)] = 1
    
  } 
  
  x1 / x2 
  
}

Psqrt = function(x){
  
  sqrt(abs(x))
  
}


## DLY: delay function
# x - vector or scalar for shift
# shift - vector or scalar of window length (zero is actual-present value)

DLY = function(x, shift){
  
  if(any(is.na(c(x, shift)))) return(NA)
  
  n = length(x)
  m = length(shift) 
  output = numeric(n)
  
  if(n == 1 | (m == 1 && shift == 0)) return(x)
  
  fOutput = .Fortran("SlidingWindow", x = as.double(x), n = as.integer(n), shift = as.double(shift), m = as.integer(m),
                     output = as.double(output), fun = as.character('DLY'))
  
  fOutput[['output']]
  
}


## SSUM: Shifted sum, sum of sliding window
# x - vector or scalar for shift
# shift - vector or scalar of window length (zero is actual-present value) 

SSUM = function(x, shift){
  
  if(any(is.na(c(x, shift)))) return(NA)
  
  n = length(x)
  m = length(shift) 
  output = numeric(n)
  
  if(n == 1 | (m == 1 && shift == 0)) return(x)
  
  fOutput = .Fortran("SlidingWindow", x = as.double(x), n = as.integer(n), shift = as.double(shift), m = as.integer(m),
                     output = as.double(output), fun = as.character('SUM'))
  
  fOutput[['output']]
  
}


## SMA: Simple Moving Average
# x - vector or scalar for shift
# shift - vector or scalar of window length (zero is actual-present value) 

SMA = function(x, shift){
  
  if(any(is.na(c(x, shift)))) return(NA)
  
  n = length(x)
  m = length(shift) 
  output = numeric(n)
  
  if(n == 1 | (m == 1 && shift == 0)) return(x)
  
  
  fOutput = .Fortran("SlidingWindow", x = as.double(x), n = as.integer(n), shift = as.double(shift), m = as.integer(m),
                     output = as.double(output), fun = as.character('SMA'))
  
  fOutput[['output']]
  
}

## RES:Reservoir function, function of simple reservoir (linear in storage response)
# x - input data for reservoir
# K - parametr nadrze
RES = function(x, K){
  
  if(any(is.na(c(x, K)))) return(NA)
  
  #   K[K == 0] = 1e-9
  K[K < 1] = 1
  
  n = length(x)
  m = length(K)
  output = numeric(n)
  
  if(n == 1) return(x)
  
  K = abs(K)
  
  if(m == 1) K = rep(K, n)
  
  fOutput = .Fortran("Reservoir", x = as.double(x), n = as.integer(n), a = as.double(K), output = as.double(output))
  
  fOutput[['output']]
  
}

# source('CHMGP_ComputeIndividual.R')

CHMGPpredict = function(model, data, testing = FALSE){
  
  if(testing){
    SFPR = SetSFPR()
    assign('SFPR', SFPR, envir = parent.env(environment())) 
    assign('DataSet', data, envir = parent.env(environment()))  # this is due to variables were set as strict
    DepArgsOfCT = DepArgsOfCTfun()
    assign('DepArgsOfCT', DepArgsOfCT, envir = parent.env(environment()))
    MeaningfullCTcombinations = MeaningfulCTcombsFun()
    assign('MeaningfullCTcombinations', MeaningfullCTcombinations, envir = parent.env(environment()))   
  }
  
  computedValues = try(eval(parse(text = model), envir = data))
  
  if(class(computedValues) == 'try-error' | any(!is.finite(computedValues))){ 
    
    computedValues = rep(NA,nrow(data))
    
  }
  
  ## Threshold of Q, for simulation of Q only!!!
  computedValues[computedValues < 0] = 0
  
  computedValues
  
}
FcriteriaPrecomputation = function(obs, sim, fitF){
  
  segmentFlows = function(obs, sim, p1, p2 = NA){
    
    seg_1 = quantile(obs, 1 - p1, na.rm = TRUE, names = FALSE)
    if(!is.na(p2)) seg_2 = quantile(obs, 1 - p2, na.rm = TRUE, names = FALSE)
    
    if(!is.na(p2)) index = which(obs >= seg_1 & obs <= seg_2) else index = which(obs >= seg_1)
    
    if(!is.na(p2)) {
      
      return(list(Qo_seg = obs[index], Qs_seg = sim[index], seg1 = seg_1, seg2 = seg_2))
      
    } else {
      
      return(list(Qo_seg = obs[index], Qs_seg = sim[index], seg1 = seg_1))
      
    }
    
  }
  
  Qoh = Qsh = Qom = Qsm = Qoi = Qsi = lQo_i1 = lQo_i2 = lQs_i1 = lQs_i2 = Qol = Qsl = lQol = lQsl = lQoL = lQsL = NA
  
  # High flows
  if(fitF == 'FHV' | fitF == 'FDC'){ 
    
    p_high = 0.02 # percentage of time that a flow is equalled or exceeded
    QsegsH = segmentFlows(obs, sim, p_high)
    Qoh = QsegsH[['Qo_seg']]
    Qsh = QsegsH[['Qs_seg']] 
    
  }
  
  # Medium flows
  if(fitF == 'FMV' | fitF == 'FDC'){ 
    
    p_medium1 = 0.2
    p_medium2 = 0.02
    QsegsM = segmentFlows(obs, sim, p_medium1, p_medium2)
    Qom = QsegsM[['Qo_seg']]
    Qsm = QsegsM[['Qs_seg']] 
    
  }
  
  # Intermediate flows  
  if(fitF == 'FMS' |fitF == 'FDC'){ 
    
    p_intermediate1 = 0.7 #m1
    p_intermediate2 = 0.2 #m2   
    QsegsI = segmentFlows(obs, sim, p_intermediate1, p_intermediate2)
    Qoi = QsegsI[['Qo_seg']]
    Qsi = QsegsI[['Qs_seg']]  
    
    if(fitF == 'FMS' |fitF == 'FDC'){
      
      Qo_i1 = quantile(obs, 1-p_intermediate1, na.rm = TRUE, names = FALSE)
      Qo_i2 = quantile(obs, 1-p_intermediate2, na.rm = TRUE, names = FALSE)
      Qs_i1 = quantile(sim, 1-p_intermediate1, na.rm = TRUE, names = FALSE)
      Qs_i2 = quantile(sim, 1-p_intermediate2, na.rm = TRUE, names = FALSE)
      
      ## Safe logarithm of Qo, Qs
      lQo_i1 = Qlog(Qo_i1)
      lQo_i2 = Qlog(Qo_i2)
      lQs_i1 = Qlog(Qs_i1)
      lQs_i2 = Qlog(Qs_i2)
      
    }
    
  }
  
  # Low flows  
  if(fitF == 'FLV' | fitF == 'FDC'){ 
    
    p_low1 = 0.9 ## WHY 1 ??? #l 
    p_low2 = 0.7 #L  
    if(p_low2 > p_low1 | (1 - p_low2) < 0.1) stop('Bad settings of FLV criteria (p_low1, p_low2)')
    QsegsL = segmentFlows(obs, sim, p_low1, p_low2)    
    Qol = QsegsL[['Qo_seg']]
    Qsl = QsegsL[['Qs_seg']]   
    
    ## Safe logarithm of Qol,L, Qsl,
    if(fitF == 'FLV' | fitF == 'FDC' ){
      
      QoL = QsegsL[['seg1']]
      QsL = quantile(sim, 1 - p_low1, na.rm = TRUE, names = FALSE)
      if(QsL == 0) QsL = quantile(sim[sim > 0], 0.1, na.rm = TRUE, names = FALSE)
      
      lQol = Qlog(Qol)
      lQsl = Qlog(Qsl)
      lQoL = Qlog(QoL)
      lQsL = Qlog(QsL)
      
    }
    
  }
  
  list(Qoh = Qoh, Qsh = Qsh, Qom = Qom, Qsm = Qsm, Qoi = Qoi, Qsi = Qsi, 
       lQo_i1 = lQo_i1, lQo_i2 = lQo_i2, lQs_i1 = lQs_i1, lQs_i2 = lQs_i2,
       Qol = Qol, Qsl = Qsl, 
       lQol = lQol, lQsl = lQsl, lQoL = lQoL, lQsL = lQsL)
  
}

## CED computation function
CEDcomputation = function(Qobs, Qsim){
  
  # empty vector
  vector.is.empty <- function(x) return(length(x) ==0 )
  
  Qobs = sort(Qobs, decreasing = TRUE)
  Qsim = sort(Qsim, decreasing = TRUE)
  
  Qsim_min=min(Qsim)
  Qsim_max=max(Qsim)
  
  Qsim2 =quantile(Qsim,0.98) # 2 percentile of simulated flows
  Qsim20=quantile(Qsim,0.80) # 20 percentile of simulated flows
  Qsim70=quantile(Qsim,0.30) # 70 percentile of simulated flows
  
  unscaled_entropy=function(start_range, end_range, Bins, tol) 
  {
    
    
    if(Bins==1)
    {
      Hu_obs = 0
      Hu_sim = 0
      
    } else {
      
      ftshistQobs=binposQobs=ftshistQsim=binposQsim=numeric()
      
      Bins_Qobs  = floor(min (Bins, (Qobs[ceiling(start_range*lenQo)]-Qobs[ceiling(end_range*lenQo)])/(2*tol)))
      Bins_Qsim  = floor(min (Bins, (Qsim[ceiling(start_range*lenQo)]-Qsim[ceiling(end_range*lenQo)])/(2*tol)))
      
      if(Bins_Qobs==0)
      {
        ftshistQobs=ftshistQobs
        binposQobs=binposQobs
      }else{
        a=min(Qobs[ceiling(start_range*lenQo):ceiling(end_range*lenQo)])
        b=max(Qobs[ceiling(start_range*lenQo):ceiling(end_range*lenQo)])
        w1=hist(Qobs[ceiling(start_range*lenQo):ceiling(end_range*lenQo)],seq(a,b,length.out=Bins_Qobs+1),plot=FALSE)
        ftshistQobs=w1[['counts']]
        binposQobs =w1[['mids']]
      }
      
      if(Bins_Qsim==0)
      {
        ftshistQsim=ftshistQsim
        binposQsim=binposQsim 
      }else{
        c=min(Qsim[ceiling(start_range*lenQo):ceiling(end_range*lenQo)])
        d=max(Qsim[ceiling(start_range*lenQo):ceiling(end_range*lenQo)])
        w2=hist(Qsim[ceiling(start_range*lenQo):ceiling(end_range*lenQo)], seq(c,d,length.out=Bins_Qsim+1),plot=FALSE)
        ftshistQsim=w2[['counts']]
        binposQsim=w2[['mids']]
      }
      
      
      Hu_obs = (-sum((ftshistQobs/sum(ftshistQobs))*log2(ftshistQobs/sum(ftshistQobs)),na.rm=TRUE)) / log2(Bins_Qobs)
      if (is.nan(Hu_obs)){Hu_obs=0}
      if (vector.is.empty(ftshistQsim)){ # In case all ftshistQsim values are 0
        Hu_sim = 0
      } else {
        Hu_sim = (-sum((ftshistQsim/sum(ftshistQsim))*log2(ftshistQsim/sum(ftshistQsim)),na.rm=TRUE)) / log2(Bins_Qsim)
      }
      if (is.nan(Hu_sim)){Hu_sim=0}
    }
    list(Hu_obs, Hu_sim)
  }
  
  
  scaled_entropy=function(quant1,quant2,quant1_sim,quant2_sim,start_range,end_range,Bins,tol) 
  {
    MAX = quant1+tol
    MIN = quant2-tol
    Bins = floor(min (Bins, (MAX-MIN)/(2*tol)))
    
    a=min(c(MIN,MAX))
    b=max(c(MIN,MAX))
    w1=hist(c(MIN,MAX),seq(a,b,length.out=Bins+1),plot=FALSE)
    ftshistQobs=w1[['counts']]
    binposQobs=w1[['mids']]
    
    binposQobsN = binposQobs - (binposQobs[2]-binposQobs[1])/2
    binposQobsN[length(binposQobsN)+1] = binposQobsN[length(binposQobsN)] + (binposQobs[2]-binposQobs[1])
    
    if (quant2_sim<MIN)
    {
      binposQobsN = c(quant2_sim,binposQobsN)
      Bins = Bins + 1
    }
    
    if (MAX < max(quant1_sim,MAX))
    {
      binposQobsN = c(binposQobsN,max(quant1_sim,MAX)+tol)
      Bins = Bins + 1
    }
    
    ftshistQobs=length(Bins)
    ftshistQsim=length(Bins)
    for (m in 1 : Bins)
    {
      ftshistQobs[m] = length(which((Qobs[ceiling(start_range*lenQo):ceiling(end_range*lenQo)] >=binposQobsN[m])&(Qobs[ceiling(start_range*lenQo):ceiling(end_range*lenQo)]<binposQobsN[m+1])))
      ftshistQsim[m] = length(which((Qsim>=binposQobsN[m])&(Qsim<binposQobsN[m+1])))
    }
    
    J=which(ftshistQobs<=0)
    if (vector.is.empty(I))
    {
      ftshistQobs = ftshistQobs
    } else {
      ftshistQobs = ftshistQobs[-J]
    }
    
    I=which(ftshistQsim<=0)
    if (vector.is.empty(I))
    {
      ftshistQsim = ftshistQsim
    } else {
      ftshistQsim  = ftshistQsim[-I]
    }
    
    Hs_obs = (-sum((ftshistQobs/sum(ftshistQobs))*log2(ftshistQobs/sum(ftshistQobs)),na.rm=TRUE)) / log2(Bins)
    if (is.nan(Hs_obs)){Hs_obs=0}
    if (vector.is.empty(ftshistQsim)) #In case all ftshistQsim values are 0
    {
      Hs_sim = 0
    } else {
      Hs_sim = (-sum((ftshistQsim/sum(ftshistQsim))*log2(ftshistQsim/sum(ftshistQsim)),na.rm=TRUE)) / log2(Bins)
    }
    if (is.nan(Hs_sim)){Hs_sim=0}
    list(Hs_obs, Hs_sim)  
  }
  
  
  # equal probable binning
  #scaled entropy
  scaled_entropy_low=function(quant1,quant2,quant1_sim,Bins,tol)
  {
    ftshistQobs=ftshistQsim=binposQobsN=as.numeric(vector())
    MAX = quant1 + tol
    MIN = max(0,quant2- tol)
    Bins = floor(min (Bins, (MAX-MIN)/(2*tol)))
    
    L_F = Qobs[ceiling(0.7*lenQo):length(Qobs)]
    SQobs = sort(L_F [L_F >=0])
    SQobs[1] = MIN
    SQobs[length(SQobs)] = MAX
    temp=floor(length(L_F)/Bins)
    temp1=length(SQobs)
    temp2=seq(1,temp1,by=temp)   
    if(Bins<=1){
      binposQobsN =binposQobsN  
    } else{
      binposQobsN = SQobs[temp2]
      binposQobsN = t(binposQobsN)
    }
    
    if (MIN > 0 && (min(Qsim)<MIN)){
      binposQobsN = c(0,binposQobsN)
      Bins = Bins + 1
    }
    
    if (MAX < max(quant1_sim,MAX)){
      binposQobsN = c(binposQobsN,max(Qsim70,MAX)+tol)
      Bins = Bins + 1
    }
    
    if (Bins+1!=length(binposQobsN)) {
      Bins
      length(binposQobsN)
    }
    
    
    for (m in 1 : Bins)
    {
      ftshistQobs[m] = length(which((L_F>=binposQobsN[m])&(L_F<binposQobsN[m+1])))
      ftshistQsim[m] = length(which((Qsim>=binposQobsN[m])&(Qsim<binposQobsN[m+1])))
    }
    
    J=which(ftshistQobs<=0)
    if (vector.is.empty(I))
    {
      ftshistQobs = ftshistQobs
    } else {
      ftshistQobs = ftshistQobs[-J]
    }
    
    I=which(ftshistQsim<=0)
    if (vector.is.empty(I))
    {
      ftshistQsim = ftshistQsim
    } else {
      ftshistQsim  = ftshistQsim[-I]
    }
    
    Hs_obs = (-sum((ftshistQobs/sum(ftshistQobs))*log2(ftshistQobs/sum(ftshistQobs)),na.rm=TRUE)) / log2(Bins)
    if (is.nan(Hs_obs)){Hs_obs=0}
    if (vector.is.empty(ftshistQsim)){ # In case all ftshistQsim values are 0
      Hs_sim = 0
    } else {
      Hs_sim = (-sum((ftshistQsim/sum(ftshistQsim))*log2(ftshistQsim/sum(ftshistQsim)),na.rm=TRUE)) / log2(Bins)
    }
    if (is.nan(Hs_sim)){Hs_sim=0}
    list(Hs_obs, Hs_sim) 
    
  }
  
  #SUSE
  SUSE = function(Hu_sim, Hu_obs, Hs_sim, Hs_obs){
    
    unscaled  = abs(Hu_sim - Hu_obs)
    scaled    = abs(Hs_sim - Hs_obs)
    criterion = max(unscaled,scaled)
    criterion
    
  }
  
  
  #computation
  #unscaled entropy
  
  Hu_obs_high=unscaled_entropy(1/lenQo,0.02,OptBins_high,0.001)[[1]]
  Hu_sim_high=unscaled_entropy(1/lenQo,0.02,OptBins_high,0.001)[[2]]
  
  Hu_obs_medium=unscaled_entropy(0.02,0.2,OptBins_medium,0.001)[[1]]
  Hu_sim_medium=unscaled_entropy(0.02,0.2,OptBins_medium,0.001)[[2]]
  
  Hu_obs_intermediate=unscaled_entropy(0.2,0.7,OptBins_intermediate,0.001)[[1]]
  Hu_sim_intermediate=unscaled_entropy(0.2,0.7,OptBins_intermediate,0.001)[[2]]
  
  Hu_obs_low=unscaled_entropy(0.7,1,OptBins_low,0.001)[[1]]
  Hu_sim_low=unscaled_entropy(0.7,1,OptBins_low,0.001)[[2]]
  
  #scaled entropy
  
  Hs_obs_high=scaled_entropy(Qobs_max,Qobs2,Qsim_max,Qsim2,1/lenQo,0.02,OptBins_high,0.001)[[1]]
  Hs_sim_high=scaled_entropy(Qobs_max,Qobs2,Qsim_max,Qsim2,1/lenQo,0.02,OptBins_high,0.001)[[2]]
  
  
  Hs_obs_medium=scaled_entropy(Qobs2,Qobs20,Qsim2,Qsim20,0.02,0.2,OptBins_medium,0.001)[[1]]
  Hs_sim_medium=scaled_entropy(Qobs2,Qobs20,Qsim2,Qsim20,0.02,0.2,OptBins_medium,0.001)[[2]]
  
  Hs_obs_intermediate=scaled_entropy(Qobs20,Qobs70,Qsim20,Qsim70,0.2,0.7,OptBins_intermediate,0.001)[[1]]
  Hs_sim_intermediate=scaled_entropy(Qobs20,Qobs70,Qsim20,Qsim70,0.2,0.7,OptBins_intermediate,0.001)[[2]]
  
  Hs_obs_low=scaled_entropy_low(Qobs70,Qobs_min,Qsim70,OptBins_low,0.001)[[1]]
  Hs_sim_low=scaled_entropy_low(Qobs70,Qobs_min,Qsim70,OptBins_low,0.001)[[2]]
  
  # SUSE
  SUSE_high         = SUSE(Hu_sim_high, Hu_obs_high, Hs_sim_high, Hs_obs_high)
  SUSE_medium       = SUSE(Hu_sim_medium, Hu_obs_medium, Hs_sim_medium, Hs_obs_medium)
  SUSE_intermediate = SUSE(Hu_sim_intermediate, Hu_obs_intermediate, Hs_sim_intermediate, Hs_obs_intermediate)
  SUSE_low          = SUSE(Hu_sim_low, Hu_obs_low, Hs_sim_low, Hs_obs_low)
  
  
  #CED
  CED = max(SUSE_high, SUSE_medium, SUSE_intermediate, SUSE_low)
  
  CED 
  
}


## CED computation function whole set
CEDcomputation_all = function(Qobs, Qsim){
  
  # empty vector
  vector.is.empty <- function(x) return(length(x) ==0 )
  
  unscaled_entropy = function(observed,simulated,NoBins){
    
    a = observed[observed >= 0]
    w1 = hist(a, seq(min(a), max(a), length.out = NoBins + 2), plot=FALSE) 
    ftshistobs = w1[['counts']]
    binposobs  = w1[['mids']]
    
    b = simulated[simulated>=0]
    w2 = hist(b, seq(min(b), max(b), length.out = NoBins + 2), plot=FALSE) 
    ftshistsim = w2[['counts']]
    binpossim  = w2[['mids']]
    
    J=which(ftshistobs<=0)
    if (vector.is.empty(I))
    {
      ftshistobs = ftshistobs
    } else {
      ftshistobs = ftshistobs[-J]
    }
    
    I=which(ftshistsim<=0)
    if (vector.is.empty(I))
    {
      ftshistsim = ftshistsim
    } else {
      ftshistsim  = ftshistsim[-I]
    }
    
    Hu_obs = (-sum((ftshistobs/sum(ftshistobs)) * log2(ftshistobs / sum(ftshistobs)),na.rm=TRUE)) / log2(NoBins)
    Hu_sim = (-sum((ftshistsim/sum(ftshistsim)) * log2(ftshistsim / sum(ftshistsim)),na.rm=TRUE)) / log2(NoBins)
    
    
    list(Hu_obs,Hu_sim)
    
  }
  
  scaled_entropy = function(observed,simulated,NoBins){
    
    ftshistobs = ftshistsim = as.numeric(vector())
    
    MAX = max(max(observed), max(simulated))
    a = observed[observed >= 0]
    b = simulated[simulated >= 0]
    MIN = min(min(a),min(b))
    
    vec = c(MIN,MAX)
    w3 = hist(vec, seq(min(vec), max(vec), length.out=NoBins+2), plot=FALSE)
    ftshistobs = w3[['counts']]
    binposobs = w3[['mids']]
    binposobsN = binposobs - (binposobs[2] - binposobs[1]) / 2
    end = length(binposobsN)
    binposobsN[end] = binposobsN[end] * 1.00001
    
    for (m in 1 : length(binposobsN)-1){
      
      ftshistobs[m] = length(which((observed >= binposobsN[m]) & (observed < binposobsN[m+1])))
      ftshistsim[m] = length(which((simulated >= binposobsN[m]) & (simulated < binposobsN[m+1])))
      
    }
    
    J=which(ftshistobs<=0)
    if (vector.is.empty(I))
    {
      ftshistobs = ftshistobs
    } else {
      ftshistobs = ftshistobs[-J]
    }
    
    I=which(ftshistsim<=0)
    if (vector.is.empty(I))
    {
      ftshistsim = ftshistsim
    } else {
      ftshistsim  = ftshistsim[-I]
    }
    
    #scaled entropy
    Hs_obs = (-sum((ftshistobs / sum(ftshistobs)) * log2(ftshistobs / sum(ftshistobs)), na.rm=TRUE)) / log2(NoBins)
    Hs_sim = (-sum((ftshistsim / sum(ftshistsim)) * log2(ftshistsim / sum(ftshistsim)), na.rm=TRUE)) / log2(NoBins)
    
    list(Hs_obs,Hs_sim)
    
  }
  
  #SUSE
  SUSE = function(Hu_sim, Hu_obs, Hs_sim, Hs_obs){
    
    unscaled  = abs(Hu_sim - Hu_obs)
    scaled    = abs(Hs_sim - Hs_obs)
    criterion = max(unscaled,scaled)
    criterion
    
  }
  
  #computation
  #unscaled entropy
  Hu_obs        = unscaled_entropy(Qobs, Qsim, OptBins)[[1]]
  Hu_sim        = unscaled_entropy(Qobs, Qsim, OptBins)[[2]]
  
  
  #scaled entropy
  Hs_obs        = scaled_entropy(Qobs, Qsim, OptBins)[[1]]
  Hs_sim        = scaled_entropy(Qobs, Qsim, OptBins)[[2]]
  
  
  CED        = SUSE(Hu_sim, Hu_obs, Hs_sim, Hs_obs)
  
  CED
  
}

#%BiasTlag

Find_Max_CCF<- function(a,b)
{
  d <- ccf(a, b, plot = FALSE)
  cor = d[['acf']][,,1]
  lag = d[['lag']][,,1]
  res = data.frame(cor,lag)
  res_max = res[which.max(res[['cor']]),]
  return(res_max)
} 

## FitnessComputation: Computation of fitness

FitnessComputation = function(computedValues, originalValues, fitnessFunction = FitnessFunction, testing = FALSE, PforTlag = NA){
  
  if(all(is.na(computedValues))) return(NA)
  
  ## FOR COMPUTATION OF WEI, ORIGINAL Q DATA HAS ZEROES IN IT!!!
  if((fitnessFunction == 'MARE' | fitnessFunction == 'Vis_2'|fitnessFunction == 'rel_d0'|fitnessFunction == 'md0') & any(originalValues == 0)){
    
    originalValues = originalValues + quantile(originalValues[originalValues > 0], 0.1)
    
  }
  
  ## Setting of Global fitness values, dependent on testing
  if(testing) FitGlobs = FitnessGlobal(originalValues, assignIt = TRUE)
  
  ## Threshold of Q, for simulation of Q only!!!
  computedValues[computedValues < 0] = 0
  
  ## Constant result modification
  if(length(computedValues) == 1) computedValues = rep(computedValues, length(originalValues))
  
  computedValuesLow = computedValues[LowValPositions]
  computedValuesHigh = computedValues[HighValPositions]
  computedValuesMedium = computedValues[MediumValPositions]
  
  difference = originalValues - computedValues
  differenceLow = OriginalValuesLow - computedValuesLow
  differenceHigh = OriginalValuesHigh - computedValuesHigh
  differenceMedium = OriginalValuesMedium - computedValuesMedium
  
  if(fitnessFunction == 'logNS0' | fitnessFunction == 'Mai0' | fitnessFunction == 'Vis_1'|fitnessFunction == 'Shafii'){
    
    # logComputedValues = Qlog(computedValues)
    # logDifference = logOriginalValues - logComputedValues
    logComputedValues = Qlog(computedValues + quantile(originalValues[originalValues > 0], 0.1))
    #     logOriginalValues = Qlog(originalValues + quantile(originalValues[originalValues > 0], 0.1))
    logDifference = logOriginalValues - logComputedValues
    
    
  }
  
  if(fitnessFunction == 'FHV' | fitnessFunction == 'FMV' | fitnessFunction == 'FMS' | 
     fitnessFunction == 'FLV' | fitnessFunction == 'FDC'){
    
    fcp = FcriteriaPrecomputation(originalValues, computedValues, fitnessFunction)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(fcp), fcp)
    
  }
  
  # P definition for Tlag computation
  if(!testing & fitnessFunction == 'Tlag') {
    
    if(!'P' %in% IndependentVariables) stop('Tlag fitness function needs a P variable (precipitation) in set of independent variables')
    P = DataSet[['P']]
    
  }
  if(testing & fitnessFunction == 'Tlag') P = PforTlag
  
  if(fitnessFunction == 'Shafii'){
    SumP=sum(DataSet[['P']]) #?????????
    #orrOriginalValues=sum(DataSet[['Q']])/SumP  #?????????????
    orrOriginalValues=sum(originalValues)/SumP
    orrlogOriginalValues=sum(logOriginalValues)/SumP  #??????????
    orrComputedValues= sum(computedValues)/SumP  #1
    orrlogComputedValues= sum(logComputedValues)/SumP  #2
    indexm1=which.min(abs(computedValues-(as.numeric(quantile(computedValues, probs = 0.2, na.rm = TRUE)))))
    computedm1=computedValues[indexm1]
    indexm2=which.min(abs(computedValues-(as.numeric(quantile(computedValues, probs = 0.7, na.rm = TRUE)))))
    computedm2=computedValues[indexm2]
    fdcmComputedValues=Qlog10(computedm2)-Qlog10(computedm1) #3
    fdchComputedValues=sum(computedValuesHigh) #4 
    computedl2=as.numeric(quantile(computedValues, probs = 0, na.rm = TRUE))
    fdclComputedValues=-1*(sum(Qlog10(computedValuesLow)-Qlog10(computedl2))) #5
    MeanComputedValues=mean(computedValues) #6
    VarComputedValues=sqrt(sum((computedValues-MeanComputedValues)^2)/(orlength-1))#7
    MedianComputedValues=median(computedValues) #8
    PeakComputedValues=max(computedValues) #9
    lag_num_Computed=numeric(orlength-1)
    for (i in 1:(orlength-1))
    {
      lag_num_Computed[i]=(computedValues[i]-MeanComputedValues)*(computedValues[i+1]-MeanComputedValues)
    }
    Lag_ac_Computed=sum(lag_num_Computed)/sum((computedValues-MeanComputedValues)^2) #10
    MeanlogComputedValues=mean(logComputedValues) #11
    VarlogComputedValues=sqrt(sum((logComputedValues-MeanlogComputedValues)^2)/(orlength-1)) #12
    Month_1=months(as.Date(DataSet[['Date']])) # ??????????????
    MaxMonthlyOriginal=max(aggregate( DataSet[['Q']] ~ Month_1 , DataSet , mean )[2]) #??????????????
    MaxMonthlyComputed=max(aggregate( computedValues ~ Month_1 , DataSet , mean )[2]) #13
    accept_tresh = .05 #acceptability treshold, choose [0.05,0.1 or 0.2] -> Shafii 2015
    sign_1   = 1-as.integer(as.logical((orrOriginalValues*(1-accept_tresh) <= orrComputedValues & orrOriginalValues*(1+accept_tresh) >= orrComputedValues)))
    sign_2   = 1-as.integer(as.logical((orrlogOriginalValues*(1-accept_tresh) <= orrlogComputedValues & orrlogOriginalValues*(1+accept_tresh) >= orrlogComputedValues)))
    sign_3   = 1-as.integer(as.logical((fdcmOriginalValues*(1-accept_tresh) <= fdcmComputedValues & fdcmOriginalValues*(1+accept_tresh) >= fdcmComputedValues)))
    sign_4   = 1-as.integer(as.logical((fdchOriginalValues*(1-accept_tresh) <= fdchComputedValues & fdchOriginalValues*(1+accept_tresh) >= fdchComputedValues)))
    sign_5   = 1-as.integer(as.logical((fdclOriginalValues*(1-accept_tresh) <= fdclComputedValues & fdclOriginalValues*(1+accept_tresh) >= fdclComputedValues)))
    sign_6   = 1-as.integer(as.logical((MeanOriginalValues*(1-accept_tresh) <= MeanComputedValues & MeanOriginalValues*(1+accept_tresh) >= MeanComputedValues)))
    sign_7   = 1-as.integer(as.logical((VarOriginalValues*(1-accept_tresh) <= VarComputedValues & VarOriginalValues*(1+accept_tresh) >= VarComputedValues)))
    sign_8   = 1-as.integer(as.logical((MedianOriginalValues*(1-accept_tresh) <= MedianComputedValues & MedianOriginalValues*(1+accept_tresh) >= MedianComputedValues)))
    sign_9   = 1-as.integer(as.logical((PeakOriginalValues*(1-accept_tresh) <= PeakComputedValues & PeakOriginalValues*(1+accept_tresh) >= PeakComputedValues)))
    sign_10  = 1-as.integer(as.logical((Lag_ac_Original*(1-accept_tresh) <= Lag_ac_Computed & Lag_ac_Original*(1+accept_tresh) >= Lag_ac_Computed)))
    sign_11  = 1-as.integer(as.logical((MeanlogOriginalValues*(1-accept_tresh) <= MeanlogComputedValues & MeanlogOriginalValues*(1+accept_tresh) >= MeanlogComputedValues)))
    sign_12  = 1-as.integer(as.logical((VarlogOriginalValues*(1-accept_tresh) <= VarlogComputedValues & VarlogOriginalValues*(1+accept_tresh) >= VarlogComputedValues)))
    sign_13  = 1-as.integer(as.logical((MaxMonthlyOriginal*(1-accept_tresh) <= MaxMonthlyComputed & MaxMonthlyOriginal*(1+accept_tresh) >= MaxMonthlyComputed)))
    
  }
  beta = mean(computedValues)/MeanOriginalValues
  #index_cor =  which(originalValues! = computedValues)
  sdo = sd(originalValues)
  sdc = sd(computedValues)
  if(is.na(sdo) | sdo == 0) sdo = .Machine[['double.xmin']] 
  if(is.na(sdc) | sdc == 0) sdc = .Machine[['double.xmin']] 
  alpha = sdc/sdo
  gamma = (sdc/mean(computedValues))/(sdo/MeanOriginalValues)
  
  r      = function()  cov(originalValues, computedValues) / (sdo * sdc)
  r0     = function()  1 - r() ## Wrong criteria - everything must go to zero. H(r) = -1 ,1 ??? !!!
  RMSE   = function()  sqrt(mean((difference)^2))
  NS0    = function()  sum(difference^2) / sum(NS0denominator^2) #14
  logNS0 = function()  sum(logDifference^2) / sum(logNS0denominator^2) #15
  RSq0   = function()  1 - r()^2
  PI0    = function()  sum((difference[2:length(originalValues)])^2) / sum((originalValues[2:length(originalValues)] - originalValues[1:(length(originalValues)-1)])^2)
  MAE    = function()  abs(mean(difference))
  MARE   = function()  ((1/length(originalValues))*(sum(abs(difference)/originalValues)))
  KG     = function(x) sqrt((r() - 1)^2 + (x - 1)^2 + (beta - 1)^2)
  VE0    = function()  sum(abs(difference)) / sum(originalValues)
  MNS0   = function() sum(abs(difference))/sum(abs(originalValues-MeanOriginalValues))
  RSD0   = function() 1-alpha
  CED0   = function() CEDcomputation(originalValues, computedValues)
  
  #yilmaz(2008): %Bias RR, FHV, FLV, FMM , FMS and Tlag
  fitnessDefinition = c(
    MARE       = 'fit = MARE()', #Mean Absolute Relative Error
    MAE        = 'fit = MAE()', # Mean Absolute Error
    RMSE       = 'fit = RMSE()', # Root Mean Square Error
    NS0        = 'fit = NS0()', #14
    logNS0     = 'fit = logNS0()', #15
    KG10       = 'fit = KG(alpha)', # KGE_2009
    KG20       = 'fit = KG(gamma)', # KGE_2012
    KGE1       = 'fit = 1-KG(alpha)', # KGE_2009
    KGE2       = 'fit = 1-KG(gamma)', # KGE_2012
    MNS0       = 'fit = MNS0()',# Modified NSE
    RSD0       = 'fit = RSD0()', # Ratio of Standard Deviation
    RSq0       = 'fit = RSq0()', # Coefficient of Determination
    PI0        = 'fit = PI0()', # Persistence Index
    VE0        = 'fit = VE0()', # Volumetric Efficiency
    r0         = 'fit = r0()', # Pearson correlation coefficient
    FHV        = 'fit = ((sum(Qsh-Qoh))/(sum(Qoh)))*100',  # remove abs : abs(((sum(Qsh-Qoh))/(sum(Qoh)))*100) %BiasFHV : Yilmaz
    FMV        = 'fit = ((sum(Qsm-Qom))/(sum(Qom)))*100', #remove abs %BiasFMV : Ley
    FMS        = 'fit = (((lQs_i1-lQs_i2)-(lQo_i1-lQo_i2))/(lQo_i1-lQo_i2))*100', #remove abs %BiasFMS : Yilmaz
    FLV        = 'fit = -1*((sum(lQsl-lQsL)-sum(lQol-lQoL))/(sum(lQol-lQoL)))*100', #remove abs %BiasFLV : Yilmaz
    FMM        = 'fit = ((log(median(computedValues))-log(median(originalValues)))/log(median(originalValues)))*100', #%BiasFMM : Yilmaz
    Tlag       = 'fit = ((Find_Max_CCF(computedValues,P)[["lag"]]-Find_Max_CCF(originalValues,P[["lag"]])/(Find_Max_CCF(originalValues,P)[["lag"]]))*100', #%BiasTlag
    RR         = 'fit = (sum(difference)/sum(originalValues))*100', #%BiasRR
    CED        = 'fit = CEDcomputation(originalValues, computedValues)',
    SUSE       = 'fit = CEDcomputation_all(originalValues, computedValues)',
    lowRMSE    = 'fit = sqrt(mean(differenceLow^2))',
    highRMSE   = 'fit = sqrt(mean(differenceHigh^2))',
    mediumRMSE = 'fit = sqrt(mean(differenceMedium^2))',
    pbias      = 'fit = (sum(difference)/sum(originalValues))*100', # Percentage Bias 
    md0        = 'fit = sum(abs(difference))/sum(abs(computedValues-MeanOriginalValues),(abs(originalValues-MeanOriginalValues)))', # modified index of agreement
    
    Borsanyi   = 'fit = sqrt((r0())^2+(NS0())^2+(VE0())^2+(KG(alpha))^2)',
    Mai0       = 'fit = sqrt(NS0()^2 + logNS0()^2)',
    Vis_1      = 'fit = sqrt(NS0()^2 + logNS0()^2 + VE0()^2)',  #Vis et al.,2015
    NSE        = 'fit = 1-NS0()',
    r          = 'fit = r()',
    Vis_3      = 'fit = sqrt((r0())^2+ (VE0())^2)', #Vis et al.,2015
    Price      = 'fit = sqrt(NS0()^2 + (MNS0())^2 + RSD0()^2)',
    CED_new    = 'fit = sqrt((CED0())^2 +(KG(alpha))^2)',
    rel_d0     = 'fit = sum((difference/originalValues)^2)/sum((sum(abs(computedValues-MeanOriginalValues),(abs(originalValues-MeanOriginalValues)))/MeanOriginalValues)^2)', #relative index of agreement
    Shafii     = 'fit = NS0()+logNS0()+sign_1+sign_2+sign_3+sign_4+sign_5+sign_6+sign_7+sign_8+sign_9+sign_10+sign_11+sign_12+sign_13' #Shafii 2015
    #AIC       = 'fit = (length(originalValues)*log(RMSE()))+(2*p)', #p: number of free parameters???
    #BIC       = 'fit = (length(originalValues)*log(RMSE()))+(p*log(length(originalValues)))',
    #C2M       = 'fit = NSE/(2-NSE)' #Mathevet et al.,2006,
  )
  
  if(fitnessFunction == 'Multi_Madsen'){
    
    fitOut = data.frame(
      lowRMSE    = eval(parse(text = fitnessDefinition['lowRMSE'])), 
      RMSE       = eval(parse(text = fitnessDefinition['RMSE'])),
      highRMSE   = eval(parse(text = fitnessDefinition['highRMSE'])), 
      MAE        = eval(parse(text = fitnessDefinition['MAE']))
    )
    
  } 
  
  else if (fitnessFunction == 'Dawson'){
    #     Dawson     = 'fit = sqrt((RMSE())^2+(RSq0())^2+(PI0())^2+(MAE())^2)', 
    fitOut = data.frame(
      RMSE       = eval(parse(text = fitnessDefinition['RMSE'])),
      RSq0       = eval(parse(text = fitnessDefinition['RSq0'])),
      PI0        = eval(parse(text = fitnessDefinition['PI0'])),
      MAE        = eval(parse(text = fitnessDefinition['MAE']))
    )
    
  } else if (fitnessFunction == 'Vis_2'){
    #     Vis_2      = 'fit = sqrt((r0())^2+(NS0())^2+ (VE0())^2 +(MARE())^2)', #Vis et al.,2015
    fitOut = data.frame(
      r0      = eval(parse(text = fitnessDefinition['r0'])),
      NS0     = eval(parse(text = fitnessDefinition['NS0'])),
      VE0     = eval(parse(text = fitnessDefinition['VE0'])),
      MARE    = eval(parse(text = fitnessDefinition['MARE']))
    )  
    
  } 
  
  else {
    
    fitOut = eval(parse(text = fitnessDefinition[fitnessFunction]))
    
  }
  
  if(fitnessFunction %in% MultiObjectiveFitnesses){
    
    fitOut = sum(scale(fitOut, center = FALSE, scale = FALSE)) ## THIS IS ONLY TEMPORARY, MUST BE PROBABLY DONE IN BETTER WAY - MULTIOBJECTIVE PROBLEM ??? !!!
    
  }
  
  fitOut
  
}

## Multi-objective criteria computation (Eucledian Distance also can be used!!!???)

MultiObjective = function(population, changeCheck = FALSE){
  
  ## BALANCED FITNESS  
  fitList = lapply(population, function(x) x[['Fitness']])
  fitDF = do.call('rbind', fitList)
  
  colMins = apply(fitDF, 2, min, na.rm = TRUE)
  
  A = max(colMins) - colMins
  
  balancedFitness = apply(fitDF, 1, function(x) sqrt(sum((x + A)^2)))
  
  balancedFitness 
  
}


## MakeEquation: Equation compilation from individual array

MakeEquation = function(individual){
  ReservoirSet <-c("WR","IR","RR","UR","FR","SR","CR")
  IA = individual[['IndArray']]
  
  # Vyreseni jedincu slozenych pouze z terminalu (nebo funkce s aritou 0 - zatim nejsou)
  if(individual[['IndLength']] == 1){
    
    return(paste0('y = ', IA[, 1]))
    
  }
  
  eq_vec = IA[, 1]
  
  for (i in individual[['IndLength']]:1){
    
    if( any(FunctionSet == IA[i, 1]) | any(ReservoirSet == IA[i, 1])){
      
      for(j in 2:ncol(IA)){
        
        if( !is.na(IA[i,j]) ){
          
          arg_row = as.numeric(IA[i,j])
          
          if( j == 2 && is.na(IA[i, (j + 1)]) ){
            
            if(!any(IA[i, 1] == SpecialFunctions)){
              
              eq_vec[i] = paste0(eq_vec[i], '(', eq_vec[arg_row])
              
            }else{
              
              eq_vec[i] = paste0(eq_vec[i], eq_vec[arg_row])
              
            }
            
          }
          
          if( j == 2 && !is.na(IA[i, (j + 1)]) ){
            
            if( any(PrefixFunctions == IA[i, 1]) ){
              
              first = paste0(eq_vec[i],'(', eq_vec[arg_row], ',')
              
            } else{
              
              first = paste0('(', eq_vec[arg_row])
              
            }
            
          }
          
          if(j > 2){
            
            if( any(PrefixFunctions == IA[i, 1]) ){
              
              eq_vec[i] = paste0(first, eq_vec[arg_row])
              
              if(j >= 3) first = paste0(eq_vec[i], ',')
              
            } else{
              
              eq_vec[i] = paste0(first, eq_vec[i], eq_vec[arg_row])
              
            }
            
          }
          
        }
        
      } # j cycle end
      
      eq_vec[i] = paste0(eq_vec[i], ')')
      
    }
    
  } # i cycle end
  
  paste0('y = ', eq_vec[1])
  
}

## ComputeIndividual: Computation of individual
# individual - individual for computation
# compVar - variables for computation

ComputeIndividual = function(individual){
  
  if(PunishOneNodeIndividuals){
    
    if(individual[['IndLength']] == 1){
      
      individual[['Fitness']] = Inf
      individual[['Changed']] = FALSE
      return(individual)
      
    }
    
  }
  
  # Transfer of individual array to equation
  individual[['Equation']] = MakeEquation(individual)
  
  # Equation evaluation
  nrDataSet = nrow(DataSet)
  computedValues = suppressWarnings(try(eval(parse(text = individual[['Equation']]), envir = DataSet), silent = TRUE))
  
  # Identification of nonsenses
  if(class(computedValues) == 'try-error' | any(!is.finite(computedValues))){ 
    
    computedValues = NA 
    
  }
  
  # Extension of constant results to length of training set
  if(length(computedValues) == 1 & !is.na(computedValues[1])){
    
    computedValues = rep(computedValues, nrDataSet)
    
  }
  individual[['Sim']] = computedValues
  # Fitness computation
  originalValues = DataSet[[DependentVariable]]
  for (i in 1:length(FitnessFunction)){
    individual[['Fitness']][i] = FitnessComputation(computedValues, originalValues,fitnessFunction = FitnessFunction[i]) 
  }
  for (i in 1:length(FitnessFunction)){
    if (is.na(individual[['Fitness']][i])){
      individual[['Fitness']]<- rep(10000,length(FitnessFunction))
    }
    if (is.infinite(individual[['Fitness']][i])){
      individual[['Fitness']] <- rep(10000,length(FitnessFunction))
    }
  }
  
  individual[['Changed']] = FALSE
  
  individual
  
}

# source('CHMGP_VariationOperators.R')

PopulationCheck = function(goodFitnesses, tournamentSize, generation){
  
  if(length(goodFitnesses) < (2 * tournamentSize)){
    
    warning(paste0('PopulationCheck: Too many NA fitnesses in population! Variability operators ineffective. Generation: ', generation))
    
  }
  
}

######## need to a change here to add crowding distance and front into cosideration of selecting the winner ############
TournamentSelection = function(goodFitnesses, tournamentSize, oneWinner = FALSE){
  
  if(oneWinner){
    
    competitors = Sample(goodFitnesses, tournamentSize)
    return(names(competitors[which.min(competitors)]))
    
  }
  
  noCompetitors = 2 * tournamentSize
  
  competitors = Sample(goodFitnesses, noCompetitors)
  
  groupA = head(competitors, tournamentSize)
  groupB = tail(competitors, tournamentSize)
  
  winners = c(groupA[which.min(groupA)], groupB[which.min(groupB)])
  
  names(winners)
  
}

NodeSelection = function(driveVector, restrictedRowsVector, terminalNodeProbability){
  
  if(length(driveVector) == 1) return(list(targetNode = 1, restriction = FALSE))
  
  # Vector of probabilities for terminal type determination (Koza's recommendation is P(function) = 0.1, P(terminal) = 0.9)
  # terminalNodeProbability = 0.9   
  selectedTerminal = as.logical(rbinom(1, 1, terminalNodeProbability))
  
  if(selectedTerminal){  
    
    targetNode = Sample(which(driveVector == 'Terminal' ), 1)
    
  } else{
    
    targetNode = Sample(which(driveVector == 'Function' ), 1)
    
  }
  
  
  if(targetNode %in% restrictedRowsVector) return(list(targetNode = targetNode, restriction = TRUE))
  
  list(targetNode = targetNode, restriction = FALSE)
  
}

RecursiveNodePointers = function(treeArray, checkRow, nodesPointers = c()){
  
  actualPointers = numeric(0)
  
  if(length(checkRow) > 0){
    
    for(i in 1:length(checkRow)){
      
      actualPointers = c(actualPointers, as.numeric(treeArray[checkRow[i], !is.na(treeArray[checkRow[i], ])]))
      
    }
    
    nodesPointers = c(nodesPointers, actualPointers)
    
    checkRow = actualPointers
    
    RecursiveNodePointers(treeArray, checkRow, nodesPointers)
    
  } else {
    
    return(nodesPointers)
    
  }
  
}

RecursiveRestrictedNodes = function(treeArray, checkRow, nodes, nodesPointers = c()){
  
  actualPointers = numeric(0)
  
  if(length(checkRow) > 0){
    
    for(i in 1:length(checkRow)){
      
      node = nodes[checkRow[i]]
      noOfNotRestrictedNodes = NonStrictArguments[node]
      
      identifiedPointers = !is.na(treeArray[ checkRow[i], ])
      
      identifiedPointers[1:noOfNotRestrictedNodes] = FALSE
      
      #       print(identifiedPointers)
      
      #       print(nodes[as.numeric(treeArray[checkRow[i], 1:noOfNotRestrictedNodes])] %in%  c('+', '-'))
      
      againRestricted = nodes[as.numeric(treeArray[checkRow[i], 1:noOfNotRestrictedNodes])] %in% StandardFunctions 
      
      identifiedPointers[againRestricted] = TRUE
      
      #       print(identifiedPointers)
      
      actualPointers = c(actualPointers, as.numeric(treeArray[checkRow[i], identifiedPointers]))
      
      #       notRestrictedNodes = NonStrictArguments[nodes[checkRow[i]]]
      #       identifiedPointers = !is.na(treeArray[ checkRow[i], ])
      #       identifiedPointers[1:notRestrictedNodes] = FALSE
      
      #       actualPointers = c(actualPointers, as.numeric(treeArray[checkRow[i], identifiedPointers]))
      
    }
    
    nodesPointers = c(nodesPointers, actualPointers)
    
    checkRow = actualPointers
    
    RecursiveNodePointers(treeArray, checkRow, nodesPointers)
    
  } else {
    
    return(nodesPointers)
    
  }
  
}

ChangeIndividual = function(newIndArray, newDriveVector, newLength, change){
  
  
  
  strictFpositions = which(newIndArray[, 1] %in% StrictFunctions)
  
  if(length(strictFpositions) > 0){
    
    restrictedRows = unique(RecursiveRestrictedNodes(newIndArray[, 2:MaxNoColumns], strictFpositions, newIndArray[, 1]))
    
    #     if(newLength > 10) browser()
    
    newIndArray[restrictedRows[newIndArray[restrictedRows, 1] %in% IndependentVariables], 1] = runif(1, ConstantRange[1], ConstantRange[2])
    
  } else {
    
    restrictedRows = NA
    
  }
  
  # Osetrit mista s promennymi tam kde nemaji byt - nahrada za konstanty
  
  #   argsOfSf = NonStrictArguments[newIndArray[strictFpositions, 1]]
  #   
  #   newIndArray[strictFpositions , (argsOfSf+1)]
  
  
  #   restrictedRows = as.numeric(as.vector(newIndArray[newIndArray[ , 1] %in% StrictFunctions, 3:ncol(newIndArray)]))
  #   if(length(restrictedRows) == 0) restrictedRows = NA
  #   if(length(restrictedRows) > 1 & min(restrictedRows, na.rm = TRUE) > 2) browser()
  
  list(IndArray = newIndArray, Equation = NA, Fitness = c(rep(NA,length(FitnessFunction))), DriveVec = newDriveVector, RestrictedRows = restrictedRows, 
       IndLength = newLength, Changed = change, Front = NA, Crowding_distance = NA,Sim =c())
  
}

AssignPointers = function(pointersArray, arrayLength){
  
  transposedPA = t(pointersArray)
  transposedPA[which(!is.na(transposedPA), arr.ind = TRUE)] = as.character(2:arrayLength)
  t(transposedPA)
  
}

SepareTrees = function(parent, cNode){
  
  if(!cNode[['variationAllowed']]){
    
    return(NA)
    
  }
  
  cn = cNode[['targetNode']]
  
  if(cn == 1){
    
    return(parent)
    
  }
  
  if(parent[['DriveVec']][cn] == 'Terminal'){
    
    newDriveVector = parent[['DriveVec']][cn] 
    
    newIndArray = array(parent[['IndArray']][cn, ], dim = c(1, MaxNoColumns))
    
    return(ChangeIndividual(newIndArray, newDriveVector, 1, TRUE))
    
  }
  
  if(parent[['DriveVec']][cn] == 'Function'){
    
    # Determination of separated nodes by pointers of node
    separatedRows = RecursiveNodePointers(parent[['IndArray']][ , 2:MaxNoColumns], cn, cn)
    
  }
  
  # Array and driving vector of separated tree
  newLength = length(separatedRows)
  
  newDriveVector = parent[['DriveVec']][separatedRows]  
  newIndArray = array(parent[['IndArray']][separatedRows, ], dim = c(newLength, MaxNoColumns))
  
  newIndArray[, 2:MaxNoColumns] = AssignPointers(newIndArray[, 2:MaxNoColumns], newLength)
  
  ChangeIndividual(newIndArray, newDriveVector, newLength, TRUE)  
  
}

TerminalToFunction = function(parent, cNode, subtree){
  
  # Determination of deleted subtree identified by pointers of deleting node
  removeRows = RecursiveNodePointers(parent[['IndArray']][ , 2:MaxNoColumns], cNode)
  
  newLength = parent[['IndLength']] - length(removeRows)
  
  newDriveVector = parent[['DriveVec']][-removeRows] 
  newDriveVector[cNode] = 'Terminal'
  
  newIndArray = array(parent[['IndArray']][-removeRows, ], dim = c(newLength, MaxNoColumns))
  newIndArray[cNode, ] = subtree[['IndArray']][1, ]
  
  newIndArray[, 2:MaxNoColumns] = AssignPointers(newIndArray[, 2:MaxNoColumns], newLength)
  
  ChangeIndividual(newIndArray, newDriveVector, newLength, TRUE)
  
}

ChangePointers = function(subtreeArray, lastPointer, startRow, stopRow){
  
  changedArray = subtreeArray[startRow:stopRow, 2:MaxNoColumns]
  
  transposedChA = t(changedArray)
  pointersPositions = which(!is.na(transposedChA), arr.ind = TRUE)
  newPointers = as.character((lastPointer + 1):(nrow(pointersPositions) + lastPointer))
  transposedChA[pointersPositions] = newPointers
  subtreeArray[startRow:stopRow, 2:MaxNoColumns] = t(transposedChA)  
  
  subtreeArray
  
}

FunctionToTerminal = function(parent, cNode, subtree){
  
  newLength = parent[['IndLength']] + subtree[['IndLength']] - 1
  
  newIndArray = array(NA, dim = c(newLength, MaxNoColumns))
  
  newDriveVector = rep(NA, newLength)
  
  # Filling new array and driving vector
  newIndArray[1:parent[['IndLength']], ] = parent[['IndArray']]
  newDriveVector[1:parent[['IndLength']]] = parent[['DriveVec']]
  
  subtreeCopy = subtree
  
  # Change of pointers
  lastMaxPointer = max(as.numeric(parent[['IndArray']][1:(cNode-1), 2:MaxNoColumns]), na.rm = TRUE)
  subtree[['IndArray']] = ChangePointers(subtree[['IndArray']], lastMaxPointer, 1, subtree[['IndLength']])
  
  # Addition of subtree root node
  newIndArray[cNode, ] = subtree[['IndArray']][1, ]
  newDriveVector[cNode] = 'Function'
  
  # Last maximal value of pointers
  lastMaxPointer = max(as.numeric(newIndArray[1:cNode, 2:MaxNoColumns]), na.rm = TRUE)
  
  # Change of original pointers by new function addition
  if(any(!is.na(newIndArray[(cNode + 1):newLength, 2:MaxNoColumns]))){
    
    newIndArray = ChangePointers(newIndArray, lastMaxPointer, (cNode + 1), newLength)
    
  }
  
  flowLength = parent[['IndLength']]
  
  for(i in 1:subtree[['IndLength']]){
    
    j = 2
    while(j <= MaxNoColumns && !is.na(subtreeCopy[['IndArray']][i, j])){
      
      subtreePointer = as.numeric(subtree[['IndArray']][i, j])
      subtreeCopyPointer = as.numeric(subtreeCopy[['IndArray']][i, j]) 
      
      if(!is.na(newDriveVector[subtreePointer])){
        
        flowLength = flowLength + 1
        
        newIndArray[(subtreePointer + 1):flowLength, ] = newIndArray[subtreePointer:(flowLength - 1), ]
        newDriveVector[(subtreePointer + 1):flowLength] = newDriveVector[subtreePointer:(flowLength - 1)]     
        
      }
      
      newIndArray[subtreePointer, ] = subtree[['IndArray']][subtreeCopyPointer, ]
      newDriveVector[subtreePointer] = subtree[['DriveVec']][subtreeCopyPointer] 
      
      if(any(subtree[['DriveVec']][subtreeCopyPointer:subtree[['IndLength']]] == 'Function', na.rm = TRUE)){
        
        lastMaxPointer = max(as.numeric(newIndArray[1:(subtreePointer - 1), 2:MaxNoColumns]), na.rm = TRUE)
        newIndArray = ChangePointers(newIndArray, lastMaxPointer, subtreePointer, newLength)
        subtree[['IndArray']] = ChangePointers(subtree[['IndArray']], lastMaxPointer, subtreeCopyPointer, subtree[['IndLength']])
        
      }
      
      j = j + 1
      
    }
    
  }
  
  ChangeIndividual(newIndArray, newDriveVector, newLength, TRUE)
  
}

AddSubtree = function(parent, cNode, subtree){
  
  # In the case, when the crossover node is not allowed, the parent is returned
  # Funguje, ale jestli radsi nezavolat mutaci na dany uzel, aby prohledavani nestagnovalo???
  if(!cNode[['variationAllowed']]){
    
    return(parent)
    
  }
  
  cn = cNode[['targetNode']]
  
  if(cn == 1){
    
    return(subtree)
    
  }
  
  parentConnectionNode = parent[['DriveVec']][cn]
  subtreeConnectionNode = subtree[['DriveVec']][1]
  
  # Case 1: crossover node in parent is terminal and subtree root is terminal
  if(parentConnectionNode == 'Terminal' & subtreeConnectionNode == 'Terminal'){
    
    parent[['IndArray']][cn, ] = subtree[['IndArray']][1, ]
    newIndividual = ChangeIndividual(parent[['IndArray']], parent[['DriveVec']], parent[['IndLength']], TRUE)
    
  }
  
  # Case 2: crossover node in parent is function and subtree root is terminal  
  if(parentConnectionNode == 'Function' & subtreeConnectionNode == 'Terminal'){
    
    newIndividual = TerminalToFunction(parent, cn, subtree)
    
  }
  
  # Case 3: crossover node in parent is terminal and subtree root is function
  if(parentConnectionNode == 'Terminal' & subtreeConnectionNode == 'Function'){
    
    newIndividual = FunctionToTerminal(parent, cn, subtree)
    
  }  
  
  # Case 4: crossover node in parent is function and subtree root is function
  if(parentConnectionNode == 'Function' & subtreeConnectionNode == 'Function'){
    
    terminalIndividual = CreateIndividual(0, 'grow', restricted = FALSE)
    newIndividual = TerminalToFunction(parent, cn, terminalIndividual)
    newIndividual = FunctionToTerminal(newIndividual, cn, subtree)
    
  } 
  
  newIndividual
  
}

NodeDepth = function(individual, nodePosition){
  
  if(nodePosition == 1) return(0)
  
  pointers = individual[['IndArray']][, 2:MaxNoColumns]
  
  oldEndPos = 1
  endPos = 1
  depth = 0
  
  while(endPos < nodePosition){
    
    endPos = max(as.numeric(pointers[1:endPos, ]), na.rm = T)
    depth = depth + 1
  }
  
  depth
  
}

ShrinkTree = function(individual, depthOfIndividual){
  
  while(depthOfIndividual > MaxDepthRun){
    
    #     shrinkNode = Sample(which(individual[['DriveVec']] == 'Function'), 1)
    shrinkNode = list(targetNode = Sample(which(individual[['DriveVec']] == 'Function'), 1), restriction = FALSE, variationAllowed = TRUE)
    
    if(shrinkNode[['targetNode']] %in% individual[['RestrictedRows']]){
      
      shrinkNode[['restriction']] = TRUE
      shrinkNode[['variationAllowed']] = FALSE
      
    }
    
    # Creation of terminal subtree with respect to restricted nodes
    if(shrinkNode[['variationAllowed']]){
      
      terminalIndividual = CreateIndividual(0, 'grow', restricted = FALSE)
      
    } else{
      
      terminalIndividual = CreateIndividual(0, 'grow', restricted = TRUE)
      #       shrinkNode[['variationAllowed']] = TRUE
      
    }
    #     terminalIndividual = CreateIndividual(0, 'grow', restricted = FALSE)
    
    individual = TerminalToFunction(individual, shrinkNode[['targetNode']], terminalIndividual)
    
    depthOfIndividual = NodeDepth(individual, individual[['IndLength']])
    
  }
  
  individual
  
}

DepthControl = function(individuals){
  
  changedIndividuals = sapply(individuals, function(x) x[['Changed']])
  
  changedIndividualsPositions = which(changedIndividuals)
  
  checkedIndividuals = individuals[changedIndividualsPositions]
  
  lengths = sapply(checkedIndividuals, function(x) x[['IndLength']])
  
  depths = mapply(function(ind, indLength) NodeDepth(ind, indLength), checkedIndividuals, lengths) 
  
  overMaxDepthPositions = which(depths > MaxDepthRun)
  
  checkedIndividuals = mapply(function(ind, depth) 
    ShrinkTree(ind, depth), 
    checkedIndividuals, depths, SIMPLIFY = FALSE)
  
  individuals[changedIndividualsPositions] = checkedIndividuals
  
  individuals
  
}

Crossover = function(parents1, parents2, terminalNodeProbability){
  # browser()
  noIndividuals = length(parents1)
  
  # Driving vectors of parents
  driveVectors1 = lapply(parents1, function(x) x[['DriveVec']])
  driveVectors2 = lapply(parents2, function(x) x[['DriveVec']])
  
  # Restricted rows vectors
  resRowVectors1 = lapply(parents1, function(x) x[['RestrictedRows']])
  resRowVectors2 = lapply(parents2, function(x) x[['RestrictedRows']])
  
  crossNodeIdentification1 = mapply(function(driveVec, resRowVec) 
    NodeSelection(driveVec, resRowVec, terminalNodeProbability), 
    driveVectors1, resRowVectors1, SIMPLIFY = FALSE)
  
  crossNodeIdentification2 = mapply(function(driveVec, resRowVec) 
    NodeSelection(driveVec, resRowVec, terminalNodeProbability), 
    driveVectors2, resRowVectors2, SIMPLIFY = FALSE) 
  
  crossNodes1 = data.frame(targetNode = sapply(crossNodeIdentification1, function(x) x[['targetNode']]),
                           restriction = sapply(crossNodeIdentification1, function(x) x[['restriction']]))  
  
  crossNodes2 = data.frame(targetNode = sapply(crossNodeIdentification2, function(x) x[['targetNode']]),
                           restriction = sapply(crossNodeIdentification2, function(x) x[['restriction']])) 
  
  #################################################
  if ("CT_new" %in% FunctionSet){
    CT_vector1 <- lapply(parents1, function(x) x[['IndArray']][1,1])
    CT_vector2 <- lapply(parents2, function(x) x[['IndArray']][1,1])
    
    CT_elements1 <- grep("CT_new", CT_vector1)
    CT_elements2 <- grep("CT_new", CT_vector2)
    
    for (i in 1:noIndividuals){
      if ((CT_vector1[[i]]=="CT_new") && (crossNodes1[i,1] %in% c(2:8)) && (CT_vector2[[i]]=="CT_new")){
        crossNodes2[i,1] <- crossNodes1[i,1]
        crossNodes2[i,2] <- TRUE
      } else if ((CT_vector1[[i]]=="CT_new") && (crossNodes1[i,1] %in% c(2:8)) && (CT_vector2[[i]] !="CT_new")){
        crossNodes2[i,1] <- crossNodes1[i,1]
        crossNodes2[i,2] <- TRUE
        parents2[[i]] <- parents2[[CT_elements2[round(runif(1,1,length(CT_elements2)))]]]
      } 
    }
    for (i in 1:noIndividuals){
      if ((CT_vector2[[i]]=="CT_new") && (crossNodes2[i,1] %in% c(2:8)) && (CT_vector1[[i]]=="CT_new")){
        crossNodes1[i,1] <- crossNodes2[i,1]
        crossNodes1[i,2] <- TRUE
      } else if ((CT_vector2[[i]]=="CT_new") && (crossNodes2[i,1] %in% c(2:8)) && (CT_vector1[[i]] !="CT_new")){
        crossNodes1[i,1] <- crossNodes2[i,1]
        crossNodes1[i,2] <- TRUE
        parents1[[i]] <- parents1[[CT_elements1[round(runif(1,1,length(CT_elements1)))]]]
      } 
    }
  } else {
    crossNodes1 <- crossNodes1
    crossNodes2 <- crossNodes2
  }
  
  #changing target node if it is 1 in marmmot
  for (u in 1:nrow(crossNodes1)){
    if(crossNodes1[u,1]==1 && parents1[[u]][['IndArray']][1,1]=="MARRMot"){
      crossNodes1[u,1]=round(runif(1,2,parents1[[u]][['IndLength']]))
      crossNodes1[u,2]=TRUE
    }
    if(crossNodes2[u,1]==1 && parents2[[u]][['IndArray']][1,1]=="MARRMot"){
      crossNodes2[u,1]=round(runif(1,2,parents2[[u]][['IndLength']]))
      crossNodes2[u,2]=TRUE
    }
  }
  
  
  variationAllowed = !xor(crossNodes1[['restriction']], crossNodes2[['restriction']])
  crossNodes1 = cbind(crossNodes1, variationAllowed)  
  crossNodes2 = cbind(crossNodes2, variationAllowed)
  
  crossNodes1 = split(crossNodes1, seq(nrow(crossNodes1))) 
  crossNodes2 = split(crossNodes2, seq(nrow(crossNodes2))) 
  
  # Separation of crossovering subtrees
  separedTrees1 = mapply(function(parent, cNode) SepareTrees(parent, cNode), parents1, crossNodes1, SIMPLIFY = FALSE)  
  separedTrees2 = mapply(function(parent, cNode) SepareTrees(parent, cNode), parents2, crossNodes2, SIMPLIFY = FALSE)  
  
  # Offsprings creation
  offsprings1 = mapply(function(parent, cNode, separedTree) AddSubtree(parent, cNode, separedTree), 
                       parents1, crossNodes1, separedTrees2, SIMPLIFY = FALSE)
  
  offsprings2 = mapply(function(parent, cNode, separedTree) AddSubtree(parent, cNode, separedTree), 
                       parents2, crossNodes2, separedTrees1, SIMPLIFY = FALSE)
  
  
  # Kontroly, smazat                     
  #   lapply(offsprings1, function(x) {
  #     a = x$IndArray
  #     if(a[1,1] %in% StrictFunctions & any(tail(a[-2,1]) %in% c('P','ET', 'CRI'))) browser() 
  #   }) 
  #   
  #   lapply(offsprings2, function(x) {
  #     a = x$IndArray
  #     if(a[1,1] %in% StrictFunctions & any(tail(a[-2,1]) %in% c('P','ET', 'CRI'))) browser() 
  #   })    
  #
  
  
  names(offsprings1) = names(offsprings2) = paste0('Off', 1:noIndividuals)
  
  # Depth control
  offsprings1 = DepthControl(offsprings1)
  offsprings2 = DepthControl(offsprings2)
  #offsprings1 <- mclapply(as.list(offsprings1), ComputeIndividual, mc.cores = n_jobs)
  #offsprings2 <- mclapply(as.list(offsprings2), ComputeIndividual, mc.cores = n_jobs)
  
  MakeEquation_new = function(individual){
    ReservoirSet <-c("WR","IR","RR","UR","FR","SR","CR")
    IA = individual[['IndArray']]
    if(individual[['IndLength']] == 1){
      return(paste0('y = ', IA[, 1]))
    }
    eq_vec = IA[, 1]
    for (i in individual[['IndLength']]:1){
      if( any(FunctionSet == IA[i, 1]) | any(ReservoirSet == IA[i, 1])){
        for(j in 2:ncol(IA)){
          if( !is.na(IA[i,j]) ){
            arg_row = as.numeric(IA[i,j])
            if( j == 2 && is.na(IA[i, (j + 1)]) ){
              if(!any(IA[i, 1] == SpecialFunctions)){
                eq_vec[i] = paste0(eq_vec[i], '(', eq_vec[arg_row])
              }else{
                eq_vec[i] = paste0(eq_vec[i], eq_vec[arg_row])
              }
            }
            if( j == 2 && !is.na(IA[i, (j + 1)]) ){
              if( any(PrefixFunctions == IA[i, 1]) ){
                first = paste0(eq_vec[i],'(', eq_vec[arg_row], ',')
              } else{
                first = paste0('(', eq_vec[arg_row])
              }
            }
            if(j > 2){
              if( any(PrefixFunctions == IA[i, 1]) ){
                eq_vec[i] = paste0(first, eq_vec[arg_row])
                if(j >= 3) first = paste0(eq_vec[i], ',')
              } else{
                eq_vec[i] = paste0(first, eq_vec[i], eq_vec[arg_row])
              }
            }
          }
        }
        eq_vec[i] = paste0(eq_vec[i], ')')
      }
    } 
    return(eq_vec)
  }
  
  #calculating offspring1
  # Ind_index <- 1:length(offsprings1)
  # Marrmot_Ind_index <-Ind_index[sapply(offsprings1, function(x) x[['IndArray']][1,1]=="MARRMot")]
  # Other_Ind_index <- Ind_index[sapply(offsprings1, function(x) x[['IndArray']][1,1]!="MARRMot")]
  # pre_index_order <-c(Marrmot_Ind_index,Other_Ind_index)
  
  Marrmot_Ind <-offsprings1[sapply(offsprings1, function(x) x[['IndArray']][1,1]=="MARRMot")]
  
  Other_Ind <-offsprings1[sapply(offsprings1, function(x) x[['IndArray']][1,1]!="MARRMot")]
  
  # if(length(Other_Ind)>0){
  #   for(vv in 1:length(Other_Ind)){
  #     Other_Ind[[vv]]<- ComputeIndividual(Other_Ind[[vv]])
  #   }
  # }
  
  if(length(Other_Ind)>0){
    Other_Ind <- mclapply(as.list(Other_Ind),ComputeIndividual, mc.cores = n_jobs)
  }
  
  Output_script_vec <- as.character(1:length(Marrmot_Ind))
  
  Matlab_code <- c()
  
  if(length(Marrmot_Ind)>0){
    
    for(gg in 1:length(Marrmot_Ind)){
      for(ml in 1:nrow(Marrmot_Ind[[gg]][['IndArray']])){
        if(Marrmot_Ind[[gg]][['IndArray']][ml,1] %in% c("P","E","T")){
          Marrmot_Ind[[gg]][['IndArray']][ml,1] <- runif(1)
        }
      }
      M_Eq <- MakeEquation_new(Marrmot_Ind[[gg]])
      Parameters <- c()
      for(hh in 1:33){
        Parameters[hh]<-max(min(eval(parse(text=noquote(M_Eq[[as.numeric(Marrmot_Ind[[gg]][['IndArray']][1,(hh+1)])]]))),1),0)
      }
      Output_script_name <- Output_script_vec[gg]
      Marrmot_Ind[[gg]][['MARRMot_Ind']]<- Output_script_name
      Model_to_use <- Code_name(Parameters[1],Parameters[2],Parameters[3],Parameters[4],Parameters[5],Parameters[6])
      Matlab_code <- c(Matlab_code,Code_generator(Model_to_use,Parameters[7:length(Parameters)],Output_script_name))
    }
    
    TO_MAT <- array(dim = c(length(Marrmot_Ind),1))
    TO_MAT[,1]<-Matlab_code
    write.csv(TO_MAT,paste(path,"/Matlab_fn/TO_MAT.csv",sep = ''))
    
    MAT_CODE <-c(paste("cd '",path,"/Matlab_fn","';",sep = ""),"A=readcell('TO_MAT.csv')","B=string();",paste("for i =2:",(length(Marrmot_Ind)+1),";B(i-1)=string(A(i,2));end",sep = ""),paste("parfor i=1:",length(Marrmot_Ind),";CommitEval(B(1,i));end",sep = ""))
    
    run_matlab_code(MAT_CODE)
    file_location <- paste(path,"/Matlab_fn",sep = '')
    
    Marrmot_Ind <- mclapply(as.list(Marrmot_Ind),Read_Marrmot, mc.cores = n_jobs)
    
    # for(x in 1:length(Marrmot_Ind)){
    #   Matlab_output <-read.csv(paste(file_location,"/Output_Q_",Output_script_vec[x],".csv",sep = ''))
    #   #unlink(paste(file_location,"/Output_Q_",Output_script_vec[x],".csv",sep = ''))
    #   Matlab_output <- t(Matlab_output)
    #   sim <- Matlab_output[,1][paste("Q_",1:nrow(DataSet),sep = "")]
    #   obs <- DataSet$Q
    #   
    #   for(mn in 1:length(FitnessFunction)){
    #     Marrmot_Ind[[x]]$Fitness[mn] <- FitnessComputation(sim,obs,FitnessFunction[mn])
    #   }
    #   Marrmot_Ind[[x]]$Sim <- sim
    #   Marrmot_Ind[[x]]$Equation <- MakeEquation(Marrmot_Ind[[x]])
    # }
    
  }
  
  offsprings1 <-c(Marrmot_Ind,Other_Ind)
  
  # Index_dataframe <-data.frame(pre_index_order,1:length(offsprings1))
  # names(Index_dataframe)<-c("pre","post")
  # Index_dataframe <-Index_dataframe[order(Index_dataframe$pre),]
  # offsprings1 <-offsprings1[Index_dataframe$post]
  
  #calculating offspring2
  # Ind_index <- 1:length(offsprings2)
  # Marrmot_Ind_index <-Ind_index[sapply(offsprings2, function(x) x[['IndArray']][1,1]=="MARRMot")]
  # Other_Ind_index <- Ind_index[sapply(offsprings2, function(x) x[['IndArray']][1,1]!="MARRMot")]
  # pre_index_order <-c(Marrmot_Ind_index,Other_Ind_index)
  
  Marrmot_Ind <-offsprings2[sapply(offsprings2, function(x) x[['IndArray']][1,1]=="MARRMot")]
  
  Other_Ind <-offsprings2[sapply(offsprings2, function(x) x[['IndArray']][1,1]!="MARRMot")]
  
  # if(length(Other_Ind)>0){
  #   for(vv in 1:length(Other_Ind)){
  #     Other_Ind[[vv]]<- ComputeIndividual(Other_Ind[[vv]])
  #   }
  # }
  
  if(length(Other_Ind)>0){
    Other_Ind <- mclapply(as.list(Other_Ind),ComputeIndividual, mc.cores = n_jobs)
  }
  
  Output_script_vec <- as.character(1:length(Marrmot_Ind))
  
  Matlab_code <- c()
  
  if(length(Marrmot_Ind)>0){
    
    for(gg in 1:length(Marrmot_Ind)){
      for(ml in 1:nrow(Marrmot_Ind[[gg]][['IndArray']])){
        if(Marrmot_Ind[[gg]][['IndArray']][ml,1] %in% c("P","E","T")){
          Marrmot_Ind[[gg]][['IndArray']][ml,1] <- runif(1)
        }
      }
      M_Eq <- MakeEquation_new(Marrmot_Ind[[gg]])
      Parameters <- c()
      for(hh in 1:33){
        Parameters[hh]<-max(min(eval(parse(text=noquote(M_Eq[[as.numeric(Marrmot_Ind[[gg]][['IndArray']][1,(hh+1)])]]))),1),0)
      }
      Output_script_name <- Output_script_vec[gg]
      Marrmot_Ind[[gg]][['MARRMot_Ind']]<- Output_script_name
      Model_to_use <- Code_name(Parameters[1],Parameters[2],Parameters[3],Parameters[4],Parameters[5],Parameters[6])
      Matlab_code <- c(Matlab_code,Code_generator(Model_to_use,Parameters[7:length(Parameters)],Output_script_name))
    }
    
    TO_MAT <- array(dim = c(length(Marrmot_Ind),1))
    TO_MAT[,1]<-Matlab_code
    write.csv(TO_MAT,paste(path,"/Matlab_fn/TO_MAT.csv",sep = ''))
    
    MAT_CODE <-c(paste("cd '",path,"/Matlab_fn","';",sep = ""),"A=readcell('TO_MAT.csv')","B=string();",paste("for i =2:",(length(Marrmot_Ind)+1),";B(i-1)=string(A(i,2));end",sep = ""),paste("parfor i=1:",length(Marrmot_Ind),";CommitEval(B(1,i));end",sep = ""))
    
    run_matlab_code(MAT_CODE)
    file_location <- paste(path,"/Matlab_fn",sep = '')
    
    Marrmot_Ind <- mclapply(as.list(Marrmot_Ind),Read_Marrmot, mc.cores = n_jobs)
    
    # for(x in 1:length(Marrmot_Ind)){
    #   Matlab_output <-read.csv(paste(file_location,"/Output_Q_",Output_script_vec[x],".csv",sep = ''))
    #   #unlink(paste(file_location,"/Output_Q_",Output_script_vec[x],".csv",sep = ''))
    #   Matlab_output <- t(Matlab_output)
    #   sim <- Matlab_output[,1][paste("Q_",1:nrow(DataSet),sep = "")]
    #   obs <- DataSet$Q
    #   
    #   for(mn in 1:length(FitnessFunction)){
    #     Marrmot_Ind[[x]]$Fitness[mn] <- FitnessComputation(sim,obs,FitnessFunction[mn])
    #   }
    #   Marrmot_Ind[[x]]$Sim <- sim
    #   Marrmot_Ind[[x]]$Equation <- MakeEquation(Marrmot_Ind[[x]])
    # }
    
  }
  
  offsprings2 <-c(Marrmot_Ind,Other_Ind)
  
  # Index_dataframe <-data.frame(pre_index_order,1:length(offsprings2))
  # names(Index_dataframe)<-c("pre","post")
  # Index_dataframe <-Index_dataframe[order(Index_dataframe$pre),]
  # offsprings2 <-offsprings2[Index_dataframe$post]
  
  
  
  # to avoid subscript out of bounds problem in some generations, let's try following code
  for (i in 1:length(offsprings1)){
    if (length(offsprings1[[i]]) <10){
      offsprings1[[i]] <- population[[round(runif(1,1,PopulationSize))]]
    } else if ((length(offsprings1[[i]]) == 10) & (length(offsprings1[[i]][['Sim']]) ==1)){
      offsprings1[[i]] <- population[[round(runif(1,1,PopulationSize))]]
    } else {
      offsprings1[[i]] <- offsprings1[[i]]
    }
  }
  
  for (i in 1:length(offsprings2)){
    if (length(offsprings2[[i]]) <10){
      offsprings2[[i]] <- population[[round(runif(1,1,PopulationSize))]]
    } else if ((length(offsprings2[[i]]) == 10) & (length(offsprings2[[i]][['Sim']]) ==1)){
      offsprings2[[i]] <- population[[round(runif(1,1,PopulationSize))]]
    } else {
      offsprings2[[i]] <- offsprings2[[i]]
    }
  }
  
  
  #offsprings1<-lapply(offsprings1,ComputeIndividual)
  #offsprings2<-lapply(offsprings2,ComputeIndividual)
  # cl <- makeCluster(n_jobs, 'FORK')
  # 
  # sub_population1 = split(offsprings1,c(1:n_jobs))
  # sub_population2 = split(offsprings2,c(1:n_jobs))
  # # Computation of fitness for both offspring groups (when the tree is changed)
  # ParComputeInd = function(sub_population) lapply(sub_population, ComputeIndividual)
  # sub_population1 = parLapply(cl, sub_population1, ParComputeInd)
  # sub_population2 = parLapply(cl, sub_population2, ParComputeInd)
  # stopCluster(cl)
  # offsprings1 = unlist(sub_population1, recursive = FALSE, use.names = FALSE)
  # offsprings2 = unlist(sub_population1, recursive = FALSE, use.names = FALSE)
  names(offsprings1) = names(offsprings2) = paste0('Off', 1:noIndividuals)
  #offsprings1 = mclapply(offsprings1, function(x) if(x[['Changed']]) ComputeIndividual(x) else x,mc.cores = getOption("mc.cores", 3L))
  #offsprings2 = mclapply(offsprings2, function(x) if(x[['Changed']]) ComputeIndividual(x) else x,mc.cores = getOption("mc.cores", 3L))
  #   browser()
  # isFitMulti = FitnessFunction %in% MultiObjectiveFitnesses
  # if(isFitMulti){
  #   
  #   parFit1 = MultiObjective(parents1)
  #   parFit2 = MultiObjective(parents2)
  #   
  #   parents1 = mapply(function(x, fit) {x[['Fitness']] = fit; x}, parents1, parFit1, SIMPLIFY = FALSE)
  #   parents2 = mapply(function(x, fit) {x[['Fitness']] = fit; x}, parents2, parFit2, SIMPLIFY = FALSE)    
  #   
  #   offFit1 = MultiObjective(offsprings1)
  #   offFit2 = MultiObjective(offsprings2)
  #   
  #   offsprings1 = mapply(function(x, fit) {x[['Fitness']] = fit; x}, offsprings1, offFit1, SIMPLIFY = FALSE)
  #   offsprings2 = mapply(function(x, fit) {x[['Fitness']] = fit; x}, offsprings2, offFit2, SIMPLIFY = FALSE)
  #   
  # }
  #
  nn <-round(runif(1,1,length(FitnessFunction)))
  betterParents = mapply(function(parent1, parent2) list(parent1, parent2)[[which.min(c(parent1[['Fitness']][nn], parent2[['Fitness']][nn]))]], 
                         parents1, parents2, SIMPLIFY = FALSE)
  
  # Determination of best individual i.e. best ind. from one row ofoffspring1, offspring2, betterParent (all together = family, without poorer parent)
  bestFamilyMembers = mapply(function(offspring1, offspring2, betterParent) {
    if(is.na(offspring1[['Fitness']][nn]) & is.na(offspring2[['Fitness']][nn])){
      return(betterParent);
    } 
    list(offspring1, offspring2)[[which.min(c(offspring1[['Fitness']][nn], offspring2[['Fitness']][nn]))]]
  },
  offsprings1, offsprings2, betterParents, SIMPLIFY = FALSE)
  
  # # The fitness must be in form of vector of all multi fitnesses, fit/rank is computed for whole population in main cycle                      
  # if(isFitMulti){
  #   
  #   bestFamilyMembers = lapply(bestFamilyMembers, function(x) ComputeIndividual(x))
  #   
  # }
  # 
  bestFamilyMembers
  
}

TreeMutation = function(individual){
  
  # Determination of mutation node 
  if (individual[['IndArray']][1,1]=="CT_new"){
    mutationNode = list(targetNode = Sample(9:individual[['IndLength']], 1), restriction = FALSE, variationAllowed = TRUE)
  } else {
    mutationNode = list(targetNode = Sample(1:individual[['IndLength']], 1), restriction = FALSE, variationAllowed = TRUE)
  }
  if(mutationNode[['targetNode']] %in% individual[['RestrictedRows']]){
    
    mutationNode[['restriction']] = TRUE
    mutationNode[['variationAllowed']] = FALSE
    
  }
  
  # Determination of mutation node depth
  mutationNodeDepth = NodeDepth(individual, mutationNode[['targetNode']])
  
  # Max. depth of new mutated subtree
  mutationTreeDepth = MaxDepthRun - mutationNodeDepth
  
  # Creation of mutation subtree
  if(mutationNode[['variationAllowed']]){
    
    mutationTree = CreateIndividual(mutationTreeDepth, 'grow', restricted = TRUE)
    #mutationTree = CreateIndividual(mutationTreeDepth, 'grow', restricted = FALSE)
    
  } else{
    
    mutationTree = CreateIndividual(mutationTreeDepth, 'grow', restricted = TRUE)
    mutationNode[['variationAllowed']] = TRUE
    
  }
  
  # Inserting of mutation subtree to mutating individual
  individual = AddSubtree(individual, mutationNode, mutationTree)
  
  # Mutant must be computed
  individual = ChangeIndividual(individual[['IndArray']], individual[['DriveVec']], individual[['IndLength']], TRUE)
  
  individual
  
}

TreeDependentMutatingIndividuals = function(mutatingIndividuals, noMutatingIndividuals, mutationType, variationProbs){
  
  if(mutationType == 'Separation'){
    
    # Individuals lengther than 1 are interesting for separation (else is just replication)
    possibleIndividuals = sapply(mutatingIndividuals, function(x) x[['IndLength']] > 1)
    P = variationProbs[['PmSeparation']]
    
  }
  
  if(mutationType == 'Constant'){
    
    possibleIndividuals = sapply(mutatingIndividuals, function(x) any(!is.na(suppressWarnings(as.numeric(x[['IndArray']][, 1])))))
    P = variationProbs[['PmConstant']]
    
  }
  
  possibleIndividuals[1] = FALSE # First position is elite individual (for sure, minimal probability, that position 1 of population will be in mutatingIndividuals)
  possibleMutants = sum(possibleIndividuals)
  noMutants = min(sum(rbinom(noMutatingIndividuals - 1, 1, P)), possibleMutants)
  MutantsPositions = Sample(which(possibleIndividuals), noMutants)
  
  MutantsPositions
  
}

SeparationMutation = function(individual){
  
  # Determination of mutation node
  if (individual[['IndArray']][1,1]=="CT_new"){
    mutationNode = c(targetNode = Sample(9:individual[['IndLength']], 1),
                     restriction = FALSE,
                     variationAllowed = TRUE)
  } else {
    mutationNode = c(targetNode = Sample(2:individual[['IndLength']], 1),
                     restriction = FALSE,
                     variationAllowed = TRUE)
  }
  
  individual = SepareTrees(individual, mutationNode)
  
  individual = ChangeIndividual(individual[['IndArray']], individual[['DriveVec']], individual[['IndLength']], TRUE)
  
  individual
  
}

NodeMutation = function(individual){
  
  usedArities = FunctionsDefinition[FunctionSet, ]
  
  if (individual[['IndArray']][1,1]=="CT_new") {
    noNodes = Sample(1:(individual[['IndLength']]-8), 1)
    mutationNodes = Sample(9:individual[['IndLength']], noNodes)
  } else {
    noNodes = Sample(1:individual[['IndLength']], 1)
    mutationNodes = Sample(1:individual[['IndLength']], noNodes)
  }
  
  functionsPositions = which(individual[['DriveVec']] == 'Function')
  mutatingFunctionsPositions = functionsPositions[functionsPositions %in% mutationNodes] 
  
  mutatingFunctions = individual[['IndArray']][mutatingFunctionsPositions, 1]  
  
  noNewFunctionNodes = length(mutatingFunctionsPositions)
  
  if(noNewFunctionNodes > 0){
    
    for(i in 1:noNewFunctionNodes){
      
      # Functions with same arities can be in node mutation process (possible problem with functions of same arity, but different strictness)
      possibleFunctions = usedArities[usedArities[ ,'Arity'] == usedArities[mutatingFunctions[i], 'Arity'], 'Functions']    
      individual[['IndArray']][mutatingFunctionsPositions[i], 1] = Sample(possibleFunctions, 1)
      
    }
    
  }
  
  terminalsPositions = which(individual[['DriveVec']] == 'Terminal')
  
  # mutating terminals positions
  mtp = terminalsPositions[terminalsPositions %in% mutationNodes] 
  
  noNewTerminalNodes = length(mtp)
  
  if(noNewTerminalNodes > 0){
    
    possibleConstants = runif(noNewTerminalNodes, ConstantRange[1], ConstantRange[2])
    
    possibleVariables = Sample(IndependentVariables, noNewTerminalNodes, replace = TRUE)
    
    variablePositions = mtp[which(individual[['IndArray']][mtp, 1]  %in% IndependentVariables)]
    constantPositions = mtp[which(!individual[['IndArray']][mtp, 1]  %in% IndependentVariables)]
    
    # number of new Terminal nodes
    nTn = length(variablePositions)
    
    # number of new Constant nodes
    nCn = length(constantPositions)
    
    if(nTn > 0){
      
      individual[['IndArray']][variablePositions, 1] = Sample(c(possibleConstants, possibleVariables), nTn)  
      
    }
    
    if(nCn > 0){
      
      individual[['IndArray']][constantPositions, 1] = Sample(possibleConstants, nCn)  
      
    }
    
  }
  
  individual = ChangeIndividual(individual[['IndArray']], individual[['DriveVec']], individual[['IndLength']], TRUE)
  
  individual
  
}

ConstantMutation = function(individual, constantMutationFactor){
  
  constantsPositions = which(individual[['DriveVec']] == 'Terminal' & !is.na(suppressWarnings(as.numeric(individual[['IndArray']][, 1]))))
  noNodes = Sample(1:length(constantsPositions), 1)
  
  # Selecton of nodes with constants. With repetitions = more mutations
  mutationNodesPositions = Sample(constantsPositions, noNodes)
  
  for(i in mutationNodesPositions){
    
    mutatingConstant = as.numeric(individual[['IndArray']][i, 1])
    actualRange = sort(c(mutatingConstant * (1 + constantMutationFactor), mutatingConstant * (1 - constantMutationFactor)))
    individual[['IndArray']][i, 1] = as.character(runif(1, actualRange[1], actualRange[2]))
    
  } 
  
  individual = ChangeIndividual(individual[['IndArray']], individual[['DriveVec']], individual[['IndLength']], TRUE)
  
  individual
  
}

Mutation = function(mutatingIndividuals, variationProbs, constantMutationFactor){
  
  noMutatingIndividuals = length(mutatingIndividuals)
  #   functionSet = computationVariables[['FunctionSet']]
  
  # Separation - individuals determination
  separationMutantsPositions = TreeDependentMutatingIndividuals(mutatingIndividuals, noMutatingIndividuals, 'Separation', variationProbs)  
  # Separation        
  mutatingIndividuals[separationMutantsPositions] = lapply(mutatingIndividuals[separationMutantsPositions], SeparationMutation)
  
  # Subtree mutation - individuals determination 
  noTreeMutants = sum(rbinom(noMutatingIndividuals, 1, variationProbs[['PmTree']]))
  treeMutantsPositions = Sample(1:noMutatingIndividuals, noTreeMutants)
  
  # Subtree mutation
  mutatingIndividuals[treeMutantsPositions] = lapply(mutatingIndividuals[treeMutantsPositions], TreeMutation)
  # Subtree mutation - control and reparation of tree depth 
  mutatingIndividuals = DepthControl(mutatingIndividuals)
  
  # Node mutation - individuals determination
  noNodeMutants = sum(rbinom(noMutatingIndividuals - 1, 1, variationProbs[['PmNode']]))
  nodeMutantsPositions = Sample(1:noMutatingIndividuals, noNodeMutants)
  
  # Node mutation
  mutatingIndividuals[nodeMutantsPositions] = lapply(mutatingIndividuals[nodeMutantsPositions], NodeMutation)
  
  # Constant mutation - individuals determination
  constantMutantsPositions = TreeDependentMutatingIndividuals(mutatingIndividuals, noMutatingIndividuals, 'Constant', variationProbs)   
  # Constant mutation
  mutatingIndividuals[constantMutantsPositions] = lapply(mutatingIndividuals[constantMutantsPositions], ConstantMutation, 
                                                         constantMutationFactor)
  
  mutantsPositions = unique(c(separationMutantsPositions, treeMutantsPositions, nodeMutantsPositions, constantMutantsPositions))
  mutantsPositions = paste0('I', mutantsPositions)
  # Computation of mutants
  delete.NULLs  <-  function(x.list){   # delele null/empty entries in a list
    x.list[unlist(lapply(x.list, length) != 0)]
  }
  names(mutatingIndividuals) = paste0('I', 1:noMutatingIndividuals)
  
  MakeEquation_new = function(individual){
    ReservoirSet <-c("WR","IR","RR","UR","FR","SR","CR")
    IA = individual[['IndArray']]
    if(individual[['IndLength']] == 1){
      return(paste0('y = ', IA[, 1]))
    }
    eq_vec = IA[, 1]
    for (i in individual[['IndLength']]:1){
      if( any(FunctionSet == IA[i, 1]) | any(ReservoirSet == IA[i, 1])){
        for(j in 2:ncol(IA)){
          if( !is.na(IA[i,j]) ){
            arg_row = as.numeric(IA[i,j])
            if( j == 2 && is.na(IA[i, (j + 1)]) ){
              if(!any(IA[i, 1] == SpecialFunctions)){
                eq_vec[i] = paste0(eq_vec[i], '(', eq_vec[arg_row])
              }else{
                eq_vec[i] = paste0(eq_vec[i], eq_vec[arg_row])
              }
            }
            if( j == 2 && !is.na(IA[i, (j + 1)]) ){
              if( any(PrefixFunctions == IA[i, 1]) ){
                first = paste0(eq_vec[i],'(', eq_vec[arg_row], ',')
              } else{
                first = paste0('(', eq_vec[arg_row])
              }
            }
            if(j > 2){
              if( any(PrefixFunctions == IA[i, 1]) ){
                eq_vec[i] = paste0(first, eq_vec[arg_row])
                if(j >= 3) first = paste0(eq_vec[i], ',')
              } else{
                eq_vec[i] = paste0(first, eq_vec[i], eq_vec[arg_row])
              }
            }
          }
        }
        eq_vec[i] = paste0(eq_vec[i], ')')
      }
    } 
    return(eq_vec)
  }
  
  #calculating mutatingIndividuals
  # Ind_index <- 1:length(mutatingIndividuals)
  # Marrmot_Ind_index <-Ind_index[sapply(mutatingIndividuals, function(x) x[['IndArray']][1,1]=="MARRMot")]
  # Other_Ind_index <- Ind_index[sapply(mutatingIndividuals, function(x) x[['IndArray']][1,1]!="MARRMot")]
  # pre_index_order <-c(Marrmot_Ind_index,Other_Ind_index)
  
  Marrmot_Ind <-mutatingIndividuals[sapply(mutatingIndividuals, function(x) x[['IndArray']][1,1]=="MARRMot")]
  
  Other_Ind <-mutatingIndividuals[sapply(mutatingIndividuals, function(x) x[['IndArray']][1,1]!="MARRMot")]
  
  # if(length(Other_Ind)>0){
  #   for(vv in 1:length(Other_Ind)){
  #     Other_Ind[[vv]]<- ComputeIndividual(Other_Ind[[vv]])
  #   }
  # }
  
  if(length(Other_Ind)>0){
    Other_Ind <- mclapply(as.list(Other_Ind),ComputeIndividual, mc.cores = n_jobs)
  }
  
  Output_script_vec <- as.character(1:length(Marrmot_Ind))
  
  Matlab_code <- c()
  
  if(length(Marrmot_Ind)>0){
    
    for(gg in 1:length(Marrmot_Ind)){
      for(ml in 1:nrow(Marrmot_Ind[[gg]][['IndArray']])){
        if(Marrmot_Ind[[gg]][['IndArray']][ml,1] %in% c("P","E","T")){
          Marrmot_Ind[[gg]][['IndArray']][ml,1] <- runif(1)
        }
      }
      M_Eq <- MakeEquation_new(Marrmot_Ind[[gg]])
      Parameters <- c()
      for(hh in 1:33){
        Parameters[hh]<-max(min(eval(parse(text=noquote(M_Eq[[as.numeric(Marrmot_Ind[[gg]][['IndArray']][1,(hh+1)])]]))),1),0)
      }
      Output_script_name <- Output_script_vec[gg]
      Marrmot_Ind[[gg]][['MARRMot_Ind']]<- Output_script_name
      Model_to_use <- Code_name(Parameters[1],Parameters[2],Parameters[3],Parameters[4],Parameters[5],Parameters[6])
      Matlab_code <- c(Matlab_code,Code_generator(Model_to_use,Parameters[7:length(Parameters)],Output_script_name))
    }
    
    TO_MAT <- array(dim = c(length(Marrmot_Ind),1))
    TO_MAT[,1]<-Matlab_code
    write.csv(TO_MAT,paste(path,"/Matlab_fn/TO_MAT.csv",sep = ''))
    
    MAT_CODE <-c(paste("cd '",path,"/Matlab_fn","';",sep = ""),"A=readcell('TO_MAT.csv')","B=string();",paste("for i =2:",(length(Marrmot_Ind)+1),";B(i-1)=string(A(i,2));end",sep = ""),paste("parfor i=1:",length(Marrmot_Ind),";CommitEval(B(1,i));end",sep = ""))
    
    run_matlab_code(MAT_CODE)
    file_location <- paste(path,"/Matlab_fn",sep = '')
    
    Marrmot_Ind <- mclapply(as.list(Marrmot_Ind),Read_Marrmot, mc.cores = n_jobs)
    
    # for(x in 1:length(Marrmot_Ind)){
    #   Matlab_output <-read.csv(paste(file_location,"/Output_Q_",Output_script_vec[x],".csv",sep = ''))
    #   #unlink(paste(file_location,"/Output_Q_",Output_script_vec[x],".csv",sep = ''))
    #   Matlab_output <- t(Matlab_output)
    #   sim <- Matlab_output[,1][paste("Q_",1:nrow(DataSet),sep = "")]
    #   obs <- DataSet$Q
    #   
    #   for(mn in 1:length(FitnessFunction)){
    #     Marrmot_Ind[[x]]$Fitness[mn] <- FitnessComputation(sim,obs,FitnessFunction[mn])
    #   }
    #   Marrmot_Ind[[x]]$Sim <- sim
    #   Marrmot_Ind[[x]]$Equation <- MakeEquation(Marrmot_Ind[[x]])
    # }
    
  }
  
  mutatingIndividuals <-c(Marrmot_Ind,Other_Ind)
  
  # Index_dataframe <-data.frame(pre_index_order,1:length(mutatingIndividuals))
  # names(Index_dataframe)<-c("pre","post")
  # Index_dataframe <-Index_dataframe[order(Index_dataframe$pre),]
  # mutatingIndividuals <-mutatingIndividuals[Index_dataframe$post]
  
  #mutatingIndividuals <- mclapply(as.list(mutatingIndividuals),ComputeIndividual,mc.cores = n_jobs)
  #mutatingIndividuals<- lapply(mutatingIndividuals,ComputeIndividual)
  # sub_mutationInd = split(mutatingIndividuals, c(1:n_jobs))
  # ParComputeMutInd = function(n) {
  #   sub_mutationInd[[n]] = delete.NULLs(sub_mutationInd[[n]][mutantsPositions])
  #   sub_mutationInd[[n]] = lapply(sub_mutationInd[[n]], ComputeIndividual)
  #   sub_mutationInd[[n]]
  # }
  # cl <- makeCluster(n_jobs,"FORK")
  # sub_mutationInd = parLapply(cl, 1:n_jobs, ParComputeMutInd)
  # stopCluster(cl)
  # 
  # #mutatingIndividuals[mutantsPositions] = mclapply(mutatingIndividuals[mutantsPositions], ComputeIndividual,mc.cores = getOption("mc.cores", 3L))
  # mutatingIndividuals[mutantsPositions] = unlist(sub_mutationInd, recursive = FALSE, use.names = FALSE)
  mutatingIndividuals   
  
}  

controlIndividual = function(x, iv, co) {
  
  a = x[['RestrictedRows']]
  if(!is.null(a) && !is.na(a[1]) && min(a > 3))  {print(x); stop(co)}
  #   if(a[1,1] %in% StrictFunctions & any(tail(a[-1,1]) %in% iv)) {print(x); stop(co)}
  
}

controlIndividualArray = function(x, iv, co) {
  
  a = x[['IndArray']]
  #   if(!is.null(a) && !is.na(a[1]) && min(a > 2))  {print(x); stop(co)}
  if(a[1,1] %in% StrictFunctions & any(tail(a[-2,1]) %in% iv)) {print(x); stop(co)}
  
}

# source('CHMGP_SuperflexFunctions.R')
ScaleParams = function(parPrcg){
  
  parPrcg[parPrcg < 0] = 0
  parPrcg[parPrcg > 1] = 1    
  mapply(function(prcg, rng){prcg * (rng[['max']] - rng[['min']]) + rng[['min']]}, parPrcg, SFPR[names(parPrcg)])
  
}

ScaleParamsManualy = function(parPrcg, functionNam){
  
  names(parPrcg) = names(formals(functionNam))
  parPrcg[parPrcg < 0] = 0
  parPrcg[parPrcg > 1] = 1    
  mapply(function(prcg, rng){prcg * (rng[['max']] - rng[['min']]) + rng[['min']]}, parPrcg, SFPR[names(parPrcg)])
  
}


## MI: Fast Reservoir 
MI = function(alpha_Qq_FR, K_Qq_FR, Ce, m_E_FR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(alpha_Qq_FR, K_Qq_FR, Ce, m_E_FR)
  names(modPars) = names(formals(MI))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MI", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR  = as.double(K_Qq_FR), Ce  = as.double(Ce), m_E_FR= as.double(m_E_FR), dT = as.double(1),
                     output = as.double(output))
  
  fOutput[['output']]
  
}


## MII: Fast Reservoir 
MII = function(Ce, Beta_Qq_UR, Smax_UR, K_Qb_UR, Beta_E_UR, SiniFr_UR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(Ce, Beta_Qq_UR, Smax_UR, K_Qb_UR, Beta_E_UR, SiniFr_UR)
  names(modPars) = names(formals(MII))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     Ce  = as.double(Ce), Beta_Qq_UR= as.double(Beta_Qq_UR), Smax_UR = as.double(Smax_UR), K_Qb_UR = as.double(K_Qb_UR), 
                     Beta_E_UR=as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), dT = as.double(1),
                     output = as.double(output))
  
  fOutput[['output']]
  
}

## M3: Unsaturated and Fast Reservoir with four parameters
MIII = function(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR)
  names(modPars) = names(formals(MIII))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MIII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR = as.double(Smax_UR), Beta_Qq_UR = as.double(Beta_Qq_UR), 
                     Beta_E_UR = as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), dT=as.double(1),
                     output = as.double(output))                  
  
  fOutput[['output']]
  
}


##M4: Unsaturated and Fast Reservoir with five parameters
MIV = function(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR)
  names(modPars) = names(formals(MIV))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MIV", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR = as.double(Smax_UR), Beta_Qq_UR = as.double(Beta_Qq_UR), 
                     Beta_E_UR = as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), dT=as.double(1),
                     output = as.double(output))
  
  fOutput[['output']]
  
}

##M5: Unsaturated and Fast Reservoir with six parameters
MV = function(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, Tlag)
  names(modPars) = names(formals(MV))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MV", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR = as.double(Smax_UR), Beta_Qq_UR = as.double(Beta_Qq_UR), 
                     Beta_E_UR = as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), dT=as.double(1), Tlag=as.double(Tlag),
                     output = as.double(output))
  
  fOutput[['output']]
  
}


##M6: Interception, Unsaturated and Fast Reservoir with eight parameters
MVI = function( alpha_Qq_FR, Beta_Qq_UR, K_Qq_FR, Ce, Smax_UR, Smax_IR, m_QE_IR, Beta_E_UR, SiniFr_UR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(alpha_Qq_FR, Beta_Qq_UR, K_Qq_FR, Ce, Smax_UR, Smax_IR, m_QE_IR, Beta_E_UR, SiniFr_UR, Tlag)
  names(modPars) = names(formals(MVI))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars) # head is due to Tlag
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MVI", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR=as.double(alpha_Qq_FR), Beta_Qq_UR=as.double(Beta_Qq_UR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR=as.double(Smax_UR), 
                     Smax_IR=as.double(Smax_IR), m_QE_IR = as.double(m_QE_IR), Beta_E_UR= as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), dT=as.double(1), Tlag=as.double(Tlag),
                     output = as.double(output))
  
  fOutput[['output']]
  
}


##M7: Unsaturated, Riparian and Fast Reservoir with seven parameters
MVII = function(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, D_R, K_Qq_RR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, D_R, K_Qq_RR, Tlag)
  names(modPars) = names(formals(MVII))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MVII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR = as.double(Smax_UR), Beta_Qq_UR = as.double(Beta_Qq_UR), 
                     Beta_E_UR = as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), D_R=as.double(D_R), K_Qq_RR=as.double(K_Qq_RR), dT=as.double(1), Tlag=as.double(Tlag),
                     output = as.double(output)) 
  #   ctrlVar1 = ctrlVar2 = numeric(n)                     
  #   fOutput = .Fortran("MVII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
  #                      alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR = as.double(Smax_UR), Beta_Qq_UR = as.double(Beta_Qq_UR), 
  #                      Beta_E_UR = as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), D_R=as.double(D_R), K_Qq_RR=as.double(K_Qq_RR), dT=as.double(1), Tlag=as.double(Tlag),
  #                      output = as.double(output), ctrlVar1 = as.double(ctrlVar1), ctrlVar2 = as.double(ctrlVar2))
  
  #   ControlVariable1 <<- fOutput[['ctrlVar1']]        
  #   ControlVariable2 <<- fOutput[['ctrlVar2']]
  
  fOutput[['output']]
  
}


##M8: Fast and slow Reservoirs
MVIII = function( K_Qq_FR, Ce, K_Qq_SR, D_S, m_E_FR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(K_Qq_FR, Ce, K_Qq_SR, D_S, m_E_FR)
  names(modPars) = names(formals(MVIII))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MVIII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), K_Qq_SR= as.double(K_Qq_SR), D_S=as.double(D_S), m_E_FR=as.double(m_E_FR),dT=as.double(1),
                     output = as.double(output))
  
  fOutput[['output']]
  
}


##M9: Unsaturated, Fast and slow Reservoirs with six parameters
MIX = function(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR)
  names(modPars) = names(formals(MIX))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MIX", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     Beta_Qq_UR=as.double(Beta_Qq_UR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), K_Qq_SR= as.double(K_Qq_SR), 
                     D_S=as.double(D_S), Smax_UR=as.double(Smax_UR), Beta_E_UR= as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), 
                     dT=as.double(1), output = as.double(output))
  
  fOutput[['output']]
  
}


##M10: Unsaturated, Fast and slow Reservoirs with five parameters
MX = function(K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, Tlag)
  names(modPars) = names(formals(MX))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MX", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), K_Qq_SR= as.double(K_Qq_SR), 
                     D_S=as.double(D_S), Smax_UR=as.double(Smax_UR), Beta_E_UR= as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), 
                     dT= as.double(1), Tlag=as.double(Tlag), output = as.double(output))
  
  fOutput[['output']]
  
}

##M11: Unsaturated, Fast and slow Reservoirs with six parameters
MXI = function(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, Tlag)
  names(modPars) = names(formals(MXI))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MXI", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     Beta_Qq_UR=as.double(Beta_Qq_UR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), K_Qq_SR= as.double(K_Qq_SR), 
                     D_S=as.double(D_S), Smax_UR=as.double(Smax_UR), Beta_E_UR= as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), 
                     dT= as.double(1), Tlag=as.double(Tlag), output = as.double(output))
  
  fOutput[['output']]
  
}


##M12: Interception, Unsaturated, Fast and slow Reservoirs with eight parameters
MXII = function(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, SiniFr_UR, Smax_UR, Smax_IR, m_QE_IR, Beta_E_UR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, SiniFr_UR, Smax_UR, Smax_IR, m_QE_IR, Beta_E_UR, Tlag)
  names(modPars) = names(formals(MXII))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MXII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     Beta_Qq_UR=as.double(Beta_Qq_UR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce),
                     K_Qq_SR= as.double(K_Qq_SR), D_S=as.double(D_S), SiniFr_UR=as.double(SiniFr_UR), Smax_UR=as.double(Smax_UR), 
                     Smax_IR=as.double(Smax_IR), m_QE_IR = as.double(m_QE_IR), Beta_E_UR= as.double(Beta_E_UR), 
                     dT=as.double(1), Tlag=as.double(Tlag), output = as.double(output))
  
  fOutput[['output']]
  
}



### combined_tank - CT
CT = function(Ce, Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR, K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, Beta_Qq_UR, Smax_UR, Beta_E_UR,
              SiniFr_UR, K_Qb_UR, mu_Qq_UR, Smax_IR, m_QE_IR, K_Qq_RR, Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, K_Qq_SR, m_E_SR, alpha_Qq_SR, P_ED_max, m_P_ED, 
              D_S, D_I, D_F, D_R, D_C, Tlag, Lag_RR, Lag_FR, Lag_SR, option_i, option_u, option_f, option_s, option_c, w_res, i_res, r_res, u_res, f_res,s_res,c_res)
{
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  x3 = DataSet[['T']]
  
  modPars = c(Ce, Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR, K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, Beta_Qq_UR, Smax_UR, Beta_E_UR,
              SiniFr_UR, K_Qb_UR, mu_Qq_UR, Smax_IR, m_QE_IR, K_Qq_RR, Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, K_Qq_SR, m_E_SR, alpha_Qq_SR, P_ED_max, m_P_ED, 
              D_S, D_I, D_F, D_R, D_C, Tlag, Lag_RR, Lag_FR, Lag_SR, option_i, option_u, option_f, option_s, option_c, w_res, i_res, r_res, u_res, f_res,s_res,c_res)
  
  names(modPars) = names(formals(CT))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  # Scaling parameters
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  # Scaling option and res parameters
  optresNam = tail(names(modPars), 15)
  sapply(optresNam, function(nm) assign(nm, round(get(nm)), envir = parent.env(environment())))
  
  ress = c(w = w_res, i = i_res, r = r_res, u = u_res, f = f_res, s = s_res,c = c_res)
  
  # Testing if function is from meaningfull combination
  #   print(0)
  testMeaningfulness = sapply(MeaningfullCTcombinations, function(x) identical(unlist(x),ress))
  if(!any(testMeaningfulness)) return(NA)
  
  # Turning on and off the arguments in dependency of res values
  nam2zero = unlist(DepArgsOfCT[!ress])
  sapply(nam2zero, function(nm) assign(nm, 0, envir = parent.env(environment())))
  #   print(1)
  
  n = length(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("combined_tank", P = as.double(x1), E = as.double(x2), T = as.double(x3), n = as.integer(n),
                     Ce = as.double(Ce),Cp_WR = as.double(Cp_WR), m_Q_WR = as.double(m_Q_WR), Kq_WR = as.double(Kq_WR), Tp_WR = as.double(Tp_WR), Tm_WR = as.double(Tm_WR), K_Qq_FR= as.double(K_Qq_FR), K_Qb_FR= as.double(K_Qb_FR), m_E_FR=as.double(m_E_FR), alpha_Qq_FR=as.double(alpha_Qq_FR), Beta_Qq_UR=as.double(Beta_Qq_UR),
                     Smax_UR=as.double(Smax_UR), Beta_E_UR= as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), K_Qb_UR= as.double(K_Qb_UR), mu_Qq_UR=as.double(mu_Qq_UR), Smax_IR=as.double(Smax_IR),
                     m_QE_IR = as.double(m_QE_IR), K_Qq_RR=as.double(K_Qq_RR),Smax_CR=as.double(Smax_CR), Umax_uCR=as.double(Umax_uCR), Smin_uCR=as.double(Smin_uCR), Beta_Qq_uCR=as.double(Beta_Qq_uCR), mu_Qq_uCR=as.double(mu_Qq_uCR), K_Qb_uCR=as.double(K_Qb_uCR), Sevmax_CR=as.double(Sevmax_CR), Beta_E_CR=as.double(Beta_E_CR), Beta_Qq_sCR=as.double(Beta_Qq_sCR), K_Qb_sCR=as.double(K_Qb_sCR), K_Qd_sCR=as.double(K_Qd_sCR), K_Qq_SR= as.double(K_Qq_SR), m_E_SR= as.double(m_E_SR), alpha_Qq_SR= as.double(alpha_Qq_SR),
                     P_ED_max= as.double(P_ED_max), m_P_ED= as.double(m_P_ED), 
                     D_S=as.double(D_S), D_I=as.double(D_I), D_F=as.double(D_F), D_R=as.double(D_R),D_C=as.double(D_C),
                     Tlag=as.double(Tlag), dT=as.double(1),Lag_RR=as.integer(Lag_RR),Lag_FR=as.integer(Lag_FR),Lag_SR=as.integer(Lag_SR),
                     option_i=as.integer(option_i),option_u=as.integer(option_u),option_f=as.integer(option_f), option_s=as.integer(option_s),option_c=as.integer(option_c), 
                     w_res=as.integer(w_res), i_res=as.integer(i_res), r_res=as.integer(r_res), u_res=as.integer(u_res), f_res=as.integer(f_res),  s_res=as.integer(s_res),c_res=as.integer(c_res),
                     output = as.double(output))
  
  fOutput[['output']]
  
}

### FUSE

FUSE = function(rferr_add,rferr_mlt, maxwatr_1 ,maxwatr_2, fracten, frchzne, fprimqb, rtfrac1, percrte, percexp,sacpmlt, sacpexp, percfrac, iflwrte, 
                baserte, qb_powr, qb_prms, qbrate_2a, qbrate_2b, sareamax, axv_bexp, loglamb, tishape, timedelay,rferr, arch1, arch2, qsurf, qperc, esoil, qintf, q_tdh)
{
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(rferr_add,rferr_mlt, maxwatr_1 ,maxwatr_2, fracten, frchzne, fprimqb, rtfrac1, percrte, percexp,sacpmlt, sacpexp, percfrac, iflwrte, 
              baserte, qb_powr, qb_prms, qbrate_2a, qbrate_2b, sareamax, axv_bexp, loglamb, tishape, timedelay,rferr, arch1, arch2, qsurf, qperc, esoil, qintf, q_tdh)
  
  names(modPars) = names(formals(FUSE))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  optresNam = tail(names(modPars), 8)
  sapply(optresNam, function(nm) assign(nm, round(get(nm)), envir = parent.env(environment())))
  
  n = length(x1)
  
  output = numeric(n)
  
  library(fuse)
  data("modlist")
  DATA <- data.frame(P=x1, E=x2)
  mid <-which(modlist[['rferr']]==rferr & modlist[['arch1']]==arch1 & modlist[['arch2']]==arch2 & modlist[['qsurf']]==qsurf & modlist[['qperc']]==qperc & modlist[['esoil']]==esoil & modlist[['qintf']]==qintf & modlist[['q_tdh']]==q_tdh)
  parameters <- data.frame(rferr_add= rferr_add,rferr_mlt=rferr_mlt, maxwatr_1=maxwatr_1 ,maxwatr_2=maxwatr_2, fracten=fracten, frchzne=frchzne, fprimqb=fprimqb, rtfrac1=rtfrac1, percrte=percrte, percexp=percexp,sacpmlt=sacpmlt, sacpexp=sacpexp, percfrac=percfrac, iflwrte=iflwrte,baserte=baserte, qb_powr=qb_powr, qb_prms=qb_prms, qbrate_2a=qbrate_2a, qbrate_2b=qbrate_2b, sareamax=sareamax, axv_bexp=axv_bexp, loglamb=loglamb, tishape=tishape, timedelay=timedelay)
  output <- fuse(DATA,mid,1,parameters)
  return(output)
  
}

####### DCT #####
DCT <- function(Ce1,Cp_WR1, m_Q_WR1, Kq_WR1, Tp_WR1, Tm_WR1, K_Qq_FR1, K_Qb_FR1,  m_E_FR1, alpha_Qq_FR1, Beta_Qq_UR1, Smax_UR1, Beta_E_UR1,
                SiniFr_UR1, K_Qb_UR1, mu_Qq_UR1, Smax_IR1, m_QE_IR1, K_Qq_RR1, Smax_CR1, Umax_uCR1, Smin_uCR1, Beta_Qq_uCR1, mu_Qq_uCR1, K_Qb_uCR1, Sevmax_CR1, Beta_E_CR1, Beta_Qq_sCR1, K_Qb_sCR1, K_Qd_sCR1, K_Qq_SR1, m_E_SR1, alpha_Qq_SR1, P_ED_max1, m_P_ED1,
                D_S1, D_I1, D_F1, D_R1, D_C1, Tlag1, Lag_RR1, Lag_FR1, Lag_SR1, option_i1, option_u1, option_f1, option_s1, option_c1,w_res1, i_res1, r_res1, u_res1, f_res1,s_res1,c_res1,Ce2,Cp_WR2, m_Q_WR2, Kq_WR2, Tp_WR2, Tm_WR2, K_Qq_FR2, K_Qb_FR2,  m_E_FR2, alpha_Qq_FR2, Beta_Qq_UR2, Smax_UR2, Beta_E_UR2,
                SiniFr_UR2, K_Qb_UR2, mu_Qq_UR2, Smax_IR2, m_QE_IR2, K_Qq_RR2, Smax_CR2, Umax_uCR2, Smin_uCR2, Beta_Qq_uCR2, mu_Qq_uCR2, K_Qb_uCR2, Sevmax_CR2, Beta_E_CR2, Beta_Qq_sCR2, K_Qb_sCR2, K_Qd_sCR2, K_Qq_SR2, m_E_SR2, alpha_Qq_SR2, P_ED_max2, m_P_ED2,
                D_S2, D_I2, D_F2, D_R2, D_C2, Tlag2, Lag_RR2, Lag_FR2, Lag_SR2, option_i2, option_u2, option_f2, option_s2, option_c2,w_res2, i_res2, r_res2, u_res2, f_res2,s_res2,c_res2,Ce3,Cp_WR3, m_Q_WR3, Kq_WR3, Tp_WR3, Tm_WR3, K_Qq_FR3, K_Qb_FR3,  m_E_FR3, alpha_Qq_FR3, Beta_Qq_UR3, Smax_UR3, Beta_E_UR3,
                SiniFr_UR3, K_Qb_UR3, mu_Qq_UR3, Smax_IR3, m_QE_IR3, K_Qq_RR3, Smax_CR3, Umax_uCR3, Smin_uCR3, Beta_Qq_uCR3, mu_Qq_uCR3, K_Qb_uCR3, Sevmax_CR3, Beta_E_CR3, Beta_Qq_sCR3, K_Qb_sCR3, K_Qd_sCR3, K_Qq_SR3, m_E_SR3, alpha_Qq_SR3, P_ED_max, m_P_ED,
                D_S3, D_I3, D_F3, D_R3, D_C3, Tlag3, Lag_RR3, Lag_FR3, Lag_SR3, option_i3, option_u3, option_f3, option_s3, option_c3,w_res3, i_res3, r_res3, u_res3, f_res3,s_res3,c_res3,R_t1,R_t2,R_t3,R_1,R_2,R_3){
  
  modPars = c(Ce1,Cp_WR1, m_Q_WR1, Kq_WR1, Tp_WR1, Tm_WR1, K_Qq_FR1, K_Qb_FR1,  m_E_FR1, alpha_Qq_FR1, Beta_Qq_UR1, Smax_UR1, Beta_E_UR1,
              SiniFr_UR1, K_Qb_UR1, mu_Qq_UR1, Smax_IR1, m_QE_IR1, K_Qq_RR1, Smax_CR1, Umax_uCR1, Smin_uCR1, Beta_Qq_uCR1, mu_Qq_uCR1, K_Qb_uCR1, Sevmax_CR1, Beta_E_CR1, Beta_Qq_sCR1, K_Qb_sCR1, K_Qd_sCR1, K_Qq_SR1, m_E_SR1, alpha_Qq_SR1, P_ED_max1, m_P_ED1,
              D_S1, D_I1, D_F1, D_R1, D_C1, Tlag1, Lag_RR1, Lag_FR1, Lag_SR1, option_i1, option_u1, option_f1, option_s1, option_c1,w_res1, i_res1, r_res1, u_res1, f_res1,s_res1,c_res1,Ce2,Cp_WR2, m_Q_WR2, Kq_WR2, Tp_WR2, Tm_WR2, K_Qq_FR2, K_Qb_FR2,  m_E_FR2, alpha_Qq_FR2, Beta_Qq_UR2, Smax_UR2, Beta_E_UR2,
              SiniFr_UR2, K_Qb_UR2, mu_Qq_UR2, Smax_IR2, m_QE_IR2, K_Qq_RR2, Smax_CR2, Umax_uCR2, Smin_uCR2, Beta_Qq_uCR2, mu_Qq_uCR2, K_Qb_uCR2, Sevmax_CR2, Beta_E_CR2, Beta_Qq_sCR2, K_Qb_sCR2, K_Qd_sCR2, K_Qq_SR2, m_E_SR2, alpha_Qq_SR2, P_ED_max2, m_P_ED2,
              D_S2, D_I2, D_F2, D_R2, D_C2, Tlag2, Lag_RR2, Lag_FR2, Lag_SR2, option_i2, option_u2, option_f2, option_s2, option_c2,w_res2, i_res2, r_res2, u_res2, f_res2,s_res2,c_res2,Ce3,Cp_WR3, m_Q_WR3, Kq_WR3, Tp_WR3, Tm_WR3, K_Qq_FR3, K_Qb_FR3,  m_E_FR3, alpha_Qq_FR3, Beta_Qq_UR3, Smax_UR3, Beta_E_UR3,
              SiniFr_UR3, K_Qb_UR3, mu_Qq_UR3, Smax_IR3, m_QE_IR3, K_Qq_RR3, Smax_CR3, Umax_uCR3, Smin_uCR3, Beta_Qq_uCR3, mu_Qq_uCR3, K_Qb_uCR3, Sevmax_CR3, Beta_E_CR3, Beta_Qq_sCR3, K_Qb_sCR3, K_Qd_sCR3, K_Qq_SR3, m_E_SR3, alpha_Qq_SR3, P_ED_max, m_P_ED,
              D_S3, D_I3, D_F3, D_R3, D_C3, Tlag3, Lag_RR3, Lag_FR3, Lag_SR3, option_i3, option_u3, option_f3, option_s3, option_c3,w_res3, i_res3, r_res3, u_res3, f_res3,s_res3,c_res3,R_t1,R_t2,R_t3,R_1,R_2,R_3)
  
  names(modPars) = names(formals(DCT))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  # Scaling parameters
  scaledPars = ScaleParams(modPars[169:174])
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  # Scaling option and res parameters
  optresNam = tail(names(modPars), 3)
  sapply(optresNam, function(nm) assign(nm, round(get(nm)), envir = parent.env(environment())))
  
  flow_1 <- CT1(Ce1,Cp_WR1, m_Q_WR1, Kq_WR1, Tp_WR1, Tm_WR1, K_Qq_FR1, K_Qb_FR1,  m_E_FR1, alpha_Qq_FR1, Beta_Qq_UR1, Smax_UR1, Beta_E_UR1,
                SiniFr_UR1, K_Qb_UR1, mu_Qq_UR1, Smax_IR1, m_QE_IR1, K_Qq_RR1, Smax_CR1, Umax_uCR1, Smin_uCR1, Beta_Qq_uCR1, mu_Qq_uCR1, K_Qb_uCR1, Sevmax_CR1, Beta_E_CR1, Beta_Qq_sCR1, K_Qb_sCR1, K_Qd_sCR1, K_Qq_SR1, m_E_SR1, alpha_Qq_SR1, P_ED_max1, m_P_ED1,
                D_S1, D_I1, D_F1, D_R1, D_C1, Tlag1, Lag_RR1, Lag_FR1, Lag_SR1, option_i1, option_u1, option_f1, option_s1, option_c1,w_res1, i_res1, r_res1, u_res1, f_res1,s_res1,c_res1)
  
  flow_2 <- CT2(Ce2,Cp_WR2, m_Q_WR2, Kq_WR2, Tp_WR2, Tm_WR2, K_Qq_FR2, K_Qb_FR2,  m_E_FR2, alpha_Qq_FR2, Beta_Qq_UR2, Smax_UR2, Beta_E_UR2,
                SiniFr_UR2, K_Qb_UR2, mu_Qq_UR2, Smax_IR2, m_QE_IR2, K_Qq_RR2, Smax_CR2, Umax_uCR2, Smin_uCR2, Beta_Qq_uCR2, mu_Qq_uCR2, K_Qb_uCR2, Sevmax_CR2, Beta_E_CR2, Beta_Qq_sCR2, K_Qb_sCR2, K_Qd_sCR2, K_Qq_SR2, m_E_SR2, alpha_Qq_SR2, P_ED_max2, m_P_ED2,
                D_S2, D_I2, D_F2, D_R2, D_C2, Tlag2, Lag_RR2, Lag_FR2, Lag_SR2, option_i2, option_u2, option_f2, option_s2, option_c2,w_res2, i_res2, r_res2, u_res2, f_res2,s_res2,c_res2)
  
  flow_3 <- CT3(Ce3,Cp_WR3, m_Q_WR3, Kq_WR3, Tp_WR3, Tm_WR3, K_Qq_FR3, K_Qb_FR3,  m_E_FR3, alpha_Qq_FR3, Beta_Qq_UR3, Smax_UR3, Beta_E_UR3,
                SiniFr_UR3, K_Qb_UR3, mu_Qq_UR3, Smax_IR3, m_QE_IR3, K_Qq_RR3, Smax_CR3, Umax_uCR3, Smin_uCR3, Beta_Qq_uCR3, mu_Qq_uCR3, K_Qb_uCR3, Sevmax_CR3, Beta_E_CR3, Beta_Qq_sCR3, K_Qb_sCR3, K_Qd_sCR3, K_Qq_SR3, m_E_SR3, alpha_Qq_SR3, P_ED_max, m_P_ED,
                D_S3, D_I3, D_F3, D_R3, D_C3, Tlag3, Lag_RR3, Lag_FR3, Lag_SR3, option_i3, option_u3, option_f3, option_s3, option_c3,w_res3, i_res3, r_res3, u_res3, f_res3,s_res3,c_res3)
  
  n = length(flow_1)
  
  output = numeric(n)
  
  routed_flow_1 <- fuserouting.sim(flow_1,R_1,1,R_t1)
  routed_flow_2 <- fuserouting.sim(flow_2,R_2,1,R_t2)
  routed_flow_3 <- fuserouting.sim(flow_3,R_3,1,R_t3)
  
  output <- routed_flow_1 + routed_flow_2 + routed_flow_3
  
  return(output)
  
}

#########################
WR <- function(Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR){
  modPars = c(Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR)
  names(modPars) = names(formals(WR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR))
}

IR <- function(Smax_IR, m_QE_IR, D_I, option_i){
  modPars = c(Smax_IR, m_QE_IR, D_I, option_i)
  names(modPars) = names(formals(IR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(Smax_IR, m_QE_IR, D_I, option_i))
}

RR <- function(K_Qq_RR, D_R, Lag_RR){
  modPars = c(K_Qq_RR, D_R, Lag_RR)
  names(modPars) = names(formals(RR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(K_Qq_RR, D_R, Lag_RR))
}

UR <- function(Beta_Qq_UR, Smax_UR, Beta_E_UR, SiniFr_UR, K_Qb_UR, mu_Qq_UR, option_u){
  modPars = c(Beta_Qq_UR, Smax_UR, Beta_E_UR, SiniFr_UR, K_Qb_UR, mu_Qq_UR, option_u)
  names(modPars) = names(formals(UR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(Beta_Qq_UR, Smax_UR, Beta_E_UR, SiniFr_UR, K_Qb_UR, mu_Qq_UR, option_u))
}

FR <- function(K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, D_F, Lag_FR, option_f){
  modPars = c(K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, D_F, Lag_FR, option_f)
  names(modPars) = names(formals(FR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, D_F, Lag_FR, option_f))
}

SR <- function(K_Qq_SR, m_E_SR, alpha_Qq_SR, Lag_SR, option_s ){
  modPars = c(K_Qq_SR, m_E_SR, alpha_Qq_SR, Lag_SR, option_s )
  names(modPars) = names(formals(SR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(K_Qq_SR, m_E_SR, alpha_Qq_SR, Lag_SR, option_s))
}

CR <- function(Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, D_C, option_c){
  modPars = c(Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, D_C, option_c)
  names(modPars) = names(formals(CR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, D_C, option_c)) 
}

CT_new <- function(WR, IR, RR, UR, FR, SR, CR, Ce, P_ED_max, m_P_ED, Tlag, D_S, w_res, i_res, r_res, u_res, f_res, s_res, c_res){
  
  modPars = c(WR, IR, RR, UR, FR, SR, CR, Ce, P_ED_max, m_P_ED, Tlag, D_S, w_res, i_res, r_res, u_res, f_res, s_res, c_res)
  names(modPars) = names(formals(CT_new))
  if(any(is.na(c(modPars)))) return(NA)
  
  WR_aug <- WR#(Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR, w_res)
  IR_aug <- IR#(Smax_IR, m_QE_IR, D_I, option_i, i_res)
  RR_aug <- RR#(K_Qq_RR, D_R, Lag_RR, r_res)
  UR_aug <- UR#(Beta_Qq_UR, Smax_UR, Beta_E_UR, SiniFr_UR, K_Qb_UR, mu_Qq_UR, option_u, u_res)
  FR_aug <- FR#(K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, D_F, Lag_FR, option_f, f_res)
  SR_aug <- SR#(K_Qq_SR, m_E_SR, alpha_Qq_SR, Lag_SR, option_s, s_res)
  CR_aug <- CR#(Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, D_C, option_c, c_res)
  
  output <- CT(Ce, WR_aug[1],WR_aug[2],WR_aug[3],WR_aug[4],WR_aug[5], FR_aug[1],FR_aug[2],FR_aug[3],FR_aug[4], UR_aug[1],UR_aug[2],UR_aug[3],UR_aug[4],UR_aug[5],UR_aug[6], IR_aug[1],IR_aug[2], RR_aug[1], CR_aug[1],
               CR_aug[2],CR_aug[3],CR_aug[4],CR_aug[5],CR_aug[6],CR_aug[7],CR_aug[8],CR_aug[9],CR_aug[10],CR_aug[11], SR_aug[1],SR_aug[2],SR_aug[3], P_ED_max, m_P_ED, D_S, IR_aug[3], FR_aug[5], RR_aug[2], 
               CR_aug[12], Tlag, RR_aug[3], FR_aug[6], SR_aug[4], IR_aug[4], UR_aug[7], FR_aug[7], SR_aug[5], CR_aug[13], w_res, i_res, r_res, u_res, f_res, s_res, c_res)
  
  return(output)
}


#################  MARRMoT Models ###########################

# Model 29 - HYMOD

# MM_29 <- function(M29_smax,M29_b,M29_a,M29_kf,M29_ks,M29_s1,M29_s2,M29_s3,M29_s4,M29_s5){
#   modPars = c(M29_smax,M29_b,M29_a,M29_kf,M29_ks,M29_s1,M29_s2,M29_s3,M29_s4,M29_s5)
#   names(modPars) = names(formals(MM_29))
#   if(any(is.na(c(modPars)))) return(NA)
#   scaledPars = ScaleParams(modPars)
#   mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
#   
#   M29_theta <- c(M29_smax,M29_b,M29_a,M29_kf,M29_ks)
#   M29_s0 <- c(M29_s1*M29_smax,M29_s2,M29_s3,M29_s4,M29_s5)
#   M29_theta_mat <- rvec_to_matlab(M29_theta)
#   M29_s0_mat <- rvec_to_matlab(M29_s0)
#   
#   file_location <- paste(path,"/Matlab_fn",sep = '')
#   output_script_name <- as.character(round(runif(1,1,1e6)))
#   
#   code <- c(paste("cd '",noquote(file_location),"';",sep = ''),"Input_Var = readmatrix('Input_Variables.csv');","input_climatology.precip = Input_Var(:,3);","input_climatology.temp   = Input_Var(:,5);",
#             "input_climatology.pet    = Input_Var(:,4);","input_climatology.delta_t  = 1;","model = 'm_29_hymod_5p_5s';",paste("input_theta =",noquote(M29_theta_mat)),paste("input_s0 =",noquote(M29_s0_mat)),
#             "input_solver.name              = 'createOdeApprox_IE';","input_solver.resnorm_tolerance = 0.1;","input_solver.resnorm_maxiter   = 6;",
#             "[output_ex,...                                                              
#  output_in,...                                                              
#  output_ss,....                                                             
#  output_waterbalance] = ...                                                 
#                     feval(model,...                                         
#                           input_climatology,...                             
#                           input_s0,...                                      
#                           input_theta,...                                   
#                           input_solver);",paste("writetable(struct2table(output_ex), '",noquote(file_location),"/output_Q_MM29_",noquote(output_script_name),".csv","')",sep = ''))
#   
#   run_matlab_code(code)
#   matlab_out <-read.csv(paste(file_location,"/output_Q_MM29_",output_script_name,".csv",sep = ''))
#   unlink(paste(file_location,"/output_Q_MM29_",output_script_name,".csv",sep = ''))
#   matlab_out <- t(matlab_out)
#   sim <- as.vector(matlab_out[((nrow(matlab_out)/2)+1):nrow(matlab_out),])
#   return(sim)
# }

#write.csv(M29_parameters,file.path(file_location,"Input_Parameters_MM29.csv"))
# if(have_matlab()){
#   get_matlab(try_defaults = TRUE, desktop = FALSE, splash = FALSE, display = FALSE, wait = TRUE)
# }
# system(paste("matlab -wait -nodesktop -nosplash -nodisplay -r \"run('",file_location,"/MM_29.m'); exit\"",sep = ''))

#M29_parameters <- cbind(M29_theta,M29_s0)
#paste("matlab -wait -nodesktop -nosplash -nodisplay -r \"run('",file_location,"/MM_29.m'); exit\"",sep = '')

# code <- c("cd 'C:/Users/E0149661/Desktop/R/Code_Server/MARMoT/Matlab_fn';","Input_Var = readmatrix('Input_Variables.csv');","input_climatology.precip = Input_Var(:,3);","input_climatology.temp   = Input_Var(:,5);",
#           "input_climatology.pet    = Input_Var(:,4);","input_climatology.delta_t  = 1;","model = 'm_29_hymod_5p_5s';","input_theta     = [ 35;3.7;0.4;0.25;0.01];","input_s0       = [15;7;3;8;22];",
#           "input_solver.name              = 'createOdeApprox_IE';","input_solver.resnorm_tolerance = 0.1;","input_solver.resnorm_maxiter   = 6;",
#           "[output_ex,...                                                              
#  output_in,...                                                              
#  output_ss,....                                                             
#  output_waterbalance] = ...                                                 
#                     feval(model,...                                         
#                           input_climatology,...                             
#                           input_s0,...                                      
#                           input_theta,...                                   
#                           input_solver);","writetable(struct2table(output_ex),'output_Q_MM29.csv')")




######## source('CHMGP_DEoptim.R')
# source('CHMGP_DDS.R')

DDS <- function(x,n) {
  r <- 0.2
  n <- 10000
  sBest <- x
  sCur <- x
  NS_In <- NSO(x)
  if (is.na(NS_In)){
    NS_In <- 100000
  } else {
    NS_In <- NS_In
  }
  CostBest <- NS_In
  dimen <- length(x)
  x_min <- 0
  x_max <- 1
  x_range <- x_max -x_min
  k <- 0
  for(i in 1:n) {
    for(j in 1:dimen) {
      if (runif(1) <(1-(log10(i)/log10(n)))) {
        k <- k +1
        sCur[j] <- sCur[j] + rnorm(1)*r*x_range
        if(sCur[j]<x_min){
          sCur[j] <- x_min + (x_min - sCur[j])
          if(sCur[j]>x_max){
            sCur[j] <-x_min
          }
        }
        else if(sCur[j]> x_max){
          sCur[j] <- x_max -(sCur[j]-x_max)
          if(sCur[j]< x_min){
            sCur[j] <- x_max
          }
        }
      }
    }
    if(k==0){
      index = round(runif(1,1,dimen))
      sCur[index] <- sCur[index] + rnorm(1)*r*x_range
      if(sCur[index]<x_min){
        sCur[index] <- x_min + (x_min - sCur[index])
        if(sCur[index]>x_max){
          sCur[index] <-x_min
        }
      }
      else if(sCur[index]> x_max){
        sCur[index] <- x_max -(sCur[index]-x_max)
        if(sCur[index]< x_min){
          sCur[index] <- x_max
        }
      }
    }
    k <- 0
    NS <- NSO(sCur)
    if (is.na(NS)){
      NS <- 100000
    } else {
      NS <- NS
    }
    if(NS < CostBest){
      sBest <- sCur
      CostBest <- NS
    } else {
      sBest <- sBest
      CostBest <- CostBest
    }
  }
  return(sBest)
}

NSO <- function(x){
  sim <- CT_new(WR(x[1],x[2],x[3],x[4],x[5]),IR(x[6],x[7],x[8],x[9]),RR(x[10],x[11],x[12]),UR(x[13],x[14],x[15],x[16],x[17],x[18],x[19]),FR(x[20],x[21],x[22],x[23],x[24],x[25],x[26]),SR(x[27],x[28],x[29],x[30],x[31]),CR(x[32],x[33],x[34],x[35],x[36],x[37],x[38],x[39],x[40],x[41],x[42],x[43],x[44]),x[45],x[46],x[47],x[48],x[49],x[50],x[51],x[52],x[53],x[54],x[55],x[56])
  obs <- DataSet[,5]
  NSO_value <- FitnessComputation(sim,obs,"NS0")
  if (is.na(NSO_value)){
    NSO_value <- 100000
  } else {
    NSO_value <- NSO_value
  }
  return(NSO_value)
}




# source('New_Marrmot_functions.R')

# make sure values between 1 and 0
Code_name <- function(x1,x2,x3,x4,x5,x6){
  Marmot <- "code_M01"
  if(x1<=0.26){
    if(x2<=0.25){
      if(x3<=0.67){
        if(x5<=0.5){
          Marmot <- "code_M44"
        } else{
          Marmot <- "code_M34"
        }
      } else{
        Marmot <- "code_M45"
      }
    } else {
      if(x3<=0.33){
        if(x5<=0.33){
          Marmot <- "code_M43"
        } else{
          if(x6<=0.5){
            Marmot <- "code_M06"
          } else{
            Marmot <- "code_M12"
          }
        }
      } else if(x3<=0.56 && x3>0.33){
        if(x6<=0.2){
          Marmot <- "code_M30"
        } else if(x6>0.2 && x6<=0.4){
          Marmot <- "code_M31"
        } else if(x6>0.4 && x6<=0.6){
          Marmot <- "code_M32"
        } else if(x6>0.6 && x6<=0.8){
          Marmot <- "code_M35"
        } else{
          Marmot <- "code_M41"
        }
      } else{
        Marmot <- "code_M37"
      }
    }
  } else{
    if(x2<=0.24){
      if(x4<=0.5){
        if(x5<=0.25){
          Marmot <- "code_M39"
        } else{
          if(x6<=0.33){
            Marmot <- "code_M18"
          } else if(x6>0.33 && x6<=0.67){
            Marmot <- "code_M22"
          } else{
            Marmot <- "code_M36"
          }
        }
      } else{
        if(x5<=0.25){
          Marmot <- "code_M42"
        } else if(x5>0.75){
          Marmot <- "code_M16"
        } else{
          if(x6<=0.5){
            Marmot <- "code_M07"
          } else{
            Marmot <- "code_M26"
          }
        }
      }
    } else{
      if(x3<=0.42){
        if(x4<=0.36){
          if(x5<=0.5){
            if(x6<=0.5){
              Marmot <- "code_M09"
            } else{
              Marmot <- "code_M11"
            }
          } else{
            if(x6<=0.5){
              Marmot <- "code_M13"
            } else{
              Marmot <- "code_M15"
            }
          }
        } else {
          if(x5<=0.14){
            Marmot <- "code_M29"
          } else if(x5>0.14 && x5<=0.28){
            Marmot <- "code_M05"
          } else if(x5>0.28 && x5<=0.86){
            if(x6<=0.25){
              Marmot <- "code_M01"
            } else if(x6>0.25 && x6<=0.5){
              Marmot <- "code_M02"
            } else if(x6>0.5 && x6<=0.75){
              Marmot <- "code_M03"
            } else {
              Marmot <- "code_M04"
            }
          } else{
            Marmot <- "code_M21"
          }
        }
      } else if(x3>0.42 && x3<=0.77){
        if(x4<=0.33){
          if(x6<=0.33){
            Marmot <- "code_M19"
          } else if (x6>0.33 && x6<=0.67){
            Marmot <- "code_M20"
          } else{
            Marmot <- "code_M23"
          }
        } else{
          if(x5<=0.5){
            if(x6<=0.33){
              Marmot <- "code_M08"
            } else if(x6>0.33 && x6<=0.67){
              Marmot <- "code_M10"
            } else{
              Marmot <- "code_M14"
            }
          } else{
            if(x6<=0.33){
              Marmot <-"code_M17"
            } else if(x6>0.33 && x6<=0.67){
              Marmot <- "code_M24"
            } else {
              Marmot <- "code_M28"
            }
          }
        }
      } else if(x3>0.77 && x3<=0.89){
        if(x5<=0.67){
          if(x6<=0.5){
            Marmot <- "code_M25"
          } else{
            Marmot <- "code_M27"
          }
        } else{
          Marmot <- "code_M46"
        }
      } else{
        if(x4<=0.33){
          Marmot <- "code_M40"
        } else{
          if(x6<=0.5){
            Marmot <- "code_M33"
          } else{
            Marmot <- "code_M38"
          }
        }
      }
    }
  }
  return(Marmot)
}
########## Code function ####################
Code_generator <- function(Model,Para_vec,Script_name){
  Code <-"empty"
  file_location <- paste(path,"/Matlab_fn",sep = '')
  
  # M01
  if(Model=="code_M01"){
    modPars <- Para_vec[1:2]
    names(modPars) <-c("M01_smax","M01_s1")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_01_collie1_1p_1s';"
    theta <-c(M01_smax)
    s0 <-c(M01_smax*M01_s1)
  }
  
  # M02
  if(Model=="code_M02"){
    modPars <- Para_vec[1:5]
    names(modPars) <-c("M02_dw","M02_betaw","M02_swmax","M02_kw","M02_s1")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_02_wetland_4p_1s';"
    theta <-c(M02_dw,M02_betaw,M02_swmax,M02_kw)
    s0 <-c(M02_swmax*M02_s1)
  }
  
  # M03
  if(Model=="code_M03"){
    modPars <- Para_vec[1:5]
    names(modPars) <-c("M03_S1max","M03_Sfc","M03_a","M03_M","M03_s1")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_03_collie2_4p_1s';"
    theta <-c(M03_S1max,M03_Sfc,M03_a,M03_M)
    s0 <-c(M03_S1max*M03_s1)
  }
  
  # M04
  if(Model=="code_M04"){
    modPars <- Para_vec[1:7]
    names(modPars) <-c("M04_S1max","M04_Sfc","M04_m","M04_a","M04_b","M04_tcbf","M04_s1")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_04_newzealand1_6p_1s';"
    theta <-c(M04_S1max,M04_Sfc,M04_m,M04_a,M04_b,M04_tcbf)
    s0 <-c(M04_S1max*M04_s1)
  }
  
  # M05
  if(Model=="code_M05"){
    modPars <- Para_vec[1:8]
    names(modPars) <-c("M05_lp","M05_d","M05_p","M05_alpha","M05_tau_q","M05_tau_s","M05_tau_d","M05_s1")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_05_ihacres_7p_1s';"
    theta <-c(M05_lp,M05_d,M05_p,M05_alpha,M05_tau_q,M05_tau_s,M05_tau_d)
    s0 <-c(M05_s1)
  }
  
  # M06
  if(Model=="code_M06"){
    modPars <- Para_vec[1:6]
    names(modPars) <-c("M06_tt","M06_ddf","M06_smax","M06_tc","M06_s1","M06_s2")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_06_alpine1_4p_2s';"
    theta <-c(M06_tt,M06_ddf,M06_smax,M06_tc)
    s0 <-c(M06_s1,M06_s2*M06_smax)
  }
  
  # M07
  if(Model=="code_M07"){
    modPars <- Para_vec[1:6]
    names(modPars) <-c("M07_x1","M07_x2","M07_x3","M07_x4","M07_s1","M07_s2")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_07_gr4j_4p_2s';"
    theta <-c(M07_x1,M07_x2,M07_x3,M07_x4)
    s0 <-c(M07_s1*M07_x1,M07_s2*M07_x3)
  }
  
  # M08
  if(Model=="code_M08"){
    modPars <- Para_vec[1:7]
    names(modPars) <-c("M08_alpha_ei","M08_m","M08_smax","M08_fc","M08_alpha_ss","M08_s1","M08_s2")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_08_us1_5p_2s';"
    theta <-c(M08_alpha_ei,M08_m,M08_smax,M08_fc,M08_alpha_ss)
    s0 <-c(M08_s1*M08_smax,M08_s2*(M08_smax-M08_s1*M08_smax))
  }
  
  # M09
  if(Model=="code_M09"){
    modPars <- Para_vec[1:8]
    names(modPars) <-c("M09_sb","M09_sfc","M09_m","M09_a","M09_b","M09_r","M09_s1","M09_s2")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_09_susannah1_6p_2s';"
    theta <-c(M09_sb,M09_sfc,M09_m,M09_a,M09_b,M09_r)
    s0 <-c(M09_s1*M09_sb,M09_s2*(M09_sb-M09_s1*M09_sb))
  }
  
  # M10
  if(Model=="code_M10"){
    modPars <- Para_vec[1:8]
    names(modPars) <-c("M10_sb","M10_phi","M10_fc","M10_r","M10_c","M10_d","M10_s1","M10_s2")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_10_susannah2_6p_2s';"
    theta <-c(M10_sb,M10_phi,M10_fc,M10_r,M10_c,M10_d)
    s0 <-c(M10_s1*M10_sb,M10_s2*(M10_sb-M10_s1*M10_sb))
  }
  
  # M11
  if(Model=="code_M11"){
    modPars <- Para_vec[1:8]
    names(modPars) <-c("M11_s1max","M11_sfc","M11_a","M11_M","M11_b","M11_lambda","M11_s1","M11_s2")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_11_collie3_6p_2s';"
    theta <-c(M11_s1max,M11_sfc,M11_a,M11_M,M11_b,M11_lambda)
    s0 <-c(M11_s1*M11_s1max,M11_s2)
  }
  
  # M12
  if(Model=="code_M12"){
    modPars <- Para_vec[1:8]
    names(modPars) <-c("M12_tt","M12_ddf","M12_smax","M12_Cfc","M12_tcin","M12_tcbf","M12_s1","M12_s2")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_12_alpine2_6p_2s';"
    theta <-c(M12_tt,M12_ddf,M12_smax,M12_Cfc,M12_tcin,M12_tcbf)
    s0 <-c(M12_s1,M12_s2*M12_smax)
  }
  
  # M13
  if(Model=="code_M13"){
    modPars <- Para_vec[1:9]
    names(modPars) <-c("M13_dw","M13_betaw","M13_swmax","M13_a","M13_th","M13_c","M13_kh","M13_s1","M13_s2")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_13_hillslope_7p_2s';"
    theta <-c(M13_dw,M13_betaw,M13_swmax,M13_a,M13_th,M13_c,M13_kh)
    s0 <-c(M13_s1*M13_swmax,M13_s2)
  }
  
  # M14
  if(Model=="code_M14"){
    modPars <- Para_vec[1:9]
    names(modPars) <-c("M14_suzmax","M14_st","M14_kd","M14_q0","M14_f","M14_chi","M14_phi","M14_s1","M14_s2")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_14_topmodel_7p_2s';"
    theta <-c(M14_suzmax,M14_st,M14_kd,M14_q0,M14_f,M14_chi,M14_phi)
    s0 <-c(M14_s1*M14_suzmax,M14_s2)
  }
  
  # M15
  if(Model=="code_M15"){
    modPars <- Para_vec[1:10]
    names(modPars) <-c("M15_fmax","M15_dp","M15_sumax","M15_lp","M15_p","M15_tp","M15_c","M15_kp","M15_s1","M15_s2")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_15_plateau_8p_2s';"
    theta <-c(M15_fmax,M15_dp,M15_sumax,M15_lp,M15_p,M15_tp,M15_c,M15_kp)
    s0 <-c(M15_s1*M15_sumax,M15_s2)
  }
  
  # M16
  if(Model=="code_M16"){
    modPars <- Para_vec[1:10]
    names(modPars) <-c("M16_s1max","M16_s2max","M16_sfc","M16_m","M16_a","M16_b","M16_tcbf","M16_d","M16_s1","M16_s2")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_16_newzealand2_8p_2s';"
    theta <-c(M16_s1max,M16_s2max,M16_sfc,M16_m,M16_a,M16_b,M16_tcbf,M16_d)
    s0 <-c(M16_s1*M16_s1max,M16_s2*M16_s2max)
  }
  
  # M17
  if(Model=="code_M17"){
    modPars <- Para_vec[1:7]
    names(modPars) <-c("M17_smax","M17_phi","M17_gam","M17_kl","M17_s1","M17_s2","M17_s3")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_17_penman_4p_3s';"
    theta <-c(M17_smax,M17_phi,M17_gam,M17_kl)
    s0 <-c(M17_s1,M17_s2,M17_s3)
  }
  
  # M18
  if(Model=="code_M18"){
    modPars <- Para_vec[1:10]
    names(modPars) <-c("M18_insc","M18_coeff","M18_sq","M18_smsc","M18_sub","M18_crak","M18_k","M18_s1","M18_s2","M18_s3")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_18_simhyd_7p_3s';"
    theta <-c(M18_insc,M18_coeff,M18_sq,M18_smsc,M18_sub,M18_crak,M18_k)
    s0 <-c(M18_s1*M18_insc,M18_s2*M18_smsc,M18_s3)
  }
  
  # M19
  if(Model=="code_M19"){
    modPars <- Para_vec[1:11]
    names(modPars) <-c("M19_sb","M19_phi","M19_fc","M19_alpha_ss","M19_beta_ss","M19_k_deep","M19_alpha_bf","M19_beta_bf","M19_s1","M19_s2","M19_s3")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_19_australia_8p_3s';"
    theta <-c(M19_sb,M19_phi,M19_fc,M19_alpha_ss,M19_beta_ss,M19_k_deep,M19_alpha_bf,M19_beta_bf)
    s0 <-c(M19_s1*M19_sb,M19_s2*(M19_sb-M19_s1*M19_sb),M19_s3)
  }
  
  # M20
  if(Model=="code_M20"){
    modPars <- Para_vec[1:11]
    names(modPars) <-c("M20_c","M20_ndc","M20_smax","M20_emax","M20_frate","M20_b","M20_dpf","M20_sdrmax","M20_s1","M20_s2","M20_s3")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_20_gsfb_8p_3s';"
    theta <-c(M20_c,M20_ndc,M20_smax,M20_emax,M20_frate,M20_b,M20_dpf,M20_sdrmax)
    s0 <-c(M20_s1,M20_s2*M20_smax,M20_s3*(M20_smax-M20_s2*M20_smax))
  }
  
  # M21
  if(Model=="code_M21"){
    modPars <- Para_vec[1:12]
    names(modPars) <-c("M21_s1max","M21_beta","M21_d","M21_percmax","M21_lp","M21_nlagf","M21_nalgs","M21_kf","M21_ks","M21_s1","M21_s2","M21_s3")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_21_flexb_9p_3s';"
    theta <-c(M21_s1max,M21_beta,M21_d,M21_percmax,M21_lp,M21_nlagf,M21_nalgs,M21_kf,M21_ks)
    s0 <-c(M21_s1*M21_s1max,M21_s2,M21_s3)
  }
  
  # M22
  if(Model=="code_M22"){
    modPars <- Para_vec[1:13]
    names(modPars) <-c("M22_ibar","M22_idelta","M22_ishift","M22_stot","M22_fsm","M22_b","M22_k1","M22_c1","M22_k2","M22_c2","M22_s1","M22_s2","M22_s3")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_22_vic_10p_3s';"
    theta <-c(M22_ibar,M22_idelta,M22_ishift,M22_stot,M22_fsm,M22_b,M22_k1,M22_c1,M22_k2,M22_c2)
    s0 <-c(M22_s1*M22_ibar,M22_s2*(M22_fsm*M22_stot),M22_s3*(1-M22_fsm)*M22_stot)
  }
  
  # M23
  if(Model=="code_M23"){
    modPars <- Para_vec[1:27]
    names(modPars) <-c("M23_af","M23_bf","M23_stot","M23_xa","M23_xf","M23_na","M23_ac","M23_bc","M23_ass","M23_bss","M23_c","M23_ag","M23_bg","M23_gf","M23_df","M23_td","M23_ab","M23_bb","M23_ga","M23_da","M23_aa","M23_ba","M23_gb","M23_db","M23_s1","M23_s2","M23_s3")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_23_lascam_24p_3s';"
    theta <-c(M23_af,M23_bf,M23_stot,M23_xa,M23_xf,M23_na,M23_ac,M23_bc,M23_ass,M23_bss,M23_c,M23_ag,M23_bg,M23_gf,M23_df,M23_td,M23_ab,M23_bb,M23_ga,M23_da,M23_aa,M23_ba,M23_gb,M23_db)
    s0 <-c(M23_s1,M23_s2*M23_stot,M23_s3*M23_s2*M23_stot)
  }
  
  # M24
  if(Model=="code_M24"){
    modPars <- Para_vec[1:9]
    names(modPars) <-c("M24_s1max","M24_tw","M24_tu","M24_se","M24_tc","M24_s1","M24_s2","M24_s3","M24_s4")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_24_mopex1_5p_4s';"
    theta <-c(M24_s1max,M24_tw,M24_tu,M24_se,M24_tc)
    s0 <-c(M24_s1*M24_s1max,M24_s2,M24_s3,M24_s4)
  }
  
  # M25
  if(Model=="code_M25"){
    modPars <- Para_vec[1:10]
    names(modPars) <-c("M25_phi","M25_rc","M25_gam","M25_k1","M25_fa","M25_k2","M25_s1","M25_s2","M25_s3","M25_s4")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_25_tcm_6p_4s';"
    theta <-c(M25_phi,M25_rc,M25_gam,M25_k1,M25_fa,M25_k2)
    s0 <-c(M25_s1,M25_s2,M25_s3,M25_s4)
  }
  
  # M26
  if(Model=="code_M26"){
    modPars <- Para_vec[1:14]
    names(modPars) <-c("M26_smax","M26_beta","M26_d","M26_percmax","M26_lp","M26_nlagf","M26_nlags","M26_kf","M26_ks","M26_imax","M26_s1","M26_s2","M26_s3","M26_s4")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_26_flexi_10p_4s';"
    theta <-c(M26_smax,M26_beta,M26_d,M26_percmax,M26_lp,M26_nlagf,M26_nlags,M26_kf,M26_ks,M26_imax)
    s0 <-c(M26_s1*M26_imax,M26_s2*M26_smax,M26_s3,M26_s4)
  }
  
  # M27
  if(Model=="code_M27"){
    modPars <- Para_vec[1:16]
    names(modPars) <-c("M27_a0","M27_b0","M27_c0","M27_a1","M27_fa","M27_fb","M27_fc","M27_fd","M27_st","M27_f2","M27_f1","M27_f3","M27_s1","M27_s2","M27_s3","M27_s4")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_27_tank_12p_4s';"
    theta <-c(M27_a0,M27_b0,M27_c0,M27_a1,M27_fa,M27_fb,M27_fc,M27_fd,M27_st,M27_f2,M27_f1,M27_f3)
    ly1 <-M27_s1*M27_st
    ly2 <-M27_s2*(M27_st-ly1)
    ly3 <-M27_s3*(M27_st-ly1-ly2)
    s0 <-c(ly1,ly2,ly3,M27_s4*(M27_st-ly1-ly2-ly3))
  }
  
  # M28
  if(Model=="code_M28"){
    modPars <- Para_vec[1:16]
    names(modPars) <-c("M28_aim","M28_a","M28_b","M28_stot","M28_fwmx","M28_flm","M28_c","M28_ex","M28_ki","M28_kg","M28_ci","M28_cg","M28_s1","M28_s2","M28_s3","M28_s4")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_28_xinanjiang_12p_4s';"
    theta <-c(M28_aim,M28_a,M28_b,M28_stot,M28_fwmx,M28_flm,M28_c,M28_ex,M28_ki,M28_kg,M28_ci,M28_cg)
    s0 <-c(M28_s1*M28_stot*M28_fwmx,M28_s2*(1-M28_fwmx)*M28_stot,M28_s3,M28_s4)
  }
  
  # M29
  if(Model=="code_M29"){
    modPars <- Para_vec[1:10]
    names(modPars) <-c("M29_smax","M29_b","M29_a","M29_kf","M29_ks","M29_s1","M29_s2","M29_s3","M29_s4","M29_s5")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_29_hymod_5p_5s';"
    theta <- c(M29_smax,M29_b,M29_a,M29_kf,M29_ks)
    s0 <- c(M29_s1*M29_smax,M29_s2,M29_s3,M29_s4,M29_s5)
  }
  
  # M30
  if(Model=="code_M30"){
    modPars <- Para_vec[1:12]
    names(modPars) <-c("M30_tcrit","M30_ddf","M30_s2max","M30_tw","M30_tu","M30_se","M30_tc","M30_s1","M30_s2","M30_s3","M30_s4","M30_s5")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_30_mopex2_7p_5s';"
    theta <- c(M30_tcrit,M30_ddf,M30_s2max,M30_tw,M30_tu,M30_se,M30_tc)
    s0 <- c(M30_s1,M30_s2*M30_s2max,M30_s3,M30_s4,M30_s5)
  }
  
  # M31
  if(Model=="code_M31"){
    modPars <- Para_vec[1:13]
    names(modPars) <-c("M31_tcrit","M31_ddf","M31_s2max","M31_tw","M31_tu","M31_se","M31_s3max","M31_tc","M31_s1","M31_s2","M31_s3","M31_s4","M31_s5")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_31_mopex3_8p_5s';"
    theta <- c(M31_tcrit,M31_ddf,M31_s2max,M31_tw,M31_tu,M31_se,M31_s3max,M31_tc)
    s0 <- c(M31_s1,M31_s2*M31_s2max,M31_s3*M31_s3max,M31_s4,M31_s5)
  }
  
  # M32
  if(Model=="code_M32"){
    modPars <- Para_vec[1:15]
    names(modPars) <-c("M32_tcrit","M32_ddf","M32_s2max","M32_tw","M32_i_alpha","M32_i_s","M32_tu","M32_se","M32_s3max","M32_tc","M32_s1","M32_s2","M32_s3","M32_s4","M32_s5")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_32_mopex4_10p_5s';"
    theta <- c(M32_tcrit,M32_ddf,M32_s2max,M32_tw,M32_i_alpha,M32_i_s,M32_tu,M32_se,M32_s3max,M32_tc)
    s0 <- c(M32_s1,M32_s2*M32_s2max,M32_s3*M32_s3max,M32_s4,M32_s5)
  }
  
  # M33
  if(Model=="code_M33"){
    modPars <- Para_vec[1:16]
    names(modPars) <-c("M33_pctim","M33_smax","M33_f1","M33_f2","M33_kuz","M33_rexp","M33_f3","M33_f4","M33_pfree","M33_klzp","M33_klzs","M33_s1","M33_s2","M33_s3","M33_s4","M33_s5")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_33_sacramento_11p_5s';"
    theta <- c(M33_pctim,M33_smax,M33_f1,M33_f2,M33_kuz,M33_rexp,M33_f3,M33_f4,M33_pfree,M33_klzp,M33_klzs)
    uztwm <- M33_smax*M33_f1
    uzfwm <- M33_f2*(M33_smax-uztwm)
    lztwm <- M33_f3*(M33_smax-uztwm-uzfwm)
    lzfwpm <- M33_f4*(M33_smax -uztwm-uzfwm-lztwm)
    lzfwsm <- (1-M33_f4)*(M33_smax -uztwm-uzfwm-lztwm)
    
    s0 <- c(M33_s1*uztwm,M33_s2*uzfwm,M33_s3*lztwm,M33_s4*lzfwpm,M33_s5*lzfwsm)
  }
  
  # M34
  if(Model=="code_M34"){
    modPars <- Para_vec[1:17]
    names(modPars) <-c("M34_smax","M34_beta","M34_d","M34_percmax","M34_lp","M34_nlagf","M34_nlags","M34_kf","M34_ks","M34_imax","M34_tt","M34_ddf","M34_s1","M34_s2","M34_s3","M34_s4","M34_s5")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_34_flexis_12p_5s';"
    theta <- c(M34_smax,M34_beta,M34_d,M34_percmax,M34_lp,M34_nlagf,M34_nlags,M34_kf,M34_ks,M34_imax,M34_tt,M34_ddf)
    s0 <- c(M34_s1,M34_s2*M34_imax,M34_s3*M34_smax,M34_s4,M34_s5)
  }
  
  # M35
  if(Model=="code_M35"){
    modPars <- Para_vec[1:17]
    names(modPars) <-c("M35_tcrit","M35_ddf","M35_s2max","M35_tw","M35_i_alpha","M35_i_s","M35_tmin","M35_trange","M35_tu","M35_se","M35_s3max","M35_tc","M35_s1","M35_s2","M35_s3","M35_s4","M35_s5")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_35_mopex5_12p_5s';"
    theta <- c(M35_tcrit,M35_ddf,M35_s2max,M35_tw,M35_i_alpha,M35_i_s,M35_tmin,M35_trange,M35_tu,M35_se,M35_s3max,M35_tc)
    s0 <- c(M35_s1,M35_s2*M35_s2max,M35_s3*M35_s3max,M35_s4,M35_s5)
  }
  
  # M36
  if(Model=="code_M36"){
    modPars <- Para_vec[1:20]
    names(modPars) <-c("M36_INSC","M36_COEFF","M36_SQ","M36_SMSC","M36_SUB","M36_CRAK","M36_EM","M36_DSC","M36_ADS","M36_MD","M36_VCOND","M36_DLEV","M36_K1","M36_K2","M36_K3","M36_s1","M36_s2","M36_s3","M36_s4","M36_s5")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_36_modhydrolog_15p_5s';"
    theta <- c(M36_INSC,M36_COEFF,M36_SQ,M36_SMSC,M36_SUB,M36_CRAK,M36_EM,M36_DSC,M36_ADS,M36_MD,M36_VCOND,M36_DLEV,M36_K1,M36_K2,M36_K3)
    s0 <- c(M36_s1*M36_INSC,M36_s2*M36_SMSC,M36_s3*M36_DSC,M36_s4,M36_s5)
  }
  
  # M37
  if(Model=="code_M37"){
    modPars <- Para_vec[1:20]
    names(modPars) <-c("M37_tt","M37_tti","M37_ttm","M37_cfr","M37_cfmax","M37_whc","M37_cflux","M37_fc","M37_lp","M37_beta","M37_k0","M37_alpha","M37_perc","M37_K1","M37_maxbas","M37_s1","M37_s2","M37_s3","M37_s4","M37_s5")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_37_hbv_15p_5s';"
    theta <- c(M37_tt,M37_tti,M37_ttm,M37_cfr,M37_cfmax,M37_whc,M37_cflux,M37_fc,M37_lp,M37_beta,M37_k0,M37_alpha,M37_perc,M37_K1,M37_maxbas)
    s0 <- c(M37_s1,M37_s2*M37_s1*M37_whc,M37_s3*M37_fc,M37_s4,M37_s5)
  }
  
  # M38
  if(Model=="code_M38"){
    modPars <- Para_vec[1:21]
    names(modPars) <-c("M38_a0","M38_b0","M38_c0","M38_a1","M38_fa","M38_fb","M38_fc","M38_fd","M38_st","M38_f2","M38_f1","M38_f3","M38_k1","M38_k2","M38_z1","M38_z2","M38_s1","M38_s2","M38_s3","M38_s4","M38_s5")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_38_tank2_16p_5s';"
    theta <- c(M38_a0,M38_b0,M38_c0,M38_a1,M38_fa,M38_fb,M38_fc,M38_fd,M38_st,M38_f2,M38_f1,M38_f3,M38_k1,M38_k2,M38_z1,M38_z2)
    ly1 <- M38_st*M38_s1
    ly2 <- M38_s2*(M38_st-ly1)
    ly3 <- M38_s3*(M38_st-ly1-ly2)
    ly4 <- M38_s4*(M38_st-ly1-ly2-ly3)
    s0 <- c(ly1,ly2,ly3,ly4,M38_s5*(M38_st-ly1-ly2-ly3-ly4))
  }
  
  # M39
  if(Model=="code_M39"){
    modPars <- Para_vec[1:21]
    names(modPars) <-c("M39_smax","M39_cmax","M39_ct","M39_c1","M39_ce","M39_dsurp","M39_kd","M39_gamd","M39_qpmax","M39_kg","M39_tau","M39_sbf","M39_kcr","M39_gamcr","M39_kor","M39_gamor","M39_s1","M39_s2","M39_s3","M39_s4","M39_s5")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_39_mcrm_16p_5s';"
    theta <- c(M39_smax,M39_cmax,M39_ct,M39_c1,M39_ce,M39_dsurp,M39_kd,M39_gamd,M39_qpmax,M39_kg,M39_tau,M39_sbf,M39_kcr,M39_gamcr,M39_kor,M39_gamor)
    s0 <- c(M39_s1*M39_smax,M39_s2*M39_dsurp,M39_s3,M39_s4*M39_sbf,M39_s5*(M39_sbf-M39_s4*M39_sbf))
  }
  
  # M40
  if(Model=="code_M40"){
    modPars <- Para_vec[1:14]
    names(modPars) <-c("M40_h","M40_y","M40_smax","M40_c","M40_g","M40_kg","M40_n","M40_nk","M40_s1","M40_s2","M40_s3","M40_s4","M40_s5","M40_s6")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_40_smar_8p_6s';"
    theta <- c(M40_h,M40_y,M40_smax,M40_c,M40_g,M40_kg,M40_n,M40_nk)
    st1 <- M40_s1*M40_smax
    st2 <- M40_s2*(M40_smax - st1)
    st3 <- M40_s3*(M40_smax - st1 - st2)
    st4 <- M40_s4*(M40_smax - st1 - st2 - st3)
    s0 <- c(st1,st2,st3,st4,M40_s5*(M40_smax - st1 - st2 - st3 - st4),M40_s6)
  }
  
  # M41
  if(Model=="code_M41"){
    modPars <- Para_vec[1:16]
    names(modPars) <-c("M41_cs","M41_cif","M41_stot","M41_cl1","M41_f1","M41_cof","M41_cl2","M41_k0","M41_k1","M41_kb","M41_s1","M41_s2","M41_s3","M41_s4","M41_s5","M41_s6")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_41_nam_10p_6s';"
    theta <- c(M41_cs,M41_cif,M41_stot,M41_cl1,M41_f1,M41_cof,M41_cl2,M41_k0,M41_k1,M41_kb)
    s0 <- c(M41_s1,M41_s2*(1-M41_f1)*M41_stot,M41_s3*M41_f1*M41_stot,M41_s4,M41_s5,M41_s6)
  }
  
  # M42
  if(Model=="code_M42"){
    modPars <- Para_vec[1:18]
    names(modPars) <-c("M42_c","M42_imax","M42_a","M42_fi2","M42_kin","M42_d50","M42_fd16","M42_sbc","M42_kb","M42_pb","M42_kh","M42_kc","M42_s1","M42_s2","M42_s3","M42_s4","M42_s5","M42_s6")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_42_hycymodel_12p_6s';"
    theta <- c(M42_c,M42_imax,M42_a,M42_fi2,M42_kin,M42_d50,M42_fd16,M42_sbc,M42_kb,M42_pb,M42_kh,M42_kc)
    s0 <- c(M42_s1*(1-M42_fi2)*M42_imax,M42_s2*M42_fi2*M42_imax,M42_s3,M42_s4,M42_s5,M42_s6)
  }
  
  # M43
  if(Model=="code_M43"){
    modPars <- Para_vec[1:18]
    names(modPars) <-c("M43_fice","M43_t0","M43_asnow","M43_tm","M43_ks","M43_aice","M43_ki","M43_a","M43_x","M43_y","M43_ksl","M43_beta","M43_s1","M43_s2","M43_s3","M43_s4","M43_s5","M43_s6")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_43_gsmsocont_12p_6s';"
    theta <- c(M43_fice,M43_t0,M43_asnow,M43_tm,M43_ks,M43_aice,M43_ki,M43_a,M43_x,M43_y,M43_ksl,M43_beta)
    s0 <- c(M43_s1,M43_s2,M43_s3,M43_s4,M43_s5*M43_a,M43_s6)
  }
  
  # M44
  if(Model=="code_M44"){
    modPars <- Para_vec[1:22]
    names(modPars) <-c("M44_rho","M44_ts","M44_tm","M44_as","M44_af","M44_gmax","M44_the","M44_phi","M44_smax","M44_fsm","M44_fsw","M44_ksat","M44_c","M44_lmax","M44_kf","M44_ks","M44_s1","M44_s2","M44_s3","M44_s4","M44_s5","M44_s6")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_44_echo_16p_6s';"
    theta <- c(M44_rho,M44_ts,M44_tm,M44_as,M44_af,M44_gmax,M44_the,M44_phi,M44_smax,M44_fsm,M44_fsw,M44_ksat,M44_c,M44_lmax,M44_kf,M44_ks)
    s0 <- c(M44_s1*M44_rho,M44_s2,M44_s3,M44_s4*M44_smax,M44_s5,M44_s6)
  }
  
  # M45
  if(Model=="code_M45"){
    modPars <- Para_vec[1:25]
    names(modPars) <-c("M45_tt","M45_ddf","M45_alpha","M45_beta","M45_stor","M45_retip","M45_fscn","M45_scx","M45_flz","M45_stot","M45_cgw","M45_resmax","M45_k1","M45_k2","M45_k3","M45_k4","M45_k5","M45_k6","M45_s1","M45_s2","M45_s3","M45_s4","M45_s5","M45_s6","M45_s7")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_45_prms_18p_7s';"
    theta <- c(M45_tt,M45_ddf,M45_alpha,M45_beta,M45_stor,M45_retip,M45_fscn,M45_scx,M45_flz,M45_stot,M45_cgw,M45_resmax,M45_k1,M45_k2,M45_k3,M45_k4,M45_k5,M45_k6)
    s0 <- c(M45_s1,M45_s2*M45_stor,M45_s3*M45_retip,M45_s4*(1-M45_flz)*M45_stot,M45_s5*M45_flz*M45_stot,M45_s6*M45_resmax,M45_s7)
  }
  
  # M46
  if(Model=="code_M46"){
    modPars <- Para_vec[1:20]
    names(modPars) <-c("M46_fap","M46_fdp","M46_dp","M46_cq","M46_d1","M46_tf","M46_fds","M46_ds","M46_d2","M46_cxq","M46_cxs","M46_cu","M46_s1","M46_s2","M46_s3","M46_s4","M46_s5","M46_s6","M46_s7","M46_s8")
    scaledPars = ScaleParams(modPars)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
    mat_model <-"model = 'm_46_classic_12p_8s';"
    theta <- c(M46_fap,M46_fdp,M46_dp,M46_cq,M46_d1,M46_tf,M46_fds,M46_ds,M46_d2,M46_cxq,M46_cxs,M46_cu)
    s0 <- c(M46_s1*M46_fdp*M46_dp,M46_s2*(1-M46_fdp)*M46_dp,M46_s3,M46_s4*M46_fds*M46_ds,M46_s5*(1-M46_fds)*M46_ds,M46_s6,M46_s7,M46_s8)
  }
  
  theta_mat <- rvec_to_matlab(theta)
  s0_mat <- rvec_to_matlab(s0)
  
  #  Code <- c(paste("cd '",noquote(file_location),"';",sep = ''),"Input_Var = readmatrix('Input_Variables.csv');","input_climatology.precip = Input_Var(:,3);","input_climatology.temp   = Input_Var(:,5);",
  #            "input_climatology.pet    = Input_Var(:,4);","input_climatology.delta_t  = 1;",mat_model,paste("input_theta =",noquote(theta_mat)),paste("input_s0 =",noquote(s0_mat)),
  #            "input_solver.name              = 'createOdeApprox_IE';","input_solver.resnorm_tolerance = 0.1;","input_solver.resnorm_maxiter   = 6;",
  #            "[output_ex,...                                                              
  # output_in,...                                                              
  # output_ss,....                                                             
  # output_waterbalance] = ...                                                 
  #                    feval(model,...                                         
  #                          input_climatology,...                             
  #                          input_s0,...                                      
  #                          input_theta,...                                   
  #                          input_solver);",paste("writetable(struct2table(output_ex), '",noquote(file_location),"/Output_Q_",noquote(Script_name),".csv","');",sep = ''))
  
  Code <-paste("cd '",noquote(file_location),"';","Input_Var = readmatrix('Input_Variables.csv');","input_climatology.precip = Input_Var(:,3);","input_climatology.temp   = Input_Var(:,5);",
               "input_climatology.pet    = Input_Var(:,4);","input_climatology.delta_t  = 1;",mat_model,"input_theta =",noquote(theta_mat),"input_s0 =",noquote(s0_mat),
               "input_solver.name              = 'createOdeApprox_IE';","input_solver.resnorm_tolerance = 0.1;","input_solver.resnorm_maxiter   = 6;",
               "[output_ex,output_in,output_ss,output_waterbalance] =feval(model,input_climatology,input_s0,input_theta,input_solver);","writetable(struct2table(output_ex), '",noquote(file_location),"/Output_Q_",noquote(Script_name),".csv","')",sep ='')
  
  # 
  # precip <-rvec_to_matlab(DataSet$P)
  # temp <-rvec_to_matlab(DataSet$T)
  # potE <- rvec_to_matlab(DataSet$E)
  # 
  # Code <-paste("input_climatology.precip =",precip,"input_climatology.temp   =",temp,
  #              "input_climatology.pet    =",potE,"input_climatology.delta_t  = 1;",mat_model,"input_theta =",noquote(theta_mat),"input_s0 =",noquote(s0_mat),
  #              "input_solver.name              = 'createOdeApprox_IE';","input_solver.resnorm_tolerance = 0.1;","input_solver.resnorm_maxiter   = 6;",
  #              "[output_ex,output_in,output_ss,output_waterbalance] =feval(model,input_climatology,input_s0,input_theta,input_solver);","writetable(struct2table(output_ex), '",noquote(file_location),"/Output_Q_",noquote(Script_name),".csv","')",sep ='')
  # 
  
  
  
  
  return(Code) 
}

########## Main GP Marmot function ##########
MARRMot <- function(o1,o2,o3,o4,o5,o6,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27){
  Parameters <-c(o1,o2,o3,o4,o5,o6,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27)
  Parameters <-unlist(lapply(Parameters, function(x) x<-max(0,x)))
  Parameters <-unlist(lapply(Parameters, function(x) x<-min(1,x)))
  
  Output_script_name <- as.character(round(runif(1,1,1e6)))
  Model_to_use <- Code_name(Parameters[1],Parameters[2],Parameters[3],Parameters[4],Parameters[5],Parameters[6])
  Matlab_code <- Code_generator(Model_to_use,Parameters[7:length(Parameters)],Output_script_name)
  
  run_matlab_code(Matlab_code)
  file_location <- paste(path,"/Matlab_fn",sep = '')
  Matlab_output <-read.csv(paste(file_location,"/Output_Q_",Output_script_name,".csv",sep = ''))
  unlink(paste(file_location,"/Output_Q_",Output_script_name,".csv",sep = ''))
  Matlab_output <- t(Matlab_output)
  sim <- Matlab_output[,1][paste("Q_",1:nrow(DataSet),sep = "")]
  return(sim)
}









#system('R CMD SHLIB CHMGP_SuperflexFunctions.f90', ignore.stdout = FALSE)
#dyn.load('CHMGP_SuperflexFunctions.so')

CHMGP = function(DependentVariable, IndependentVariables, DataSet, Weights = NULL,
                 FunctionSet = c('+', '-', '*', '/', 'sqrt','^'), ConstantRange = c(-1, 1), StrictFunctionsArgumentsRange = c(lower = 0, upper = 1),
                 PopulationSize = 500, NumberOfGenerations = 50, FitnessFunction = c('RMSE','RMSE','RMSE','RMSE','RMSE'), 
                 MaxDepthIni = 2, MaxDepthRun = 3, TournamentSize = 4, RoundingFactor = 6, 
                 VariationProbs = c(Pc = 0.7, PrSelected = 0.05, PmTree = 0.5, PmSeparation = 0.2, PmNode = 0.5, PmConstant = 0.7),
                 TerminalNodeProbs = 0.7, PunishOneNodeIndividuals = FALSE, SimpleOutput = TRUE, YacasSimplification = FALSE, 
                 DEoptimization = FALSE, DEiterMax = 200, DEfitness = 'SameAsGPfit',
                 DynamicPopulation = FALSE, CTtypesOccurence = FALSE, n_jobs=4
){
  
  if(DEoptimization) library('DEoptim')
  
  if(YacasSimplification){
    
    yacTest = system('which yacas', intern = TRUE)
    if(length(yacTest) == 0) stop('Yacas is not installed!')
    
  }
  library('parallel')
  library('nsga2R')
  library('fuse')
  
  P_E_funs = c("TANK", "MI", "MII", "MIII", "MIV", "MV", "MVI", "MVII", "MVIII", "MIX", "MX","MXI","MXII", "CT","FUSE","DCT","CT_new","MM_29")
  
  SetGlobalVariables(FunctionSet, IndependentVariables, ConstantRange, StrictFunctionsArgumentsRange, RoundingFactor, DataSet, DependentVariable, 
                     FitnessFunction, Weights, PunishOneNodeIndividuals, MaxDepthRun)
  
  constantMutationFactors = seq(2, 0.2, length = NumberOfGenerations)
  
  #######################
  #create individuals in first population
  population = FirstGeneration(MaxDepthIni, PopulationSize)
  
  #######################
  #Apply DDS algorithm to part of initial population
  DDS_indArray = array(c(NA), dim = c(64, 20))
  DDS_indArray[1,1] <- "CT_new"
  DDS_indArray[2,1] <- "WR"
  DDS_indArray[3,1] <- "IR"
  DDS_indArray[4,1] <- "RR"
  DDS_indArray[5,1] <- "UR"
  DDS_indArray[6,1] <- "FR"
  DDS_indArray[7,1] <- "SR"
  DDS_indArray[8,1] <- "CR"
  DDS_indArray[1,2:20] <-as.character(c(2:20))
  DDS_indArray[2,2:6] <- as.character(c(21:25))
  DDS_indArray[3,2:5] <- as.character(c(26:29))
  DDS_indArray[4,2:4] <- as.character(c(30:32))
  DDS_indArray[5,2:8] <- as.character(c(33:39))
  DDS_indArray[6,2:8] <- as.character(c(40:46))
  DDS_indArray[7,2:6] <- as.character(c(47:51))
  DDS_indArray[8,2:14] <- as.character(c(52:64))
  DDS_indArray[9:64,1] <- as.character(c(1:56))
  
  DDS_IND <- list(IndArray = DDS_indArray, Equation = NA, Fitness = c(rep(NA,length(FitnessFunction))), DriveVec = c(rep("Function",8),rep("Terminal",56)),
                  RestrictedRows = c(2:64), IndLength = 64, Changed = TRUE, Front = NA, Crowding_distance = NA,Sim =c())
  
  no_DDS_IND <- n_jobs*2 
  DDS_IND_list <- list()
  for (i in 1:no_DDS_IND){
    DDS_IND_list[[i]] <-runif(56,0,1)
  }
  
  Sep_DDS_IND_list_1 <- mclapply(as.list(DDS_IND_list),DDS,mc.cores = n_jobs)
  #Sep_DDS_IND_list_1 <- lapply(DDS_IND_list,DDS)
  for (i in 1:no_DDS_IND){
    DDS_IND[['IndArray']][9:64,1] <-as.character(Sep_DDS_IND_list_1[[i]])
    DDS_IND_list[[i]] <- DDS_IND
  }
  
  population[1:no_DDS_IND] <- DDS_IND_list
  ###################
  #Remove unsuitable CT_new individuals
  
  for (r in 1:PopulationSize) {
    if ((population[[r]][['IndArray']][1,1]=="CT_new") && (population[[r]][['IndLength']] <= 20)){
      population[[r]] <- CreateIndividual(2,"full",F)
    } else {
      population[[r]]<- population[[r]]
    }
  }
  
  ####################
  #calculate fitness for 1st generation
  #population = lapply(population, ComputeIndividual)
  population <- mclapply(as.list(population),ComputeIndividual, mc.cores = n_jobs)
  #population <- lapply(population,ComputeIndividual)
  names(population) = paste0('I', 1:PopulationSize)
  
  ###################
  #need to calculate front and crowding distances
  # calculating front 
  front_matrix <-t(sapply(population, function(x) x[['Fitness']]))
  
  front_list <-fastNonDominatedSorting(front_matrix)
  
  for (i in 1: length(front_list)){
    for (k in 1:length(front_list[[i]])){
      population[[front_list[[i]][k]]][['Front']] <- i 
    }
  }
  ##################################
  # calculating crowding distance
  
  rnkIndex <- integer(PopulationSize)
     i <- 1
     while (i <= length(front_list)) {
        rnkIndex[front_list[[i]]] <- i
        i <- i + 1
      }
  front_matrix <- cbind(front_matrix,rnkIndex)
  
  
  objrange <- apply(front_matrix[,1:length(FitnessFunction),drop=F],2,max) - apply(front_matrix[,1:length(FitnessFunction),drop=F],2,min)
  cr_dis <- crowdingDist4frnt(front_matrix,front_list,objrange)
  front_matrix <-cbind(front_matrix,apply(cr_dis,1,sum))
  
  for (i in 1:PopulationSize){
    population[[i]][['Crowding_distance']] <- sum(cr_dis[i,])
  }
  
  #####################
  isFitMulti = FitnessFunction %in% MultiObjectiveFitnesses
  if(isFitMulti){
    
    generationFitnesses = MultiObjective(population)
    
  } else{
    
    generationFitnesses = sapply(population, function(x) x[['Fitness']])
    generationFitnesses = generationFitnesses[1,]
  }
  
  ##########################
  #printing & saving results 
  
  runResults = vector("list",PopulationSize)
  printGenerationResults = c(0, floor(seq(0.1, 0.9, 0.1) * NumberOfGenerations), NumberOfGenerations)
  cat('Generation:', printGenerationResults[1], 'Completed,  ')

  ########################
  
  tenthDEopt = FALSE
  if(tenthDEopt){
    
    tenthSeq = seq(10, NumberOfGenerations, 10) + 1
    tenthDEoptResults = data.frame(Generation = tenthSeq, Fit = numeric(length(tenthSeq)), OptModel = character(length(tenthSeq)), stringsAsFactors = FALSE)
    
  }
  
  if(CTtypesOccurence) {
    
    # List for storage of equations from each generation
    eqStorage = vector(mode = 'list', NumberOfGenerations + 1)
    eqStorage = lapply(eqStorage, function(x){character(PopulationSize)})
    
    eqStorage[[1]] = sapply(population, '[[', 'Equation')
    
  }
  ################# 
  #main cycle
  
  for(q in 1:NumberOfGenerations){
    
    goodFitnesses = generationFitnesses[!is.na(generationFitnesses)]
    
    PopulationCheck(goodFitnesses, TournamentSize, q)
    
    newPopulation = vector('list', PopulationSize)
    
    # creating a mating pool
    front_matrix <- cbind(1:PopulationSize,front_matrix)
    matingPool <- tournamentSelection(front_matrix,PopulationSize,TournamentSize)
    
    noOfNewIndividuals = min(sum(rbinom(PopulationSize, 1, VariationProbs[['Pc']])), PopulationSize)
    noOfNewIndividuals = ceiling(noOfNewIndividuals/n_jobs)*n_jobs
    
    if(noOfNewIndividuals > 0){
      parents = matrix(NA, nrow = noOfNewIndividuals, ncol = 2)
      colnames(parents) = c('Parent1', 'Parent2')
      for (i in 1:noOfNewIndividuals){
        ran_ind_cr <- round(runif(2,1,PopulationSize))
        parents[i,1] <- matingPool[ran_ind_cr[1],][1]
        parents[i,2] <- matingPool[ran_ind_cr[2],][1]
      }
      
      positionRange = 1:noOfNewIndividuals
      newPopulation[positionRange] = Crossover(population[parents[,1]], population[parents[,2]], TerminalNodeProbs)
      
    }
    
    noOfMutants = PopulationSize - noOfNewIndividuals
    Mutating_ind <- c(rep(NA,noOfMutants))
    for (i in 1:noOfMutants){
      ran_ind_mu <- round(runif(1,1,PopulationSize))
      Mutating_ind[i] <- matingPool[ran_ind_mu,1]
    }
    if(noOfNewIndividuals < PopulationSize){
      
      rangeOfMutatingPopulation = (noOfNewIndividuals + 1):PopulationSize
      newPopulation[rangeOfMutatingPopulation] = Mutation(population[Mutating_ind], VariationProbs, constantMutationFactors[q])  
    }  
    
    ####################
    
    accPopulation <- c(population,newPopulation)
    
    names(accPopulation) = paste0('I', 1:length(accPopulation))
    
    # for (i in 1:length(accPopulation)){
    #   cat(i,'-',length(accPopulation[[i]][['Fitness']]),',')
    # }
    
    #creating front matrix
    front_matrix <-t(sapply(accPopulation, function(x) x[['Fitness']]))
    
    front_list <-fastNonDominatedSorting(front_matrix)
    for (i in 1: length(front_list)){
      for (k in 1:length(front_list[[i]])){
        accPopulation[[front_list[[i]][k]]][['Front']] <- i 
      }
    }
    
    #crowding_distance
    rnkIndex <- integer(PopulationSize*2)
    i <- 1
    while (i <= length(front_list)) {
      rnkIndex[front_list[[i]]] <- i
      i <- i + 1
    }
    front_matrix <- cbind(front_matrix,rnkIndex)
    
    
    objrange <- apply(front_matrix[,1:length(FitnessFunction),drop=F],2,max) - apply(front_matrix[,1:length(FitnessFunction),drop=F],2,min)
    cr_dis <- crowdingDist4frnt(front_matrix,front_list,objrange)
    front_matrix <-cbind(front_matrix,apply(cr_dis,1,sum))
    
    for (i in seq_along(accPopulation)){
      accPopulation[[i]][['Crowding_distance']] <- sum(cr_dis[i,])
    }
    
    #############
    #NSGA selection
    
    m <- 1
    Ind_cunt <- 0
    if (length(front_list[[1]])>= PopulationSize){
      margin_front <- 1
    } else {
      while(Ind_cunt +length(front_list[[m]]) < PopulationSize){
        Ind_cunt <- Ind_cunt + length(front_list[[m]])
        m <- m + 1
        margin_front <- m
      }
    }
    
    nextPopulation = list()
    
    if (margin_front ==1){
      margin_list <- accPopulation[front_list[[margin_front]]]
      crw_dis <- c(rep(0,length(front_list[[margin_front]])))
      for (i in 1:length(front_list[[margin_front]])){
        crw_dis[i] <- margin_list[[i]][['Crowding_distance']]
      }
      remain_ind <- data.frame(index = c(1:length(front_list[[margin_front]])),crw_dis = crw_dis)
      sort_remain_ind <- remain_ind[order(remain_ind[['crw_dis']],decreasing = T),]
      nextPopulation <- accPopulation[sort_remain_ind[1:PopulationSize,1]]
      population <- nextPopulation
    } else {
      len_front <- 0
      for (i in 1:(margin_front-1)){
        nextPopulation <- c(nextPopulation,accPopulation[front_list[[i]]])
        len_front <- len_front + length(front_list[[i]])
      }
      if (length(nextPopulation) == PopulationSize){
        population = nextPopulation
      } else {
        margin_list <- accPopulation[front_list[[margin_front]]]
        crw_dis <- c(rep(0,length(front_list[[margin_front]])))
        for (i in 1:length(front_list[[margin_front]])){
          crw_dis[i] <- margin_list[[i]][['Crowding_distance']]
        }
        
        remain_ind <- data.frame(index = c(1:length(front_list[[margin_front]])),crw_dis = crw_dis)
        sort_remain_ind <- remain_ind[order(remain_ind[['crw_dis']],decreasing = T),]
        nextPopulation <- c(nextPopulation,accPopulation[sort_remain_ind[1:(PopulationSize- Ind_cunt),1]])
        population <- nextPopulation
      }
    }
    
    names(population) = paste0('I', 1:PopulationSize)
    
    # need to calculate front and crowding distances
    # calculating front
    front_matrix <-t(sapply(population, function(x) x[['Fitness']]))
    front_list <-fastNonDominatedSorting(front_matrix)
    for (i in 1: length(front_list)){
      for (k in 1:length(front_list[[i]])){
        population[[front_list[[i]][k]]][['Front']] <- i 
      }
    }
    
    #################################
    # calculating crowding distance
    rnkIndex <- integer(PopulationSize)
    i <- 1
    while (i <= length(front_list)) {
      rnkIndex[front_list[[i]]] <- i
      i <- i + 1
    }
    front_matrix <- cbind(front_matrix,rnkIndex)
    
    
    objrange <- apply(front_matrix[,1:length(FitnessFunction),drop=F],2,max) - apply(front_matrix[,1:length(FitnessFunction),drop=F],2,min)
    #objrange <-c((max(front_matrix[,1])-min(front_matrix[,1])),(max(front_matrix[,2])-min(front_matrix[,2])),(max(front_matrix[,3])-min(front_matrix[,3])))
    cr_dis <- crowdingDist4frnt(front_matrix,front_list,objrange)
    front_matrix <-cbind(front_matrix,apply(cr_dis,1,sum))
    
    for (i in 1:PopulationSize){
      population[[i]][['Crowding_distance']] <- sum(cr_dis[i,])
    }
    
    ########################
    if(isFitMulti){
      
      generationFitnesses = MultiObjective(population)
      
    } else{
      
      generationFitnesses = sapply(population, function(x) x[['Fitness']])
      
    }
    
    runResults = SaveRunResults(population, front_list, runResults, q, NumberOfGenerations)
    if (q == NumberOfGenerations){
      runResults <- runResults[1:length(front_list[[1]])]
    } else {
      runResults <- runResults
    }
    if(any(q == printGenerationResults)) cat('Generation:', q, 'Completed,  ')    
  } 
  runResults
}








setwd('../')

# Import relevant libraries 
library('parallel')
library('nsga2R')
library('fuse')
library('matlabr')


###############
# Setting CHMGP parameters

IndependentRuns <- 10    #Running the entire GP process multiple times from scratch, 
#each time starting with a different random population.
PopulationSize <- 500
NumberOfGenerations <- 50   #higher number of generations allows the population to evolve for a longer time
n_jobs <- 1  #24
specFunct = c('MARRMot')
criterion = c('VE0','KG10','NS0','logNS0')

constantMutationFactors = seq(2, 0.2, length = NumberOfGenerations) # strength of mutation: Early generations: Mutations are stronger (value 2).
# Later generations: Mutations are weaker (value 0.2).
VariationProbs = c(Pc = 0.7, PrSelected = 0.05, PmTree = 0.5, PmSeparation = 0.3, PmNode = 0.3, PmConstant = 0.7)

# Pc = 0.7 (Crossover 70%) Like parents creating offspring. Two models mix their "genes" to create a new model.
# PrSelected = 0.05 (Selection 5%) Chance that a model is chosen to pass its genes.
# PmTree = 0.5 (Big Mutation 50%) 50% chance of changing the whole structure of a model.
# PmSeparation = 0.3 (Small Mutation 30%) A smaller subtree gets changed instead of the whole structure.
# PmNode = 0.3 (Node Mutation 30%) A single node (small part of the model) is changed.
# PmConstant = 0.7 (Constant Mutation 70%) If a model has fixed values (like numbers in formulas), they have a 70% chance to change.

TournamentSize <- 4  #We need 500 winners to form/create the next generation; Each tournament picks 4 random equations (TournamentSize = 4).
# We run 500 tournaments (one per needed equation); same equations can win multiple times, unless replace = TRUE to prevent dupllciate

TerminalNodeProbs <- 0.9  ##  in 90% of cases, the algorithm will choose a "terminal node" (e.g. Constants (e.g., 3,5) or Variables (x, Q)) instead of a function/operator (+/-, sin())when building or mutating an equation.
dataNameTraining = 'luxHue_train'   #to train the model--- here string used, to dynamically change name
toSearch = 'Q'       # 'Q' column of dataframe 'luxHue_train' 
depVari = paste0(toSearch)
dataSetTrain = get(dataNameTraining)   #getting the actual DF from the string

setwd('./Matlab_fn')
write.csv(dataSetTrain,"Input_Variables.csv")
setwd('../')

Date=as.Date(dataSetTrain[['Date']])  #Date=strptime(dataSetTrain[['Date']],format = '%d/%m/%Y %H:%M')
DependentVariable = depVari
IndependentVariables = c('P','E','T')
FunctionSet = c(specFunct, '+', '-', '*', '/')
ConstantRange = c(0, 1)
StrictFunctionsArgumentsRange = c(lower = 0, upper = 1)
FitnessFunction = criterion
MaxDepthIni <- 1
MaxDepthRun = 3
RoundingFactor = 3
results = vector(mode = 'list', IndependentRuns)
namesRes = paste0('run_', 1:IndependentRuns)
names(results) = namesRes
DataSet <- dataSetTrain
Weights = NULL
PunishOneNodeIndividuals = FALSE
SimpleOutput = TRUE
YacasSimplification = FALSE
DEoptimization = FALSE
DEiterMax = 200
DEfitness = 'SameAsGPfit'
DynamicPopulation = FALSE
CTtypesOccurence = FALSE

# Link Global variables
SetGlobalVariables(FunctionSet, IndependentVariables, ConstantRange, StrictFunctionsArgumentsRange, RoundingFactor, DataSet, DependentVariable, 
                   FitnessFunction, Weights, PunishOneNodeIndividuals, MaxDepthRun)

################
# Saving console content

sink("Blackwater_Marrmot_Lumped_1", append=FALSE, split=FALSE)

################
# Running several Independent runs

ptm = proc.time()

for(p in namesRes){
  
  print(p)
  
  ################
  # Initialization
  
  population = FirstGeneration(MaxDepthIni, PopulationSize)
  
  # Ind_index <- 1:length(population)
  # Marrmot_Ind_index <-Ind_index[sapply(population, function(x) x[['IndArray']][1,1]=="MARRMot")]
  # Other_Ind_index <- Ind_index[sapply(population, function(x) x[['IndArray']][1,1]!="MARRMot")]
  # pre_index_order <-c(Marrmot_Ind_index,Other_Ind_index)
  
  Marrmot_Ind <-population[sapply(population, function(x) x[['IndArray']][1,1]=="MARRMot")]
  
  Other_Ind <-population[sapply(population, function(x) x[['IndArray']][1,1]!="MARRMot")]
  
  # if(length(Other_Ind)>0){
  #   for(vv in 1:length(Other_Ind)){
  #     Other_Ind[[vv]]<- ComputeIndividual(Other_Ind[[vv]])
  #   }
  # }
  
  if(length(Other_Ind)>0){
    Other_Ind <- mclapply(as.list(Other_Ind),ComputeIndividual, mc.cores = n_jobs)
  }
  
  Output_script_vec <- as.character(1:length(Marrmot_Ind))
  
  MakeEquation_new = function(individual){
    ReservoirSet <-c("WR","IR","RR","UR","FR","SR","CR")
    IA = individual[['IndArray']]
    if(individual[['IndLength']] == 1){
      return(paste0('y = ', IA[, 1]))
    }
    eq_vec = IA[, 1]
    for (i in individual[['IndLength']]:1){
      if( any(FunctionSet == IA[i, 1]) | any(ReservoirSet == IA[i, 1])){
        for(j in 2:ncol(IA)){
          if( !is.na(IA[i,j]) ){
            arg_row = as.numeric(IA[i,j])
            if( j == 2 && is.na(IA[i, (j + 1)]) ){
              if(!any(IA[i, 1] == SpecialFunctions)){
                eq_vec[i] = paste0(eq_vec[i], '(', eq_vec[arg_row])
              }else{
                eq_vec[i] = paste0(eq_vec[i], eq_vec[arg_row])
              }
            }
            if( j == 2 && !is.na(IA[i, (j + 1)]) ){
              if( any(PrefixFunctions == IA[i, 1]) ){
                first = paste0(eq_vec[i],'(', eq_vec[arg_row], ',')
              } else{
                first = paste0('(', eq_vec[arg_row])
              }
            }
            if(j > 2){
              if( any(PrefixFunctions == IA[i, 1]) ){
                eq_vec[i] = paste0(first, eq_vec[arg_row])
                if(j >= 3) first = paste0(eq_vec[i], ',')
              } else{
                eq_vec[i] = paste0(first, eq_vec[i], eq_vec[arg_row])
              }
            }
          }
        }
        eq_vec[i] = paste0(eq_vec[i], ')')
      }
    } 
    return(eq_vec)
  }
  
  Matlab_code <- c()
  
  if(length(Marrmot_Ind)>0){
    
    for(gg in 1:length(Marrmot_Ind)){
      for(ml in 1:nrow(Marrmot_Ind[[gg]][['IndArray']])){
        if(Marrmot_Ind[[gg]][['IndArray']][ml,1] %in% c("P","E","T")){
          Marrmot_Ind[[gg]][['IndArray']][ml,1] <- runif(1)
        }
      }
      M_Eq <- MakeEquation_new(Marrmot_Ind[[gg]])
      Parameters <- c()
      for(hh in 1:33){
        Parameters[hh]<-max(min(eval(parse(text=noquote(M_Eq[[as.numeric(Marrmot_Ind[[gg]][['IndArray']][1,(hh+1)])]]))),1),0)
      }
      Output_script_name <- Output_script_vec[gg]
      Marrmot_Ind[[gg]][['MARRMot_Ind']]<- Output_script_name
      Model_to_use <- Code_name(Parameters[1],Parameters[2],Parameters[3],Parameters[4],Parameters[5],Parameters[6])
      Matlab_code <- c(Matlab_code,Code_generator(Model_to_use,Parameters[7:length(Parameters)],Output_script_name))
    }
    
    TO_MAT <- array(dim = c(length(Marrmot_Ind),1))
    TO_MAT[,1]<-Matlab_code
    write.csv(TO_MAT,paste(path,"/Matlab_fn/TO_MAT.csv",sep = ''))
    
    MAT_CODE <-c(paste("cd '",path,"/Matlab_fn","';",sep = ""),"A=readcell('TO_MAT.csv')","B=string();",paste("for i =2:",(length(Marrmot_Ind)+1),";B(i-1)=string(A(i,2));end",sep = ""),paste("parfor i=1:",length(Marrmot_Ind),";CommitEval(B(1,i));end",sep = ""))
    
    run_matlab_code(MAT_CODE)
    file_location <- paste(path,"/Matlab_fn",sep = '')
    
    Read_Marrmot <- function(Individual){
      Matlab_output <-read.csv(paste(file_location,"/Output_Q_",Individual[['MARRMot_Ind']],".csv",sep = ''))
      Matlab_output <- t(Matlab_output)
      sim <- Matlab_output[,1][paste("Q_",1:nrow(DataSet),sep = "")]
      obs <- DataSet$Q
      for(mn in 1:length(FitnessFunction)){
        Individual$Fitness[mn] <- FitnessComputation(sim,obs,FitnessFunction[mn])
      }
      Individual$Sim <- sim
      Individual$Equation <- MakeEquation(Individual)
      Individual <-Individual[1:(length(Individual)-1)]
      
      Individual
      
    }
    
    Marrmot_Ind <- mclapply(as.list(Marrmot_Ind),Read_Marrmot, mc.cores = n_jobs)
    
    # for(x in 1:length(Marrmot_Ind)){
    #   Matlab_output <-read.csv(paste(file_location,"/Output_Q_",Output_script_vec[x],".csv",sep = ''))
    #   #unlink(paste(file_location,"/Output_Q_",Output_script_vec[x],".csv",sep = ''))
    #   Matlab_output <- t(Matlab_output)
    #   sim <- Matlab_output[,1][paste("Q_",1:nrow(DataSet),sep = "")]
    #   obs <- DataSet$Q
    #   
    #   for(mn in 1:length(FitnessFunction)){
    #     Marrmot_Ind[[x]]$Fitness[mn] <- FitnessComputation(sim,obs,FitnessFunction[mn])
    #   }
    #   Marrmot_Ind[[x]]$Sim <- sim
    #   Marrmot_Ind[[x]]$Equation <- MakeEquation(Marrmot_Ind[[x]])
    # }
    
  }
  
  population <-c(Marrmot_Ind,Other_Ind)
  
  # Index_dataframe <-data.frame(pre_index_order,1:length(population))
  # names(Index_dataframe)<-c("pre","post")
  # Index_dataframe <-Index_dataframe[order(Index_dataframe$pre),]
  # population <-population[Index_dataframe$post]
  
  #population <- mclapply(as.list(population),ComputeIndividual, mc.cores = n_jobs)
  
  names(population) = paste0('I', 1:PopulationSize)
  
  # calculating fronts
  front_matrix <-t(sapply(population, function(x) x[['Fitness']]))
  
  front_list <-fastNonDominatedSorting(front_matrix)
  for (i in 1: length(front_list)){
    for (k in 1:length(front_list[[i]])){
      population[[front_list[[i]][k]]][['Front']] <- i
    }
  }
  
  # calculating crowding distance
  rnkIndex <- integer(PopulationSize)
  i <- 1
  while (i <= length(front_list)) {
    rnkIndex[front_list[[i]]] <- i
    i <- i + 1
  }
  front_matrix <- cbind(front_matrix,rnkIndex)
  
  objrange <- apply(front_matrix[,1:length(FitnessFunction),drop=F],2,max) - apply(front_matrix[,1:length(FitnessFunction),drop=F],2,min)
  cr_dis <- crowdingDist4frnt(front_matrix,front_list,objrange)
  front_matrix <-cbind(front_matrix,apply(cr_dis,1,sum))
  
  for (i in 1:PopulationSize){
    population[[i]][['Crowding_distance']] <- sum(cr_dis[i,])
  }
  
  printGenerationResults = c(0, floor(seq(0.1, 0.9, 0.1) * NumberOfGenerations), NumberOfGenerations)
  cat('Generation:', printGenerationResults[1], 'Completed,  ')
  
  ##################
  # Main cycle
  
  for (q in 1:NumberOfGenerations){
    
    newPopulation = vector('list', PopulationSize)
    
    # creating a mating pool
    front_matrix <- cbind(1:PopulationSize,front_matrix)
    matingPool <- tournamentSelection(front_matrix,PopulationSize,TournamentSize)
    
    noOfNewIndividuals = min(sum(rbinom(PopulationSize, 1, VariationProbs[['Pc']])), PopulationSize)
    noOfNewIndividuals = ceiling(noOfNewIndividuals/n_jobs)*n_jobs
    
    # Crossover
    if(noOfNewIndividuals > 0){
      parents = matrix(NA, nrow = noOfNewIndividuals, ncol = 2)
      colnames(parents) = c('Parent1', 'Parent2')
      for (i in 1:noOfNewIndividuals){
        ran_ind_cr <- round(runif(2,1,PopulationSize))
        parents[i,1] <- matingPool[ran_ind_cr[1],][1]
        parents[i,2] <- matingPool[ran_ind_cr[2],][1]
      }
      
      positionRange = 1:noOfNewIndividuals
      newPopulation[positionRange] = Crossover(population[parents[,1]], population[parents[,2]], TerminalNodeProbs)
      
    }
    
    # Mutation
    noOfMutants = PopulationSize - noOfNewIndividuals
    Mutating_ind <- c(rep(NA,noOfMutants))
    for (i in 1:noOfMutants){
      ran_ind_mu <- round(runif(1,1,PopulationSize))
      Mutating_ind[i] <- matingPool[ran_ind_mu,1]
    }
    if(noOfNewIndividuals < PopulationSize){
      
      rangeOfMutatingPopulation = (noOfNewIndividuals + 1):PopulationSize
      newPopulation[rangeOfMutatingPopulation] = Mutation(population[Mutating_ind], VariationProbs, constantMutationFactors[q])  
    }
    
    # to avoid subscript out of bounds problem in some generations, let's try following code
    for (i in 1:length(newPopulation)){
      if (length(newPopulation[[i]]) <10){
        newPopulation[[i]] <- population[[round(runif(1,1,PopulationSize))]]
      } else if ((length(newPopulation[[i]]) == 10) & (length(newPopulation[[i]][['Sim']]) ==1)){
        newPopulation[[i]] <- population[[round(runif(1,1,PopulationSize))]]
      } else {
        newPopulation[[i]] <- newPopulation[[i]]
      }
    }
    
    # Creating Parent + Child Population
    accPopulation <- c(population,newPopulation)
    
    names(accPopulation) = paste0('I', 1:length(accPopulation))
    
    # Creating front matrix for Parent + Child Population
    front_matrix <-t(sapply(accPopulation, function(x) x[['Fitness']]))
    
    front_list <-fastNonDominatedSorting(front_matrix)
    
    for (i in 1: length(front_list)){
      for (k in 1:length(front_list[[i]])){
        accPopulation[[front_list[[i]][k]]][['Front']] <- i
      }
    }
    
    # Crowding_distance for Parent + Child Population
    rnkIndex <- integer(PopulationSize*2)
    i <- 1
    while (i <= length(front_list)) {
      rnkIndex[front_list[[i]]] <- i
      i <- i + 1
    }
    front_matrix <- cbind(front_matrix,rnkIndex)
    
    objrange <- apply(front_matrix[,1:length(FitnessFunction),drop=F],2,max) - apply(front_matrix[,1:length(FitnessFunction),drop=F],2,min)
    cr_dis <- crowdingDist4frnt(front_matrix,front_list,objrange)
    front_matrix <-cbind(front_matrix,apply(cr_dis,1,sum))
    
    for (i in seq_along(accPopulation)){
      accPopulation[[i]][['Crowding_distance']] <- sum(cr_dis[i,])
    }
    
    # NSGA Selection
    m <- 1
    Ind_cunt <- 0
    if (length(front_list[[1]])>= PopulationSize){
      margin_front <- 1
    } else {
      while(Ind_cunt +length(front_list[[m]]) < PopulationSize){
        Ind_cunt <- Ind_cunt + length(front_list[[m]])
        m <- m + 1
        margin_front <- m
      }
    }
    
    nextPopulation = list()
    
    if (margin_front ==1){
      margin_list <- accPopulation[front_list[[margin_front]]]
      crw_dis <- c(rep(0,length(front_list[[margin_front]])))
      for (i in 1:length(front_list[[margin_front]])){
        crw_dis[i] <- margin_list[[i]][['Crowding_distance']]
      }
      remain_ind <- data.frame(index = c(1:length(front_list[[margin_front]])),crw_dis = crw_dis)
      sort_remain_ind <- remain_ind[order(remain_ind[['crw_dis']],decreasing = T),]
      nextPopulation <- accPopulation[sort_remain_ind[1:PopulationSize,1]]
      population <- nextPopulation
    } else {
      len_front <- 0
      for (i in 1:(margin_front-1)){
        nextPopulation <- c(nextPopulation,accPopulation[front_list[[i]]])
        len_front <- len_front + length(front_list[[i]])
      }
      if (length(nextPopulation) == PopulationSize){
        population = nextPopulation
      } else {
        margin_list <- accPopulation[front_list[[margin_front]]]
        crw_dis <- c(rep(0,length(front_list[[margin_front]])))
        for (i in 1:length(front_list[[margin_front]])){
          crw_dis[i] <- margin_list[[i]][['Crowding_distance']]
        }
        
        remain_ind <- data.frame(index = c(1:length(front_list[[margin_front]])),crw_dis = crw_dis)
        sort_remain_ind <- remain_ind[order(remain_ind[['crw_dis']],decreasing = T),]
        nextPopulation <- c(nextPopulation,accPopulation[sort_remain_ind[1:(PopulationSize- Ind_cunt),1]])
        population <- nextPopulation
      }
    }
    
    names(population) = paste0('I', 1:PopulationSize)
    
    # Calculating front for selected population
    front_matrix <-t(sapply(population, function(x) x[['Fitness']]))
    front_list <-fastNonDominatedSorting(front_matrix)
    
    for (i in 1: length(front_list)){
      for (k in 1:length(front_list[[i]])){
        population[[front_list[[i]][k]]][['Front']] <- i
      }
    }
    
    # Calculating crowding distance for selected population
    rnkIndex <- integer(PopulationSize)
    i <- 1
    while (i <= length(front_list)) {
      rnkIndex[front_list[[i]]] <- i
      i <- i + 1
    }
    front_matrix <- cbind(front_matrix,rnkIndex)
    
    objrange <- apply(front_matrix[,1:length(FitnessFunction),drop=F],2,max) - apply(front_matrix[,1:length(FitnessFunction),drop=F],2,min)
    cr_dis <- crowdingDist4frnt(front_matrix,front_list,objrange)
    front_matrix <-cbind(front_matrix,apply(cr_dis,1,sum))
    
    for (i in 1:PopulationSize){
      population[[i]][['Crowding_distance']] <- sum(cr_dis[i,])
    }
    if(any(q == printGenerationResults)) cat('Generation:', q, 'Completed,  ')
  }
  
  # Running Time
  etmh = (proc.time() - ptm)['elapsed']/3600
  etmd = (proc.time() - ptm)['elapsed']/86400
  message(paste(p, 'Elapsed time:', etmh, 'hours, (ie.',etmd, 'days)'))
  
  cat('Elapsed time:',p, etmh,'hours')
  
  ################
  # Saving Results
  
  runResults = vector("list",PopulationSize)
  runResults = SaveRunResults(population, front_list, runResults, q, NumberOfGenerations)
  runResults <- runResults[front_list[[1]]]
  results[[p]] <- runResults
  saveRDS(results, paste0('Backup_',dataNameTraining, '_Blackwater_Marrmot_Lumped_1.rds'))
}
saveRDS(results, paste0(dataNameTraining, '_Blackwater_Marrmot_Lumped_1.rds'))
sink()