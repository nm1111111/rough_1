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






