#alfalfa
#install.packages("dplyr")
library(univariateML)
library(VineCopula)

library(dplyr)
library(copula)
library(copBasic)

geoid_stctid = read.csv("/home/uallxr001/clim_shock/geoid_stctid.csv")
spi_fin = read.csv("/home/uallxr001/clim_shock/growing_per_spi_alfalfa.csv",check.names = FALSE)       #change
colnames(spi_fin)[1] <- "GEOID"

spi_fin_stct = merge(geoid_stctid, spi_fin, by = "GEOID")
spi_fin_stct = subset(spi_fin_stct, select = -GEOID)
spi_fin_stct = spi_fin_stct %>%
  filter(!duplicated(STCTID))
rownames(spi_fin_stct) <- spi_fin_stct$STCTID
spi_fin_stct$STCTID <-NULL
colnames(spi_fin_stct)<- c(1961:2019)


alfalfa_yld_ct_std_pivot = read.csv("/home/uallxr001/clim_shock/alfalfa_ct_std_yld_pivot_full.csv", check.names = FALSE, row.names='STCTID')

alfalfa_yld_ct_std_pivot1 = alfalfa_yld_ct_std_pivot[, which(as.numeric(colnames(alfalfa_yld_ct_std_pivot))>=1971 & as.numeric(colnames(alfalfa_yld_ct_std_pivot))<=2000)]
alfalfa_yld_ct_std_pivot2 = alfalfa_yld_ct_std_pivot[,  which(as.numeric(colnames(alfalfa_yld_ct_std_pivot))>=1991 & as.numeric(colnames(alfalfa_yld_ct_std_pivot))<=2020)]


spi_fin1 = t(read.csv("/home/uallxr001/clim_shock/growing_per_spi_new_prec_alfalfa.csv",check.names = FALSE, row.names = "Date")  )     #change
spi_fin1 = spi_fin1[order(row.names(spi_fin1)), ]
spi_fin_stct1 = spi_fin1[, which(as.numeric(colnames(spi_fin1))>=1971 & as.numeric(colnames(spi_fin1))<=2000)]        #change 
spi_fin_stct2 = spi_fin1[, which(as.numeric(colnames(spi_fin1))>=1991 & as.numeric(colnames(spi_fin1))<=2020)]       #change 
spi_fin_stct_tot = spi_fin1[, which(as.numeric(colnames(spi_fin1))>=1971 & as.numeric(colnames(spi_fin1))<=2020)]

calculate_mse <- function(true, simulated, n, k) {
  mse <- sum((true - simulated)^2)
  mse <- mse/(n-k)   #dividing by n-k to get the unbiased estimates of squared dif
  return(mse)
}

calculate_aic <- function(mse, n, k) {
  aic = 2*k + n*log(mse) + 2*k*(k+1)/(n-k-1)   #Poulomi HESS
  return(aic)
}

copula_main <- function(yld_pivot, spi_pivot, col_nam1 = "SPI", col_nam2 = "YLD"){
  copula_ct_all_list = data.frame()
  cat("STCTID", "Family", "Par", "Par2", "Avg Yld","cop08", "cop99", "cop130", "cop149", "cop199", 
      "ProbDr08", "ProbDr99", "ProbDr130", "ProbDr149", "ProbDr199",
      "cdf08", "cdf99", "cdf130", "cdf149", "cdf199", "cdf_Yldavg",'tau',
      'p-val',"pval_cvm", "pval_ks", "stat_cvm", "stat_ks", "len")
  
  for(i in 1:nrow(yld_pivot)){
    temp_stctid = row.names(yld_pivot)[i]
    
    tryCatch({
      temp_x1 =  spi_pivot[temp_stctid,, drop = FALSE]
      temp_x2 = yld_pivot[temp_stctid,]
      m = rbind(temp_x1, temp_x2)
      row.names(m) = c("SPI", "YLD")
      m = na.omit(t(m))
      
      bestdist_spi = model_select(na.omit(t(temp_x1)), models = univariateML_models)
      bestdist_yld = model_select(na.omit(t(temp_x2)), models = univariateML_models)
      
      set.seed(1234)
      temp1 = 1:5
      u = pobs(m)[,1]
      v = pobs(m)[,2]
      
      cop = BiCopSelect(u, v, familyset = temp1, indeptest = TRUE, level =0.1)
      family = cop$family
      
      mean_yld = mean(as.matrix(alfalfa_yld_ct_std_pivot[temp_stctid, which(as.numeric(colnames(alfalfa_yld_ct_std_pivot))>=1971 & as.numeric(colnames(alfalfa_yld_ct_std_pivot))<=2020)]), na.rm = TRUE)
      
      
      if(family>0){
        gof_pars = BiCopGofTest(u,v, cop,method = "kendall", B = 500)
        pval_cvm  =gof_pars$p.value.CvM
        pval_ks = gof_pars$p.value.KS
        stat_cvm = gof_pars$statistic.CvM
        stat_ks = gof_pars$statistic.KS
      }
      if(family==0){

        cop_gauss = BiCopSelect(u, v, familyset = 1, indeptest = FALSE)
        cop_t = BiCopSelect(u, v, familyset = 2, indeptest = FALSE)
        cop_clay = BiCopSelect(u, v, familyset = 3, indeptest = FALSE)
        cop_gumb = BiCopSelect(u, v, familyset = 4, indeptest = FALSE)
        cop_fran = BiCopSelect(u, v, familyset = 5, indeptest = FALSE)
        
        plack_par = summary(fitCopula(copula = plackettCopula(), data = cbind(u,v)))$coefficient[1]
        plack_cop = plackettCopula(plack_par)    
        
        cdf_emp = C.n(pobs(m),X = pobs(m))
        cdf_ind = unname(u*v)
        cdf_plack = unname(pCopula(pobs(m), plack_cop))
        cdf_gauss = BiCopCDF(u,v,cop_gauss)
        cdf_t = BiCopCDF(u,v,cop_t)
        cdf_clay = BiCopCDF(u,v,cop_clay)
        cdf_gumb = BiCopCDF(u,v,cop_gumb)
        cdf_fran = BiCopCDF(u,v,cop_fran)
        
        mse_ind = calculate_mse(true = cdf_emp, simulated = cdf_ind, n = nrow(m), k = 0)
        mse_plack = calculate_mse(true = cdf_emp, simulated = cdf_ind, n = nrow(m), k = 1)
        mse_gauss = calculate_mse(true = cdf_emp, simulated = cdf_gauss, n = nrow(m), k = 1)
        mse_t = calculate_mse(true = cdf_emp, simulated = cdf_t, n = nrow(m), k = 2)
        mse_clay = calculate_mse(true = cdf_emp, simulated = cdf_clay, n = nrow(m), k = 1)
        mse_gumb = calculate_mse(true = cdf_emp, simulated = cdf_gumb, n = nrow(m), k = 1)
        mse_fran = calculate_mse(true = cdf_emp, simulated = cdf_fran, n = nrow(m), k = 1)
        
        aic_ind = calculate_aic(mse_ind, n = nrow(m), k = 0)
        aic_plack = calculate_aic(mse_plack, n = nrow(m), k = 1)
        aic_gauss = calculate_aic(mse_gauss, n = nrow(m), k = 1)
        aic_t = calculate_aic(mse_t, n = nrow(m), k = 2)
        aic_clay = calculate_aic(mse_clay, n = nrow(m), k = 1)
        aic_gumb = calculate_aic(mse_gumb, n = nrow(m), k = 1)
        aic_fran = calculate_aic(mse_fran, n = nrow(m), k = 1)     
        
        min_mse = min(mse_plack, mse_gauss, mse_t, mse_clay,
                      mse_gumb, mse_fran)
        
        min_aic = min(aic_plack, aic_gauss, aic_t, aic_clay,
                      aic_gumb, aic_fran)
        
      if (aic_ind < min_aic) {
          pval_cvm  =-1
          pval_ks = -1
          stat_cvm = -1
          stat_ks = -1
        family = 0
      }
        
        if (aic_ind >= min_aic) {
        cop = BiCopSelect(u, v, familyset = temp1, indeptest = FALSE)
        family = cop$family
        
        gof_pars = BiCopGofTest(u,v, cop,method = "kendall", B = 500)
        pval_cvm  =gof_pars$p.value.CvM
        pval_ks = gof_pars$p.value.KS
        stat_cvm = gof_pars$statistic.CvM
        stat_ks = gof_pars$statistic.KS
        
        }
          
      }
      
      
      tau = cor(u,v, method = "kendall")
      
      cop_prob_sev_08dr = BiCopCDF(pml(-0.8,bestdist_spi),pml(mean_yld, bestdist_yld), cop)
      cop_prob_sev_99dr = BiCopCDF(pml(-0.99,bestdist_spi),pml(mean_yld, bestdist_yld), cop)
      cop_prob_sev_130dr = BiCopCDF(pml(-1.3,bestdist_spi),pml(mean_yld, bestdist_yld), cop)
      cop_prob_sev_149dr = BiCopCDF(pml(-1.49,bestdist_spi),pml(mean_yld, bestdist_yld), cop)
      cop_prob_sev_199dr = BiCopCDF(pml(-1.99,bestdist_spi),pml(mean_yld, bestdist_yld), cop)
      
      
      cond_prob_yld_sev_08dr = cop_prob_sev_08dr/pml(-0.8,bestdist_spi)
      cond_prob_yld_sev_99dr = cop_prob_sev_99dr/pml(-0.99,bestdist_spi)
      cond_prob_yld_sev_130dr = cop_prob_sev_130dr/pml(-1.30,bestdist_spi)
      cond_prob_yld_sev_149dr = cop_prob_sev_149dr/pml(-1.49,bestdist_spi)
      cond_prob_yld_sev_199dr = cop_prob_sev_199dr/pml(-1.99,bestdist_spi)
      
      
      if(pval_cvm>=0 & pval_cvm<0.05){
        set.seed(123)
        plack_par = summary(fitCopula(copula = plackettCopula(), data = cbind(u,v)))$coefficient[1]
        plack_cop = plackettCopula(plack_par)
        plack_gof = gofCopula(plack_cop, x = cbind(u,v), N=500)
        pval_cvm = plack_gof$p.value
        stat_cvm = plack_gof$statistic
        
        
        
        if(pval_cvm>=0.05){
          family = 1000
          cop_prob_sev_08dr =pCopula(c(pml(-0.8,bestdist_spi),pml(mean_yld, bestdist_yld)), plack_cop)
          cop_prob_sev_99dr =pCopula(c(pml(-0.99,bestdist_spi),pml(mean_yld, bestdist_yld)), plack_cop)
          cop_prob_sev_130dr =pCopula(c(pml(-1.30,bestdist_spi),pml(mean_yld, bestdist_yld)), plack_cop)
          cop_prob_sev_149dr =pCopula(c(pml(-1.49,bestdist_spi),pml(mean_yld, bestdist_yld)), plack_cop)
          cop_prob_sev_199dr =pCopula(c(pml(-1.99,bestdist_spi),pml(mean_yld, bestdist_yld)), plack_cop)
          
          cond_prob_yld_sev_08dr = cop_prob_sev_08dr/pml(-0.8,bestdist_spi)
          cond_prob_yld_sev_99dr = cop_prob_sev_99dr/pml(-0.99,bestdist_spi)
          cond_prob_yld_sev_130dr = cop_prob_sev_130dr/pml(-1.30,bestdist_spi)
          cond_prob_yld_sev_149dr = cop_prob_sev_149dr/pml(-1.49,bestdist_spi)
          cond_prob_yld_sev_199dr = cop_prob_sev_199dr/pml(-1.99,bestdist_spi)
          
          
        }
      }      
      
      cat(temp_stctid, family, round(cop$par,3), round(cop$par2,3), round(mean_yld,3),
          round(cop_prob_sev_08dr,2), round(cop_prob_sev_99dr,2), round(cop_prob_sev_130dr,2), round(cop_prob_sev_149dr,2), round(cop_prob_sev_199dr,2),
          round(cond_prob_yld_sev_08dr,2), round(cond_prob_yld_sev_99dr,2), round(cond_prob_yld_sev_130dr,2), round(cond_prob_yld_sev_149dr,2),round(cond_prob_yld_sev_199dr,2),
          pml(-0.8,bestdist_spi), pml(-0.99,bestdist_spi), pml(-1.3,bestdist_spi), pml(-1.49,bestdist_spi),pml(-1.99,bestdist_spi),
          pml(mean_yld, bestdist_yld), tau, BiCopIndTest(u, v)$p.value, pval_cvm, pval_ks, stat_cvm, stat_ks ,nrow(m), '\n')
      
      temp_info = c(temp_stctid, family, cop$par, cop$par2, mean_yld,
                    cop_prob_sev_08dr, cop_prob_sev_99dr, cop_prob_sev_130dr, cop_prob_sev_149dr,cop_prob_sev_199dr,
                    cond_prob_yld_sev_08dr, cond_prob_yld_sev_99dr, cond_prob_yld_sev_130dr, cond_prob_yld_sev_149dr,cond_prob_yld_sev_199dr,
                    pml(-0.8,bestdist_spi), pml(-0.99,bestdist_spi), pml(-1.3,bestdist_spi), pml(-1.49,bestdist_spi),pml(-1.99,bestdist_spi),
                    pml(mean_yld, bestdist_yld), tau, BiCopIndTest(u, v)$p.value, pval_cvm, pval_ks, stat_cvm, stat_ks, nrow(m))
      
      copula_ct_all_list =rbind(copula_ct_all_list, temp_info)
    },
    error = function(cond) {an.error<<-TRUE})
    
    
  }
  colnames(copula_ct_all_list) <- c("STCTID", "Family", "Par", "Par2", "Avg Yld","cop08", "cop99", "cop130", "cop149", "cop199", 
                                    "ProbDr08", "ProbDr99", "ProbDr130", "ProbDr149", "ProbDr199",
                                    "cdf08", "cdf99", "cdf130", "cdf149", "cdf199", "cdf_Yldavg",'tau',
                                    'p-val',"pval_cvm", "pval_ks", "stat_cvm", "stat_ks", "len")
  
  return(copula_ct_all_list)
  
}


cond_prob_spi3_stdyld = copula_main(yld_pivot=alfalfa_yld_ct_std_pivot2, spi_pivot=spi_fin_stct2)

write.csv(cond_prob_spi3_stdyld, "/home/uallxr001/clim_shock/r_cop_output/alfalfa_cond_prob_and_GOF_and_aic_spi_growper_R_1991_20.csv", row.names=FALSE)
