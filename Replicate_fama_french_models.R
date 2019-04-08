# Dissecting Anomalies with a Five-Factor Model , Fama and French, 2016
# Original paper:   doi:10.1093/rfs/hhv043
# - Data Replication
# By Sherry He


#### Note that some discrepancies from the 2016 paper might be due to 
#     rounding, or the changes in data source (CRSP Data), as stated on
#     http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/changes_crsp.html

library(dplyr)
library(GRS.test)

########## Step 0 ##########
## Define functions for data processing
process_input_data <- function(input_df){
    colnames(input_df)[1] <- "YYYYMM"
    n_col = ncol(input_df)
    for (i in 1:n_col) {
      input_df[,i] <- as.numeric(as.character(input_df[,i]))}
    return(input_df)
}

filter_by_date <- function(input_df){
  output_df <- NULL
  output_df <- filter(input_df,as.numeric(YYYYMM)>=196307 & as.numeric(YYYYMM)<=201412)
  return(output_df)
}

########## Table 1 ##########
## Averages, standard deviations, and t-statistics for monthly factor returns, 
## July 1963–December 2014 (618 months)
#############################

# Read and process factors: RM-RF SMB HML RMW CMA
three_factors_m <- read.csv("F-F_Research_Data_5_Factors_2x3.CSV",header=TRUE, skip=3)
three_factors_m <- process_input_data(three_factors_m)
three_factors_m_filtered <- filter_by_date(three_factors_m)

# Read and process factor MOM
factor_mom <- read.csv("F-F_Momentum_Factor.CSV",header=TRUE, skip=13)
factor_mom <- process_input_data(factor_mom)
factor_mom_filtered <- filter_by_date(factor_mom)

# Formulate Table 1
three_factors_m_filtered$MOM <- factor_mom_filtered$Mom
Table_1 <- rbind (round(apply(three_factors_m_filtered[c(2:6,8)],2,mean,na.rm=TRUE),digits = 2),
                  round(apply(three_factors_m_filtered[c(2:6,8)],2,sd,na.rm=TRUE),digits = 2),
                  round(apply(three_factors_m_filtered[c(2:6,8)],2,function(x) t.test(x,mu=0)$statistic),digits = 2))
rownames(Table_1) = c('Mean','SD','t-statistic')
print(Table_1)

########## Table 2 ##########
## Summary statistics for tests of three-, four-, and five-factor models,
## July 1963–December 2014 (618 months)
## Using Average Value Weighted Returns -- Monthly
## "Mkt.RF" , "SMB" , "HML" , "RMW" , "CMA" 
#############################

RF <- three_factors_m_filtered$RF

T2_colnames<- c("GRS",
                "GRS_p_val",
                "mean_abs_a",
                "mean_ai_mean_ri",
                "mean_a2_mean_r2",
                "mean_a_var_mean_a2",
                "mean_adj_R2" )

# Function to generate one row for a model
generate_t2_row <- function(input_ret, factor_string){
  return_mat <- input_ret
  no_of_factor <- length(factor_string)
  factor_mat = data.matrix(three_factors_m_filtered[,factor_string])
  result <- GRS.test(return_mat,factor_mat)
  # calculate alpha
  reg.func <- function (y, m) {
    mod <- lm(y ~ m)
    mod$coefficients["(Intercept)"]
  }
  alphas <- as.vector(apply(factor_mat, 2, reg.func, return_mat))
  # calculate mean(ai^2)/mean(ri^2)
  R_i <- colMeans(return_mat) - mean(three_factors_m_filtered$Mkt.RF) 
  reg.func_error <- function (y, m) {
    mod <- lm(y ~ m)
    summary(mod)$coefficients[1,2]
  }
  alphas_error <- as.vector(apply(return_mat, 2, reg.func_error, factor_mat))
  # calculate adjustedR squared
  reg.func_adj_R2 <- function (y, m) {
    mod <- lm(y ~ m)
    summary(mod)$adj.r.squared
  }
  alphas_adj_R2 <- as.vector(apply(return_mat, 2, reg.func_adj_R2, factor_mat))
  
  output_row <- c(result$GRS.stat,
                  result$GRS.pval,
                  mean(abs(alphas)),
                  mean(abs(alphas))/mean(abs(R_i)),
                  mean(alphas^2)/mean(R_i^2),
                  mean(alphas_error^2)/mean(alphas^2),
                  mean(alphas_adj_R2)
                  )
}

################################
##### 25 Size-β portfolios ##### 
size_beta_25 <- read.csv("25_Portfolios_ME_BETA_5x5.CSV",header=TRUE, skip=16)
size_beta_25 <- size_beta_25[1:668,] 
size_beta_25 <- process_input_data(size_beta_25)
size_beta_25_filtered <- filter_by_date(size_beta_25)

# Form GRS test factor and return matrix
return_mat_1 <- data.matrix(size_beta_25_filtered[,2:26])
return_mat_1 <- sweep(return_mat,2,RF, FUN='-')

size_beta_25_Mkt <- generate_t2_row(return_mat_1,c("Mkt.RF")) 
size_beta_25_Mkt_SMB_HML <- generate_t2_row(return_mat_1,c("Mkt.RF" , "SMB" , "HML" )) 
size_beta_25_Mkt_SMB_HML_RMW <- generate_t2_row(return_mat_1,c("Mkt.RF" , "SMB" , "HML" , "RMW" )) 
size_beta_25_Mkt_SMB_HML_CMA <- generate_t2_row(return_mat_1,c("Mkt.RF" , "SMB" , "HML" , "CMA" )) 
size_beta_25_Mkt_SMB_RMW_CMA <- generate_t2_row(return_mat_1,c("Mkt.RF" , "SMB" , "RMW" , "CMA" )) 
size_beta_25_Mkt_SMB_HML_RMW_CMA <- generate_t2_row(return_mat_1,c("Mkt.RF" , "SMB" , "HML", "RMW" , "CMA" )) 

#################################
##### 35 Size-NI portfolios ##### 
size_NI_35 <- read.csv("25_Portfolios_ME_NI_5x5.CSV",header=TRUE, skip=18)
size_NI_35 <- size_NI_35[1:668,] 
size_NI_35 <- process_input_data(size_NI_35)
size_NI_35_filtered <- filter_by_date(size_NI_35)

# Form GRS test factor and return matrix
return_mat_2 <- data.matrix(size_NI_35_filtered[,2:36])
return_mat_2 <- sweep(return_mat_2,2,RF, FUN='-')

size_NI_35_Mkt_SMB_HML <- generate_t2_row(return_mat_2,c("Mkt.RF" , "SMB" , "HML" )) 
size_NI_35_Mkt_SMB_HML_RMW <- generate_t2_row(return_mat_2,c("Mkt.RF" , "SMB" , "HML" , "RMW" )) 
size_NI_35_Mkt_SMB_HML_CMA <- generate_t2_row(return_mat_2,c("Mkt.RF" , "SMB" , "HML" , "CMA" )) 
size_NI_35_Mkt_SMB_RMW_CMA <- generate_t2_row(return_mat_2,c("Mkt.RF" , "SMB" , "RMW" , "CMA" )) 
size_NI_35_Mkt_SMB_HML_RMW_CMA <- generate_t2_row(return_mat_2,c("Mkt.RF" , "SMB" , "HML", "RMW" , "CMA" )) 

##################################
##### 25 Size-Var portfolios ##### 
size_var_25 <- read.csv("25_Portfolios_ME_VAR_5x5.csv",header=TRUE, skip=19)
size_var_25 <- size_var_25[1:668,] 
size_var_25 <- process_input_data(size_var_25)
size_var_25_filtered <- filter_by_date(size_var_25)

# Form GRS test factor and return matrix
return_mat_3 <- data.matrix(size_var_25_filtered[,2:26])
return_mat_3 <- sweep(return_mat_3,2,RF, FUN='-')

size_var_25_Mkt_SMB_HML <- generate_t2_row(return_mat_3,c("Mkt.RF" , "SMB" , "HML" )) 
size_var_25_Mkt_SMB_HML_RMW <- generate_t2_row(return_mat_3,c("Mkt.RF" , "SMB" , "HML" , "RMW" )) 
size_var_25_Mkt_SMB_HML_CMA <- generate_t2_row(return_mat_3,c("Mkt.RF" , "SMB" , "HML" , "CMA" )) 
size_var_25_Mkt_SMB_RMW_CMA <- generate_t2_row(return_mat_3,c("Mkt.RF" , "SMB" , "RMW" , "CMA" )) 
size_var_25_Mkt_SMB_HML_RMW_CMA <- generate_t2_row(return_mat_3,c("Mkt.RF" , "SMB" , "HML", "RMW" , "CMA" )) 

###########################################
##### 25 Size-Residual Var portfolios ##### 
size_Rvar_25 <- read.csv("25_Portfolios_ME_RESVAR_5x5.csv",header=TRUE, skip=19)
size_Rvar_25 <- size_Rvar_25[1:668,] 
size_Rvar_25 <- process_input_data(size_Rvar_25)
size_Rvar_25_filtered <- filter_by_date(size_Rvar_25)

# Form GRS test factor and return matrix
return_mat_4 <- data.matrix(size_Rvar_25_filtered[,2:26])
return_mat_4 <- sweep(return_mat_4,2,RF, FUN='-')

size_Rvar_25_Mkt_SMB_HML <- generate_t2_row(return_mat_4,c("Mkt.RF" , "SMB" , "HML" )) 
size_Rvar_25_Mkt_SMB_HML_RMW <- generate_t2_row(return_mat_4,c("Mkt.RF" , "SMB" , "HML" , "RMW" )) 
size_Rvar_25_Mkt_SMB_HML_CMA <- generate_t2_row(return_mat_4,c("Mkt.RF" , "SMB" , "HML" , "CMA" )) 
size_Rvar_25_Mkt_SMB_RMW_CMA <- generate_t2_row(return_mat_4,c("Mkt.RF" , "SMB" , "RMW" , "CMA" )) 
size_Rvar_25_Mkt_SMB_HML_RMW_CMA <- generate_t2_row(return_mat_4,c("Mkt.RF" , "SMB" , "HML", "RMW" , "CMA" )) 

#################################
##### 25 Size-AC portfolios ##### 
size_AC_25 <- read.csv("25_Portfolios_ME_AC_5x5.csv",header=TRUE, skip=18)
size_AC_25 <- size_AC_25[1:668,] 
size_AC_25 <- process_input_data(size_AC_25)
size_AC_25_filtered <- filter_by_date(size_AC_25)

# Form GRS test factor and return matrix
return_mat_5 <- data.matrix(size_AC_25_filtered[,2:26])
return_mat_5 <- sweep(return_mat_5,2,RF, FUN='-')

size_AC_25_Mkt_SMB_HML <- generate_t2_row(return_mat_5,c("Mkt.RF" , "SMB" , "HML" )) 
size_AC_25_Mkt_SMB_HML_RMW <- generate_t2_row(return_mat_5,c("Mkt.RF" , "SMB" , "HML" , "RMW" )) 
size_AC_25_Mkt_SMB_HML_CMA <- generate_t2_row(return_mat_5,c("Mkt.RF" , "SMB" , "HML" , "CMA" )) 
size_AC_25_Mkt_SMB_RMW_CMA <- generate_t2_row(return_mat_5,c("Mkt.RF" , "SMB" , "RMW" , "CMA" )) 
size_AC_25_Mkt_SMB_HML_RMW_CMA <- generate_t2_row(return_mat_5,c("Mkt.RF" , "SMB" , "HML", "RMW" , "CMA" )) 

#########################################
##### 25 Size-Prior 2–12 portfolios ##### 
# (the sum of a stock’s monthly returns from t−12 to t−2))
size_Prior12_2_25 <- read.csv("25_Portfolios_ME_Prior_12_2.csv",header=TRUE, skip=11)
size_Prior12_2_25 <- size_Prior12_2_25[1:1106,] 
size_Prior12_2_25 <- process_input_data(size_Prior12_2_25)
size_Prior12_2_25_filtered <- filter_by_date(size_Prior12_2_25)

# Form GRS test factor and return matrix
return_mat_6 <- data.matrix(size_Prior12_2_25_filtered[,2:26])
return_mat_6 <- sweep(return_mat_6,2,RF, FUN='-')

size_Prior12_2_25_Mkt_SMB_HML <- generate_t2_row(return_mat_6,c("Mkt.RF" , "SMB" , "HML" )) 
size_Prior12_2_25_Mkt_SMB_HML_RMW <- generate_t2_row(return_mat_6,c("Mkt.RF" , "SMB" , "HML" , "RMW" )) 
size_Prior12_2_25_Mkt_SMB_HML_CMA <- generate_t2_row(return_mat_6,c("Mkt.RF" , "SMB" , "HML" , "CMA" )) 
size_Prior12_2_25_Mkt_SMB_RMW_CMA <- generate_t2_row(return_mat_6,c("Mkt.RF" , "SMB" , "RMW" , "CMA" )) 
size_Prior12_2_25_Mkt_SMB_HML_RMW_CMA <- generate_t2_row(return_mat_6,c("Mkt.RF" , "SMB" , "HML", "RMW" , "CMA" )) 

size_Prior12_2_25_Mkt_SMB_HML_MOM <- generate_t2_row(return_mat_6,c("Mkt.RF" , "SMB" , "HML" , "MOM" )) 
size_Prior12_2_25_Mkt_SMB_HML_RMW_MOM <- generate_t2_row(return_mat_6,c("Mkt.RF" , "SMB" , "HML" , "RMW" , "MOM" )) 
size_Prior12_2_25_Mkt_SMB_HML_CMA_MOM <- generate_t2_row(return_mat_6,c("Mkt.RF" , "SMB" , "HML" , "CMA" , "MOM" )) 
size_Prior12_2_25_Mkt_SMB_RMW_CMA_MOM <- generate_t2_row(return_mat_6,c("Mkt.RF" , "SMB" , "RMW" , "CMA" , "MOM" )) 
size_Prior12_2_25_Mkt_SMB_HML_RMW_CMA_MOM <- generate_t2_row(return_mat_6,c("Mkt.RF" , "SMB" , "HML", "RMW" , "CMA" , "MOM" )) 

#########################################
Table_2<-data.frame(rbind(
  size_beta_25_Mkt,
  size_beta_25_Mkt_SMB_HML,
  size_beta_25_Mkt_SMB_HML_RMW,
  size_beta_25_Mkt_SMB_HML_CMA,
  size_beta_25_Mkt_SMB_RMW_CMA,
  size_beta_25_Mkt_SMB_HML_RMW_CMA,
  
  size_NI_35_Mkt_SMB_HML,
  size_NI_35_Mkt_SMB_HML_RMW, 
  size_NI_35_Mkt_SMB_HML_CMA,
  size_NI_35_Mkt_SMB_RMW_CMA,
  size_NI_35_Mkt_SMB_HML_RMW_CMA,
  
  size_var_25_Mkt_SMB_HML,
  size_var_25_Mkt_SMB_HML_RMW,
  size_var_25_Mkt_SMB_HML_CMA,
  size_var_25_Mkt_SMB_RMW_CMA,
  size_var_25_Mkt_SMB_HML_RMW_CMA,
  
  size_Rvar_25_Mkt_SMB_HML,
  size_Rvar_25_Mkt_SMB_HML_RMW,
  size_Rvar_25_Mkt_SMB_HML_CMA,
  size_Rvar_25_Mkt_SMB_RMW_CMA,
  size_Rvar_25_Mkt_SMB_HML_RMW_CMA,
  
  size_AC_25_Mkt_SMB_HML,
  size_AC_25_Mkt_SMB_HML_RMW,
  size_AC_25_Mkt_SMB_HML_CMA,
  size_AC_25_Mkt_SMB_RMW_CMA,
  size_AC_25_Mkt_SMB_HML_RMW_CMA,
  
  size_Prior12_2_25_Mkt_SMB_HML,
  size_Prior12_2_25_Mkt_SMB_HML_RMW ,
  size_Prior12_2_25_Mkt_SMB_HML_CMA ,
  size_Prior12_2_25_Mkt_SMB_RMW_CMA ,
  size_Prior12_2_25_Mkt_SMB_HML_RMW_CMA ,
  
  size_Prior12_2_25_Mkt_SMB_HML_MOM ,
  size_Prior12_2_25_Mkt_SMB_HML_RMW_MOM ,
  size_Prior12_2_25_Mkt_SMB_HML_CMA_MOM ,
  size_Prior12_2_25_Mkt_SMB_RMW_CMA_MOM ,
  size_Prior12_2_25_Mkt_SMB_HML_RMW_CMA_MOM))

colnames(Table_2) = T2_colnames

print(Table_2) 

########## Table 3 ##########
## Average excess returns and characteristics of stocks in the 25 Size-β portfolios,
## July 1963–December 2014 (618 months)

#############################




data(data)
factor.mat = data[1:342,2:5]            # Fama-French 3-factor model
ret.mat = data[1:342,8:ncol(data)]      # 25 size-BM portfolio returns
GRS.test(ret.mat,factor.mat)$GRS.stat   # See Table 9C of Fama-French (1993)


