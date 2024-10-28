# Dr. Zhenming Su
# Institute for Fisheries Research 
# Michigan Department of Natural Resources 
#      and University of Michigan  
# ANN ARBOR, MI 48100

# Contact Zhenming Su (suz@michigan.gov) for any questions
# Revised for "Evaluation of Bus-route and Aerial-access Methods for Great Lakes Recreational Fisheries Surveys"
#  by Zhenming Su, 2/9/2024 

# Revised: Aug 8, 2021
# Revised: 2/9/2024

# Software built for the research paper:
# Zhenming Su & David Clapp (2013): Evaluation of Sample Design and Estimation Methods for Great Lakes
#   Angler Surveys, Transactions of the American Fisheries Society, 142:1, 234-246
# To link to this article: http://dx.doi.org/10.1080/00028487.2012.728167

f_sum_stat <- function(est, v_est, true_value)
{
		Bias <- est - true_value
		RB   <- 100 * Bias / true_value 
		RAB  <- 100 * abs(Bias) / true_value
        #cat("est: ", est, "v_est: ", v_est,  "\n")

		CVRG <- 0
	    if (!is.na(v_est)) {
			# 95% interval
			#LB <- est - 2 * sqrt(v_est)
			#UB <- est + 2 * sqrt(v_est)
            # 80% interval: qz = 1.28
			LB <- est - 1.28 * sqrt(v_est)
			UB <- est + 1.28 * sqrt(v_est)
			if ((true_value >= LB) & (true_value < UB))
				CVRG <- 1
			else
				CVRG <- 0
        }
        # Error square for calculation RMSE
		# Jan 4, 2022
		Err2 <- (est - true_value)^2
		
    c(est, sqrt(v_est)/est, RB, RAB, Err2, CVRG)
}

f_angler_survey_simu <- function(Species, n_simu = 100, dailyEst = TRUE, est_by_daytype = TRUE, true_value, 
                                 AerialProportion, PropCorr, NumDays, FValue, strat_by_wkday= TRUE, 
                                 RovingCounting, RovingInterviewing, ndays_s = 20, ncnts = 2, 
                                 progressive_count_interval, max_nints = 30, is_simu_pop = TRUE, 
                                 simu_Year, simu_Month, angler_trip_pop, n_shifts = 1, shift_times = c(7, 19), 
                                 n_shifts_int = 2, shift_times_int = c(4, 12, 22), 
                                 shifts_non_overlapping = FALSE, opt_int_sampling, 
                                 exclude_method, n_days_excluded, imputation_method, n_impute)
    
{
	res.mat <- matrix(NA, nrow = n_simu, ncol = 18, byrow = T)
	res.mat <- data.frame(res.mat)
	names(res.mat) <- c("BH", "cv_BH", "RB_BH", "RAB_BH", "MSE_BH", "COVG_BH", "Effort", "cv_effort", "RB_effort", "RAB_effort", "MSE_BH", "COVG_effort", "Catch", "cv_catch", "RB_catch", "RAB_catch", "MSE_BH", "COVG_catch")
	
	for (i in 1:n_simu)
	{
		sample_data <- make_sample(is_simu_pop, angler_trip_pop, RovingCounting, RovingInterviewing, ndays_s, 
		                n_shifts, shift_times, n_shifts_int, shift_times_int, shifts_non_overlapping, AerialProportion, ncnts, 
						        progressive_count_interval, max_nints, strat_by_wkday, opt_int_sampling)

		if (i==1){
			sample_data$interviews$Waterbody <- 1
			interviews_simu <- sample_data$interviews
			sample_data$counts$Waterbody <- 1
			counts_simu <- sample_data$counts
		}else
		{
			sample_data$interviews$Waterbody <- i
			interviews_simu <- rbind(interviews_simu, sample_data$interviews)
			sample_data$counts$Waterbody <- i
			counts_simu <- rbind(counts_simu, sample_data$counts)
		}

	   catch_effort <- f_creel_estimation(simu_Year, simu_Month, shift_times, dailyEst, est_by_daytype, Species, 
	                                      exclude_method, n_days_excluded, sample_data$interviews, sample_data$counts, 
	                                      imputation_method, n_impute, RovingCounts, AerialProportion, PropCorr, NumDays, FValue)

     #print(c(catch_effort$v_angler_hours, catch_effort$v_catch, catch_effort$v_boat_hours), digits = 2)
		 res.mat[i, 1:6] <- f_sum_stat(catch_effort$boat_hours, catch_effort$v_boat_hours, true_value$boat_hours)
		 res.mat[i, 7:12] <- f_sum_stat(catch_effort$angler_hours, catch_effort$v_angler_hours, true_value$effort)
		 res.mat[i, 13:18] <- f_sum_stat(catch_effort$catch, catch_effort$v_catch, true_value$catch)
     print(res.mat[i, c(1, 2, 3, 6, 7, 8, 9, 12, 13, 14, 15, 18)], digits = 2, sep=":\t")
	}
	
	#Save interviews and counts
	#write.csv(counts_simu, "simu_counts.csv", row.names = FALSE)
	#write.csv(interviews_simu, "simu_interviews.csv", row.names = TRUE)

	res.mat
}


f_summaries_table <- function(x, True_Pop_Values)
{
	# 1-6:   "BH", "CV_BH", "RB_BH", "RAB_BH", "Err2_BH", "COVG_BH"
    # 7-12:  "effort", "CV_effort", "RB_effort", "RAB_effort", "Err2_effort", "COVG_effort"
    # 13-18: "catch", "CV_catch", "RB_catch", "RAB_catch", "Err2_catch", "COVG_catch"

	# Mean Est, Mean CV, Mean Bias, MAPE
	summ <- apply(x, 2, mean, na.rm = TRUE)

	# Variance
    VAR <- apply(x[,c(1,7,13)],2, var, na.rm = TRUE)

	VEST <- (summ[c(2,8,14)]*summ[c(1,7,13)])^2
	
	# RMSE
  summ[c(5,11,17)] <- sqrt(summ[c(5,11,17)])
	names(summ) <- c("BH", "CV_BH", "RB_BH", "MAPE_BH", "RMSE_BH", "COVG_BH",
					 "E", "CV_E", "RB_E", "MAPE_E", "RMSE_E", "COVG_E", 
					 "C", "CV_C", "RB_C", "MAPE_C", "RMSE_C", "COVG_C")

	#matplot(x[,c(1,7,13)], type="l", pch = c(1,2,3))
	#abline(h = True_Pop_Values)

	res <- data.frame(BH = c(summ[1:6], sqrt(VAR[1])/True_Pop_Values$boat_hours, sqrt(VEST[1])/summ[1], VEST[1]/VAR[1]), 
	                  E = c(summ[7:12], sqrt(VAR[2])/True_Pop_Values$effort, sqrt(VEST[2])/summ[7],VEST[2]/VAR[2]), 
	                  C = c(summ[13:18], sqrt(VAR[3])/True_Pop_Values$catch, sqrt(VEST[3])/summ[13], VEST[3]/VAR[3]))
	rownames(res) <- c("EST", "cv", "RB", "MAPE", "RMSE", "COVG", "CV", "cv_EST","VEST2VAR")
	res
}


f_summaries <- function(x, True_Pop_Values)
{

#"BH", "V_BH", "RB_BH", "RAB_BH", "RMSE_BH", "COVG_BH", "effort", "V_effort", "RB_effort", "RAB_effort", "RMSE_effort", "COVG_effort", "catch", "V_catch", "RB_catch", "RAB_catch", "RMSE_catch", "COVG_catch")
np = length(names(x))
summ <- numeric(np)
summ <- apply(x, 2, mean, na.rm = TRUE)
names(summ) <- names(x)
# 
#  par(mfrow=c(1,3))
#  plot(x[,1], main = "BoatHour estimates", type = "l")
#  abline(h = True_Pop_Values$boat_hours, lwd = 2, col = "blue")
#  plot(x[,7], main = "Effort estimates", type = "l")
#  abline(h = True_Pop_Values$effort, lwd = 2, col = "blue")
#  plot(x[,13], main = "Catch estimates", type = "l")
#  abline(h = True_Pop_Values$catch, lwd = 2, col = "blue")

res <- summ
}

f_simu_creel <- function(R_NM, len_shift, n_simu = 200)
{
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # sampling options
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  simu_cond <- f_simu_conditions()

  for (EST_METHOD in simu_cond$EstMethods)
  {
    if (EST_METHOD== "DailyEst") 
       dailyEst <- TRUE
    else
       dailyEst <- FALSE
    

    for (opt_int_sampling in simu_cond$OptIntSampling)
    {
      SUMM_FN <- paste("SUMM_", R_NM, "_", EST_METHOD, "_", opt_int_sampling, sep="")

       nr <-  length(simu_cond$NDAYS) * length(simu_cond$NSHIFTS)  * length(simu_cond$NCNTS)   #24

      if (opt_int_sampling== "threshhold") 
         nr <- nr * length(simu_cond$MaxNInt)

      res <- matrix(0, nrow= nr, ncol= 45)
      res <- data.frame(res)

      colnames(res) <- c("NDAY", "MaxNInt", "NSHIFTS", "NCNTS", "ESTIMATOR","INT_SAMPLING", 
          "E","TSE_E",  "SE_E", "RB_E", "MAE_E", "RMSE_E", "COVG_E", 
          "C", "TSE_C", "SE_C", "RB_C", "MAE_C",  "RMSE_C", "COVG_C", 
          "B", "TSE_B", "SE_B", "RB_B", "MAE_B",  "RMSE_B", "COVG_B", 
          "E_med", "RB_E_med", "C_med", "RB_C_med", "B_med", "RB_B_med",
          "E_q025", "RB_E_q025", "C_q025", "RB_C_q025", "B_q025", "RB_B_q025",
          "E_q975", "RB_E_q975", "C_q975", "RB_C_q975",  "B_q975", "RB_B_q975")  

      isimu <- 1
      for (dd in simu_cond$NDAYS)
      {
        # number of days to be sampled
        ndays_s <- dd; 
        for (n_shifts in simu_cond$NSHIFTS)
        {
        if (n_shifts > 3) stop("n_shifts > 3 is not handled!")

          for (ncnts in simu_cond$NCNTS)
          {
            for (mint in simu_cond$MaxNInt)
            {
              n_cnts <-  ncnts
			        # No subsampling fishing day, ENTIRE
              if (n_shifts == 1) 
              {
                max_nints <- 1 * mint
                n_cnts <- 1 * ncnts
              } 

              if (n_shifts == 2) max_nints <- mint
              if (n_shifts == 3) max_nints <- mint #round(2 * mint/3)

              res[isimu,1:6] <- c(dd, max_nints, n_shifts, n_cnts, EST_METHOD, opt_int_sampling)

              temp_res <- f_angler_survey_simu(isimu, n_simu, dailyEst, do_daytype_est, UseTwoStageVar,
                                               true_value, NumDays, FValue, strat_by_wkday, ndays_s, 
                                               n_cnts, max_nints,is_simu_pop, True_Pop, n_shifts,len_shift, 
                                               opt_int_sampling, exclude_method, n_days_excluded, 
                                               imputation_method, n_impute)
              
              write.table(temp_res, paste(SUMM_FN, "_mid.txt", sep=""), col.names = NA, sep ="\t")

              summ <- round(f_summaries("mean", temp_res), 2)
              res[isimu, 7:27] <- summ

              print(isimu)
              print(res[isimu, 1:27])

              summ <- round(f_summaries("median", temp_res), 2)
              res[isimu, 28:33] <- summ
              summ <- round(f_summaries("quart025", temp_res), 2)
              res[isimu, 34:39] <- summ
              summ <- round(f_summaries("quart975", temp_res), 2)
              res[isimu, 40:45] <- summ

              write.table(res, paste(SUMM_FN, "_res.txt", sep=""), col.names = NA, sep ="\t")
              rm(temp_res)
              rm(summ)

              isimu <- isimu + 1

              if (opt_int_sampling== "fix_prop") break
            }
          }
       }
      }
      res <- as.data.frame(res)

      SUMM_FN <- paste(SUMM_FN, ".txt", sep="")

      print("*************************************")
      print(format(Sys.time(), "%a %b %d %H:%M:%S %Y"))
      print(SUMM_FN)
      print("*************************************")
      print(res)

      write.table(res, SUMM_FN, col.names = NA, sep ="\t")
    }
  }
}

f_simu_conditions <- function()
{
EM = c("MonthlyEst","DailyEst")
ND = c(31, 20, 10)
NS = c(1, 2)
NC = c(1, 2, 4, 8)
NI = c(10, 50, 100, 4000)

OptIntSampling = c("threshhold", "fix_prop")

#EM = c("DailyEst")
ND = c(31)
NS = c(3)
NC = c(18)
NI = c(10)

OptIntSampling = c("threshhold")

return(list(EstMethods = EM, OptIntSampling = OptIntSampling, NDAYS = ND, NSHIFTS = NS, NCNTS = NC, MaxNInt = NI))
}

