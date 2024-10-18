# Revised: Aug 8, 2021

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
		# Jan 3, 2022
		Err2 <- (est - true_value)^2
		
    c(est, sqrt(v_est)/est, RB, RAB, Err2, CVRG)
}

f_angler_survey_bs_simu <- function(n_simu = 100, Species, est_by_daytype = TRUE, true_value, 
                                    DailyBHTrue, NumDays, duration_route, route_proto, 
                                    strat_by_wkday= TRUE, ndays_s = 20, is_simu_pop = TRUE, 
                                    angler_trip_pop, n_shifts = 2, shift_times = c(4, 12, 22), 
                                    shifts_non_overlapping = FALSE, adjust_overlap_prob, 
                                    expansion_method, opt_int_sampling, max_nints, fpc)
{
	res.mat <- matrix(NA, nrow = n_simu, ncol = 18, byrow = T)
	res.mat <- data.frame(res.mat)
	res.VWD <- rep(0, n_simu)
	res.V_2stage <- rep(0, n_simu)

	names(res.mat) <- c("BH", "CV_BH", "RB_BH", "RAB_BH", "Err2_BH", "COVG_BH", 
	                    "effort", "CV_effort", "RB_effort", "RAB_effort", "Err2_effort", "COVG_effort", "catch", "CV_catch", "RB_catch", "RAB_catch", "Err2_catch", "COVG_catch")
	
	for (i in 1:n_simu)
	{
    sample_data <- make_sample_bus_route(is_simu_pop, angler_trip_pop, route_proto, duration_route, ndays_s, 
                                         n_shifts, shift_times, shifts_non_overlapping, adjust_overlap_prob, expansion_method,
                                         strat_by_wkday, opt_int_sampling, max_nints)
		
		n_repl <- length(sample_data$interviews_repl)
		if (n_repl==1){
			r <- 1
			if (i==1){
				sample_data$interviews_repl[[r]]$Waterbody <- 1
				interviews_simu <- sample_data$interviews_repl[[r]]
				counts_simu <- sample_data$counts_repl[[r]]
			}else{
				sample_data$interviews_repl[[r]]$Waterbody <- i
				interviews_simu <- rbind(interviews_simu, sample_data$interviews_repl[[r]])
				counts_simu <- rbind(counts_simu, sample_data$counts_repl[[r]])
			}
    }
    
		if ((!adjust_overlap_prob)&& (expansion_method)){  
		  total_route_time = shift_times[length(shift_times)]-shift_times[1]
		}
		
	  catch_effort <- f_bus_route_creel_estimation(Species, est_by_daytype, sample_data$interviews_repl, sample_data$counts_repl, 
	                                               total_route_time, route_proto, NumDays, DailyBHTrue, V_BD, fpc)
	  
		res.VWD[i] <- catch_effort$V_WD
		res.V_2stage[i] <- catch_effort$V_BH_2stage

    #print(c(catch_effort$v_angler_hours, catch_effort$v_catch, catch_effort$v_boat_hours), digits = 2)
		res.mat[i, 1:6] <- f_sum_stat(catch_effort$boat_hours, catch_effort$v_boat_hours, true_value$boat_hours)
		res.mat[i, 7:12] <- f_sum_stat(catch_effort$angler_hours, catch_effort$v_angler_hours, true_value$effort)
		res.mat[i, 13:18] <- f_sum_stat(catch_effort$catch, catch_effort$v_catch, true_value$catch)
    cat("V_2stage" , res.V_2stage[i], "\n")
    print(res.mat[i, c(1, 2, 3, 6, 7, 8, 9, 12, 13, 14, 15, 18)], digits = 2, sep=":\t")
	}
	
	#Save interviews and counts
	#write.csv(counts_simu, "simu_counts.csv", row.names = FALSE)
	#write.csv(interviews_simu, "simu_interviews.csv", row.names = TRUE)

	list(Estimates = res.mat, VWD = res.VWD, V_2stage = res.V_2stage)
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

	matplot(x[,c(1,7,13)], type="l", pch = c(1,2,3))
	abline(h = True_Pop_Values)
	
	res <- data.frame(BH = c(summ[1:6],VAR[1],VEST[1],VEST[1]/VAR[1], sqrt(VAR[1])/summ[1]),
	                  E = c(summ[7:12], VAR[2], VEST[2],VEST[2]/VAR[2], sqrt(VAR[2])/summ[7]),
	                  C = c(summ[13:18], VAR[3], VEST[3], VEST[3]/VAR[3], sqrt(VAR[3])/summ[13]))
	rownames(res) <- c("EST", "CVest", "RB", "MAPE", "RMSE", "COVG", "VAR", "VEST","VEST2VAR","CV_TRUE")
	res
}


f_summaries <- function(x, True_Pop_Values)
{

 #"BH", "CV_BH", "RB_BH", "RAB_BH", "RMSE_BH", "COVG_BH", "effort", "CV_effort", "RB_effort", "RAB_effort", "RMSE_effort", "COVG_effort",  "catch", "CV_catch", "RB_catch", "RAB_catch", "RMSE_catch", "COVG_catch")
 np = length(names(x))
 summ <- numeric(np)
 summ <- apply(x, 2, mean, na.rm = TRUE)
 names(summ) <- names(x)

 par(mfrow=c(1,3))
 plot(x[,1], main = "BoatHour estimates", type = "l")
 abline(h = True_Pop_Values$boat_hours, lwd = 2, col = "blue")
 plot(x[,7], main = "Effort estimates", type = "l")
 abline(h = True_Pop_Values$effort, lwd = 2, col = "blue")
 plot(x[,13], main = "Catch estimates", type = "l")
 abline(h = True_Pop_Values$catch, lwd = 2, col = "blue")

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

              temp_res <- f_angler_survey_simu(isimu, n_simu, dailyEst, do_daytype_est, UseTwoStageVar, true_value, NumDays, FValue, strat_by_wkday, ndays_s, n_cnts, max_nints,is_simu_pop, True_Pop, n_shifts,len_shift,opt_int_sampling, exclude_method, n_days_excluded, imputation_method, n_impute)
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

