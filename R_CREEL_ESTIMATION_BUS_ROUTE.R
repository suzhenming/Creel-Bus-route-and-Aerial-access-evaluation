#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Bus-route creel survey estimation
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dr. Zhenming Su
# Institute for Fisheries Research 
# Michigan Department of Natural Resources 
#      and University of Michigan  
# ANN ARBOR, MI 48100

# Contact Zhenming Su (suz@michigan.gov) for any questions

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Built: Sept 1, 2021
#        Sept 5, 2021 worked  
# Jan/04/2021
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f_var_mean_bootstrap <- function(x){
	B = 5000
	n <- length(x)
	x.mean <- mean(x)
	x.boot <- matrix(sample(x, B*n, repl=TRUE), nrow = B)
	mean.boot <- rowMeans(x.boot)
	bias <- mean(mean.boot - x.mean)

	var(mean.boot)
}


f_catch_rate_bus_route <- function(CdHd)
{
	# ---------------------------------------------------------
	# Ratio of means estimation of catch rate and variance 
	#   Used for completed trips
	# ---------------------------------------------------------
	Cd <- CdHd[,1]
    Hd <- CdHd[,2]
	
	# ratio of means
	Rd  <- mean(Cd)/mean(Hd)

    SSQ_C   <- crossprod(Cd); #sum of squares
    SSQ_H   <- crossprod(Hd);
    SPD_C_H <- crossprod(Cd, Hd);

    # calculate variance of Rd
	n_int <- length(Hd);
    
	mHd <- mean(Hd)
    term1 <- 1/(n_int * mHd * mHd);
    term2 <- 2 * Rd * SPD_C_H;
    term3 <- Rd * Rd * SSQ_H;

    # Eq. 2
    V_R <- term1 * (SSQ_C - term2 + term3)/(n_int - 1);

	c(Rd, V_R)
}

bus_route_site_estimation <-  function(BoatCounts_sites, total_route_time, wait_time, trip_interviews)
{
  # Make daily estimates of access sites on the bus route 
	
	# Total boat counts made at each access site of the daily bus-route
    # BoatCounts_sites <- aggregate(list(FishBoatCounts = boat_counts$FishingBoatCounts), list(Day = as.factor(boat_counts$Day), AccessSite = boat_counts$AccessSite), sum)	

	# Averages of site-specific weighted trip lengths or catch for each sampled day
	#  Q can be the trip length or catch of a species
  Avg_trip_Q <- aggregate(list(m_trip_Q = trip_interviews$Q), list(Day = trip_interviews$Day, AccessSite = trip_interviews$AccessSite), mean)	
  # Expand the average Q at each site and day by the corresponding total number of trips counted
  route_daily_est <- merge(BoatCounts_sites, Avg_trip_Q, by=c("Day","AccessSite"), all.x = TRUE)
	route_daily_est$m_trip_Q[is.na(route_daily_est$m_trip_Q)] <- 0
	route_daily_est$site_total_Q <- route_daily_est$count * route_daily_est$m_trip_Q

	# Creel times at each access site
  route_daily_est$wait_time <- 0
	for (s in wait_time$SiteName){
	  route_daily_est[route_daily_est$AccessSite == s,]$wait_time <- wait_time[wait_time$SiteName == s,]$CreelTime_Min/60
  }
  # Creel estimates per unit of wait time
	route_daily_est$site_Q_per_wait_time <- route_daily_est$site_total_Q/route_daily_est$wait_time 
	
	# total daily site estimates
  route_daily_est$daily_site_Q <- total_route_time * route_daily_est$site_Q_per_wait_time
  route_daily_est
}

week_days <- function(date_str) {
    # use strftime function with "%u" format parameter to return the weekday
	weekdays <- as.numeric(strftime(as.Date(date_str, "%m/%d/%Y"), "%u"))
	weekdays
}


f_bus_route_total_boat_hours <- function(est_by_daytype, ints, boat_counts, total_route_time, 
                                         wait_times, NumDays, DailyBHTrue, fpc)
{
	# ----------------------------------------------
	# monthly effort estimation  
	#   only for one month 
	# ----------------------------------------------
    # Added within-day variance calculation on 12/30/2021

	# Calculate trip lengths in hours
	ints$trip_length <- ints$ET - ints$ST #as.numeric(difftime(ints$ET, ints$ST, units = "hours"))
	
	# trip lengths must be greater than zero
	if (any(ints$trip_length <= 0)){
	  #cat("trip_lengths are <= 0: ", "\n")
	  #print(ints[ints$trip_length <= 0, ])

	  ints <- ints[ints$trip_length > 0, ]
	}
	#stopifnot(ints$trip_length >= 0)

	# boat hours
	TripLenW <- ints$trip_length/ints$prob_sampling	

	trip_lengths <- data.frame(Q = TripLenW, AccessSite = ints$AccessSite, Day = ints$DAY)
	route_daily_site_BH_est <- bus_route_site_estimation(boat_counts, total_route_time, wait_times, trip_lengths)
	
	daily_BH_est <- aggregate(list(BH_day = route_daily_site_BH_est$daily_site_Q),list(Day = as.factor(route_daily_site_BH_est$Day)),sum)

	if (est_by_daytype){
	   date_str <- paste(ints$MONTH[1], "/", daily_BH_est$Day,"/", ints$YEAR[1], sep="")

	   daily_BH_est$weekday <- week_days(date_str)
	   daily_BH_est$DayType <- ifelse(daily_BH_est$weekday<6, "WD", "WE")

	   n_day <- as.numeric(tapply(daily_BH_est$BH_day, list(DayType = as.factor(daily_BH_est$DayType)), length))

	   psu_BH_daytype <- split(daily_BH_est$BH_day, as.factor(daily_BH_est$DayType));

	   mean_monthly_boat_hours <- as.numeric(tapply(daily_BH_est$BH_day, list(DayType = as.factor(daily_BH_est$DayType)), mean))
	   V_BH_PSU <- as.numeric(tapply(daily_BH_est$BH_day, list(DayType = as.factor(daily_BH_est$DayType)), var))
	   monthly_boat_hours <- mean_monthly_boat_hours * NumDays
	   
	   if (!fpc){
	     v_monthly_boat_hours <- NumDays^2  * V_BH_PSU / n_day  
	   }else{
	     v_monthly_boat_hours <- NumDays^2  * V_BH_PSU / n_day * (1-n_day/NumDays)
	   }
	   DailyBHTrue_daytype <- split(DailyBHTrue[daily_BH_est$Day,]$BH, as.factor(daily_BH_est$DayType));
     V_BD <- unlist(lapply(DailyBHTrue_daytype, var))
	   
	   # PSU term of two-stage sampling variance
	   v_psu_term <- (NumDays^2/n_day) * (1 - n_day/NumDays) * V_BH_PSU 
	   names(v_psu_term)=c("WD","WE")
	   
	   # Daily sampling error 
     ErrWD_BH_d_WD <- (psu_BH_daytype$WD - DailyBHTrue_daytype$WD)
     ErrWD_BH_d_WE <- (psu_BH_daytype$WE - DailyBHTrue_daytype$WE)
     ErrWD_BH_d <- list(WD = ErrWD_BH_d_WD, WE = ErrWD_BH_d_WE)

	   VBp_Boot <- NumDays^2 * (n_day/(n_day-1))* (1 - n_day/NumDays) * unlist(lapply(psu_BH_daytype, f_var_mean_bootstrap))
	   #VBp_Boot <- NumDays^2 * (n_day/(n_day-1)) * 2* (NumDays - n_day)/(NumDays-1) * unlist(lapply(psu_BH_daytype, f_var_mean_bootstrap))
	   #VBp_Boot <- NumDays^2 * (n_day/(n_day-1)) * (NumDays - n_day)/(NumDays-1) * unlist(lapply(psu_BH_daytype, f_var_mean_bootstrap))
	}
	else
	{
	   n_day <- nrow(daily_BH_est)
	   mean_monthly_boat_hours <- mean(daily_BH_est$BH_day)
	   psu_BH_daytype <- daily_BH_est$BH_day

	   v_mean_monthly_boat_hours <- var(daily_BH_est$BH_day)/n_day
	   monthly_boat_hours <- mean_monthly_boat_hours * sum(NumDays)
	   v_monthly_boat_hours <- v_mean_monthly_boat_hours * (sum(NumDays))^2
	   
	   VBp_Boot <- NumDays^2 * unlist(lapply(psu_BH_daytype, f_var_mean_bootstrap))
	   #VBp_Boot <- NumDays^2 * (n_day/(n_day-1)) * 2 * (NumDays - n_day)/(NumDays-1) * unlist(lapply(psu_BH_daytype, f_var_mean_bootstrap))
	   #VBp_Boot <- NumDays^2 * (n_day/(n_day-1)) * (NumDays - n_day)/(NumDays-1) * unlist(lapply(psu_BH_daytype, f_var_mean_bootstrap))
	   
	   # Daily sampling error^2 
     ErrWD_BH_d <- (psu_BH_daytype - DailyBHTrue[daily_BH_est$Day,]$BH)^2

     V_BD <- var(DailyBHTrue[daily_BH_est$Day,]$BH)
	}

	list(boat_hours = sum(monthly_boat_hours), 
	     v_boat_hours_boot = sum(VBp_Boot), 
	     v_boat_hours_approx = sum(v_monthly_boat_hours), 
	     daily_BH_est = psu_BH_daytype, 
	     ErrWD_BH_d=ErrWD_BH_d, 
	     V_BD = V_BD, 
	     v_psu_term = v_psu_term)
}

f_bus_route_total_angler_hours <- function(est_by_daytype, ints, boat_counts, total_route_time, 
                                             wait_times, NumDays, fpc)
{
  # ----------------------------------------------
  # monthly effort estimation  
  #   only for one month 
  # ----------------------------------------------
  
  # Convert time strings,such as "10:50 AM" to the POSIXct format
  #ints$st <- as.POSIXct(paste(ints$SDAY, " ", ints$STIME, sep=""), format = "%m/%d/%Y %I:%M %p") 
  #ints$et <- as.POSIXct(paste(ints$EDAY, " ", ints$ETIME, sep=""), format = "%m/%d/%Y %I:%M %p") 
  # Calculate trip lengths in hours
  ints$trip_length <- ints$ET - ints$ST # as.numeric(difftime(ints$ET, ints$ST, units = "hours"))
  
  # trip lengths must be greater than zero
  # trip lengths must be greater than zero
	if (any(ints$trip_length <= 0)){
	  #cat("trip_lengths are <= 0: ", "\n")
      #print(ints[ints$trip_length <= 0, ])

	  ints <- ints[ints$trip_length > 0, ]
	}
	#stopifnot(ints$trip_length >= 0)
	
	# Total party angler hours
	ints$party_ang_hours <- ints$ANGCNT * ints$trip_length 
	ints$party_ang_hours_w <- ints$party_ang_hours/ints$prob_sampling
	
	# angler hours
	party_ang_hours <- data.frame(Q = ints$party_ang_hours_w, AccessSite = ints$AccessSite, Day = ints$DAY)
	route_daily_site_AH_est <- bus_route_site_estimation(boat_counts, total_route_time, wait_times, party_ang_hours)
	
	daily_AH_est <- aggregate(list(AH_day = route_daily_site_AH_est$daily_site_Q),  list(Day = as.factor(route_daily_site_AH_est$Day)), sum)

  if (est_by_daytype){
      date_str <- paste(ints$MONTH[1], "/", daily_AH_est$Day,"/", ints$YEAR[1], sep="")
      
      daily_AH_est$weekday <- week_days(date_str)
      daily_AH_est$DayType <- ifelse(daily_AH_est$weekday<6, "WD", "WE")
      
      n_day <- as.numeric(tapply(daily_AH_est$AH_day, list(DayType = as.factor(daily_AH_est$DayType)), length))
      mean_monthly_angler_hours <- as.numeric(tapply(daily_AH_est$AH_day, list(DayType = as.factor(daily_AH_est$DayType)), mean))
      v_mean_monthly_angler_hours <- as.numeric(tapply(daily_AH_est$AH_day, list(DayType = as.factor(daily_AH_est$DayType)), var))
      v_mean_monthly_angler_hours <-  v_mean_monthly_angler_hours / n_day
      
      monthly_angler_hours <- mean_monthly_angler_hours * NumDays
      if (!fpc){
         v_monthly_angler_hours <- v_mean_monthly_angler_hours * NumDays^2
      }else{
         v_monthly_angler_hours <- v_mean_monthly_angler_hours * NumDays^2 * (1-n_day/NumDays)
      }
	}
	else
	{
	  n_day <- nrow(daily_AH_est)
	  monthly_angler_hours <- mean(daily_AH_est$AH_day) * sum(NumDays)
	  v_mean_monthly_angler_hours <- var(daily_AH_est$AH_day)/n_day
	  v_monthly_angler_hours <- v_mean_monthly_angler_hours * (sum(NumDays))^2
	}

	list(angler_hours = sum(monthly_angler_hours), v_angler_hours = sum(v_monthly_angler_hours))
}

f_bus_route_total_catch <- function(species, ints, boat_counts, total_route_time, wait_times, NumDays, fpc)
{
  # Total catch for a multiple day period
  
  #species = c("WAE","NOP","LMB","SMB","MUS","WHB","YEP","LWF","LHR","BLG","CCF","CWS","RWF","RKB","PSF","BKT","WHP","SPL","BCR","TMU","DRU","PKS","ATS","CAR","SMT","STN","OTH","COS_REL","CHS_REL","RBT_REL","BNT_REL","LAT_REL","FAT_REL","LAT.FAT.UNK_REL","WAE_REL","NOP_REL","LMB_REL","SMB_REL","MUS_REL","WHB_REL","YEP_REL","LWF_REL","LHR_REL","BLG_REL","CCF_REL","CWS_REL","RWF_REL","RKB_REL","PSF_REL","BKT_REL","WHP_REL","SPL_REL","BCR_REL","TMU_REL","DRU_REL","PKS_REL","ATS_REL","SMT_REL","STN_REL","OTH_REL","COS_NLEG","CHS_NLEG","RBT_NLEG","BNT_NLEG","LAT_NLEG","FAT_NLEG","LAT.FAT.UNK_NLEG","WAE_NLEG","NOP_NLEG","LMB_NLEG","SMB_NLEG","MUS_NLEG","WHB_NLEG")
  
  #species = c("YEP")
  
  monthly_C <- rep(0, length(species))
  names(monthly_C) <- species
  v_monthly_C <- rep(0, length(species))
  names(v_monthly_C) <- species
	i <- 1
	for (s in species) {
	  catch <- ints[, s]
      catch[is.na(catch)] <- 0
	  if (sum(catch) > 0){
	    weighted_catch <- catch/ints$prob_sampling
			
	    Cd_df <- data.frame(Q = weighted_catch, AccessSite = ints$AccessSite, Day = ints$DAY)
        route_daily_est <- bus_route_site_estimation(boat_counts, total_route_time, wait_times, Cd_df)

	    daily_Catch_est <- aggregate(list(C_day = route_daily_est$daily_site_Q),  list(Day = as.factor(route_daily_est$Day)), sum)	

      if (est_by_daytype){
          date_str <- paste(ints$MONTH[1], "/", daily_Catch_est$Day,"/", ints$YEAR[1], sep="")
          daily_Catch_est$weekday <- week_days(date_str)
          daily_Catch_est$DayType <- ifelse(daily_Catch_est$weekday<6, "WD", "WE")
          
          n_day <- as.numeric(tapply(daily_Catch_est$C_day, list(DayType = as.factor(daily_Catch_est$DayType)), length))
          mean_monthly_catch <- as.numeric(tapply(daily_Catch_est$C_day, list(DayType = as.factor(daily_Catch_est$DayType)), mean))
          v_mean_monthly_catch <- as.numeric(tapply(daily_Catch_est$C_day, list(DayType = as.factor(daily_Catch_est$DayType)), var))
          v_mean_monthly_catch <-  v_mean_monthly_catch / n_day
          
          monthly_C[i] <- sum(mean_monthly_catch * NumDays)
          if (!fpc){
             v_monthly_C[i] <- sum(v_mean_monthly_catch * NumDays^2)
          }else{
             v_monthly_C[i] <- sum(v_mean_monthly_catch * NumDays^2 * (1-n_day/NumDays))
          }
	    }
	    else
	    {		
	      mean_monthly_catch <- mean(daily_Catch_est$C_day)
	      n_day <- nrow(daily_Catch_est)
	      v_mean_monthly_catch <- var(daily_Catch_est$C_day)/n_day
	      
	      monthly_C[i] <- mean_monthly_catch * sum(NumDays)
	      v_monthly_C[i] <- v_mean_monthly_catch * (sum(NumDays))^2
	    }
	  }
	  i <- i + 1
    }

    data.frame(catch = monthly_C[monthly_C!=0], v_catch = v_monthly_C[monthly_C!=0])
}

f_bus_route_creel_estimation <- function(Species, est_by_daytype, interviews_repl, BoatCounts_repl, total_route_time, wait_times, 
                                         NumDays, DailyBHTrue, V_BD, fpc)
{
	if ((!est_by_daytype))
	{
		NumDays <- sum(NumDays)
	}
	
	n_repl <- length(interviews_repl)
	ints <- interviews_repl[[1]]
	nd <- length(unique(ints$DAY))

	daily_BH_est_WD <- vector(mode="list", length=n_repl) 
	daily_BH_est_WE <- vector(mode="list", length=n_repl) 
	
	ErrWD_BH_d_WD <- vector(mode="list", length=n_repl) 
	ErrWD_BH_d_WE <- vector(mode="list", length=n_repl) 
	
	v_psu_term_WD <- rep(0, n_repl) 
	v_psu_term_WE <-  rep(0, n_repl)
	
	boat_hours <- rep(0, n_repl) 
	v_boat_hours <-  rep(0, n_repl)
	angler_hours <-  rep(0, n_repl)
	v_angler_hours <-  rep(0, n_repl)
	catch <-  rep(0, n_repl)
	v_catch <-  rep(0, n_repl)

	for (r in 1:n_repl){
		interviews_repl[[r]]$DAYTYPE <- ifelse(interviews_repl[[r]]$DOW < 6, 1, 2)
	}
	
	for (r in 1:n_repl){
		# est_by_daytype is only applied to daily effort estimation, not catch rate estimation
		# ---------------------------------------------------
		# psu is (1) day for the daily estimator, and 
		#        (2) a period of multiple days of a daytype for the multiple-day estimator
		#        (3) a month for the multiple-day estimator if daytype is not considered  
		#psu_stats <- f_psu_stat(est_by_daytype, interviews, simu_Month)
	   
		tot_boat_hours <- f_bus_route_total_boat_hours(est_by_daytype, interviews_repl[[r]], 
		                                               BoatCounts_repl[[r]], total_route_time, wait_times, 
		                                               NumDays, DailyBHTrue, fpc)
		
		daily_BH_est_WD[[r]] <- tot_boat_hours$daily_BH_est$WD
		names(daily_BH_est_WD)[r] = paste("R", r, sep="")
		daily_BH_est_WE[[r]] <- tot_boat_hours$daily_BH_est$WE
		names(daily_BH_est_WE)[r] = paste("R", r, sep="")

		v_psu_term_WD[r] <- tot_boat_hours$v_psu_term[1] #$WD
		v_psu_term_WE[r] <- tot_boat_hours$v_psu_term[2] #$WE

		ErrWD_BH_d_WD[[r]] <- tot_boat_hours$ErrWD_BH_d$WD
        names(ErrWD_BH_d_WD)[r] = paste("R", r, sep="")
		ErrWD_BH_d_WE[[r]] <- tot_boat_hours$ErrWD_BH_d$WE
        names(ErrWD_BH_d_WE)[r] = paste("R", r, sep="")

		tot_angler_hours <- f_bus_route_total_angler_hours(est_by_daytype, interviews_repl[[r]], BoatCounts_repl[[r]], 
		                                                   total_route_time, wait_times, NumDays, fpc)	
		total_catch <- f_bus_route_total_catch(Species, interviews_repl[[r]], BoatCounts_repl[[r]], 
		                                        total_route_time, wait_times, NumDays, fpc)

    boat_hours[r] <- tot_boat_hours$boat_hours  
    v_boat_hours[r] <- tot_boat_hours$v_boat_hours_approx  
    angler_hours[r] <- tot_angler_hours$angler_hours
    v_angler_hours[r] <- tot_angler_hours$v_angler_hours
    catch[r] <- total_catch$catch[1] 
    v_catch[r] <- total_catch$v_catch[1]
  }

	# Two stage variance
	BH_d_WD <- t(as.data.frame(daily_BH_est_WD))
	v_bh_wd <- mean(v_psu_term_WD) + NumDays[1] * mean(apply(BH_d_WD, 2, var))

	BH_d_WE <- t(as.data.frame(daily_BH_est_WE))
	v_bh_we <- mean(v_psu_term_WE) + NumDays[2] * mean(apply(BH_d_WE, 2, var))

	V_BH_2stage <- v_bh_wd + v_bh_we
	
  # Mean of daily mean err_d^2
	# Within day variance for the stratum
	# daily BH error: Err_BH_d
  ErrWD_WD <- t(as.data.frame(ErrWD_BH_d_WD))
	#V_WD_WD = var(ErrWD_WD)
  ErrWD_WE <- t(as.data.frame(ErrWD_BH_d_WE))
	#V_WD_WE = var(ErrWD_WE)
	V_WD = var(c(as.numeric(ErrWD_WD), as.numeric(ErrWD_WE)))

  m_boat_hours = mean(boat_hours)
  m_v_boat_hours = mean(v_boat_hours)

  m_angler_hours = mean(angler_hours)
  m_v_angler_hours = mean(v_angler_hours)
	
  m_catch = mean(catch)
  m_v_catch = mean(v_catch)

  list(angler_hours = m_angler_hours, v_angler_hours = m_v_angler_hours, 
      boat_hours = m_boat_hours, v_boat_hours = m_v_boat_hours, 
      catch = m_catch, v_catch = m_v_catch, 
      V_BH_2stage = V_BH_2stage, V_WD = V_WD)
}



