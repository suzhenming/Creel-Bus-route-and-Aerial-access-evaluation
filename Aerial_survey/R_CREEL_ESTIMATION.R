#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Daily and Multiple-day estimation
#    for Michigan Great Lakes
#    creel surveys
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dr. Zhenming Su
# Institute for Fisheries Research 
# Michigan Department of Natural Resources 
#      and University of Michigan  
# ANN ARBOR, MI 48100

# Contact Zhenming Su (suz@michigan.gov) for any questions
# Revised for "the Evaluation of Bus-route and Aerial-access Methods for Great Lakes Recreational Fisheries Surveys"
#  by Zhenming Su, 2/9/2024 

# Software built for the research paper:
# Zhenming Su & David Clapp (2013): Evaluation of Sample Design and Estimation Methods for Great Lakes
#   Angler Surveys, Transactions of the American Fisheries Society, 142:1, 234-246
# To link to this article: http://dx.doi.org/10.1080/00028487.2012.728167

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Built: Oct 2, 2009
# Revised: Aug 8, 2021
#   Daytype estimates for catch rate and catch were added.
#   Provided same estimates as MiCreel 
#
# Added calculation for incomp trip interviews 0n Aug 28, 2021
# Aug 25: fixed problems with numang 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


f_approximateBayesianBootstrap <- function(nc, yobs, v_yobs)
{
	# ------------------------------------------------------------------------------------
	# ABB: 

	# ------------------------------------------------------------------------------------
	# (1) first draws nobs values randomly with replacement from yobs
	# ------------------------------------------------------------------------------------
	nobs <- length(yobs)

	ind_star <- sample(1:nobs, size = nobs, replace=TRUE)

	yobs_star <- yobs[ind_star]
	v_yobs_star <- v_yobs[ind_star]

	# ------------------------------------------------------------------------------------
    # (2) then draw nmis = n - nobs components of ymis  
	#   randomly with replacement from yobs_star
	# ------------------------------------------------------------------------------------
	nmis = nc - nobs

	ind_mis <- sample(1:nobs, size = nmis, replace=TRUE)

	ymis <- yobs_star[ind_mis]
	v_ymis <- v_yobs_star[ind_mis]	

	list(ymis = ymis, v_ymis = v_ymis)
}

f_interview_exclude_days <- function(exclude_method, n_days_excluded, interviews, month)
{
  # ------------------------------------------------------------------------------------
  # Exclude n_days_excluded days of interviews
  # ------------------------------------------------------------------------------------
  sampled_days <- interviews$DAY
  
  ints_day <- sort(unique(sampled_days))
  n_days_ints <- length(ints_day) 
  nd_kept <- n_days_ints - n_days_excluded
  switch(
    exclude_method,
    
    # randomly exclude n_days_excluded of days
    ints_day <- sample(ints_day, size = nd_kept),
    
    # exclude first nd_excl days
    ints_day <- ints_day[-n_days_excluded],
    
    # exclude last nd_excl days
    ints_day <- ints_day[1:nd_kept]
  )	
  
  # get the interviews for kept days
  int <- interviews[interviews$DAY %in% ints_day, ]
}

f_catch_rate_ration_of_means <- function(CdiHdi)
{
  # ---------------------------------------------------------
  # Ratio-of-means method is used for estimating 
  #   catch rate and its estimated variance 
  #   for completed trips
  # ---------------------------------------------------------
  Cdi <- CdiHdi[,1]  # party catches
  Hdi <- CdiHdi[,2]  # party angler-hours
  # ratio of means
  Rd  <- mean(Cdi)/mean(Hdi)
  
  SSQ_C   <- crossprod(Cdi); #sum of squares
  SSQ_H   <- crossprod(Hdi);
  SPD_C_H <- crossprod(Cdi, Hdi);
  
  # Variance of Rd
  n_int <- length(Hdi);
  mHd <- mean(Hdi)
  term1 <- 1/(n_int * mHd * mHd);
  term2 <- 2 * Rd * SPD_C_H;
  term3 <- Rd * Rd * SSQ_H;
  
  # Eq. 2
  V_R <- term1 * (SSQ_C - term2 + term3)/(n_int - 1);
  
  c(Rd, V_R)
}

f_catch_rate_mean_of_ratios <- function(CdiHdi)
{
  # ---------------------------------------------------------
  # Mean of Ratios method is used for estimating 
  #   catch rate and its estimated variance 
  #   for im-completed trips	
  # ---------------------------------------------------------
  Cdi <- CdiHdi[,1]
  Hdi <- CdiHdi[,2]
  
  Ratios <- Cdi/Hdi
  
  # mean of ratios 
  Rd  <- mean(Ratios)
  
  SUM_Ratios <- sum(Ratios); #sum of ratios
  SSQ_Ratios <- sum(Ratios*Ratios);
  
  # variance of Rd
  k_d <- length(Hdi);
  
  # Eq. xx
  V_R <- (SSQ_Ratios - SUM_Ratios * SUM_Ratios / k_d) /(k_d * (k_d - 1));
  
  c(Rd, V_R)
}

f_psu_interview_stat <- function(Species, use_daily_estimator, est_by_daytype, 
                                 exclude_method, n_days_excluded=0, interviews, month)
{
  # ---------------------------------------------------
  # A psu is a day for the daily estimator, and 
  #          a multiple-day period of  
  #            for the multiple-day estimator 
  
  # ----------------------------------------------
  # psu catch rate estimation and 
  #   psu level summary information
  # ----------------------------------------------
  make_MultiDay_estimates = !use_daily_estimator
  
  # Testing several methods for missing interview data in some days
  ints <- f_interview_exclude_days(exclude_method, n_days_excluded, interviews, month)
  n_ints <- nrow(ints)
  
  # Daily estimator
  psu_ints <- ints$DAY
  
  # Multiple-day estimator
  if (make_MultiDay_estimates)
  {
    if (est_by_daytype)
      # daytype as psu
      psu_ints <- ifelse(ints$DOW < 6, 1, 2)
    else
      # one psu
      psu_ints <- rep(1, n_ints)
  }
  ints$psu <- psu_ints
  
  # Calculate psu summaries 
  #   get unique rows containing the following columns:
  #   unique(ints[, c("Waterbody", "FSITE", "YEAR", "MONTH",  "psu")])
  psu_summ = ints[!duplicated(ints[,c("Waterbody",  "YEAR", "MONTH",  "psu")]),]	
  psu_summ = psu_summ[,c("Waterbody",  "YEAR", "MONTH", "DOW", "DAY", "psu")]
  
  # Numbers of interviews per psu (e.g., day or daytype)
  psu_summ$NINT    <- as.vector(table(psu_ints))
  psu_summ$DAYTYPE <- ifelse(psu_summ$DOW < 6, 1, 2)	
  if ((make_MultiDay_estimates) & (!est_by_daytype)){
    psu_summ$DAYTYPE <- 0
  }
  
  # Numbers of sampled days
  psus_DAYs = ints[!duplicated(ints[,c("Waterbody",  "YEAR", "MONTH", "DAY",  "psu")]), ]
  psu_summ$NDAYS <- tapply(psus_DAYs$DAY, list(as.factor(psus_DAYs$psu)), length)
  
  # Number of completed trips or incompleted trips per psu
  n_trips_psu_comptrip = tapply(ints$CTRIP, list(CTRIP = as.factor(ints$CTRIP), psu = as.factor(ints$psu)), length)
  
  if(length(n_trips_psu_comptrip[rownames(n_trips_psu_comptrip)=="1",])>0){
    psu_summ$NCOMPTRIPS <- n_trips_psu_comptrip[rownames(n_trips_psu_comptrip)=="1",]
  } else {
    psu_summ$NCOMPTRIPS <- rep(0, length(n_trips_psu_comptrip))
  }
  if(length(n_trips_psu_comptrip[rownames(n_trips_psu_comptrip)=="2",])>0){
    psu_summ$NINCOMPTRIPS <- n_trips_psu_comptrip[rownames(n_trips_psu_comptrip)=="2",]
  } else {
    psu_summ$NINCOMPTRIPS <- rep(0, length(n_trips_psu_comptrip))
  }
  # % completed trips
  psu_summ$PercentComptrips <- psu_summ$NCOMPTRIPS/(psu_summ$NCOMPTRIPS+psu_summ$NINCOMPTRIPS)
  
  # TripHours
  st <- as.POSIXct(paste(ints$SDAY, " ", ints$STIME, sep=""),format = "%Y-%m-%d %H:%M")  # format = "%m/%d/%Y %H:%M") 
  et <- as.POSIXct(paste(ints$EDAY, " ", ints$ETIME, sep=""), format = "%Y-%m-%d %H:%M") 
  
  TRIPHOURS <- difftime(et, st, units = "hours")
  ints$TRIPHOURS <- TRIPHOURS
  
  # Trip Angler-hours
  #ints$ANGCNT <- ifelse(ints$individualInterview == "Yes", 1, ints$ANGCNT)
  ints$party_ang_hours <- ints$ANGCNT * ints$TRIPHOURS 
  
  n_psu <- nrow(psu_summ)
  psu_f <- as.factor(ints$psu)
  
  meanTripHours <- tapply(TRIPHOURS, list(psu_f), mean)
  v_TripHours   <- tapply(TRIPHOURS, list(psu_f), var)/psu_summ$NINT
  psu_summ$TripHours <- meanTripHours
  psu_summ$V_TripHours <- v_TripHours
  
  # Mean party size: anglers/party
  if (mean(psu_summ$PercentComptrips) >= 0.8){
    # completed trips
    mean_party_size <- tapply(ints$ANGCNT, list(psu_f), mean, na.rm = TRUE)
    v_party_size    <- tapply(ints$ANGCNT, list(psu_f), var, na.rm = TRUE)/psu_summ$NINT
    PartyAngHours   <- ints$ANGCNT * ints$TRIPHOURS
    meanPartyAngHours <- tapply(PartyAngHours, list(psu_f), mean)
    v_PartyAngHours   <- tapply(PartyAngHours, list(psu_f), var)/psu_summ$NINT
    
  } else {
    # imcompleted trips
    NumAng_psu <- unique(ints[, c("DAY", "party_no", "psu", "ANGCNT", "TRIPHOURS")])
    mean_party_size <- tapply(NumAng_psu$ANGCNT, list(NumAng_psu$psu), mean, na.rm = TRUE)
    v_party_size    <- tapply(NumAng_psu$ANGCNT, list(NumAng_psu$psu), var, na.rm = TRUE)/psu_summ$NINT	   
    PartyAngHours <- NumAng_psu$ANGCNT * NumAng_psu$TRIPHOURS
    meanPartyAngHours <- tapply(PartyAngHours, list(NumAng_psu$psu), mean)
    v_PartyAngHours   <- tapply(PartyAngHours, list(NumAng_psu$psu), var)/psu_summ$NINT
  }
  
  # Mean party size
  psu_summ$PARTY_SIZE <- mean_party_size
  psu_summ$V_PARTY_SIZE <- v_party_size
  
  psu_summ$PartyAngHours <- meanPartyAngHours
  psu_summ$V_PartyAngHours <- v_PartyAngHours
  
  # psu catch rates 
  if (Species == "WAE"){
    Cdi <- ints$WAE 
  }else{
    Cdi <- ints$YEP
  }
  
  Hdi <- as.numeric(ints$party_ang_hours)
  
  C_H <- data.frame(C = Cdi, H = Hdi, psu = psu_f)
  C_H <- C_H[C_H$H > 0,]
  if (mean(psu_summ$PercentComptrips) >= 0.8)
    catch_rate <- by(C_H, list(C_H$psu), f_catch_rate_ration_of_means, simplify = TRUE)
  else
    catch_rate <- by(C_H, list(C_H$psu), f_catch_rate_mean_of_ratios, simplify = TRUE)
  
  Rd <- matrix(NA, nrow = n_psu, byrow = T, ncol =2)
  for (i in 1:n_psu) 
    Rd[i,] <- catch_rate[[i]] 
  
  psu_summ$catch_rate <- Rd[,1]
  psu_summ$v_catch_rate <- Rd[,2]
  
  psu_summ <- na.omit(psu_summ)
}


f_psu_boat_instant_hours <- function(use_daily_estimator, est_by_daytype, counts, month, NumDays, FValue, FreqPropNight)
{
  # ---------------------------------------------------
  # A psu is (1) a day for the daily estimator, and 
  #          (2) a multiple-day period of a daytype 
  #                for the multiple-day estimator 
  # 
  # Daily boat-hours estimates from instant counts
  # ---------------------------------------------------
  make_Multiple_day_estimates = !use_daily_estimator
  
  sampled_days <- counts$DAY
  days_cnt <- sort(unique(sampled_days))
  
  ndays_cnt <- length(days_cnt)
  s_months <- rep(month, ndays_cnt)
  sampled_days <- as.factor(sampled_days)
  
  wkday <- tapply(counts$DOW, list(sampled_days), mean)
  daytype <- ifelse(wkday < 6, 1, 2)
  
  # Calculate boat-hours from instant counts based on fishable-hours
  boat_hours_cnts <- counts$COUNT * FValue / (1 - FreqPropNight)
  psu_boat_hours <- tapply(boat_hours_cnts, list(sampled_days), mean)
  
  if (make_Multiple_day_estimates)
  {
    # Multiple-day estimates
    if (est_by_daytype)
    {
      psu_boat_hours[daytype==1] <- psu_boat_hours[daytype==1] * NumDays[1]
      psu_boat_hours[daytype==2] <- psu_boat_hours[daytype==2] * NumDays[2]
    } else
    {
      psu_boat_hours <- psu_boat_hours * NumDays
    }
  }
  
  # Variance of psu BH estimates
  n_cnts_daily <- tapply(boat_hours_cnts, list(sampled_days), length)
  v_psu_boat_hours <- tapply(boat_hours_cnts, list(sampled_days), var)/n_cnts_daily
  
  if (make_Multiple_day_estimates)
  {
    # Variance of Multiple-day estimates
    if (est_by_daytype)
    {
      v_psu_boat_hours[daytype==1] <- v_psu_boat_hours[daytype==1] * (NumDays[1]^2)
      v_psu_boat_hours[daytype==2] <- v_psu_boat_hours[daytype==2] * NumDays[2]^2
    } else
    {
      v_psu_boat_hours <- v_psu_boat_hours * NumDays^2
    }
  }
  
  BH <- data.frame(DAY = days_cnt, MONTH = s_months, DAYTYPE = daytype, NCNT = n_cnts_daily, 
                   PARTY_EFFORT = psu_boat_hours, V_PARTY_EFFORT = v_psu_boat_hours)
  
  if (use_daily_estimator & (trunc(mean(n_cnts_daily))!= 1))
  {
    BH <- na.omit(BH)
  }
  else
    BH
}

f_psu_boat_hours_Aerial_Proportion <- function(use_daily_estimator, est_by_daytype, Freq_Boats_Present, counts, month, NumDays, FValue)
{
  # ---------------------------------------------------
  # A psu is a day for the daily estimator, and 
  #      multiple days of a daytype 
  #       for the multiple-day estimator 
  # 
  # Daily boat-hours estimates from instant counts
  #   from proportion method, Lockwood 2001 
  # ---------------------------------------------------
  make_Multiple_day_estimates = !use_daily_estimator
  
  sampled_days <- counts$DAY
  days_cnt <- sort(unique(sampled_days))
  
  ndays_cnt <- length(days_cnt)
  s_months <- rep(month, ndays_cnt)
  sampled_days <- as.factor(sampled_days)
  
  wkday <- tapply(counts$DOW, list(sampled_days), mean)
  daytype <- ifelse(wkday < 6, 1, 2)
  
  # Calculate boat-hours from instant counts
  # Expansion values based on hourly proportion of trip-hours
  exp_bh = 1/Freq_Boats_Present[floor(counts$BEGHR)]
  
  # Expanding to daily boat-hours
  boat_hours_cnts <- counts$COUNT * exp_bh 
  
  # Mean daily boat-hours
  psu_boat_hours <- tapply(boat_hours_cnts, list(sampled_days), mean)
  
  if (make_Multiple_day_estimates)
  {
    # Multiple-day estimates
    if (est_by_daytype)
    {
      psu_boat_hours[daytype==1] <- psu_boat_hours[daytype==1] * NumDays[1]
      psu_boat_hours[daytype==2] <- psu_boat_hours[daytype==2] * NumDays[2]
    } else
    {
      psu_boat_hours <- psu_boat_hours * NumDays
    }
  }
  
  n_cnts_daily <- tapply(boat_hours_cnts, list(sampled_days), length)
  v_psu_boat_hours <- tapply(boat_hours_cnts, list(sampled_days), var)/n_cnts_daily
  if (make_Multiple_day_estimates)
  {
    # Multiple-day estimates
    if (est_by_daytype)
    {
      v_psu_boat_hours[daytype==1] <- v_psu_boat_hours[daytype==1] * (NumDays[1]^2)
      v_psu_boat_hours[daytype==2] <- v_psu_boat_hours[daytype==2] * NumDays[2]^2
    } else
    {
      v_psu_boat_hours <- v_psu_boat_hours * NumDays^2
    }
  }
  
  BH <- data.frame(DAY = days_cnt, MONTH = s_months, DAYTYPE = daytype, NCNT = n_cnts_daily, 
                   PARTY_EFFORT = psu_boat_hours,  
                   V_PARTY_EFFORT = v_psu_boat_hours)
  
  if (use_daily_estimator & (trunc(mean(n_cnts_daily))!= 1))
  {
    BH <- na.omit(BH)
  }
  else
    BH
}

#f_daily_multiple_day_period_total_BH <- function(psu_boat_hours, v_psu_boat_hours, NumDays)
#{
#   # Number of days in the count sample
#   nSDays <- length(psu_boat_hours)
#
#   Bp <- NumDays * mean(psu_boat_hours)  # boat hours
#   VBp <- ((NumDays^2) * (1-nSDays/NumDays) * var(psu_boat_hours)/nSDays +(NumDays/nSDays) * sum(v_psu_boat_hours));
#
#   list(Bp = Bp, VBp = VBp)	
#}

f_var_angler_hours <- function(x)
{
  BoatHours <- x[1]
  VBoatHours <- x[2]
  MeanAngCnts <- x[3]
  VMeanAngCnts <- x[4]
  
  VAnglerHours <- (BoatHours^2 * VMeanAngCnts + MeanAngCnts^2 * VBoatHours - VMeanAngCnts * VBoatHours); #Eq. 24
  VAnglerHours
}

f_dailyEstimator_impute <- function(imputation_method, year, month, 
                                    int_days, cnt_days, 
                                    psu_boat_hours_d, v_psu_boat_hours_d, 
                                    party_size_d, v_party_size_d)
{
  # ------------------------------------------------------------------------------------
  # Expand sampled days to all those days with interviews and counts
  #  and fill in data for missing days
  # Only used for the Daily-estimator
  # ------------------------------------------------------------------------------------
  
  # Merge sampled days in count and interview
  all_sampled_days <- sort(unique(c(int_days, cnt_days)))
  ns_days <- length(all_sampled_days)
  
  wkdys <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  dates <- paste(year, "-", month, "-", all_sampled_days, sep="")
  wkday <- match(weekdays(as.Date(dates), abbreviate = T), wkdys)	
  daytype <- ifelse(wkday < 6, 1, 2)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # For party effort
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ------------------------------------------------------------------------------------
  # Fillin effort for missing days with one of Imputation_Methods
  # ------------------------------------------------------------------------------------
  #Imputation_Methods <- c("complete_case", "mean", "approximate Bayesian bootstrap")
  switch(imputation_method,
         # complete-case
         { 
           psu_boat_hours_d_ex <- rep(NA, ns_days)
           v_psu_boat_hours_d_ex <- rep(NA, ns_days)
           
           # copy the party effort of non-missing days to the expanded vector
           psu_boat_hours_d_ex[all_sampled_days %in% cnt_days] <- psu_boat_hours_d
           v_psu_boat_hours_d_ex[all_sampled_days %in% cnt_days] <- v_psu_boat_hours_d
         },  
         # mean impute
         {
           mean_effort  <- mean(psu_boat_hours_d)
           v_effort_obs <- var(psu_boat_hours_d)
           #v_effort_obs <- mean(v_psu_boat_hours_d)
           psu_boat_hours_d_ex <- rep(mean_effort, ns_days)
           v_psu_boat_hours_d_ex <- rep(v_effort_obs, ns_days)
           
           # copy the party effort of non-missing days to the expanded vector
           psu_boat_hours_d_ex[all_sampled_days %in% cnt_days] <- psu_boat_hours_d
           v_psu_boat_hours_d_ex[all_sampled_days %in% cnt_days] <- v_psu_boat_hours_d
           
         },  
         # approximate Bayesian bootstrap
         {
           psu_boat_hours_d_ex <- rep(NA, ns_days)
           v_psu_boat_hours_d_ex <- rep(NA, ns_days)
           
           psu_boat_hours_d_mis <- f_approximateBayesianBootstrap(ns_days, psu_boat_hours_d, v_psu_boat_hours_d)
           
           # observed days
           psu_boat_hours_d_ex[all_sampled_days %in% cnt_days] <- psu_boat_hours_d
           v_psu_boat_hours_d_ex[all_sampled_days %in% cnt_days] <- v_psu_boat_hours_d
           if (length(psu_boat_hours_d_mis$ymis) > 0)
           {
             psu_boat_hours_d_ex[is.na(psu_boat_hours_d_ex)] <- psu_boat_hours_d_mis$ymis
             v_psu_boat_hours_d_ex[is.na(v_psu_boat_hours_d_ex)] <- psu_boat_hours_d_mis$v_ymis
           }
         }
  )    
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # For party_size
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ------------------------------------------------------------------------------------
  # Fillin party_size of missing days with mean party_size	
  # ------------------------------------------------------------------------------------
  #Imputation_Methods <- c("complete_case", "mean", "approximate Bayesian bootstrap")
  switch(imputation_method,
         # complete-case
         {
           party_size_ex <- rep(NA, ns_days)
           v_party_size_ex <- rep(NA, ns_days)
           
           # copy the party_size of non-missing days to the expanded vector
           party_size_ex[all_sampled_days %in% int_days] <- party_size_d
           v_party_size_ex[all_sampled_days %in% int_days] <- v_party_size_d
           
         },  
         # mean impute
         {
           mean_party_size <- mean(party_size_d)
           v_mean_party_size <- var(party_size_d)
           #v_mean_party_size <-  mean(v_party_size_d)
           party_size_ex <- rep(mean_party_size, ns_days)
           v_party_size_ex <- rep(v_mean_party_size, ns_days)
           
           # copy the party effort of non-missing days to the expanded vector
           party_size_ex[all_sampled_days %in% int_days] <- party_size_d
           v_party_size_ex[all_sampled_days %in% int_days] <- v_party_size_d
         },
         # approximate Bayesian bootstrap  
         {
           party_size_ex <- rep(NA, ns_days)
           v_party_size_ex <- rep(NA, ns_days)
           
           party_size_mis <- f_approximateBayesianBootstrap(ns_days, party_size_d, v_party_size_d)
           party_size_ex[all_sampled_days %in% int_days] <- party_size_d
           v_party_size_ex[all_sampled_days %in% int_days] <- v_party_size_d
           
           if (length(party_size_mis$ymis) > 0)
           {
             party_size_ex[is.na(party_size_ex)] <- party_size_mis$ymis
             v_party_size_ex[is.na(v_party_size_ex)] <- party_size_mis$v_ymis
           }
         }
  )
  
  angler_effort <- party_size_ex * psu_boat_hours_d_ex
  
  df <- data.frame(B = psu_boat_hours_d_ex, VB = v_psu_boat_hours_d_ex, A = party_size_ex, VA = v_party_size_ex)
  v_angler_effort <- by(df, list(all_sampled_days), f_var_angler_hours, simplify = TRUE)
  
  V_E <- matrix(NA, nrow = ns_days, byrow = T, ncol = 1)
  for (i in 1:ns_days) V_E[i] <- v_angler_effort[[i]] 
  
  s_months <- rep(month, ns_days)
  
  data.frame(DAY = all_sampled_days, MONTH = s_months, DAYTYPE = daytype, 
             PARTY_SIZE = party_size_ex, V_PARTY_SIZE = v_party_size_ex,
             PARTY_EFFORT = psu_boat_hours_d_ex, V_PARTY_EFFORT = v_psu_boat_hours_d_ex, 
             angler_effort = angler_effort, v_angler_effort = V_E)
}

var_mean_bootstrap <- function(x){
	B = 1000
	n <- length(x)
	x.mean <- mean(x)
	x.boot <- matrix(sample(x, B*n, repl=TRUE), nrow = B)
	mean.boot <- rowMeans(x.boot)
	bias <- mean(mean.boot - x.mean)

	var(mean.boot)
}


f_dailyEstimator_total_effort <- function(daytype, ncnts_per_day, imputation_method, psu_boat_hours, v_psu_boat_hours, angler_effort, v_angler_effort, NumDays)
{
  # total effort in period
  if (imputation_method == 1) # for complete-case
  {
    ind <- which(is.na(psu_boat_hours))
    psu_boat_hours[ind]   <- 0
    v_psu_boat_hours[ind] <- 0	  
  }
  
  psu_BH_daytype <- split(psu_boat_hours, daytype)
  v_psu_BH_daytype <- split(v_psu_boat_hours, daytype)
  
  nSDays <- unlist(lapply(psu_BH_daytype, length))
  
  Bp <- NumDays * unlist(lapply(psu_BH_daytype, mean))  # boat hours
  
  if (ncnts_per_day > 1)
  {
    vp <-  unlist(lapply(psu_BH_daytype, var))
    sv <-  unlist(lapply(v_psu_BH_daytype, sum))
    VBp <- ((NumDays^2) * (1-nSDays/NumDays) * vp/nSDays +(NumDays/nSDays) * sv)
  }
  else
  {	   
    VBp <- (NumDays^2/nSDays) * unlist(lapply(psu_BH_daytype, var))
    #VBp_Boot <- NumDays^2 * unlist(lapply(psu_BH_daytype, var_mean_bootstrap))
  }
  
  if (imputation_method == 1) # for complete-case 
  {
    ind <- which(is.na(angler_effort))
    angler_effort[ind] <- 0
    v_angler_effort[ind] <- 0
  }
  
  a_eff <- split(angler_effort, daytype)
  v_a_eff <- split(v_angler_effort, daytype)
  
  nSDays <- unlist(lapply(a_eff, length))
  
  Ep <- NumDays * unlist(lapply(a_eff, mean))
  if (ncnts_per_day > 1)
    #VEp <- (NumDays^2 * (1-nSDays/NumDays) * var(angler_effort)/nSDays + (NumDays/nSDays)*sum(v_angler_effort))
  {
    vp <-  unlist(lapply(a_eff, var))
    sv <-  unlist(lapply(v_a_eff, sum))
    VEp <- ((NumDays^2) * (1-nSDays/NumDays) * vp/nSDays +(NumDays/nSDays) * sv)
  }
  else {
    VEp <- (NumDays^2/nSDays) * unlist(lapply(a_eff, var))
    #VEp_Boot <-  NumDays^2 * unlist(lapply(a_eff, var_mean_bootstrap))	
  }
  
  list(Bp = sum(Bp), VBp = sum(VBp), Ep = sum(Ep), VEp = sum(VEp))	
}


f_multipleDayEst_total_effort <- function(ncnts_per_day, psu_boat_hours, v_psu_boat_hours, meanAngCnt, VmeanAngCnt, NumDays)
{
  # psu is month here
  psu_boat_hours <- na.omit(psu_boat_hours)
  
  nSDays <- length(psu_boat_hours)
  
  Bp <- mean(psu_boat_hours) # boat hours
  
  if (ncnts_per_day > 1)
  {
    VBp <- ((1-nSDays/NumDays) * var(psu_boat_hours)/nSDays +(1/(nSDays*NumDays)) * sum(v_psu_boat_hours, na.rm = TRUE))
  }
  else
  { 
    # mjack <- apply(t(1:nSDays), 2, function(i) mean(psu_boat_hours[-i]))
    # VBp   <- (nSDays - 1) * crossprod((mjack - rep(mean(psu_boat_hours),nSDays)))/ nSDays
    
    VBp <- var(psu_boat_hours)/nSDays
  }
  
  Ep <- meanAngCnt * Bp
  x <- c(Bp, VBp, meanAngCnt, VmeanAngCnt)
  VEp <- f_var_angler_hours(x)
  
  list(Bp = Bp, VBp = VBp, Ep = Ep, VEp = VEp)	
}


f_multipleDayEst_total_effort_daytype <- function(daytype, ncnts_per_day, psu_boat_hours, v_psu_boat_hours, 
                                                  meanAngCnt, VmeanAngCnt, NumDays)
{
  # psu can contain days belong to a daytype (weekdays or weekend days)
  psu_boat_hours <- na.omit(psu_boat_hours)
  psu_BH_daytype <- split(psu_boat_hours, daytype)
  v_psu_BH_daytype <- split(v_psu_boat_hours, daytype)
  
  nSDays <- unlist(lapply(psu_BH_daytype, length))
  Bp <- unlist(lapply(psu_BH_daytype, mean))  # boat hours
  if (ncnts_per_day > 1)
  {
    #  VBp <- ((1-nSDays/NumDays) * var(psu_boat_hours)/nSDays +(1/(nSDays*NumDays)) * sum(v_psu_boat_hours, na.rm = TRUE))
    vp <-  unlist(lapply(psu_BH_daytype, var))
    sv <-  unlist(lapply(v_psu_BH_daytype, sum))
    VBp <- ((1-nSDays/NumDays) * vp /nSDays +(1/(nSDays*NumDays)) * sv)
  } else  { 
    # mjack <- apply(t(1:nSDays), 2, function(i) mean(psu_boat_hours[-i]))
    # VBp   <- (nSDays - 1) * crossprod((mjack - rep(mean(psu_boat_hours),nSDays)))/ nSDays
    
    VBp <- var(psu_boat_hours)/nSDays
  }
  
  # Angler-hours
  Ep <- meanAngCnt * Bp
  
  VEp <- c(0,0)
  x <- c(Bp[1], VBp[1], meanAngCnt[1], VmeanAngCnt[1])
  VEp[1] <- f_var_angler_hours(x)
  x <- c(Bp[2], VBp[2], meanAngCnt[2], VmeanAngCnt[2])
  VEp[2] <- f_var_angler_hours(x)
  
  list(Bp = (Bp), VBp = (VBp), Ep = (Ep), VEp = (VEp))	
}

f_dailyEstimator_impute_catch_rate <- function(imputation_method, n_impute, 
                                               interview_days, all_days_sampled, catch_rate, v_catch_rate)
{
  # Fill in catch rates for days without interviews
  
  # ------------------------------------------------------------------------------------
  # Fill-in catch_rate of missing days with one of Imputation_Methods	
  #   Imputation_Methods <- c("complete_case", "mean", "approximate Bayesian bootstrap")
  # ------------------------------------------------------------------------------------
  ns_days <- length(all_days_sampled)
  
  # ------------------------------------------------------------------------------------
  switch(imputation_method,
         # complete-case
         {
           catch_rate_ex <- rep(NA, ns_days)
           v_catch_rate_ex <- rep(NA, ns_days)
           
           catch_rate_ex[all_days_sampled %in% interview_days] <- catch_rate
           
           v_catch_rate_ex[all_days_sampled %in% interview_days] <- v_catch_rate
         },  
         # mean impute
         {
           mean_catch_rate <- mean(catch_rate)
           v_catch_rate <- var(catch_rate)
           #v_catch_rate <- mean(v_catch_rate)
           
           catch_rate_ex <- rep(mean_catch_rate ,ns_days)
           v_catch_rate_ex <- rep(v_catch_rate ,ns_days)
           
           catch_rate_ex[all_days_sampled %in% interview_days] <- catch_rate
           v_catch_rate_ex[all_days_sampled %in% interview_days] <- v_catch_rate
           
         },
         # approximate Bayesian bootstrap
         {
           catch_rate_ex <- rep(NA, ns_days)
           v_catch_rate_ex <- rep(NA, ns_days)
           
           catch_rate_mis <- f_approximateBayesianBootstrap(ns_days, catch_rate, v_catch_rate)
           catch_rate_ex[all_days_sampled %in% interview_days] <- catch_rate
           v_catch_rate_ex[all_days_sampled %in% interview_days] <- v_catch_rate
           
           catch_rate_ex[is.na(catch_rate_ex)] <- catch_rate_mis$ymis
           v_catch_rate_ex[is.na(v_catch_rate_ex)] <- catch_rate_mis$v_ymis
         }
  )
  
  list(catch_rate = catch_rate_ex, v_catch_rate = v_catch_rate_ex)
}

f_dailyEstimator_Catch <-  function(daily_angler_effort, v_daily_angler_effort, daily_catch_rate, v_daily_catch_rate)
{
  
  daily_catch <- daily_angler_effort * daily_catch_rate

  v_daily_catch <- (daily_catch_rate^2 * v_daily_angler_effort 
              + daily_angler_effort^2 * v_daily_catch_rate 
              - v_daily_catch_rate * v_daily_angler_effort); # eq. 75

  list(daily_catch = daily_catch, v_daily_catch = v_daily_catch)
}

f_dailyEstimator_totalCatch <- function(ncnts_per_day, imputation_method, dailycatch, v_dailycatch, NumDays)
{
  # Total catch for a stratum: multiple day period
  if (imputation_method == 1) #only for complete case
  {
    ind <- which(is.na(v_dailycatch))
    v_dailycatch[ind] <- 0	  	   
    
    ind <- which(is.na(dailycatch))
    dailycatch[ind] <- 0	  	   
  }
  
  df <- data.frame(C = dailycatch, VC = v_dailycatch)
  
  # Number of sampled days
  nSDays <- length(df$C)
  
  # Stratum catch estimator and its estimated variance
  catch <- NumDays * mean(df$C) #, na.rm =TRUE)
  if (ncnts_per_day > 1)
    VC <- (NumDays^2 * (1 - nSDays / NumDays) * var(df$C)/nSDays + (NumDays/nSDays) * sum(df$VC)) # Eq. 48
  else
  {
    VC <- NumDays^2 * var(df$C)/nSDays
  }
  
  c(catch, VC)	
}

f_dailyEstimator_totalCatch_Daytype <- function(daytype, ncnts_per_day, imputation_method, dailycatch, v_dailycatch, NumDays)
{
  # total catch in period
  if (imputation_method == 1) #only for complete case
  {
    ind <- which(is.na(v_dailycatch))
    v_dailycatch[ind] <- 0	  	   
    
    ind <- which(is.na(dailycatch))
    dailycatch[ind] <- 0	  	   
  }
  
  C_daytype <- split(dailycatch, daytype)
  V_C_daytype <- split(v_dailycatch, daytype)
  
  # Number of sampled days
  nSDays <- unlist(lapply(C_daytype, length))
  
  # Stratum catch estimator and its estimated variance
  catch <- NumDays * unlist(lapply(C_daytype, mean))  # catch by daytype
  
  if (ncnts_per_day > 1)
  {
    vcatch <-  unlist(lapply(C_daytype, var))
    sv <-  unlist(lapply(V_C_daytype, sum))
    VC <- ((NumDays^2) * (1-nSDays/NumDays) * vcatch/nSDays +(NumDays/nSDays) * sv)
  }
  else
  {	   
    VC <- NumDays * NumDays * unlist(lapply(C_daytype, var))/nSDays
  }

  c(sum(catch), sum(VC))	
}

f_multipleDayEst_total_catch <- function(est_by_daytype, Ep, VEp, R, VR)
{
  # total catch in period
  Cp <- Ep * R
  df = data.frame(Ep, VEp, R, VR)
  VC <- f_var_angler_hours(df)
  as.data.frame(cbind(Cp = Cp, VC=VC$Ep))
}

freq <- function(Hr, ST, ET) {
	ifelse ((ST <= Hr) & (Hr < ET), 1, 0)
}

f_Freq_Boats_Present <- function(Interviews, ShiftCnt){
  # Cal frequency of boats present at each hourly interval of a day 
  #    over a monthly period
  Hr <- seq(1, 23, by = 1)
  ST <- Interviews$ST
  ET <- Interviews$ET
  n <- nrow(Interviews)
  
  Freq <- matrix(0, nrow = n, ncol = length(Hr))
  
  for (i in 1:n){
    Freq[i,] <- freq(Hr, ST[i], ET[i])
  }
  Freq <- apply(Freq, 2, sum)
  Freq <- Freq/sum(Freq)
  Freq <- Freq+0.000001  # avoid zero prob
  
  sumFreq <- sum(Freq)
  FreqNight <- (sumFreq - sum(Freq[ShiftCnt[1]:ShiftCnt[2]]))/sumFreq
  FreqProp <- Freq/sumFreq
  
  names(FreqProp) <- as.character(Hr)
  
  list(FreqPropNight = FreqNight, FreqProp = FreqProp)
}

f_creel_estimation <- function(simu_Year, simu_Month, shift_times_cnts, use_daily_estimator, est_by_daytype, Species, 
                               exclude_method, n_days_excluded, interviews, counts, imputation_method, n_impute, 
                               RovingCounts, AerialProportion, PropCorr, NumDays, FValue)
{
  # Use multiple-day estimator or daily estimator?
  use_multi_day_estimator = !use_daily_estimator
  
  # Imputing missing daily catch rate
  if (imputation_method != 3) n_impute <- 1
  if (use_multi_day_estimator) n_impute <- 1
  
  if ((!est_by_daytype))
  {
    NumDays <- sum(NumDays)
  }
  
  # Used for the fixed count-time aerial survey method or 
  #  the method of adjusting day-time effort estimate by the proportion of day-time effort when PropCorr = TRUE
  #list(FreqPropNight = FreqNight, FreqProp = FreqProp)
  BH_freq <- f_Freq_Boats_Present(interviews, shift_times_cnts)
  
  TE <- matrix(NA, nrow = n_impute, ncol = 2, byrow = TRUE) # total effort
  TB <- matrix(NA, nrow = n_impute, ncol = 2, byrow = TRUE) # total boat hours
  TC <- matrix(NA, nrow = n_impute, ncol = 2, byrow = TRUE) # total catch
  for (i_impute in 1:n_impute)
  {
    # est_by_daytype is only applied to daily effort estimation, not catch rate estimation
    # ---------------------------------------------------
    # psu is (1) a day for the daily estimator, and 
    #        (2) a period of multiple days of a daytype for the multiple-day estimator
    #        (3) a month for the multiple-day estimator if daytype is not considered  
    psuCatchRates <- f_psu_interview_stat(Species, use_daily_estimator, est_by_daytype, exclude_method, n_days_excluded, interviews, simu_Month)
    
    if (AerialProportion)
      # fixed count-time aerial survey
      psuBoatHours <- f_psu_boat_hours_Aerial_Proportion(use_daily_estimator, 
                                                         est_by_daytype, BH_freq$FreqProp, counts, simu_Month, NumDays, FValue)
    else 
    {
      if (!PropCorr) #  the method of adjusting day-time effort estimate by the proportion of day-time effort when PropCorr = TRUE
        BH_freq$FreqPropNight = 0
      
      psuBoatHours <- f_psu_boat_instant_hours(use_daily_estimator, 
                                               est_by_daytype, counts, simu_Month, NumDays, FValue, BH_freq$FreqPropNight)
    }
    
    # Effort estimation
    if (use_daily_estimator)
    {	
      # Fill-in daily statistics for days with no interview or count
      daily_E <- f_dailyEstimator_impute(imputation_method, simu_Year, simu_Month,  
                                         psuCatchRates$DAY, psuBoatHours$DAY, 
                                         psuBoatHours$PARTY_EFFORT, psuBoatHours$V_PARTY_EFFORT,
                                         psuCatchRates$PARTY_SIZE, psuCatchRates$V_PARTY_SIZE)
    }
    
    if (use_multi_day_estimator)
    {
      # Multiple-day estimation
      if (!est_by_daytype)
      {
        tot_effort <- f_multipleDayEst_total_effort(mean(psuBoatHours$NCNT), psuBoatHours$PARTY_EFFORT, psuBoatHours$V_PARTY_EFFORT,  psuCatchRates$PARTY_SIZE,  psuCatchRates$V_PARTY_SIZE, NumDays)
      } else
      {
        tot_effort <- f_multipleDayEst_total_effort_daytype(psuBoatHours$DAYTYPE, mean(psuBoatHours$NCNT), psuBoatHours$PARTY_EFFORT, psuBoatHours$V_PARTY_EFFORT, psuCatchRates$PARTY_SIZE, psuCatchRates$V_PARTY_SIZE, NumDays)
      }
      
    } else
    {
      # Daily estimation
      if (!est_by_daytype)
      {
        daytype <- rep(1, length(daily_E$DAYTYPE))
      }
      else
        daytype <- daily_E$DAYTYPE
      
      tot_effort <- f_dailyEstimator_total_effort(daytype, mean(psuBoatHours$NCNT), imputation_method, daily_E$PARTY_EFFORT,  daily_E$V_PARTY_EFFORT, daily_E$angler_effort, daily_E$v_angler_effort, NumDays)
    }
    TB[i_impute, 1] <- sum(tot_effort$Bp)
    TB[i_impute, 2] <- sum(tot_effort$VBp)  #VAR
    
    TE[i_impute, 1] <- sum(tot_effort$Ep)
    TE[i_impute, 2] <- sum(tot_effort$VEp)  #VAR
    
    # Catch rate estimation
    if (use_daily_estimator)
    {
      # impute daily catch rate
      daily_catch_rate <- f_dailyEstimator_impute_catch_rate(imputation_method, n_impute, 
                                                       psuCatchRates$DAY, daily_E$DAY, 
                                                       psuCatchRates$catch_rate, psuCatchRates$v_catch_rate)
      # daily catch estimation
      dailycatch <- f_dailyEstimator_Catch(daily_E$angler_effort, daily_E$v_angler_effort, 
                                           daily_catch_rate$catch_rate, daily_catch_rate$v_catch_rate)
    }
    
    # Stratum catch estimation
    if (use_multi_day_estimator)
    {
       totcatch <- f_multipleDayEst_total_catch(est_by_daytype, tot_effort$Ep, tot_effort$VEp, 
                                                psuCatchRates$catch_rate,  psuCatchRates$v_catch_rate)
    }
    else
    {
      # Daily estimation
      if (!est_by_daytype)
      {
        NDays <- sum(NumDays)
        totcatch <- f_dailyEstimator_totalCatch(mean(psuBoatHours$NCNT), imputation_method, 
                                                dailycatch$daily_catch, dailycatch$v_daily_catch, NDays)
      }
      else
        totcatch <- f_dailyEstimator_totalCatch_Daytype(daytype, mean(psuBoatHours$NCNT), imputation_method, 
                                                        dailycatch$daily_catch, dailycatch$v_daily_catch, NumDays)			
    }
    
    TC[i_impute, 1] <- sum(totcatch[1])  #Est
    TC[i_impute, 2] <- sum(totcatch[2])  #VAR
  }
  
  #total_boat_hours <- mean(TB[,1]) 
  if (imputation_method == 3 & use_daily_estimator)	
  {
    total_BH <- mean(TB[,1]) 
    W <- mean(TB[,2]) 
    B <- var(TB[,1])
    v_total_BH <- W + (1 + 1/n_impute) * B
    
    total_effort <- mean(TE[,1]) 
    W <- mean(TE[,2]) 
    B <- var(TE[,1])
    v_total_effort <- W + (1 + 1/n_impute) * B
    
    total_catch <- mean(TC[,1])
    W <- mean(TC[,2]) 
    B <- var(TC[,1])
    v_total_catch <- W + (1 + 1/n_impute) * B
  }
  else
  {
    total_BH <- (TB[,1]) 
    v_total_BH <- (TB[,2]) 
    
    total_effort <- (TE[,1]) 
    v_total_effort <- (TE[,2]) 
    
    total_catch <- (TC[,1])
    v_total_catch <- (TC[,2]) 
  }
  
  list(angler_hours = total_effort, v_angler_hours = v_total_effort, boat_hours = total_BH, v_boat_hours = v_total_BH, catch = total_catch, v_catch = v_total_catch)
}


