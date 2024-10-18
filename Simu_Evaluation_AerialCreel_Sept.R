#########################################
# Simulation of inland creel
#   Aug 4, 2021
# Revised for aerial survey simulation:
#   Dec 15, 2021 
#########################################

setwd("C:/Users/suz/OneDrive - State of Michigan DTMB/Documents/MiCreel/ResearchStudy/BusRouteStudy/Simulations/RGenerate_Simu_Aerial")

# True simulated population: Simu_2021_Belleville_Pop.txt
# TP_NM = "../True_Pop/Simu_BusRoute_2021_Pop.csv"  # July
TP_NM = "../True_Pop/Simu_2021_LK_ERIE_bus_route_Pop_9.csv"  # Sept

#True_Pop <- read.table(TP_NM, header =TRUE, fill = TRUE, sep ="\t", stringsAsFactors = FALSE)
True_Pop <- read.csv(TP_NM, header =TRUE, fill = TRUE, stringsAsFactors = FALSE)

Species = "YEP"
if (Species == "WAE") {
  True_catch = sum(True_Pop$WAE)
}else{
  True_catch = sum(True_Pop$YEP)
}

True_Pop_Values <- list(boat_hours = sum(True_Pop$TRIPHOURS), 
                        effort = sum(True_Pop$PARTYHOURS), 
                        catch = True_catch)
True_Pop_Values

set.seed(9413335);
True_Pop$position = runif(nrow(True_Pop), 0, 1)

# Survey project settings

(Simu_Year = True_Pop$YEAR[1])
(Simu_Month = True_Pop$MONTH[1])

m2d1 = paste(Simu_Year, "-", Simu_Month+1, "-", 1, sep="")
m1d1 = paste(Simu_Year, "-", Simu_Month, "-", 1, sep="")
DaysInMon = difftime(as.Date(m2d1), as.Date(m1d1) )

Nweekdays <- function(a, b, holidays, weekend) { 
  possible_days <- seq(a, b, "days")
  # Count all days that are not weekend and
  # are not holidays
  sum(!weekdays(possible_days) %in% weekend & !possible_days %in% holidays)
}

weekend <-  c("Saturday", "Sunday")
holidays <- as.Date(c(paste(Simu_Year, "-", 9, "-", 6, sep="")))

m1d2 = paste(Simu_Year, "-", Simu_Month, "-", DaysInMon, sep="")
nwds = Nweekdays(as.Date(m1d1), as.Date(m1d2), holidays, weekend)
nwes = DaysInMon - nwds

#NumDays <- c(20, 11); #July 2005
#NumDays <- c(21, 10); #July 2007 July 4
(NumDays <- c(nwds, nwes));

#Great Lakes Fishable Hours per Day per Month from Svboda AUG 1995:  
#(12, 12, 13, 15, 16, 18, 18, 17, 16, 15, 14, 13) 
# for months  (1   2   3   4   5   6   7   8   9  10  11  12).

source("R_CODE/R_SURVEY_SAMPLING.R")    # Schedule sampling
source("R_CODE/R_CREEL_ESTIMATION.R")   # Estimation
source("R_CODE/R_SIMU_MI_ANG_SURVEY.R") # Simulation

library("compiler")
enableJIT(3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sampling design options
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Instantaneous or progressive counting
progressive_Counting = TRUE   #TRUE
# Time interval for creating count starttimes within a sampling day
progressive_count_interval = 1  # one hour

RovingInterviewing = FALSE
AerialProportion =  FALSE
PropCorr =  FALSE #  TRUE

strat_sampling_by_wkday = T
is_completed_trip_interviews = T

max_nints <- 1000
#OptIntSampling = "threshhold"
#OptIntSampling = "fix_prop"
OptIntSampling = "all"

# excluding some days for interviewing to mimic actual survey situation
n_days_excluded <- 0; 
exclude_method <- 1; 
# imputation options
imputation_methods <- c("complete_case", "mean_inpute", "ABB")
imputation_method <- 1; 
n_impute <- 0

# simulation options

is_simu_pop <- T;
#shifts_non_overlapping = T
est_by_daytype = T

# Estimation options
dailyEst = T

# Count shifts
n_shifts_airplane = 1

# Sept
# Airplane flight times
if (n_shifts_airplane == 1)
{
  shift_times_airplane = c(7, 17)  # daytime c(7, 17) # all-day  #c(7, 21) same length as interview shift 
} else {
	n_shifts_airplane = 1
	shift_times_airplane = c(7, 17)  # daytime c(7, 17)
}

# Interview shifts
n_shifts_int = 2
shifts_non_overlapping =  FALSE 
if (n_shifts_int == 1)
{
  shift_times_int = c(3, 23)  # 
} else if (n_shifts_int == 2) 
{
	if (shifts_non_overlapping) {
	  # two shifts
	  shift_times_int = c(3, 13, 23)
	}else{
	  # Two overlapping shifts
      # Shift A: 7:00 AM~3:30 PM  Shift B: 2:30 PM~11:00 PM  
	  shift_times_int =  c(7, 15, 13, 21) # c(5, 15, 11, 21)  #c(7, 15, 13, 21)
	}
} else if (n_shifts_int == 3) 
{
  #Three shifts, two are sampled
  shift_times_int = c(4, 10, 16, 21)
} else
	{
	n_shifts_int = 1
	shift_times_int = c(5, 21)
}

#Sampling
res_scenarios <- data.frame()

set.seed(53335);
set.seed(33599);
n_simu <- 500
for (ncnts in c(1, 2)){  #c(1,2)
  for (ndays_s in c(10, 15, 20, 25)) {  #c(10,15,20,25)
    for (FValue in c(10)) {  #FishableHoursAerial <- 10  # 16  #for Sept  c(10, 11, 12, 13, 16)
      
      #ncnts=2
      #ndays_s=20
      #FValue = 11 # shift_times_airplane[2]-shift_times_airplane[1]

      if (dailyEst) {
        if (est_by_daytype)
          FN <- paste("res_xx_PropCorr_nD_", ndays_s, "_SHint", n_shifts_int, "_SHc", n_shifts_airplane, "_F" , FValue, "_ncnts_" , ncnts, sep="")
        else
          FN <- paste("res_DailyEst_sday_prop_shift7-19_", ndays_s, "_shifts_int", n_shifts_int, "_c", n_shifts_airplane,"_F" , FValue,  "_ncnts_" , ncnts, sep="")
      }else{
        if (est_by_daytype)
          FN <- paste("res_MultiDayEst_by_daytype_sday", ndays_s, "_shifts_int", n_shifts_int, "_c", n_shifts_airplane,"_F" , FValue,  "_ncnts_" , ncnts, sep="")
        else
          FN <- paste("res_MultiDayEst_sday", ndays_s, "_shifts_int", n_shifts_int, "_c", n_shifts_airplane,"_F" , FValue,  "_ncnts_" , ncnts, sep="")
      }
      
      scenario_df <- data.frame(n_simu = n_simu, Simu_Month = Simu_Month, Species = Species, ndays_s = ndays_s, n_shifts_int = n_shifts_int, 
                                  n_shifts_airplane = n_shifts_airplane, ncnts = ncnts,
                                  FValue = FValue)
      
      options(scipen=100)
      #debug(f_angler_survey_simu)
      #debug(f_psu_catch_rate)
      #debug(make_roving_interviews)
      #debug(f_creel_estimation)
      #debug(make_sample)
      
      #write.table(counts, "simu_counts_progressive.txt", col.names = NA, sep ="\t")
      #write.table(interviews, "simu_interviews_incomp.txt", col.names = NA, sep ="\t")
      
      #debug(f_psu_interview_stat)
      #Rprof()
      
      assign(FN, f_angler_survey_simu(Species, n_simu, dailyEst, est_by_daytype, True_Pop_Values, 
                                      AerialProportion, PropCorr, NumDays, FValue, 
                                      strat_sampling_by_wkday, progressive_Counting, 
                                      RovingInterviewing, ndays_s, ncnts, progressive_count_interval, 
                                      max_nints, is_simu_pop, Simu_Year, Simu_Month,True_Pop, 
                                      n_shifts_airplane, shift_times_airplane, 
                                      n_shifts_int, shift_times_int, shifts_non_overlapping, 
                                      OptIntSampling, exclude_method, n_days_excluded, imputation_method, n_impute))
      
      
      #Rprof(NULL)
      
      #debug(f_summaries)
      simu_res <- get(FN)
      summ <- round(f_summaries(simu_res, True_Pop_Values), 2)
      
      cat("\n" )
      cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" )
      
      cat("ndays_s", ndays_s, "ncnts", ncnts, "\n" )
      print(FN)
      print(True_Pop_Values)
      summ
      
      summ2 = round(f_summaries_table(simu_res, True_Pop_Values), 2)
      print(summ2)
      write.table(summ2, paste(".\\Results_", Simu_Month ,"\\summ2_", FN, ".txt", sep = ""), sep ="\t")
      
      simu_res$Scenario = FN
      
      summf = t(c(Scenario=FN, summ))
      print(summf)
      
      res_scen = cbind(scenario_df, summf)
      res_scenarios <- rbind(res_scenarios, res_scen)
        
      write.table(summf, paste(".\\Results_", Simu_Month ,"\\summ_", FN, ".txt", sep = ""), sep ="\t", row.names = F)
      
      write.table(simu_res, paste(".\\Results_", Simu_Month ,"\\", FN, ".txt", sep = ""), sep ="\t", col.names = NA)
    }
  }
}

write.table(res_scenarios, paste(".\\Results_", Simu_Month ,"\\Results_scenarios_PropCorr", PropCorr,  "_", FN, ".txt", sep = ""), sep ="\t", row.names = F)


