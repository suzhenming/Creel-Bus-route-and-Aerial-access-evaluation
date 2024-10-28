####################################################
# Simulation evaluation of bus-route creel method
#  by Dr. Zhenming Su
#  Started Sept 5, 2021
#
####################################################

setwd("../Simu_bus_route")

library(lubridate)

# True simulated population: for September
TP_NM = "../True_Pop/Simu_2021_LK_ERIE_bus_route_Pop_9.csv"

True_Pop <- read.csv(TP_NM, header =TRUE, fill = TRUE, stringsAsFactors = FALSE)
#(True_Pop_Values <- list(boat_hours = sum(True_Pop$TRIPHOURS), effort = sum(True_Pop$PARTYHOURS), catch = sum(True_Pop$WAE)))

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

DailyBHTrue <- aggregate(True_Pop$TRIPHOURS, list(True_Pop$DAY), sum)
names(DailyBHTrue) <- c("DAY","BH")
V_BD <- var(DailyBHTrue$BH)

# Survey project settings
(Simu_Year = True_Pop$YEAR[1])
(Simu_Month = True_Pop$MONTH[1])

#m2d1 = paste(Simu_Year, "-", Simu_Month+1, "-", 1, sep="")
m1d1 = paste(Simu_Year, "-", Simu_Month, "-", 1, sep="")
#(DaysInMon = as.numeric(difftime(as.Date(m2d1), as.Date(m1d1), unit = "day")))
SurveyYearMon = ymd(m1d1)
(DaysInMon = as.numeric(days_in_month(SurveyYearMon)))

nWeekdays <- function(d1, d2, holidays, weekend) { 
  # d1, d2 are date objects
  days <- seq(d1, d2, "days")
  # Count all days that are not weekend and
  # are not holidays
  sum(!weekdays(days) %in% weekend & !days %in% holidays)
}
weekend <-  c("Saturday", "Sunday")
holidays <- 1  # as.Date(c(paste(Simu_Year, "-", 5, "-", 29, sep="")))
m1d2 = paste(Simu_Year, "-", Simu_Month, "-", DaysInMon, sep="")
nwds = nWeekdays(as.Date(m1d1), as.Date(m1d2), holidays, weekend)
nwes = DaysInMon - nwds

(NumDays <- c(nwds, nwes));

# Bus-route sampling and estimation
library("compiler")
enableJIT(3)
source("R_CODE_BUS_ROUTE/R_SURVEY_SAMPLING_BUS_ROUTE.R")    # Schedule sampling
source("R_CODE_BUS_ROUTE/R_CREEL_ESTIMATION_BUS_ROUTE.R")   # Estimation
source("R_CODE_BUS_ROUTE/R_SIMU_MI_ANG_SURVEY_BUS_ROUTE.R") # Simulation

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sampling design options
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Day and shift sampling
#ndays_s <- 25

# Stratification of days by weekday or weekend
strat_sampling_by_wkday = T

#Day and shift sampling
n_shifts = 2
shift_std_len = 8 #for two shifts

AStart <- 6 
BEnd <- 22

expansion_method = F  # Don't what is this for
if (n_shifts == 1)
{
  # One shift: 18 hours
  # shift_times = c(4, 22)
  shift_times = c(AStart, BEnd)
  shift_length = shift_times[2]-shift_times[1]

  shifts_overlapping  = F
  adjust_overlap_prob = F
} else if (n_shifts == 2) 
{
  shifts_overlapping = T
  if (!shifts_overlapping) {
    # Two non-overlapping shifts: 8 hours
    # shift_times = c(6, 14, 22)
    AEnd <- AStart+shift_std_len
    BEnd = AEnd+shift_std_len
    shift_times = c(AStart, AEnd, BEnd)
    shift_length = shift_std_len
    
    adjust_overlap_prob = F
  }else{
    # Two overlapping shifts: 8 hours
    overlapping_time = 1
    AEnd   <- AStart + shift_std_len + overlapping_time # shift_times[2] #16
    BStart   <- BEnd - shift_std_len -  overlapping_time   #shift_times[4] #22  
    #shift_times = c(7, 15, 14, 22)
    
    shift_times = c(AStart, AEnd, BStart, BEnd)  
    shift_length = (shift_times[4]-shift_times[1])/2+overlapping_time
    
    adjust_overlap_prob = T 
    expansion_method = F
  }
} else if (n_shifts == 3) {
  #Three shifts, two are sampled: 6 hours
  shift_times = c(4, 10, 16, 22)
  shift_length = (shift_times[4]-shift_times[1])/3

  shifts_overlapping = F  
  adjust_overlap_prob = F  
} else {
  n_shifts = 1
  shift_times = c(4, 22)
  shift_length = shift_times[2]-shift_times[1]
  
  shifts_overlapping = F  
  adjust_overlap_prob = F  
}

# Bus-route settings
# 2 sites on a route
Travel_time = 20 #minutes
creel_time_total = 60*shift_length - 2*Travel_time
(creel_time_total/60)

if (n_shifts==1){
   #One shifts
   # total_route_time = 16 hours
   route_proto = structure(list(
		site_num = 1:2, 
		SiteName = c("Sterling", "Bolles"), 
		# CreelTime_Min = c(615, 425),  #shift_times = c(4, 22)
		CreelTime_Min = creel_time_total * c(0.6, 0.4),
		TravelTime_Min = c(Travel_time, Travel_time)
   ), 
   class = "data.frame", row.names = c("1", "2"))
} else {
   #Two shifts
   # total_route_time = 8 hours
   route_proto = structure(list(
		site_num = 1:2, 
		SiteName = c("Sterling", "Bolles"), 
		#CreelTime_Min = c(260, 180), 
		# CreelTime_Min = c(350, 90),
		# CreelTime_Min = c(90, 350),
		#TravelTime_Min = c(20, 20)
		CreelTime_Min = creel_time_total * c(0.6, 0.4),
		TravelTime_Min = c(Travel_time, Travel_time)		
   ), 
   class = "data.frame", row.names = c("1", "2"))
}
(total_route_time <- (sum(route_proto$CreelTime_Min) + sum(route_proto$TravelTime_Min))/60)

# Interview sampling options
max_nints <- 50
OptIntSampling = "threshhold"
OptIntSampling = "all"

# Estimation options
est_by_daytype = T
is_simu_pop <- T;

#Sampling variance
fpc = TRUE

# simulation options
options(scipen=100)

set.seed(5335);
n_simu <- 5000
res_scenarios <- data.frame()

for (ndays_s in c(10,15,20,25)) 
{  
  FN <- paste("bus_route_fpc_", fpc, "_species_", Species, "_sday", ndays_s, "_shifts_", n_shifts, '_OL_', shifts_overlapping, "_Adj_OL-", adjust_overlap_prob, sep="")
  scenario_df <- data.frame(n_simu = n_simu, Simu_Month = Simu_Month, 
                            Species = Species, ndays_s = ndays_s, 
                            n_shifts = n_shifts, 
                            shifts_overlapping = shifts_overlapping,
                            adjust_overlap_prob = adjust_overlap_prob
  )
  
  assign(FN, f_angler_survey_bs_simu(n_simu, Species, est_by_daytype, True_Pop_Values, DailyBHTrue, NumDays, 
                                     total_route_time, route_proto, 
                                     strat_sampling_by_wkday, ndays_s, is_simu_pop, True_Pop, 
                                     n_shifts, shift_times, !shifts_overlapping, adjust_overlap_prob, expansion_method,
                                     OptIntSampling, max_nints, fpc))
  
  
  simu_res <- get(FN)
  summ <- round(f_summaries(simu_res$Estimates, True_Pop_Values), 2)
  FN
  summ
  # Variance component 
  V_WD <- (mean(simu_res$VWD))
  
  cat("Prop V_BD = ", V_BD/(V_WD+V_BD), "\n")
  cat("V_2stage = ", mean(simu_res$V_2stage), "\n")
  
  sum_tab = round(f_summaries_table(simu_res$Estimates, True_Pop_Values), 4)
  sum_tab=t(sum_tab)
  
  BH_stat = sum_tab[1,]
  names(BH_stat) = paste(names(sum_tab[1,]),"BH", sep="_")
  
  E_stat = sum_tab[2,]
  names(E_stat) = paste(names(sum_tab[2,]),"E", sep="_")
  
  C_stat = sum_tab[3,]
  names(C_stat) = paste(names(sum_tab[3,]),"H", sep="_")
  
  res_scen = cbind(scenario_df, t(BH_stat), t(E_stat), t(C_stat))
  res_scenarios <- rbind(res_scenarios, res_scen)
  
  simu_res$Scenario = FN
  summf = t(c(Scenario=FN, summ))
  summf

  write.table(sum_tab, paste(".\\Results_Bus_route\\Sept\\summ_tab_", FN, ".txt", sep = ""), sep ="\t", row.names = T)
  write.table(summf, paste(".\\Results_Bus_route\\Sept\\summ_", FN, ".txt", sep = ""), sep ="\t", row.names = F)
  write.table(simu_res, paste(".\\Results_Bus_route\\Sept\\", FN, ".txt", sep = ""), col.names = NA, sep ="\t")
}

write.table(res_scenarios, paste(".\\Results_Bus_route\\Sept\\Results_scenarios_", FN, ".txt", sep = ""), sep ="\t", row.names = F)


