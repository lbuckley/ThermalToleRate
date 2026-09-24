## Script for calculating daily injury under natural temperature fluctuations

install.packages("tidyverse",repos = "http://cran.us.r-project.org")
install.packages("cowplot",repos = "http://cran.us.r-project.org")
install.packages("rootSolve",repos = "http://cran.us.r-project.org")
install.packages("zoo",repos = "http://cran.us.r-project.org")
install.packages("readxl",repos = "http://cran.us.r-project.org")
install.packages("tictoc",repos = "http://cran.us.r-project.org")
install.packages("lubridate",repos = "http://cran.us.r-project.org")
library(tidyverse)
library(cowplot)###plot settings
library(scales)
library(rootSolve)
library(zoo)
library(readxl)##read excel files
library(tictoc)##calculate the time running a programme
library(lubridate)##deal with date and time


theme_set(theme_cowplot())
Sys.setlocale("LC_ALL", "English")
Colors <- data.frame("Red" = c("#F8766D"), "Green" = c("#00BA38"), "Blue" = c("#619CFF"))

## input TDT Parameters
setwd("/Users/lbuckley/Library/CloudStorage/GoogleDrive-lbuckley@uw.edu/My\ Drive/Buckley/Work/AphidThermalTol/data/aphidTDT/Thermal\ injury\ calculation")

TDT <- read.csv("tdt.csv")
TDT_high <- TDT %>%
  filter(type == "warm")
{
  MD_6_high <- TDT_high %>%
    filter(species == "MD", age == 6)
  
  RP_6_high <- TDT_high %>%
    filter(species == "RP", age == 6)
  
  SA_6_high <- TDT_high %>%
    filter(species == "SA", age == 6)
  }

TDT_low <- TDT %>%
  filter(type == "cold")

{
  MD_6_low <- TDT_low %>%
    filter(species == "MD", age == 6)
  
  RP_6_low <- TDT_low %>%
    filter(species == "RP", age == 6)
  
  SA_6_low <- TDT_low %>%
    filter(species == "SA", age == 6)
  }
### Calculated Tc, here we define Tc as the Temperature causing 50% mortality after one week(7*24*60=10080min)'s exposure.
{
  Tc_MD = round((log(10080,10)-MD_6_high$intercept)/(MD_6_high$slope),2)
  
  Tc_SA = round((log(10080,10)-SA_6_high$intercept)/(SA_6_high$slope),2)
  
  Tc_RP = round((log(10080,10)-RP_6_high$intercept)/(RP_6_high$slope),2)
  
  Tc_cold_MD = round((log(10080,10)-MD_6_low$intercept)/(MD_6_low$slope),2)
  
  Tc_cold_SA = round((log(10080,10)-SA_6_low$intercept)/(SA_6_low$slope),2)
  
  Tc_cold_RP = round((log(10080,10)-RP_6_low$intercept)/(RP_6_low$slope),2)
}

## Importing temperature data and convert to minutely recordings
###summer temperature data for three sites
Wuhan_summer_temps <- read.csv("wuhan 2016-05-01 to 2016-08-31.csv")

####  Creating rows with "linearized" temperature change (1-min resolution instead of 1-hour)  ####

#Selecting relevant columns and fixing date/time as POSIXct-values (to understand that it is time, not characters)
Wuhan_summer_temps <- Wuhan_summer_temps %>% 
  dplyr::select(name, datetime, temp) %>%
  dplyr::mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d T %H:%M:%S", tz = "GMT"))

#Using zoo (Z's ordered observations) to order time and creating 1-min resolution (60 s)
## summer 
Wuhan_summer_seq <- zoo(order.by = (as.POSIXct( seq (min(Wuhan_summer_temps$datetime), max(Wuhan_summer_temps$datetime), by = 60))))

#Merging temperatures with the time-sequence (NAs produced for all the "new" minutes)
Wuhan_summer_merge <- merge(zoo(Wuhan_summer_temps[3], order.by = Wuhan_summer_temps$datetime), Wuhan_summer_seq)

#Using na.approx to interpolate temp where NA is present (in all the "new" minutes) - formatting as data frame and creating columns with time and date
Wuhan_summer <- na.approx(Wuhan_summer_merge) %>%
  as.data.frame() %>%
  rownames_to_column("datetime") %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "GMT"), 
         date = as.POSIXct(format(datetime, format = "%Y-%m-%d")),
         city = "Wuhan") 

Wuhan_summer_data <- Wuhan_summer

#Checking that all original time points are present and that the temperature is identical (returns 0 rows if TRUE)
Wuhan_summer_temps %>%
  full_join(Wuhan_summer_data, by = "datetime") %>%
  filter(temp.x != temp.y)   

#Quick plotting of the interpolated (red line) and original data (black points)
ggplot(Wuhan_summer_data, aes(datetime, temp))+
  geom_line(col = "red", size = 1)+
  geom_point(data = Wuhan_summer_temps)+
  scale_x_datetime(breaks = date_breaks("1 week"))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

## Calculation of injury accumulation
### Wuhan summer

#The data frame CITY_SEASON_data can now be used for calculation of injury accumulation
####  From Scientific Reports script  ####
#Creating new columns for the time/temperature-specific injury accumulation rate

Wuhan_summer_data$add_inj_MD_6 = NA
Wuhan_summer_data$add_inj_RP_6 = NA
Wuhan_summer_data$add_inj_SA_6 = NA

#Additive injury is determined per time interval (in minutes)

tic()
for(time in 1:nrow(Wuhan_summer_data))
{
  Wuhan_summer_data$add_inj_MD_6[time] = ifelse(any(Wuhan_summer_data$temp[time+1]>=Tc_MD,Wuhan_summer_data$temp[time]>=Tc_MD),
                                                (100*(Wuhan_summer_data$datetime[time+1]-Wuhan_summer_data$datetime[time]))/(10^(MD_6_high[[4]]*max(Wuhan_summer_data$temp[time+1],Wuhan_summer_data$temp[time])+MD_6_high[[9]])),0)
  
  Wuhan_summer_data$add_inj_RP_6[time] = ifelse(any(Wuhan_summer_data$temp[time+1]>=Tc_RP,Wuhan_summer_data$temp[time]>=Tc_RP),
                                                (100*(Wuhan_summer_data$datetime[time+1]-Wuhan_summer_data$datetime[time]))/(10^(RP_6_high[[4]]*max(Wuhan_summer_data$temp[time+1],Wuhan_summer_data$temp[time])+RP_6_high[[9]])),0)
  
  Wuhan_summer_data$add_inj_SA_6[time] = ifelse(any(Wuhan_summer_data$temp[time+1]>=Tc_SA,Wuhan_summer_data$temp[time]>=Tc_SA),
                                                (100*(Wuhan_summer_data$datetime[time+1]-Wuhan_summer_data$datetime[time]))/(10^(SA_6_high[[4]]*max(Wuhan_summer_data$temp[time+1],Wuhan_summer_data$temp[time])+SA_6_high[[9]])),0)
  
}
toc()

#Summing the accumulated injury in new columns

Wuhan_summer_data$cum_inj_MD_6 = cumsum(Wuhan_summer_data$add_inj_MD_6)
Wuhan_summer_data$cum_inj_RP_6 = cumsum(Wuhan_summer_data$add_inj_RP_6)
Wuhan_summer_data$cum_inj_SA_6 = cumsum(Wuhan_summer_data$add_inj_SA_6)


#Checking that injury was calculated for all time points (not possible to calculate for the last time point)
summary(Wuhan_summer_data)

#Saving under another name to avoid messing with the data that took a long time to calculate
Wuhan_summer_injury <- Wuhan_summer_data
#write.csv(Wuhan_summer_data, "Wuhan_summer_injury.csv")


#Quick overview of accumulated injury over the full period (temperature adjusted to allow showing all data in the same plot)
ggplot(Wuhan_summer_injury)+
  geom_line(aes(datetime, cum_inj_MD_6/1000), color=Colors$Red, size = 1.1)+
  geom_line(aes(datetime, cum_inj_RP_6/1000), color=Colors$Green, size = 1.1)+
  geom_line(aes(datetime, cum_inj_SA_6/1000), color=Colors$Blue, size = 1.1)+
  geom_line(aes(datetime, temp), color = "darkgrey", alpha = .8)+
  geom_hline(yintercept =  30, color = "purple")+
  scale_x_datetime(breaks = date_breaks("1 week"))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

#Summarizing injury accumulation within a day (midnight-midnight)
Wuhan_summer_daily_injury <- Wuhan_summer_injury %>%
  mutate(date = format(date, "%Y-%m-%d")) %>%
  group_by(date, city) %>%
  summarise(
    
    daily_injury_MD_6 = sum(add_inj_MD_6),
    daily_injury_RP_6 = sum(add_inj_RP_6),
    daily_injury_SA_6 = sum(add_inj_SA_6),
  ) %>%filter(date != "2016-09-01")


summary(Wuhan_summer_daily_injury)

Wuhan_summer_hourly_temps <- Wuhan_summer_temps %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "GMT"), 
         date = as.POSIXct(format(datetime, format = "%Y-%m-%d")))


Wuhan_summer_daily_injury_long <- Wuhan_summer_daily_injury %>%
  select(date, city, daily_injury_MD_6, daily_injury_RP_6, daily_injury_SA_6) %>%
  pivot_longer(!c(date, city), names_to = "Species", values_to = "Daily_injury") %>%
  mutate(Daily_injury = Daily_injury/100)


Wuhan_summer_without_zero <- Wuhan_summer_daily_injury_long %>%
  filter(Daily_injury != 0)

