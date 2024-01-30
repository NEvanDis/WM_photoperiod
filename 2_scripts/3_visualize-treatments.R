# Analysis of winter moth egg photoperiod experiments 1999-2000 and 2013-2014 ####
# Manipulated 1) photoperiod and 2) photoperiod+temperature the eggs received during development
# Visualize photoperiod and temperature treatments used

# Open R project in main folder

# Load packages
#-----------------------------------
library(tidyverse)
library(readxl)
library(cowplot)
theme_set(theme_cowplot()) #white background instead of grey -> don't load if want grey grid


# Load data 
#-----------------------------------
temp <- list( cold=read.csv("1_data/settemps_cold-cor.csv") %>% mutate(Treatment="Cold"),
                   warm=read.csv("1_data/settemps_warm-cor.csv") %>% mutate(Treatment="Warm") )
temp <- lapply(temp, function(X){ 
  X <- X %>% mutate(Date=as.Date(Date))
  X_h1 <- X %>% mutate(temp=h1, Temp_cycle="day_min") %>% select(Treatment, Date, temp, Temp_cycle) 
  X_h7 <- X %>% mutate(temp=h7, Temp_cycle="day_max") %>% select(Treatment, Date, temp, Temp_cycle)
  X_h19 <- X %>% mutate(temp=h19, Temp_cycle="day_average") %>% select(Treatment, Date, temp, Temp_cycle)
  
  X <- rbind(X_h1, X_h7, X_h19) # formatted for plotting
  })
lapply(temp, tail)
temp <- data.table::rbindlist(temp, use.names=T)


photo <- list( control=read_xls("1_data/photoperiodic schedule.xls", sheet=3, skip=1),
                 vearly=read_xls("1_data/photoperiodic schedule.xls", sheet=1, skip=1),
                 early=read_xls("1_data/photoperiodic schedule.xls", sheet=2, skip=1),
                 late=read_xls("1_data/photoperiodic schedule.xls", sheet=4, skip=1),
                 vlate=read_xls("1_data/photoperiodic schedule.xls", sheet=5, skip=1) )
photo <- lapply(photo, function(X){
  colnames(X) <- c("date", "day", "actual_date", "sunrise", "sunset", "remarks")
  X <- X %>% mutate(
    date=as.character(openxlsx::convertToDate(date)) %>% as.Date, #backtransform dates from numerical to date format
    actual_date=as.Date(actual_date)) %>%
    mutate(
    start_light=( as.integer(gsub("\\d{4}-\\d{2}-\\d{2} (\\d{2}):\\d{2}:\\d{2}", "\\1", sunrise)) +
                        as.integer(gsub("\\d{4}-\\d{2}-\\d{2} \\d{2}:(\\d{2}):\\d{2}", "\\1", sunrise))/60 ),
    end_light=( as.integer(gsub("\\d{4}-\\d{2}-\\d{2} (\\d{2}):\\d{2}:\\d{2}", "\\1", sunset)) +
                    as.integer(gsub("\\d{4}-\\d{2}-\\d{2} \\d{2}:(\\d{2}):\\d{2}", "\\1", sunset))/60 ),
    Treatment=ifelse(date==actual_date, "  0 weeks", 
                     ifelse(date-actual_date==-28, " -4 weeks", 
                            ifelse(date-actual_date==-14, " -2 weeks",
                                   ifelse(date-actual_date==14, "+2 weeks", "+4 weeks")))) # spaces in names to make legend look nicer aligned
  ) 
  Y <- X %>% mutate(day_h=end_light - start_light, Treatment=ifelse(is.na(Treatment), X$Treatment[1], Treatment) )
  
  if(X$Treatment[1]=="+4 weeks"){ # cut-off for plotting purposes
    Y <- filter(Y, actual_date<=as.Date("2014-04-02"))
  }
  
  return(Y)
})
lapply(photo, tail)
photo <- data.table::rbindlist(photo, use.names=T) %>% 
  mutate(Treatment=factor(Treatment, levels=c(" -4 weeks", " -2 weeks", "  0 weeks", "+2 weeks", "+4 weeks")))
# NB: note that the photoperiods are from the photoperiod experiment, which started on January 15 (~1 month after egg laying)
# For the photoperiod-temperature experiment, the same [-2 weeks] and [+2 weeks] treatments were used, but starting earlier on December 12 (just after egg laying)


# Visualize
#-----------------------------------

# -- Photoperiod treatments
phot_treat <- ggplot(data=photo, aes(x=actual_date, y=day_h, col=Treatment))+
  geom_smooth() +
  scale_x_date(limits=c(as.Date("2014-01-15"), as.Date("2014-06-05")), expand=c(0,0))+
  labs(y="Photoperiod treatment (light:dark)")+
  scale_y_continuous(breaks=seq(8, 15, by=1), labels=c("8:16", " ", "10:14", " ", "12:12", " ", "14:10", " "))+
  theme(axis.title.y=element_text(size=16, vjust=0.5), axis.title.x=element_blank(),
        axis.text.y=element_text(size=14), axis.text.x=element_text(size=18), 
        legend.title = element_text(size=16), legend.text = element_text(size=14))+
  background_grid(major = "xy") #get grid lines

phot_treat1 <- phot_treat + theme(legend.position="none")
phot_treat1

ggsave(filename="methods/phot_treats_wlegend.png", plot=phot_treat, device="png", width=180, height=100, units="mm", dpi="print")
ggsave(filename="methods/phot_treats.png", plot=phot_treat1, device="png", width=180, height=100, units="mm", dpi="print")


# -- Temperature treatments
temp_averages <- data.frame(temp=c(3.935, 5.29), Treatment=c("Cold", "Warm")) # taken from Supplemental table S3

temp_treat <- ggplot(data=temp, aes(x=Date))+
  scale_color_manual(values=c("purple","red","blue"))+
  facet_wrap(~ Treatment, ncol=1)+
  geom_point(aes(y=temp, col=Temp_cycle))+
  geom_smooth(aes(y=temp, col=Temp_cycle))+
  geom_hline(data=temp_averages, aes(yintercept=temp), linewidth=1.5, linetype="twodash", col="aquamarine2")+
  ylab(expression("Temperature treatment"*~degree*C))+
  coord_cartesian(ylim = c(-8, 17))+ scale_y_continuous(breaks=seq(-5,15,by=5), labels=c("  -5", "   0", "   5", "   10", "  15")) +
  scale_x_date(limits=c(as.Date("2000-01-15"), as.Date("2000-06-05")), expand=c(0,0))+
  theme(axis.title.y=element_text(size=18, vjust=-0.5), axis.title.x=element_blank(),
        axis.text.y=element_text(size=14), axis.text.x=element_text(size=18), strip.text = element_text(size=18),
        legend.title = element_text(size=16), legend.text = element_text(size=14))+
  background_grid(major = "xy") #get grid lines

temp_treat1 <- temp_treat + theme(legend.position="none")
temp_treat1

ggsave(filename="methods/temp_treats_wlegend.png", plot=temp_treat, device="png", width=180, height=150, units="mm", dpi="print")
ggsave(filename="methods/temp_treats.png", plot=temp_treat1, device="png", width=180, height=150, units="mm", dpi="print")
