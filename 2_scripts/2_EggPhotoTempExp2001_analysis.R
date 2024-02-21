# Analysis of winter moth egg temperature-photoperiod experiment 2000-2001 ####
# Manipulated photoperiod and temperature the eggs received during development, kept in climate cabinets
# RQ: what is the relative contribution of photoperiod and temp on the timing of egg hatching?

# Open R project in main folder

# Load packages
#-----------------------------------
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot()) #white background instead of grey -> don't load if want grey grid
library(lme4)
library(lmerTest)


# Load data ####
#-----------------------------------
d_2001 <- read.csv("1_data/db_EggPhotTemp2001_d50.csv") # 2001 experiment d50 egg hatch dates per subclutch
head(d_2001)

table(is.na(d_2001$D50Calc)) # for 1 subclutches no D50 available, N=79 subclutches
length(unique(d_2001$TubeID)) # for 20 females
table(d_2001$AreaShortName) # all from OH area
table(d_2001$NovemberDate) # Catch dates
table(d_2001$Treatment) # in total 4 treatments, KK = koud kort (cold year 1973 + short day length=early season photoperiod), WL = warm lang (warm year 1999 + long day length=late season photoperiod)
summary(d_2001$Eggs) # subclutches of at least 14 eggs
aggregate(Eggs~TubeID, d_2001, sum) # range of clutch size
aggregate(ClutchID~TubeID+Treatment, d_2001, length) # number of subclutches per female per treatment
filter(d_2001, is.na(D50Calc)) # female and treatment for which subclutch was excluded because <10 eggs hatched
unique(d_2001$ClutchID) %>% length() # total number of subclutches


# -- Prep data for analysis
d_2001 <- d_2001 %>% mutate(TempTreat=as.factor(ifelse(Treatment=="KK" | Treatment=="KL", "Cold", "Warm")),
                            PhotTreat=as.factor(ifelse(Treatment=="KK" | Treatment=="WK", "-2 weeks", "+2 weeks")), TubeID=as.factor(TubeID)) %>% 
  mutate(TempTreat=factor(TempTreat, levels=c("Cold", "Warm")),
         PhotTreat=factor(PhotTreat, levels=c("-2 weeks", "+2 weeks"))) %>%
  mutate(eggs_devtime= as.numeric(as.Date(D50Calc, origin=as.Date("2001-03-31"))-as.Date(NovemberDate, origin=as.Date("2000-10-31")))) %>%
  filter(!is.na(D50Calc))
head(d_2001)


# CHECK TEMP LOGGER DATA ####
#----------------------------------
d_temp2001 <- read.csv("1_data/db_EggPhotTemp2001_temploggers.csv") # load temp logger data
d_temp2001 <- d_temp2001 %>% mutate(Treatment=as.factor(LocationName), DateTime=as.POSIXct(paste0(Year, "-", Month, "-", Day, " ", Hour, ":", Minute))) %>% 
  select(Treatment, DateTime, LoggerName, Temperature)
head(d_temp2001)


# -- Use temp average from Dec 12 (experiment start) to median egg hatching date per treatment
medianhatchtreatment <- aggregate(D50Calc~Treatment, data=d_2001, median) # get median egg hatching date per treatment
medianhatchtreatment <- medianhatchtreatment %>% mutate(Date=as.Date(D50Calc, origin=as.Date("2001-03-31"))) #ignores the decimals

# --- Get temperature data for that time period
temps2001 <- d_temp2001 %>% group_split(Treatment)
names(temps2001) <- lapply(temps2001, function(X) X$Treatment[1]) %>% unlist()
lapply(temps2001, function(X) summary(X$DateTime)) # all same minimum date of Dec 12
lapply(temps2001, nrow)

temps2001 <- lapply(temps2001, function(X) filter(X, as.Date(DateTime)<medianhatchtreatment$Date[medianhatchtreatment$Treatment==X$Treatment[1]])) # subset treatment on median egg hatching date
lapply(temps2001, nrow)
lapply(temps2001, function(X) summary(X$DateTime))

Cold_temp <- mean(rbind(temps2001$KK, temps2001$KL)$Temperature)
Warm_temp <- mean(rbind(temps2001$WK, temps2001$WL)$Temperature)
Warm_temp - Cold_temp # mean temp difference between warm and cold treatments
  
# --- Calculate mean temps per incubator
meantemps2001 <- aggregate(Temperature~Treatment, data=data.table::rbindlist(temps2001), mean)
colnames(meantemps2001)[2] <- "MeanTemp"
meantemps2001$sd <- aggregate(Temperature~Treatment, data=data.table::rbindlist(temps2001), sd)$Temperature
# write.csv(meantemps2001, file="results/meantemps_2001.csv", row.names=F)
rm(temps2001, medianhatchtreatment, d_temp2001) # clean up R environment



#-----------------------------------
# Experiment analysis ####
#-----------------------------------

#-- Visualize hatching dates ####
devtime_avg <- Rmisc::summarySE(filter(d_2001, !is.na(eggs_devtime)), measurevar="eggs_devtime", groupvars=c("TempTreat","PhotTreat")) # average D50 hatch dates per photoperiod treatment
devtime_avg

raw_devtime <- ggplot(data=devtime_avg, aes(x=PhotTreat, y=eggs_devtime))+
  scale_colour_manual(values=c("grey27", "orangered2"))+ #"dodgerblue4"
  facet_grid(~TempTreat)+
  geom_jitter(data=d_2001, aes(col=PhotTreat), alpha=0.3, size=2.5, height=0, width=0.25)+ # jitter in width but not height
  geom_point(size=5, col="black") +
  geom_errorbar(aes(ymax = eggs_devtime+se, ymin=eggs_devtime-se), width=0.1, col="black") +
  labs(y="Egg development time (days)", x="Photoperiod treatment")+
  scale_y_continuous(lim=c(125, 165), breaks=seq(-40,200, by=5))+
  theme(axis.title.y=element_text(size=18, vjust=2), axis.title.x=element_text(size=18, vjust=-0.5),
        axis.text=element_text(size=16), legend.position="none", strip.text = element_text(size=18))
raw_devtime
# ggsave(filename="results/eggdevtime_raw_2001.png", plot=raw_devtime , device="png", width=200, height=150, units="mm", dpi="print")


#-- Linear model ####
lm1 <- lmer(eggs_devtime ~ TempTreat*PhotTreat + (1|TubeID), data=d_2001) #random effect: correct for split-brood design, not fitting slope because <4 treatments with replicates
summary(lm1)
anova1 <- anova(lm1) # find interaction effect BUT check if can be explained by temp logger data instead! ####
anova1$mod <- "lm1"

ranova1 <- ranova(lm1) # random intercept significant
ranova1$mod <- "lm1"

lm_res1 <- summary(lm1)$coefficients %>% as.data.frame

#--- check residuals
plot(lm1) #equal variance? ok, one major outlier! +-2SD from mean
qqnorm(resid(lm1)) #normally distributed? okish
qqline(resid(lm1))

  #- check outlier effect
d_2001$resid <- residuals(lm1)
d_2001$fitted <- predict(lm1)
filter(d_2001, abs(resid)> (abs(mean(d_2001$resid)) + (2*sd(d_2001$resid))) ) # outlier 2SD above the mean
filter(d_2001, abs(resid)> (abs(mean(d_2001$resid)) + (3*sd(d_2001$resid))) ) # outlier even 3SD above the mean

outlier_plot <- ggplot(data=d_2001, aes(x=fitted, y=resid))+
  geom_point(shape=1, size=2, col="blue")+
  geom_hline(yintercept=0)+
  labs(y="Residuals", x="Fitted")+
  scale_y_continuous(breaks=seq(-10,10, by=2))+ #scale_x_continuous(lim=c(3.6, 5.9), breaks=seq(3,7, by=0.5))+
  theme(axis.title.y=element_text(size=18, vjust=2), axis.title.x=element_text(size=18, vjust=-0.5),
        axis.text=element_text(size=16))+
  background_grid(major="xy")
outlier_plot
# ggsave(filename="results/outlier_2001.png", plot=outlier_plot, device="png", width=200, height=150, units="mm", dpi="print")


# -- Exclude outlier and rerun model
lm2 <- lmer(eggs_devtime ~ TempTreat*PhotTreat + (1|TubeID), data=filter(d_2001, abs(resid)< (abs(mean(d_2001$resid)) + (3*sd(d_2001$resid))) )) 
summary(lm2)
anova2 <- anova(lm2) #excl outlier doesn't change the result
anova2$mod <- "lm2"

ranova2 <- ranova(lm2) # random intercept significant
ranova2$mod <- "lm2"

lm_res2 <- summary(lm2)$coefficients %>% as.data.frame

plot(lm2) # residual variance much better now
qqnorm(resid(lm2))
qqline(resid(lm2))


#-- save results
# write.csv(rbind(lm_res1, lm_res2), file="results/output_lmer2001.csv", row.names=T)
# write.csv(rbind(anova1, anova2), file="results/anova_lmer2001.csv", row.names=T)
# write.csv(rbind(ranova1, ranova2), file="results/ranova_lmer2001.csv", row.names=T)



# Post-hoc test ####
#----------------------------------
summary(lm2) # final model = without outlier ####

posthoc_res1 <- emmeans::emmeans(lm2, pairwise ~ PhotTreat|TempTreat) # compare photoperiod treatments within temperature treatment
posthoc_res2 <- emmeans::emmeans(lm2, pairwise ~ TempTreat|PhotTreat) # compare temp treatments within photoperiod treatment

# write.csv(rbind(posthoc_res1$contrasts, posthoc_res2$contrasts), file="results/posthoc2001_contrasts.csv", row.names=F)
# write.csv(posthoc_res1$emmeans, file="results/posthoc2001_emmeans.csv", row.names=F)


# Plot random intercept variation (raw) ####
#----------------------------------
tube_devtime_avg <- Rmisc::summarySE(filter(d_2001, !is.na(eggs_devtime)), measurevar="eggs_devtime", groupvars=c("TubeID"))
tube_devtime_avg <- tube_devtime_avg %>% arrange(eggs_devtime) %>% mutate(Female=1:nrow(tube_devtime_avg))
summary(tube_devtime_avg$eggs_devtime)

d_2001 <- merge(d_2001, tube_devtime_avg[,c("TubeID", "Female")], by=c("TubeID"), all.x=TRUE)

tube_raw_devtime <- ggplot(data=tube_devtime_avg, aes(x=Female, y=eggs_devtime))+
  scale_colour_discrete(guide = "none")+ scale_shape_manual(values=c(17, 15))+
  geom_jitter(data=d_2001, aes(col=TubeID, shape=TempTreat), alpha=0.3, size=3, height=0, width=0.3)+ # jitter in width but not height
  geom_errorbar(aes(ymax = eggs_devtime+se, ymin=eggs_devtime-se), width=0.1, col="black") +
  geom_point(size=5, col="black", alpha=0.9) +
  labs(y="Egg development time (days)", x="Female")+
  scale_y_continuous(breaks=seq(-40,200, by=5))+ scale_x_continuous(breaks=seq(1,20, by=1))+ 
  theme(axis.title.y=element_text(size=18, vjust=2), axis.title.x=element_text(size=18, vjust=-0.5),
        axis.text.y=element_text(size=16), axis.text.x=element_text(size=12),
        legend.title = element_text(size=16), legend.text=element_text(size=15), legend.position = c(.85, .15))
tube_raw_devtime
# ggsave(filename="results/eggdevtime_raw_2001_perfemale.png", plot=tube_raw_devtime , device="png", width=230, height=150, units="mm", dpi="print")
