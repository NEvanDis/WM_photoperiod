# Analysis of winter moth egg photeriod experiment 2013-2014 ####
# Manipulated photoperiod the eggs received during development, kept in climate cabinets at a constant temp
# RQ: does photoperiod affect the timing of egg hatching?

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
d_2014 <- read.csv("1_data/db_EggPhot2014_d50.csv") # 2014 experiment d50 egg hatch dates per subclutch
head(d_2014)

table(is.na(d_2014$D50Calc)) # for 3 subclutches no D50 available, N=462 subclutches
length(unique(d_2014$TubeID)) # for 31 females
table(d_2014$AreaShortName) # all from OH area
table(d_2014$NovemberDate) # Catch dates
table(d_2014$Treatment) # in total 5 photoperiod treatments, 0 weeks = control, e.g. +2 weeks means photoperiod was manipulated to be 2 weeks later than the control treatment
summary(d_2014$Eggs) # subclutches of at least 12 eggs
aggregate(Eggs~TubeID, d_2014, sum) # range of clutch size


# -- Prep data for analysis
d_2014 <- d_2014 %>% mutate(Treatment=as.factor(Treatment), TubeID=as.factor(TubeID)) %>% # Treatment and TubeID(=MotherID) as factors
  mutate(Treatment=factor(Treatment, levels=c("-4 weeks", "-2 weeks", "0 weeks", "+2 weeks", "+4 weeks"))) %>%
  mutate(eggs_devtime= as.numeric(as.Date(D50Calc, origin=as.Date("2014-03-31"))-as.Date(NovemberDate, origin=as.Date("2013-10-31")))) # calculate egg development time as days from catch date of the female (=proxy of timing of egg laying)
head(d_2014)



# CHECK TEMP LOGGER DATA ####
#----------------------------------
boxes <- read.csv("1_data/db_EggPhot2014_boxes.csv") # data 2014
boxes

d_temp1 <- read.csv("1_data/db_EggPhot2014_temploggers_box1-5.csv")
d_temp2 <- read.csv("1_data/db_EggPhot2014_temploggers_box6-10.csv")
d_temp3 <- read.csv("1_data/db_EggPhot2014_temploggers_box11-15.csv")
d_temp <- rbind(d_temp1, d_temp2, d_temp3)
rm(d_temp1, d_temp2, d_temp3)
gc()

# -- Prep data
str(d_temp)
table(is.na(d_temp$Remarks))
d_temp <- d_temp %>% mutate(Box=as.factor(LocationName), DateTime=as.POSIXct(paste0(Year, "-", Month, "-", Day, " ", Hour, ":", Minute)))
d_temp <- merge(d_temp, boxes[,c("Location", "Treatment")], by.x="Box", by.y="Location", all.x=TRUE) %>%
  mutate(box=gsub("CT10-(box)(\\d+)","\\1\\2",Box)) %>% # Extract box info
  mutate(Box=ifelse(box==Box, gsub("CT10-(box).(\\d)","\\1\\2", Box), box)) %>%
  select(Box, Treatment, DateTime, LoggerName, Temperature) %>%
  filter(DateTime>as.Date("2014-1-15")) # start of the experiment
head(d_temp)
table(d_temp$Box)

# --- Use temp average from experiment start (Jan15) to median egg hatching date per treatment
medianhatchtreatment <- aggregate(D50Calc~Treatment, data=d_2014, median) # get median egg hatching date per treatment
medianhatchtreatment <- medianhatchtreatment %>% mutate(Date=as.Date(D50Calc, origin=as.Date("2014-03-31"))) #ignores the decimals

d_temp <- d_temp %>% group_split(Box)
d_temp <- lapply(d_temp, function(X) filter(X, as.Date(DateTime)<medianhatchtreatment$Date[medianhatchtreatment$Treatment==X$Treatment[1]])) # subset treatment on median egg hatching date
lapply(d_temp, nrow)
d_temp <- data.table::rbindlist(d_temp)

summary(d_temp$DateTime)


# -- Plot
temp <- ggplot(data=d_temp, aes(x=DateTime, y=Temperature)) +
  facet_wrap(~Box)+
  geom_line(aes(col=Treatment))+
  geom_smooth() + theme(axis.text.x=element_text(size=7))
temp
# ggsave(filename="results/temps2014.png", plot=temp , device="png", width=200, height=150, units="mm", dpi="print")


# -- Get averages
meantemps <- aggregate(Temperature~Box, data=d_temp, mean)
meantemps$sd <- aggregate(Temperature~Box, data=d_temp, sd)$Temperature
meantemps <- merge(meantemps, unique(d_temp[,c("Box", "Treatment")]), by="Box", all.x=TRUE)
meantemps1 <- aggregate(Temperature~Treatment, data=d_temp, mean)
meantemps1$sd <- aggregate(Temperature~Treatment, data=d_temp, sd)$Temperature
meantemps1$Box <- "combined"
meantemps <- rbind(meantemps, meantemps1) %>% arrange(Treatment, Box)
meantemps
rm(meantemps1)
# write.csv(meantemps, file="results/meantemps_2014.csv", row.names=FALSE)



#-----------------------------------
# Experiment analysis ####
#-----------------------------------

# Visualize hatching dates ####
#----------------------------------
devtime_avg <- Rmisc::summarySE(filter(d_2014, !is.na(eggs_devtime)), measurevar="eggs_devtime", groupvars=c("Treatment")) # average D50 hatch dates per photoperiod treatment
devtime_avg

raw_devtime <- ggplot(data=devtime_avg, aes(x=Treatment, y=eggs_devtime))+
  geom_jitter(data=d_2014, aes(col=Treatment), alpha=0.3, size=2.5, height=0, width=0.25)+ # jitter in width but not height
  geom_point(size=4, col="black") +
  geom_errorbar(aes(ymax = eggs_devtime+se, ymin=eggs_devtime-se), width=0.1, col="black") +
  labs(y="Egg development time (days)", x="Photoperiod treatment")+
  scale_y_continuous(breaks=seq(-40,200, by=5))+
  theme(axis.title.y=element_text(size=18, vjust=2), axis.title.x=element_text(size=18, vjust=-0.5),
        axis.text=element_text(size=16), legend.position="none")
raw_devtime
# ggsave(filename="results/eggdevtime_raw_2014.png", plot=raw_devtime , device="png", width=200, height=150, units="mm", dpi="print")


# -- Check out pattern: shorter dev times missing in -4 weeks variance
filter(d_2014, is.na(eggs_devtime)) # NAs are from different Treatments and different females, so not all from -4 weeks
min(filter(d_2014, Treatment=="-4 weeks")$eggs_devtime, na.rm=TRUE)

females <- filter(d_2014, eggs_devtime<93.5)$TubeID %>% unique() #check which females in other treatments have early egg hatch dates that are missing in -4 weeks Treatment
filter(d_2014, TubeID %in% females & Treatment=="-4 weeks") # 8 females in total

fems <- ggplot(data=filter(d_2014, TubeID %in% females), aes(x=Treatment, y=eggs_devtime))+
  geom_hline(yintercept = 93.4, alpha=0.6, linetype="dashed", linewidth=1)+
  geom_jitter(aes(col=Treatment), alpha=0.3, size=2.5, height=0, width=0.25)+ # jitter in width but not height
  labs(y="Egg development time (days)", x="Photoperiod treatment")+
  scale_y_continuous(lim=c(90,115), breaks=seq(-40,200, by=5))+
  theme(axis.title.y=element_text(size=18, vjust=2), axis.title.x=element_text(size=18, vjust=-0.5),
        axis.text=element_text(size=16), legend.position="none")
fems
# ggsave(filename="results/eggdevtime_raw_2014_earlyfemales.png", plot=fems , device="png", width=200, height=150, units="mm", dpi="print")



# Linear model ####
#----------------------------------
lm1 <- lmer(eggs_devtime ~ Treatment + (Treatment|TubeID), data=d_2014) #random effect: correct for split-brood design and allow for random main effect slopes
summary(lm1) # get singular boundary error, probably because random slope estimate is close to 0 (i.e. there is no such effect) -> correlations between random slope treatments are all 1 or -1
anova1 <- anova(lm1)
anova1$mod <- "lm1"

ranova1 <- ranova(lm1) # random slope not significant
ranova1$mod <- "lm1"
#MuMIn::r.squaredGLMM(lm1) 

lm2 <- lmer(eggs_devtime ~ Treatment + (1|TubeID), data=d_2014) # drop random slopes
summary(lm2)
anova2 <- anova(lm2) # do find an effect of photoperiod treatment
anova2$mod <- "lm2"

ranova2 <- ranova(lm2) # and random intercept effect
ranova2$mod <- "lm2"


# -- Final model
lm_final <- lm2
summary(lm_final)
lm_res <- summary(lm_final)$coefficients %>% as.data.frame

plot(lm_final) #equal variance? ok
qqnorm(resid(lm_final)) #normally distributed? okish
qqline(resid(lm_final))
 
# write.csv(lm_res, file="results/output_lmer2014.csv", row.names=T)
# write.csv(rbind(anova1, anova2), file="results/anova_lmer2014.csv", row.names=T)
# write.csv(rbind(ranova1, ranova2), file="results/ranova_lmer2014.csv", row.names=T)


# Plot random intercept variation (raw) ####
#----------------------------------
tube_devtime_avg <- Rmisc::summarySE(filter(d_2014, !is.na(eggs_devtime)), measurevar="eggs_devtime", groupvars=c("TubeID")) # get average egg development time per female
tube_devtime_avg <- tube_devtime_avg %>% arrange(eggs_devtime) %>% mutate(Female=1:nrow(tube_devtime_avg))
summary(tube_devtime_avg$eggs_devtime)

d_2014 <- merge(d_2014, tube_devtime_avg[,c("TubeID", "Female")], by=c("TubeID"), all.x=TRUE) # add averages to sub-clutch level data
d_2014 <- d_2014 %>% mutate(Photoperiod=ifelse(Treatment=="-4 weeks" | Treatment=="+4 weeks", as.character(Treatment),"")) # highlight most extreme photoperiod treatments in the figure

tube_raw_devtime <- ggplot(data=tube_devtime_avg, aes(x=Female, y=eggs_devtime))+
  scale_colour_discrete(guide="none")+ scale_shape_manual(values=c(3, 17, 15))+ scale_size_manual(values=c(2,3.5,3.5))+
  geom_jitter(data=d_2014, aes(col=TubeID, shape=Photoperiod, size=Photoperiod), alpha=0.3, height=0, width=0.25)+ # jitter in width but not height; size=2
  geom_errorbar(aes(ymax = eggs_devtime+se, ymin=eggs_devtime-se), width=0.1, col="black") +
  geom_point(size=5.5, col="black", alpha=0.9) + #size=4
  labs(y="Egg development time (days)", x="Female")+
  scale_y_continuous(breaks=seq(-40,200, by=5))+ scale_x_continuous(breaks=seq(1,31, by=1))+ 
  theme(axis.title.y=element_text(size=18, vjust=2), axis.title.x=element_text(size=18, vjust=-0.5),
        axis.text.y=element_text(size=16), axis.text.x=element_text(size=12),
        legend.title = element_text(size=16), legend.text=element_text(size=15), legend.position = c(.8, .15))
tube_raw_devtime
# ggsave(filename="results/eggdevtime_raw_2014_perfemale_wshape.png", plot=tube_raw_devtime , device="png", width=230, height=150, units="mm", dpi="print")


# Post-hoc test ####
#----------------------------------
summary(lm_final)

emm_treat <- emmeans::emmeans(lm_final, "Treatment") # compare between treatments
posthoc_res <- pairs(emm_treat)
# write.csv(posthoc_res, file="results/posthoc2014.csv", row.names=T)


# Ordered heterogeneity test ####
#----------------------------------

# Step 1: Get the estimates from the analysis summary panel
lm_res
lm_coef <- data.frame(fixEff="(Intercept)", coef=0) # set the reference level to zero = -4 weeks
lm_coef <- rbind(lm_coef, data.frame(fixEff=rownames(lm_res)[2:5], coef=lm_res$Estimate[2:5])) %>%
  mutate(Treatment=gsub("Treatment(.+)","\\1",fixEff)) %>% 
  mutate(Treatment=as.factor(ifelse(Treatment=="(Intercept)", "-4 weeks", Treatment))) %>%
  mutate(Treatment=factor(Treatment, levels=c("-4 weeks", "-2 weeks", "0 weeks", "+2 weeks", "+4 weeks"))) %>%
           arrange(Treatment)
rownames(lm_coef) <- 1:nrow(lm_coef)
lm_coef$order <- seq(1,5, by=1) # Order you want to test ####
lm_coef 


# Step 2: get the p-value from the ANOVA
pval <- anova2$`Pr(>F)` # use Treatment effect p-value


# Step 3: Calculation combining the value of the model comparison and the estimates  					

# -- Visualize:
plot <- ggplot(data=lm_coef, aes(x=order, y=coef))+ 
  geom_point(size=4)+
  labs(y="LM coefficient", x="Tested order")+
  scale_x_continuous(breaks=seq(1,5, by=1), labels=lm_coef$Treatment)+
  theme(axis.title = element_text(size=18), axis.text=element_text(size=16))
plot # NB: If ordered, should lay on a nice correlation line
# ggsave(filename="results/OH_test2014.png", plot=plot , device="png", width=200, height=150, units="mm", dpi="print")


# -- Calculate ordered heterogeneity test statistic: RsPc= [Rs(spearman rank correlation)*(1-pvalue obtained by the non-ordered test)]	
ordcor <- cor.test(lm_coef$order, lm_coef$coef, method="spearman")

OHtest <- ordcor$estimate*(1-pval)
names(OHtest) <- "OH test"
OHtest


# Step 4: look up the OHtest value into Table1 of Rice and Gaines 1994 (absolute value), and get the corresponding p value 
# (for a two tailed test you multiply by 2)
# Reference: Rice & Gaines (1994). The Ordered-Heterogeneity Family of Tests, Biometrics 50(3), https://doi.org/10.2307/2532788 ####
alpha <- c(0.001, 0.010, 0.025, 0.050, 0.1, 0.2, 0.3, 0.4, 0.5) # alphas taken from table in paper for Number of groups = 10
1-alpha # OHtest is negative, so alpha=1-alpha (i.e. left tail of the distribution)

# So for one-tailed test that we would observe this order by chance 5% of the time (i.e. P=0.05), critical value for OH test = 0.509
# and critical value for alpha 1% OH=0.369, so result = P>0.1
