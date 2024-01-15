# clean workspace
rm(list=ls())
gc()
graphics.off()

# load libraries
library(dplyr);library(lubridate);library(stringr)
library(ggplot2);library(lme4);library(lmerTest);library(effects)
library(sjPlot);library(tinytex);library(tidyverse);library(ggthemes)
library(bayestestR); library(tidybayes)

Sys.setenv(TZ='UTC') # change time zone for session
options(scipen = 999)

# Read in data and filter
#setwd() #set working directory

load("Patch_Exploration_DayNight.RData")

#summary stats
x <- df3 %>% group_by(tag) %>% summarise(n_per_ind = n())
summary(x$n_per_ind)
sd(x$n_per_ind)

y <- df3 %>% group_by(tod_start, tideID) %>% dplyr::summarise(sum_duration = (mean(duration_s)/60/60))
y %>% filter(tod_start == "day") %>% summary(sum_duration)
y %>% filter(tod_start == "night") %>% summary(sum_duration)


#multivariate mixed effect models ####
library(lme4);library(arm);library(MuMIn);library(tidyverse)
library(plyr);library(broom);library(coda);library(grid)
library(gridExtra);library(brms); library(broom.mixed); library(merTools);
library(tidybayes);library(parallel)
my.cores <- detectCores()

df3$duration_min_log <- log10(df3$duration_min)
df3$dist_in_patch_m_log <- log10(df3$dist_in_patch_m)
df3$disp_in_patch_m_log <- log10(1+df3$disp_in_patch_m)
df3$year <- as.factor(df3$year)

brm_1 <- bf(mvbind(scale(dist_in_patch_m_log),scale(disp_in_patch_m_log),scale(duration_min_log)) ~
               explor_Log10_mean_speed * tod_start +
               night_day_ratio +
               scale_space +  inout + 
               (1|p|tag) + (1|tideID) + (1|year)) + set_rescor(TRUE)

# brm <- brm(brm_1, data = df3, 
#         family = "gaussian",
#         control = list(adapt_delta=0.99),
#         warmup = 500,
#         iter = 1000, thin=2,
#         chains = 1, #inits = "random",
#         seed = 12345,
#         cores = my.cores) 

#save(brm, file= "brm_2022-07-07-1chain.Rdata")
load(file= "brm_2022-07-07-1chain.Rdata") #result of the above model
#as.mcmc(brm)
#brm <- add_criterion(brm, "loo")
summary(brm)
tab_model(brm)
bayes_R2(brm)

#check fit
p1<- pp_check(brm, resp = "scaledistinpatchmlog")
p2<- pp_check(brm, resp = "scaledispinpatchmlog")
p3<- pp_check(brm, resp = "scaledurationminlog")

conditional_effects(brm, "explor_Log10_mean_speed:tod_start", resp = "scaledistinpatchmlog")
conditional_effects(brm, "explor_Log10_mean_speed:tod_start", resp = "scaledispinpatchmlog")
conditional_effects(brm, "explor_Log10_mean_speed:tod_start", resp = "scaledurationminlog")

library(emmeans)
emmeans(brm,  ~ explor_Log10_mean_speed * tod_start, interaction = c("pairwise","pairwise"))

emtrends(brm2, pairwise ~ tod_start, var = "explor_Log10_mean_speed")

emmip(brm, tod_start ~ explor_Log10_mean_speed, cov.reduce = range)
emmip(brm2, tod_start ~ explor_Log10_mean_speed, cov.reduce = range)
emmip(brm3, tod_start ~ explor_Log10_mean_speed, cov.reduce = range)


conditional_effects(brm, 
                    effects = "explor_Log10_mean_speed:tod_start", 
                    resp = "scaledistinpatchmlog",
                    spaghetti = T, nsamples = 150) %>% 
  plot(points = T,
       point_args = c(alpha = 1/3, size = 1), 
       mean = T, 
       rug = T, 
       theme= theme_classic())

conditional_effects(brm, 
                    effects = "explor_Log10_mean_speed:tod_start", 
                    resp = "scaledispinpatchmlog",
                    spaghetti = T, nsamples = 150) %>% 
  plot(points = T,
       point_args = c(alpha = 1/3, size = 1), 
       mean = T, 
       rug = T, 
       theme= theme_classic())

conditional_effects(brm, 
                    effects = "explor_Log10_mean_speed:tod_start", 
                    resp = "scaledurationminlog",
                    spaghetti = T, nsamples = 150) %>% 
  plot(points = T,
       point_args = c(alpha = 1/3, size = 1), 
       mean = T, 
       rug = T, 
       theme= theme_classic())

conditional_effects(brm2, 
                    effects = "explor_Log10_mean_speed:tod_start", 
                    resp = "n_log",
                    spaghetti = T, nsamples = 150) %>% 
  plot(points = T,
       point_args = c(alpha = 1/3, size = 1), 
       mean = T, 
       rug = T, 
       theme= theme_classic())

conditional_effects(brm3, 
                    effects = "explor_Log10_mean_speed:tod_start", 
                    resp = "dist_bw_patch2_log",
                    spaghetti = T, nsamples = 150) %>% 
  plot(points = T,
       point_args = c(alpha = 1/3, size = 1), 
       mean = T, 
       rug = T, 
       theme= theme_classic())


#plot aggregated speed ####
#plot speed 
df3$speed_in_patch_mm = df3$dist_in_patch_m/df3$duration_min

df <- df3 %>% 
  mutate(expl = case_when(explor_Log10_mean_speed < 0.4 ~ "Slower", 
                          explor_Log10_mean_speed >= 0.4 & explor_Log10_mean_speed <= 0.7 ~ "middle",
                          explor_Log10_mean_speed > 0.7 ~ "Faster")) 
df2 <- df %>% group_by(tide_cont,expl) %>% summarise(speed_in_patch_ms= mean(speed_in_patch_ms),
                                                     sum_duration_m = sum(sum_duration))

df2 <- df2 %>% mutate(period = case_when(tide_cont <= 445 ~ "NA",
                                        tide_cont > 445 & tide_cont <= 560 ~ "Before day 560",
                                        tide_cont > 560 ~ "After day 560"))
#discussion figure
df2 %>% 
  filter(expl != "middle" & period != "NA") %>% 
  mutate(across(period, factor, levels= c("Before day 560", "After day 560"))) %>% 
  ggplot(aes(y=sum_duration_m, x=expl, fill=expl, color=expl)) +
  geom_point(position="jitter", alpha=0.4, aes(shape=expl)) + geom_boxplot(alpha=0.6) + 
  theme_classic() + 
  facet_wrap(~period) +
  scale_color_manual(values=c("#EDCB64","#7496D2")) +
  scale_fill_manual(values=c("#EDCB64","#7496D2")) +
  ylab("Summed duration tracked in residence patches (h)") + 
  xlab("Exploration speed measured in captivity") +
  theme(axis.title.y = element_text(size=14), axis.title.x = element_text(size=14), 
        legend.background = element_rect(), legend.position= "none") +
  theme(plot.background = element_rect(colour = "black", fill=NA, size=1))

#discussion figure
ggplot(data= df, aes(x = tide_cont, y = (speed_in_patch_ms))) + 
  geom_point(data = subset(df, expl != "middle"), aes(col=expl),alpha=0.3, size=.05) +
  geom_line(data= subset(df2, expl != "middle"), 
            aes(group=as.factor(expl),
                col=as.factor(expl)),alpha=1, size=1) +
  theme_classic() + 
  ylab("Speed in a patch (m/s)")+ 
  scale_color_manual(values=c("#EDCB64","#7496D2")) +
  xlab("Tidal time") +
  theme(legend.position= "none") +
  xlim(445,600) + 
  ylim(0, 40)

#numbers ####
#length(unique(df4$tag))
#length(unique(df2$tideID))

length(unique(df3$tag))
length(unique(df3$tideID))

df3 %>%  group_by(year, month) %>%
  summarise(sum_duration = sum(duration_s)/60/60,
            n=n(),
            ntag = length(unique(tag)))

#day
df_tide_day <- df3 %>% filter(tod_start == "day") 
length(unique(df_tide_day$tag))
length(unique(df_tide_day$tideID))

#night 
df_tide_night <- df3 %>% filter(tod_start == "night") 
length(unique(df_tide_night$tag))
length(unique(df_tide_night$tideID))
df_tide_night %>% group_by(year) %>%
  summarise(sum_duration = sum(duration_s)/60/60,
            n=n(),
            ntag = length(unique(tag)))

#duration
df3 %>% ggplot() +
  geom_histogram(aes(x=duration_min), bins= 100) 
summary(df3$duration_min)

duration <- df3 %>% 
  group_by(tag, year, explor_Log10_mean_speed, tideID, tod_start) %>%
  summarise(sum_duration = sum(duration_s)/60/60)

duration %>% 
  ggplot(aes(x=explor_Log10_mean_speed, y=sum_duration)) + 
  facet_wrap(.~tod_start + year) +
  geom_point(alpha=0.2) + theme_bw() + 
  ylab("Sum of residence patch durations per tide (h)") + 
  xlab("Exploration speed (log10 cm/s)")

summary(lmer(data= duration, 
             sum_duration ~ explor_Log10_mean_speed*tod_start + 
               (1|tag) + (1|tideID) + (1|year)))
df3 %>% 
  group_by(year, tod_start) %>%
  summarise(sum_duration = sum(duration_s)/60/60,
            ntag = length(unique(tag)),
            ntide = length(unique(tideID)))
#

# tracked birds histogram ####
h <- read.csv("Hist.csv")

#Suppl. fig.
ggplot(h, aes(x=Score, group=Type, fill=Type)) +
   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins=50) +
   theme_classic() +
   xlab("Exploration speed (log10 cm/s)") +
   ylab("Number of individuals") +
   theme(legend.position= c(0.23, 0.95)) +
   scale_fill_manual(values=c("#999999", "#E69F00"),
                     name="",
                     breaks=c("all", "tracked"),
                     labels=c("All knots", "Knots tracked on Griend mudflats")) +
   theme(axis.title.y = element_text(size=16), axis.title.x = element_text(size=15))+
   theme(legend.title = element_blank(),
         panel.grid = element_blank(),
         panel.border =element_rect(colour = "black", fill=NA, size=0.5),
         legend.key = element_blank(),
         strip.text.x = element_text(size=17),
         strip.background =element_blank(),
         plot.title = element_text(hjust = 0.5, size= 20))

# space by water level ####
wl <- read.csv("spacebywaterlevel.csv")

ggplot(wl, aes(x = waterlevel, y = area_km2))+
  geom_point(col = "blue")+
  scale_x_continuous(expand = c(0,0), limits = c(-150,150), n.breaks = 6)+
  scale_y_continuous(expand = c(0,0),limits = c(0,6), n.breaks = 7)+
  theme_classic() +  
  ylab("Area (km2)") +
  xlab("Water level (m)") + 
  theme(plot.background = element_rect(colour = "black", fill=NA, size=1)) 
  


df2 <- left_join(df2,wl, by=c("waterlevel_start_cm"="waterlevel"))
str(df2)
# repeatability
## Do exploratory birds switch patches more often? ####
duration <- df3 %>% 
  group_by(tag, year, night_day_ratio, explor_Log10_mean_speed, tideID, tod_start) %>%
  dplyr::summarise(sum_duration = sum(duration_s)/60/60)

df_3h <-  merge(df3, duration[,c("tag","tideID", "night_day_ratio", "tod_start", "sum_duration")], by=c("tag", "tideID", "tod_start"), all.x=T, all.y=T) #%>% 
  #filter(sum_duration >= 3)

patches_per_tide <- df_3h %>%
  #select(tag, tideID, tod_start, patch, sum_duration) %>% 
  group_by(tag, explor_Log10_mean_speed, tideID, tod_start)%>% # make sure unique id and explor score are the same otherwise will have to merge separately.
  dplyr::summarise(n = n()) 

patches_per_tide <-  merge(patches_per_tide, duration[,c("tag","tideID", "night_day_ratio", "tod_start", "sum_duration")], 
                by=c("tag", "tideID", "tod_start"), all.x=T, all.y=T) 

patches_per_tide <- patches_per_tide %>%
  mutate(year = factor(substr(tideID, 1, 4)),
         tide_cont = as.numeric(substr(tideID, 5, 7)),
         n_per_tide = n/sum_duration)

summary(factor(patches_per_tide$tod_start))
hist((patches_per_tide$n))

#patches_per_tide <- patches_per_tide %>%  filter(sum_duration > 1)

#count
# patches_per_tide2 <- patches_per_tide %>% filter(explor_Log10_mean_speed <1,n>1, sum_duration > 1)
# 
# m_npatches <- glmer(data=patches_per_tide2, 
#                    n ~ 
#                      explor_Log10_mean_speed * tod_start + #tide_cont +
#                      (1|tag) + (1|tideID)+ (1|year), 
#                    family = "poisson") 

patches_per_tide$n_log <- log10(patches_per_tide$n_per_tide) ### Residuals don't look great. Use log. 

brm_2 <- bf(n_log ~
              explor_Log10_mean_speed * tod_start + 
              night_day_ratio +
              (1|tag) + (1|tideID)+ (1|year))

brm2 <- brm(brm_2, data = patches_per_tide, 
           family = "gaussian",
           control = list(adapt_delta=0.99),
           warmup = 500,
           iter = 1000, thin=2,
           chains = 2, #inits = "random",
           seed = 12345,
           cores = my.cores) 
summary(brm2)

#save(brm2, file= "C:/Users/sersoy/OneDrive - NIOZ/PhD/PhD_chapters/Chapter 3- Personality and patch movement/Data and R/Final_2022-06-01/brm_PATCHVISITS_2022-07-07-2chain.Rdata")
load(file= "brm_PATCHVISITS_2022-07-07-2chain.Rdata")

p5 <- pp_check(brm2, resp = "n_log")
# m_npatches <- lmer(data=patches_per_tide, 
#                    n_log ~ 
#                      explor_Log10_mean_speed * tod_start + 
#                      night_day_ratio +
#                      (1|tag) + (1|tideID)+ (1|year) ) 

#hist(residuals(m_npatches))
plot(m_npatches)
tab_model(m_npatches)
plot_model(m_npatches, type = "pred", terms= c("explor_Log10_mean_speed", "tod_start")) +
  theme_bw()

library(emmeans)
library(rstanarm)
emmeans(brm2, pairwise ~ explor_Log10_mean_speed * tod_start )
#emmeans(m_npatches,  ~ explor_Log10_mean_speed * tod_start, interaction = c("pairwise","pairwise"))

interact_plot(m_npatches, pred = explor_Log10_mean_speed, modx = tod_start, 
              partial.residuals = T,
              interval = TRUE, int.width = 0.95, 
              x.label = "Exploration speed (log10 cm/s)",
              y.label = "Number of patch visits per tide (log10)",
              legend.main = "Day-Night",
              modx.labels = c("Day","Night"),
              plot.points = T,
               point.size = 0.8,
               point.alpha = 0.1,
              line.thickness = 1.1,
              # pred.point.size = T,
              #rug= T,
              colors = c("orange","black")) +
  theme_classic() + 
  theme(legend.position= "none") +
  ylim(-0.05, 1.5) 

plot(predictorEffect("explor_Log10_mean_speed", m_npatches, residuals=TRUE))

e1.mm1 <- predictorEffect("explor_Log10_mean_speed", m_npatches)
plot(e1.mm1, lines=list(multiline=TRUE, col=c("orange","black")), 
     confint=list(style="auto"),
     axes=list(y=list(type="response", 
                      lab="Number of patch visits per tide (log10)")),
     xlab="Exploration speed (log10 cm/s)",
     type="l", main="")

plot(effect(c("explor_Log10_mean_speed", "tod_start"), m_npatches))

newdf <- expand.grid(explor_Log10_mean_speed=seq(-0.15,1.25,by=0.2),
                     tod_start = levels(factor(patches_per_tide$tod_start)),
                     tag=levels(factor(patches_per_tide$tag)),
                     tideID=levels(factor(patches_per_tide$tideID)),
                     night_day_ratio = mean(patches_per_tide$night_day_ratio),
                     year=levels(factor(patches_per_tide$year)))

newdf$n_log <- predict(m_npatches,newdata=newdf,type="response")

tidal2 <- patches_per_tide

tidal2$expl_round <- format(round(tidal2$explor_Log10_mean_speed, 2), nsmall = 2)

tidal <- tidal2 %>% 
  group_by(expl_round, tod_start) %>% 
  summarise(se = sd(n_log)/sqrt(length(n_log)),
            n_log = mean(n_log)) #%>% 
            #expl_round = format(round(explor_Log10_mean_speed, 2), nsmall = 2),
            #explor_Log10_mean_speed=unique(expl_round)) %>% 
  #mutate(explor_Log10_mean_speed = as.numeric(explor_Log10_mean_speed)) 

#tidal <-tidal[!duplicated(tidal), ]

colnames(tidal)[1] <- "explor_Log10_mean_speed"
tidal$explor_Log10_mean_speed <- as.numeric(tidal$explor_Log10_mean_speed)

library(ggthemes)
ggplot(patches_per_tide, aes(x=explor_Log10_mean_speed, y=n_log, color=as.factor(tod_start))) +
  geom_errorbar(aes(ymin=n_log-se, ymax=n_log+se),
                data=tidal, size=0.1, width=0, alpha=0.5) +
  geom_point(data=tidal, aes(x=explor_Log10_mean_speed, y=n_log),
             size=1, alpha= 0.5) +
  scale_colour_manual(values=c("orange","black")) +
  ylab("Number of patch visits per tide (log10)")+ 
  geom_smooth(method="lm", se=TRUE, data= newdf, 
              aes(x=explor_Log10_mean_speed, y=n_log),
              size=1) +
  theme_base() + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + #, limits = c(0.1,0.7)) +
  xlab("Exploration speed (log10 cm/s)") +
  theme(legend.position= "none") 

tidal2 <- tidal2 %>% 
  mutate(expl_round = case_when(
    explor_Log10_mean_speed < -0.15 ~ -0.15,
    explor_Log10_mean_speed >= -0.15 & explor_Log10_mean_speed < -0.175 ~ -0.175,
    explor_Log10_mean_speed >= -0.175 & explor_Log10_mean_speed < -0.1 ~ -0.1,
    explor_Log10_mean_speed >= -0.1 & explor_Log10_mean_speed < 0 ~ 0,
    explor_Log10_mean_speed >= 0 & explor_Log10_mean_speed < 0.025 ~ 0.025,
    explor_Log10_mean_speed >= 0.025 & explor_Log10_mean_speed < 0.05 ~ 0.05,
    explor_Log10_mean_speed >= 0.05 & explor_Log10_mean_speed < 0.075 ~ 0.075,
    explor_Log10_mean_speed >= 0.075 & explor_Log10_mean_speed < 0.1 ~ 0.1,
    explor_Log10_mean_speed >= 0.1 & explor_Log10_mean_speed < 0.125 ~ 0.125,
    explor_Log10_mean_speed >= 0.125 & explor_Log10_mean_speed < 0.15 ~ 0.15,
    explor_Log10_mean_speed >= 0.15 & explor_Log10_mean_speed < 0.175 ~ 0.175,
    explor_Log10_mean_speed >= 0.175 & explor_Log10_mean_speed < 0.2 ~ 0.2,
    explor_Log10_mean_speed >= 0.2 & explor_Log10_mean_speed < 0.225 ~ 0.225,
    explor_Log10_mean_speed >= 0.225 & explor_Log10_mean_speed < 0.25 ~ 0.25,
    explor_Log10_mean_speed >= 0.25 & explor_Log10_mean_speed < 0.275 ~ 0.275,
    explor_Log10_mean_speed >= 0.275 & explor_Log10_mean_speed < 0.3 ~ 0.3,
    explor_Log10_mean_speed >= 0.3 & explor_Log10_mean_speed < 0.325 ~ 0.325,
    explor_Log10_mean_speed >= 0.325 & explor_Log10_mean_speed < 0.35 ~ 0.35,
    explor_Log10_mean_speed >= 0.35 & explor_Log10_mean_speed < 0.375 ~ 0.375,
    explor_Log10_mean_speed >= 0.375 & explor_Log10_mean_speed < 0.4 ~ 0.4,
    explor_Log10_mean_speed >= 0.4 & explor_Log10_mean_speed < 0.425 ~ 0.425,
    explor_Log10_mean_speed >= 0.425 & explor_Log10_mean_speed < 0.45 ~ 0.45,
    explor_Log10_mean_speed >= 0.45 & explor_Log10_mean_speed < 0.475 ~ 0.475,
    explor_Log10_mean_speed >= 0.475 & explor_Log10_mean_speed < 0.5 ~ 0.5,
    explor_Log10_mean_speed >= 0.5 & explor_Log10_mean_speed < 0.525 ~ 0.525,
    explor_Log10_mean_speed >= 0.525 & explor_Log10_mean_speed < 0.55 ~ 0.55,
    explor_Log10_mean_speed >= 0.55 & explor_Log10_mean_speed < 0.575 ~ 0.575,
    explor_Log10_mean_speed >= 0.575 & explor_Log10_mean_speed < 0.6 ~ 0.6,
    explor_Log10_mean_speed >= 0.6 & explor_Log10_mean_speed < 0.625 ~ 0.625,
    explor_Log10_mean_speed >= 0.625 & explor_Log10_mean_speed < 0.65 ~ 0.65,
    explor_Log10_mean_speed >= 0.65 & explor_Log10_mean_speed < 0.675 ~ 0.675,
    explor_Log10_mean_speed >= 0.675 & explor_Log10_mean_speed < 0.7 ~ 0.7,
    explor_Log10_mean_speed >= 0.7 & explor_Log10_mean_speed < 0.725 ~ 0.725,
    explor_Log10_mean_speed >= 0.725 & explor_Log10_mean_speed < 0.75 ~ 0.75,
    explor_Log10_mean_speed >= 0.75 & explor_Log10_mean_speed < 0.775 ~ 0.775,
    explor_Log10_mean_speed >= 0.775 & explor_Log10_mean_speed < 0.8 ~ 0.8,
    explor_Log10_mean_speed >= 0.8 & explor_Log10_mean_speed < 0.825 ~ 0.825,
    explor_Log10_mean_speed >= 0.825 & explor_Log10_mean_speed < 0.85 ~ 0.85,
    explor_Log10_mean_speed >= 0.85 & explor_Log10_mean_speed < 0.875 ~ 0.875,
    explor_Log10_mean_speed >= 0.875 & explor_Log10_mean_speed < 0.9 ~ 0.9,
    explor_Log10_mean_speed >= 0.9 & explor_Log10_mean_speed < 0.925 ~ 0.925,
    explor_Log10_mean_speed >= 0.925 & explor_Log10_mean_speed < 0.95 ~ 0.95,
    explor_Log10_mean_speed >= 0.95 & explor_Log10_mean_speed < 0.975 ~ 0.975,
    explor_Log10_mean_speed >= 0.975 & explor_Log10_mean_speed < 1 ~ 1,
    explor_Log10_mean_speed >= 1 & explor_Log10_mean_speed < 1.025 ~ 1.025,
    explor_Log10_mean_speed >= 1.025 & explor_Log10_mean_speed < 1.05 ~ 1.05,
    explor_Log10_mean_speed >= 1.05 & explor_Log10_mean_speed < 1.075 ~ 1.075,
    explor_Log10_mean_speed >= 1.075 & explor_Log10_mean_speed < 1.1 ~ 1.1,
    explor_Log10_mean_speed >= 1.1 & explor_Log10_mean_speed < 1.125 ~ 1.125,
    explor_Log10_mean_speed >= 1.125 & explor_Log10_mean_speed < 1.15 ~ 1.15,
    explor_Log10_mean_speed >= 1.15 & explor_Log10_mean_speed < 1.175 ~ 1.175,
    explor_Log10_mean_speed >= 1.175 & explor_Log10_mean_speed < 1.2 ~ 1.2,
    explor_Log10_mean_speed <= 1.25 ~ 1.25
  ))

tidal <- tidal2 %>% 
  group_by(expl_round, tod_start) %>% 
  summarise(se = sd(n_log)/sqrt(length(n_log)),
            n_log = mean(n_log)) #%>% 
#expl_round = format(round(explor_Log10_mean_speed, 2), nsmall = 2),
#explor_Log10_mean_speed=unique(expl_round)) %>% 
#mutate(explor_Log10_mean_speed = as.numeric(explor_Log10_mean_speed)) 

#tidal <-tidal[!duplicated(tidal), ]

colnames(tidal)[1] <- "explor_Log10_mean_speed"
tidal$explor_Log10_mean_speed <- as.numeric(tidal$explor_Log10_mean_speed)


#newdf1 <- newdf %>% mutate(expl= case_when(explor_Log10_mean_speed < 0 ~ "Slower-explorers",
#                                          explor_Log10_mean_speed > 0.9 ~ "Faster-explorers"))
#newdf1 <- na.omit(newdf1)

ggplot(data=newdf1, aes(y= n_log, x=tod_start, col=expl)) + 
  geom_boxplot(notch =T, outlier.alpha=0.006) +
  theme_classic() + 
  #facet_wrap(~inout) +
  scale_color_manual(values = c("seagreen4", "orange2"),
                     name="") +
  ylab("Number of patch visits (log10)") + 
  theme(axis.title.y = element_text(size=13),
        legend.title = element_blank(),
        #legend.position = "none",
        panel.grid = element_blank(),
        panel.border =element_rect(colour = "black", fill=NA, size=0.5),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=12),
        strip.background =element_blank()) +
  scale_x_discrete(labels= labels)

## Do exploratory birds faster in a patch? ####
#ggplot(data=df3, aes(col=dist_in_patch_m, x=nfixes, y=duration_min)) + 
#  geom_point(alpha= 0.2) +
#  scale_color_viridis_c() +theme_bw()

#spd <- df3 %>% filter( speed_in_patch_ms < 0.7)
hist((df3$speed_median_ms))

df3$speed_ms_log <- log10(df3$speed_in_patch_ms)

m_spd <- lmer(data=df3, 
              speed_ms_log ~ 
                     explor_Log10_mean_speed * tod_start +  
                     scale_space + inout +
                     night_day_ratio + 
                     tide_cont + 
                     scale(nfixes) + 
                     (1|tag) + (1|tideID) + (1|year),
                     weight = sum_duration)
#plot(m_spd)
tab_model(m_spd)
plot_model(m_spd, type = "pred", terms= c("explor_Log10_mean_speed", "tod_start")) 

interact_plot(m_spd, pred = explor_Log10_mean_speed, modx = tod_start, 
              partial.residuals = T,
              interval = TRUE, int.width = 0.95, 
              x.label = "Exploration speed (log10 cm/s)",
              y.label = "Speed in a residence patch (log10 m/s)",
              legend.main = "Day-Night",
              modx.labels = c("Day","Night"),
              plot.points = T,
              point.size = 0.3,
              point.alpha = 0.1,
              pred.point.size = T,
              #rug= T,
              colors = c("orange","black")) +
  theme_classic()

newdf <- expand.grid(explor_Log10_mean_speed=seq(-0.15,1.25,by=0.2),
                     tod_start = levels(factor(df3$tod_start)),
                     tag=levels(factor(df3$tag)),
                     tideID=levels(factor(df3$tideID)),
                     night_day_ratio = mean(df3$night_day_ratio),
                     scale_space= -0.02371, #mean
                     inout = levels(factor(df3$inout)),
                     year=levels(factor(df3$year)))

newdf$speed_ms_log <- predict(m_spd,newdata=newdf,type="response")

tidal2 <- df3

tidal2 <- tidal2 %>% 
  mutate(expl_round = case_when(
    explor_Log10_mean_speed < -0.15 ~ -0.15,
    explor_Log10_mean_speed >= -0.15 & explor_Log10_mean_speed < -0.175 ~ -0.175,
    explor_Log10_mean_speed >= -0.175 & explor_Log10_mean_speed < -0.1 ~ -0.1,
    explor_Log10_mean_speed >= -0.1 & explor_Log10_mean_speed < 0 ~ 0,
    explor_Log10_mean_speed >= 0 & explor_Log10_mean_speed < 0.025 ~ 0.025,
    explor_Log10_mean_speed >= 0.025 & explor_Log10_mean_speed < 0.05 ~ 0.05,
    explor_Log10_mean_speed >= 0.05 & explor_Log10_mean_speed < 0.075 ~ 0.075,
    explor_Log10_mean_speed >= 0.075 & explor_Log10_mean_speed < 0.1 ~ 0.1,
    explor_Log10_mean_speed >= 0.1 & explor_Log10_mean_speed < 0.125 ~ 0.125,
    explor_Log10_mean_speed >= 0.125 & explor_Log10_mean_speed < 0.15 ~ 0.15,
    explor_Log10_mean_speed >= 0.15 & explor_Log10_mean_speed < 0.175 ~ 0.175,
    explor_Log10_mean_speed >= 0.175 & explor_Log10_mean_speed < 0.2 ~ 0.2,
    explor_Log10_mean_speed >= 0.2 & explor_Log10_mean_speed < 0.225 ~ 0.225,
    explor_Log10_mean_speed >= 0.225 & explor_Log10_mean_speed < 0.25 ~ 0.25,
    explor_Log10_mean_speed >= 0.25 & explor_Log10_mean_speed < 0.275 ~ 0.275,
    explor_Log10_mean_speed >= 0.275 & explor_Log10_mean_speed < 0.3 ~ 0.3,
    explor_Log10_mean_speed >= 0.3 & explor_Log10_mean_speed < 0.325 ~ 0.325,
    explor_Log10_mean_speed >= 0.325 & explor_Log10_mean_speed < 0.35 ~ 0.35,
    explor_Log10_mean_speed >= 0.35 & explor_Log10_mean_speed < 0.375 ~ 0.375,
    explor_Log10_mean_speed >= 0.375 & explor_Log10_mean_speed < 0.4 ~ 0.4,
    explor_Log10_mean_speed >= 0.4 & explor_Log10_mean_speed < 0.425 ~ 0.425,
    explor_Log10_mean_speed >= 0.425 & explor_Log10_mean_speed < 0.45 ~ 0.45,
    explor_Log10_mean_speed >= 0.45 & explor_Log10_mean_speed < 0.475 ~ 0.475,
    explor_Log10_mean_speed >= 0.475 & explor_Log10_mean_speed < 0.5 ~ 0.5,
    explor_Log10_mean_speed >= 0.5 & explor_Log10_mean_speed < 0.525 ~ 0.525,
    explor_Log10_mean_speed >= 0.525 & explor_Log10_mean_speed < 0.55 ~ 0.55,
    explor_Log10_mean_speed >= 0.55 & explor_Log10_mean_speed < 0.575 ~ 0.575,
    explor_Log10_mean_speed >= 0.575 & explor_Log10_mean_speed < 0.6 ~ 0.6,
    explor_Log10_mean_speed >= 0.6 & explor_Log10_mean_speed < 0.625 ~ 0.625,
    explor_Log10_mean_speed >= 0.625 & explor_Log10_mean_speed < 0.65 ~ 0.65,
    explor_Log10_mean_speed >= 0.65 & explor_Log10_mean_speed < 0.675 ~ 0.675,
    explor_Log10_mean_speed >= 0.675 & explor_Log10_mean_speed < 0.7 ~ 0.7,
    explor_Log10_mean_speed >= 0.7 & explor_Log10_mean_speed < 0.725 ~ 0.725,
    explor_Log10_mean_speed >= 0.725 & explor_Log10_mean_speed < 0.75 ~ 0.75,
    explor_Log10_mean_speed >= 0.75 & explor_Log10_mean_speed < 0.775 ~ 0.775,
    explor_Log10_mean_speed >= 0.775 & explor_Log10_mean_speed < 0.8 ~ 0.8,
    explor_Log10_mean_speed >= 0.8 & explor_Log10_mean_speed < 0.825 ~ 0.825,
    explor_Log10_mean_speed >= 0.825 & explor_Log10_mean_speed < 0.85 ~ 0.85,
    explor_Log10_mean_speed >= 0.85 & explor_Log10_mean_speed < 0.875 ~ 0.875,
    explor_Log10_mean_speed >= 0.875 & explor_Log10_mean_speed < 0.9 ~ 0.9,
    explor_Log10_mean_speed >= 0.9 & explor_Log10_mean_speed < 0.925 ~ 0.925,
    explor_Log10_mean_speed >= 0.925 & explor_Log10_mean_speed < 0.95 ~ 0.95,
    explor_Log10_mean_speed >= 0.95 & explor_Log10_mean_speed < 0.975 ~ 0.975,
    explor_Log10_mean_speed >= 0.975 & explor_Log10_mean_speed < 1 ~ 1,
    explor_Log10_mean_speed >= 1 & explor_Log10_mean_speed < 1.025 ~ 1.025,
    explor_Log10_mean_speed >= 1.025 & explor_Log10_mean_speed < 1.05 ~ 1.05,
    explor_Log10_mean_speed >= 1.05 & explor_Log10_mean_speed < 1.075 ~ 1.075,
    explor_Log10_mean_speed >= 1.075 & explor_Log10_mean_speed < 1.1 ~ 1.1,
    explor_Log10_mean_speed >= 1.1 & explor_Log10_mean_speed < 1.125 ~ 1.125,
    explor_Log10_mean_speed >= 1.125 & explor_Log10_mean_speed < 1.15 ~ 1.15,
    explor_Log10_mean_speed >= 1.15 & explor_Log10_mean_speed < 1.175 ~ 1.175,
    explor_Log10_mean_speed >= 1.175 & explor_Log10_mean_speed < 1.2 ~ 1.2,
    explor_Log10_mean_speed <= 1.25 ~ 1.25
  ))

tidal <- tidal2 %>% 
  group_by(expl_round, tod_start) %>% 
  summarise(se = sd(speed_ms_log)/sqrt(length(speed_ms_log)),
            speed_ms_log = mean(speed_ms_log)) #%>% 

colnames(tidal)[1] <- "explor_Log10_mean_speed"
tidal$explor_Log10_mean_speed <- as.numeric(tidal$explor_Log10_mean_speed)

library(ggthemes)
ggplot(df3, aes(x=explor_Log10_mean_speed, y=speed_ms_log, color=as.factor(inout))) +
  geom_errorbar(aes(ymin=speed_ms_log-se, ymax=speed_ms_log+se),
                data=tidal, size=0.1, width=0, alpha=0.5) +
  geom_point(data=tidal, aes(x=explor_Log10_mean_speed, y=speed_ms_log),
             size=1, alpha= 0.5) +
  #scale_colour_manual(values=c("orange","black")) +
  scale_colour_manual(values=c("#999999","#0072B2")) +
  facet_wrap(~tod_start) +
  ylab("Speed in patch (log10 m/s)")+ 
  geom_smooth(method="lm", se=TRUE, data= newdf, 
              aes(x=explor_Log10_mean_speed, y=speed_ms_log),
              size=1) +
  theme_base() + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(-0.8,-0.6)) +
  xlab("Exploration speed (log10 cm/s)") +
  theme(legend.position= "none") 

ggplot(df3, aes(x=explor_Log10_mean_speed, y=speed_ms_log, color=as.factor(tod_start))) +
  geom_errorbar(aes(ymin=speed_ms_log-se, ymax=speed_ms_log+se),
                data=tidal, size=0.1, width=0, alpha=0.5) +
  geom_point(data=tidal, aes(x=explor_Log10_mean_speed, y=speed_ms_log),
             size=1, alpha= 0.5) +
  scale_colour_manual(values=c("black", "orange")) +
  #scale_colour_manual(values=c("#999999","#0072B2")) +
  #facet_wrap(~tod_start) +
  ylab("Speed in patch (log10 m/s)")+ 
  geom_smooth(method="lm", se=TRUE, data= newdf, 
              aes(x=explor_Log10_mean_speed, y=speed_ms_log),
              size=1) +
  theme_base() + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(-0.8,-0.6)) +
  xlab("Exploration speed (log10 cm/s)") +
  theme(legend.position= "none") 




## Do exploratory birds stay less time in a patch? ####
hist(log10(df3$duration_min))
df3$duration_min_log <- log10(df3$duration_min)

m_duration <- lmer(data=df3, 
                   duration_min_log ~ 
                     explor_Log10_mean_speed * tod_start +
                     night_day_ratio +
                     scale_space +  inout +
                     (1|tag) + (1|tideID) + (1|year))
plot(m_duration)
tab_model(m_duration)
plot_model(m_duration, type = "pred", terms= c("explor_Log10_mean_speed", "tod_start")) 

emmeans(m_duration,  ~ explor_Log10_mean_speed * tod_start, interaction = c("pairwise","pairwise"))

e1.mm2 <- predictorEffect("explor_Log10_mean_speed", m_duration)
plot(e1.mm2, lines=list(multiline=TRUE, col=c("orange","black")), 
     confint=list(style="auto"),
     axes=list(y=list(type="response", 
                      lab="Duration in patch (log10 min)")),
     xlab="Exploration speed (log10 cm/s)",
     type="l", main="")

newdf <- expand.grid(explor_Log10_mean_speed=seq(-0.15,1.25,by=0.2),
                     tod_start = levels(factor(df3$tod_start)),
                     scale_space= -0.02371, #mean
                     tag=levels(factor(df3$tag)),
                     tideID=levels(factor(df3$tideID)),
                     inout = levels(factor(df3$inout)),
                     night_day_ratio = mean(df3$night_day_ratio),
                     year=levels(factor(df3$year)))

newdf$duration_min_log <- predict(m_duration,newdata=newdf,type="response")

tidal <- df3 %>% 
  group_by(tag, tod_start) %>% 
  summarise(se = sd(duration_min_log)/sqrt(length(duration_min_log)),
            duration_min_log = mean(duration_min_log),
            expl_round = format(round(explor_Log10_mean_speed, 2), nsmall = 2),
            explor_Log10_mean_speed=unique(expl_round)) %>% 
  mutate(explor_Log10_mean_speed = as.numeric(explor_Log10_mean_speed))

tidal <-tidal[!duplicated(tidal), ]
summary(df3$duration_min_log)

ggplot(df3, aes(x=explor_Log10_mean_speed, y=duration_min_log, color=as.factor(tod_start))) +
  geom_errorbar(aes(ymin=duration_min_log-se, ymax=duration_min_log+se),
                data=tidal, size=0.1, width=0, alpha=0.5) +
  geom_point(data=tidal, aes(x=explor_Log10_mean_speed, y=duration_min_log),
             size=1, alpha= 0.5) +
  scale_colour_manual(values=c("orange","black")) +
  ylab("Duration in patch (log10 min)")+ 
  geom_smooth(method="lm", data= newdf, aes(x=explor_Log10_mean_speed,y=duration_min_log),
              size=1) +
  theme_base() + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + #, limits = c(0.8, 1.8)) +
  xlab("Exploration speed (log10 cm/s)") +
  theme(legend.position= "none")

tidal2 <- df3

tidal2 <- tidal2 %>% 
  mutate(expl_round = case_when(
    explor_Log10_mean_speed < -0.15 ~ -0.15,
    explor_Log10_mean_speed >= -0.15 & explor_Log10_mean_speed < -0.175 ~ -0.175,
    explor_Log10_mean_speed >= -0.175 & explor_Log10_mean_speed < -0.1 ~ -0.1,
    explor_Log10_mean_speed >= -0.1 & explor_Log10_mean_speed < 0 ~ 0,
    explor_Log10_mean_speed >= 0 & explor_Log10_mean_speed < 0.025 ~ 0.025,
    explor_Log10_mean_speed >= 0.025 & explor_Log10_mean_speed < 0.05 ~ 0.05,
    explor_Log10_mean_speed >= 0.05 & explor_Log10_mean_speed < 0.075 ~ 0.075,
    explor_Log10_mean_speed >= 0.075 & explor_Log10_mean_speed < 0.1 ~ 0.1,
    explor_Log10_mean_speed >= 0.1 & explor_Log10_mean_speed < 0.125 ~ 0.125,
    explor_Log10_mean_speed >= 0.125 & explor_Log10_mean_speed < 0.15 ~ 0.15,
    explor_Log10_mean_speed >= 0.15 & explor_Log10_mean_speed < 0.175 ~ 0.175,
    explor_Log10_mean_speed >= 0.175 & explor_Log10_mean_speed < 0.2 ~ 0.2,
    explor_Log10_mean_speed >= 0.2 & explor_Log10_mean_speed < 0.225 ~ 0.225,
    explor_Log10_mean_speed >= 0.225 & explor_Log10_mean_speed < 0.25 ~ 0.25,
    explor_Log10_mean_speed >= 0.25 & explor_Log10_mean_speed < 0.275 ~ 0.275,
    explor_Log10_mean_speed >= 0.275 & explor_Log10_mean_speed < 0.3 ~ 0.3,
    explor_Log10_mean_speed >= 0.3 & explor_Log10_mean_speed < 0.325 ~ 0.325,
    explor_Log10_mean_speed >= 0.325 & explor_Log10_mean_speed < 0.35 ~ 0.35,
    explor_Log10_mean_speed >= 0.35 & explor_Log10_mean_speed < 0.375 ~ 0.375,
    explor_Log10_mean_speed >= 0.375 & explor_Log10_mean_speed < 0.4 ~ 0.4,
    explor_Log10_mean_speed >= 0.4 & explor_Log10_mean_speed < 0.425 ~ 0.425,
    explor_Log10_mean_speed >= 0.425 & explor_Log10_mean_speed < 0.45 ~ 0.45,
    explor_Log10_mean_speed >= 0.45 & explor_Log10_mean_speed < 0.475 ~ 0.475,
    explor_Log10_mean_speed >= 0.475 & explor_Log10_mean_speed < 0.5 ~ 0.5,
    explor_Log10_mean_speed >= 0.5 & explor_Log10_mean_speed < 0.525 ~ 0.525,
    explor_Log10_mean_speed >= 0.525 & explor_Log10_mean_speed < 0.55 ~ 0.55,
    explor_Log10_mean_speed >= 0.55 & explor_Log10_mean_speed < 0.575 ~ 0.575,
    explor_Log10_mean_speed >= 0.575 & explor_Log10_mean_speed < 0.6 ~ 0.6,
    explor_Log10_mean_speed >= 0.6 & explor_Log10_mean_speed < 0.625 ~ 0.625,
    explor_Log10_mean_speed >= 0.625 & explor_Log10_mean_speed < 0.65 ~ 0.65,
    explor_Log10_mean_speed >= 0.65 & explor_Log10_mean_speed < 0.675 ~ 0.675,
    explor_Log10_mean_speed >= 0.675 & explor_Log10_mean_speed < 0.7 ~ 0.7,
    explor_Log10_mean_speed >= 0.7 & explor_Log10_mean_speed < 0.725 ~ 0.725,
    explor_Log10_mean_speed >= 0.725 & explor_Log10_mean_speed < 0.75 ~ 0.75,
    explor_Log10_mean_speed >= 0.75 & explor_Log10_mean_speed < 0.775 ~ 0.775,
    explor_Log10_mean_speed >= 0.775 & explor_Log10_mean_speed < 0.8 ~ 0.8,
    explor_Log10_mean_speed >= 0.8 & explor_Log10_mean_speed < 0.825 ~ 0.825,
    explor_Log10_mean_speed >= 0.825 & explor_Log10_mean_speed < 0.85 ~ 0.85,
    explor_Log10_mean_speed >= 0.85 & explor_Log10_mean_speed < 0.875 ~ 0.875,
    explor_Log10_mean_speed >= 0.875 & explor_Log10_mean_speed < 0.9 ~ 0.9,
    explor_Log10_mean_speed >= 0.9 & explor_Log10_mean_speed < 0.925 ~ 0.925,
    explor_Log10_mean_speed >= 0.925 & explor_Log10_mean_speed < 0.95 ~ 0.95,
    explor_Log10_mean_speed >= 0.95 & explor_Log10_mean_speed < 0.975 ~ 0.975,
    explor_Log10_mean_speed >= 0.975 & explor_Log10_mean_speed < 1 ~ 1,
    explor_Log10_mean_speed >= 1 & explor_Log10_mean_speed < 1.025 ~ 1.025,
    explor_Log10_mean_speed >= 1.025 & explor_Log10_mean_speed < 1.05 ~ 1.05,
    explor_Log10_mean_speed >= 1.05 & explor_Log10_mean_speed < 1.075 ~ 1.075,
    explor_Log10_mean_speed >= 1.075 & explor_Log10_mean_speed < 1.1 ~ 1.1,
    explor_Log10_mean_speed >= 1.1 & explor_Log10_mean_speed < 1.125 ~ 1.125,
    explor_Log10_mean_speed >= 1.125 & explor_Log10_mean_speed < 1.15 ~ 1.15,
    explor_Log10_mean_speed >= 1.15 & explor_Log10_mean_speed < 1.175 ~ 1.175,
    explor_Log10_mean_speed >= 1.175 & explor_Log10_mean_speed < 1.2 ~ 1.2,
    explor_Log10_mean_speed <= 1.25 ~ 1.25
  ))

tidal <- tidal2 %>% 
  group_by(expl_round, tod_start) %>% 
  summarise(se = sd(duration_min_log)/sqrt(length(duration_min_log)),
            duration_min_log = mean(duration_min_log)) #%>% 
#expl_round = format(round(explor_Log10_mean_speed, 2), nsmall = 2),
#explor_Log10_mean_speed=unique(expl_round)) %>% 
#mutate(explor_Log10_mean_speed = as.numeric(explor_Log10_mean_speed)) 

#tidal <-tidal[!duplicated(tidal), ]

colnames(tidal)[1] <- "explor_Log10_mean_speed"
tidal$explor_Log10_mean_speed <- as.numeric(tidal$explor_Log10_mean_speed)


## Do exploratory birds move more in a patch? ####
summary(df3$dist_in_patch_m)
df3$dist_in_patch_m_log <- log10(df3$dist_in_patch_m)
hist(df3$dist_in_patch_m_log)

m_distance <- lmer( data=df3,
                    dist_in_patch_m_log ~
                      explor_Log10_mean_speed  * tod_start + 
                      #scale(nfixes) + 
                      tide_cont + 
                      scale_space + 
                      inout + 
                      night_day_ratio +
                      (1|tag) + (1|tideID) + (1|year))
plot(m_distance)
tab_model(m_distance)
plot_model(m_distance, type = "pred", terms= c("tod_start","explor_Log10_mean_speed [0.1,0.9]")) 

emmeans(m_distance,  ~ explor_Log10_mean_speed * tod_start, interaction = c("pairwise","pairwise"))

e1.mm3 <- predictorEffect("explor_Log10_mean_speed", m_distance)
plot(e1.mm3, lines=list(multiline=TRUE, col=c("orange","black")), 
     confint=list(style="auto"),
     axes=list(y=list(type="response", 
                      lab="Distance in patch (log10 m)")),
     xlab="Exploration speed (log10 cm/s)",
     type="l", main="")

newdf <- expand.grid(explor_Log10_mean_speed=seq(-0.15,1.25,by=0.2),
                     tod_start = levels(factor(df3$tod_start)),
                     scale_space= -0.02371, #mean
                     tag=levels(factor(df3$tag)),
                     tideID=levels(factor(df3$tideID)),
                     inout = levels(factor(df3$inout)),
                     night_day_ratio = mean(df3$night_day_ratio),
                     year=levels(factor(df3$year)))

newdf$dist_in_patch_m_log <- predict(m_distance,newdata=newdf,type="response")

tidal <- df3 %>% 
  group_by(tag, tod_start) %>% 
  summarise(se = sd(dist_in_patch_m_log)/sqrt(length(dist_in_patch_m_log)),
            dist_in_patch_m_log = mean(dist_in_patch_m_log),
            expl_round = format(round(explor_Log10_mean_speed, 2), nsmall = 2),
            explor_Log10_mean_speed=unique(expl_round)) %>% 
  mutate(explor_Log10_mean_speed = as.numeric(explor_Log10_mean_speed))

tidal <-tidal[!duplicated(tidal), ]

ggplot(df3, aes(x=explor_Log10_mean_speed, y=dist_in_patch_m_log, color=as.factor(tod_start))) +
  geom_errorbar(aes(ymin=dist_in_patch_m_log-se, ymax=dist_in_patch_m_log+se),
                data=tidal, size=0.1, width=0, alpha=0.5) +
  geom_point(data=tidal, aes(x=explor_Log10_mean_speed, y=dist_in_patch_m_log),
             size=1, alpha= 0.5) +
  scale_colour_manual(values=c("orange","black")) +
  ylab("Distance in patch (log10 m)")+ 
  geom_smooth(method="lm", data= newdf, aes(x=explor_Log10_mean_speed,y=dist_in_patch_m_log),
              size=1) +
  theme_base() + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + #, limits = c(1.9,2.8)) +
  xlab("Exploration speed (log10 cm/s)") +
  theme(legend.position= "none")

tidal2 <- df3

tidal2 <- tidal2 %>% 
  mutate(expl_round = case_when(
    explor_Log10_mean_speed < -0.15 ~ -0.15,
    explor_Log10_mean_speed >= -0.15 & explor_Log10_mean_speed < -0.175 ~ -0.175,
    explor_Log10_mean_speed >= -0.175 & explor_Log10_mean_speed < -0.1 ~ -0.1,
    explor_Log10_mean_speed >= -0.1 & explor_Log10_mean_speed < 0 ~ 0,
    explor_Log10_mean_speed >= 0 & explor_Log10_mean_speed < 0.025 ~ 0.025,
    explor_Log10_mean_speed >= 0.025 & explor_Log10_mean_speed < 0.05 ~ 0.05,
    explor_Log10_mean_speed >= 0.05 & explor_Log10_mean_speed < 0.075 ~ 0.075,
    explor_Log10_mean_speed >= 0.075 & explor_Log10_mean_speed < 0.1 ~ 0.1,
    explor_Log10_mean_speed >= 0.1 & explor_Log10_mean_speed < 0.125 ~ 0.125,
    explor_Log10_mean_speed >= 0.125 & explor_Log10_mean_speed < 0.15 ~ 0.15,
    explor_Log10_mean_speed >= 0.15 & explor_Log10_mean_speed < 0.175 ~ 0.175,
    explor_Log10_mean_speed >= 0.175 & explor_Log10_mean_speed < 0.2 ~ 0.2,
    explor_Log10_mean_speed >= 0.2 & explor_Log10_mean_speed < 0.225 ~ 0.225,
    explor_Log10_mean_speed >= 0.225 & explor_Log10_mean_speed < 0.25 ~ 0.25,
    explor_Log10_mean_speed >= 0.25 & explor_Log10_mean_speed < 0.275 ~ 0.275,
    explor_Log10_mean_speed >= 0.275 & explor_Log10_mean_speed < 0.3 ~ 0.3,
    explor_Log10_mean_speed >= 0.3 & explor_Log10_mean_speed < 0.325 ~ 0.325,
    explor_Log10_mean_speed >= 0.325 & explor_Log10_mean_speed < 0.35 ~ 0.35,
    explor_Log10_mean_speed >= 0.35 & explor_Log10_mean_speed < 0.375 ~ 0.375,
    explor_Log10_mean_speed >= 0.375 & explor_Log10_mean_speed < 0.4 ~ 0.4,
    explor_Log10_mean_speed >= 0.4 & explor_Log10_mean_speed < 0.425 ~ 0.425,
    explor_Log10_mean_speed >= 0.425 & explor_Log10_mean_speed < 0.45 ~ 0.45,
    explor_Log10_mean_speed >= 0.45 & explor_Log10_mean_speed < 0.475 ~ 0.475,
    explor_Log10_mean_speed >= 0.475 & explor_Log10_mean_speed < 0.5 ~ 0.5,
    explor_Log10_mean_speed >= 0.5 & explor_Log10_mean_speed < 0.525 ~ 0.525,
    explor_Log10_mean_speed >= 0.525 & explor_Log10_mean_speed < 0.55 ~ 0.55,
    explor_Log10_mean_speed >= 0.55 & explor_Log10_mean_speed < 0.575 ~ 0.575,
    explor_Log10_mean_speed >= 0.575 & explor_Log10_mean_speed < 0.6 ~ 0.6,
    explor_Log10_mean_speed >= 0.6 & explor_Log10_mean_speed < 0.625 ~ 0.625,
    explor_Log10_mean_speed >= 0.625 & explor_Log10_mean_speed < 0.65 ~ 0.65,
    explor_Log10_mean_speed >= 0.65 & explor_Log10_mean_speed < 0.675 ~ 0.675,
    explor_Log10_mean_speed >= 0.675 & explor_Log10_mean_speed < 0.7 ~ 0.7,
    explor_Log10_mean_speed >= 0.7 & explor_Log10_mean_speed < 0.725 ~ 0.725,
    explor_Log10_mean_speed >= 0.725 & explor_Log10_mean_speed < 0.75 ~ 0.75,
    explor_Log10_mean_speed >= 0.75 & explor_Log10_mean_speed < 0.775 ~ 0.775,
    explor_Log10_mean_speed >= 0.775 & explor_Log10_mean_speed < 0.8 ~ 0.8,
    explor_Log10_mean_speed >= 0.8 & explor_Log10_mean_speed < 0.825 ~ 0.825,
    explor_Log10_mean_speed >= 0.825 & explor_Log10_mean_speed < 0.85 ~ 0.85,
    explor_Log10_mean_speed >= 0.85 & explor_Log10_mean_speed < 0.875 ~ 0.875,
    explor_Log10_mean_speed >= 0.875 & explor_Log10_mean_speed < 0.9 ~ 0.9,
    explor_Log10_mean_speed >= 0.9 & explor_Log10_mean_speed < 0.925 ~ 0.925,
    explor_Log10_mean_speed >= 0.925 & explor_Log10_mean_speed < 0.95 ~ 0.95,
    explor_Log10_mean_speed >= 0.95 & explor_Log10_mean_speed < 0.975 ~ 0.975,
    explor_Log10_mean_speed >= 0.975 & explor_Log10_mean_speed < 1 ~ 1,
    explor_Log10_mean_speed >= 1 & explor_Log10_mean_speed < 1.025 ~ 1.025,
    explor_Log10_mean_speed >= 1.025 & explor_Log10_mean_speed < 1.05 ~ 1.05,
    explor_Log10_mean_speed >= 1.05 & explor_Log10_mean_speed < 1.075 ~ 1.075,
    explor_Log10_mean_speed >= 1.075 & explor_Log10_mean_speed < 1.1 ~ 1.1,
    explor_Log10_mean_speed >= 1.1 & explor_Log10_mean_speed < 1.125 ~ 1.125,
    explor_Log10_mean_speed >= 1.125 & explor_Log10_mean_speed < 1.15 ~ 1.15,
    explor_Log10_mean_speed >= 1.15 & explor_Log10_mean_speed < 1.175 ~ 1.175,
    explor_Log10_mean_speed >= 1.175 & explor_Log10_mean_speed < 1.2 ~ 1.2,
    explor_Log10_mean_speed <= 1.25 ~ 1.25
  ))

tidal <- tidal2 %>% 
  group_by(expl_round, tod_start) %>% 
  summarise(se = sd(dist_in_patch_m_log)/sqrt(length(dist_in_patch_m_log)),
            dist_in_patch_m_log = mean(dist_in_patch_m_log)) #%>% 

colnames(tidal)[1] <- "explor_Log10_mean_speed"
tidal$explor_Log10_mean_speed <- as.numeric(tidal$explor_Log10_mean_speed)



## Do exploratory birds move fast within patches? ####
#df2$speed_in_patch_ms <- df2$dist_in_patch_cm / df2$duration_s
#df2 <- df2 %>% filter(speed_in_patch_ms >0.4)

df3$speed_patch_log <- log10(df3$speed_in_patch_ms)
df3_finite <- df3 %>%  filter_all(all_vars(!is.infinite(.)))

#hist(df3_finite$speed_patch_log)
#hist(df2$speed_in_patch_ms)

m_speed <- lmer(speed_patch_log ~ 
                  explor_Log10_mean_speed * tod_start + 
                  explor_Log10_mean_speed * inout +
                  explor_Log10_mean_speed * scale_space +
                  (1|tag) + (1|tideID) + (1|year),
                data= df3_finite)
plot(m_speed)
tab_model(m_speed)
plot_model(m_speed, type = "pred", terms= c("explor_Log10_mean_speed", "inout", "tod_start")) 
plot_model(m_speed, type = "pred", terms= c("explor_Log10_mean_speed", "scale_space")) 
plot_model(m_speed, type = "pred", terms= c("explor_Log10_mean_speed", "inout")) 

newdf <- expand.grid(explor_Log10_mean_speed=seq(-0.15,1.25,by=0.2),
                     tod_start = levels(factor(df2$tod_start)),
                     inout = levels(factor(df2$inout)),
                     scale_space=mean(df2$scale_space),
                     tag=levels(factor(df2$tag)),
                     tideID=levels(factor(df2$tideID)),
                     year=levels(factor(df2$year)))

newdf$speed_patch_log <- predict(m_speed,newdata=newdf,type="response")

tidal <- df2 %>% 
  group_by(tag, inout) %>% 
  summarise(se = sd(speed_patch_log)/sqrt(length(speed_patch_log)),
            speed_patch_log = mean(speed_patch_log),
            explor_Log10_mean_speed=unique(expl_round)) %>% 
  mutate(explor_Log10_mean_speed = as.numeric(explor_Log10_mean_speed))

ggplot(df2, aes(x=explor_Log10_mean_speed, y=speed_patch_log, color=as.factor(inout))) +
  geom_errorbar(aes(ymin=speed_patch_log-se, ymax=speed_patch_log+se),
                data=tidal, size=0.1, width=0, alpha=0.5) +
  geom_point(data=tidal, aes(x=explor_Log10_mean_speed, y=speed_patch_log),
             size=1, alpha= 0.5) +
  scale_colour_manual(values=c("orange","black")) +
  ylab("Speed in patch (m/s)")+ 
  geom_smooth(method="lm", data= newdf, aes(x=explor_Log10_mean_speed,y=speed_patch_log, 
                                            lintype=as.factor(inout)),
              size=1) +
  theme_base() + 
  scale_linetype_manual(values=c("dashed","solid")) +
  facet_wrap(~tod_start) + 
  xlab("Exploration speed (log10 cm/s)") +
  theme(legend.position= "none") + 
  ylim(1.16,1.66)


#boxplot
newdf <- newdf %>% mutate(expl= case_when(explor_Log10_mean_speed < 0 ~ "Slower-explorers",
                                          explor_Log10_mean_speed > 0.9 ~ "Faster-explorers"))
newdf1 <- na.omit(newdf)

ggplot(data=newdf1, aes(y= speed_patch_unlog, x=inout, col=expl)) + 
  geom_boxplot(notch =T, outlier.alpha=0.006) +
  theme_classic() + facet_wrap(~tod_start) +
  scale_color_manual(values = c("seagreen4", "orange2"),
                     name="") +
  ylab("Speed in a patch (min)") + 
  theme(axis.title.y = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.15,0.92),
        panel.grid = element_blank(),
        panel.border =element_rect(colour = "black", fill=NA, size=0.5),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=12),
        strip.background =element_blank()) +
  scale_x_discrete(labels= labels) 

## Do exploratory birds displace more in a patch? ####
hist(log10(1+df3$disp_in_patch_m))
df3$disp_in_patch_m_log <- log10(1+df3$disp_in_patch_m)
df3_finite <- df3 %>%  filter_all(all_vars(!is.infinite(.)))
hist(df3_finite$disp_in_patch_m_log)

df3_finite <- df3_finite %>% filter(disp_in_patch_m_log > -0.5)

m_disp <- lmer( data=df3_finite,
                disp_in_patch_m_log ~
                      explor_Log10_mean_speed  * tod_start +
                      scale_space +  inout + night_day_ratio +
                      (1|tag) + (1|tideID) + (1|year))
tab_model(m_disp)
plot(m_disp)
plot_model(m_disp, type = "pred", terms= c("explor_Log10_mean_speed", "tod_start")) 

emmeans(m_disp,  ~ explor_Log10_mean_speed * tod_start, interaction = c("pairwise","pairwise"))

e1.mm4 <- predictorEffect("explor_Log10_mean_speed", m_disp)
plot(e1.mm4, lines=list(multiline=TRUE, col=c("orange","black")), 
     confint=list(style="auto"),
     axes=list(y=list(type="response", 
                      lab="Displacement in patch (log10 m)")),
     xlab="Exploration speed (log10 cm/s)",
     type="l", main="")


newdf <- expand.grid(explor_Log10_mean_speed=seq(-0.15,1.25,by=0.2),
                     tod_start = levels(factor(df3_finite$tod_start)),
                     scale_space= -0.02371, #mean
                     tag=levels(factor(df3_finite$tag)),
                     tideID=levels(factor(df3_finite$tideID)),
                     inout = levels(factor(df3_finite$inout)),
                     night_day_ratio = mean(df3_finite$night_day_ratio),
                     year=levels(factor(df3_finite$year)))

newdf$disp_in_patch_m_log <- predict(m_disp,newdata=newdf,type="response")

tidal <- df3_finite %>% 
  group_by(tag, tod_start) %>% 
  summarise(se = sd(disp_in_patch_m_log)/sqrt(length(disp_in_patch_m_log)),
            disp_in_patch_m_log = mean(disp_in_patch_m_log),
            expl_round = format(round(explor_Log10_mean_speed, 2), nsmall = 2),
            explor_Log10_mean_speed=unique(expl_round)) %>% 
  mutate(explor_Log10_mean_speed = as.numeric(explor_Log10_mean_speed))

tidal <-tidal[!duplicated(tidal), ]

ggplot(df3_finite, aes(x=explor_Log10_mean_speed, y=disp_in_patch_m_log, color=as.factor(tod_start))) +
  geom_errorbar(aes(ymin=disp_in_patch_m_log-se, ymax=disp_in_patch_m_log+se),
                data=tidal, size=0.1, width=0, alpha=0.5) +
  geom_point(data=tidal, aes(x=explor_Log10_mean_speed, y=disp_in_patch_m_log),
             size=1, alpha= 0.5) +
  scale_colour_manual(values=c("orange","black")) +
  ylab("Displacement in patch (log10 m)")+ 
  geom_smooth(method="lm", data= newdf, aes(x=explor_Log10_mean_speed,y=disp_in_patch_m_log),
              size=1) +
  theme_base() + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(1.6,2.1)) +
  xlab("Exploration speed (log10 cm/s)") +
  theme(legend.position= "none")

tidal2 <- df3_finite

tidal2 <- tidal2 %>% 
  mutate(expl_round = case_when(
    explor_Log10_mean_speed < -0.15 ~ -0.15,
    explor_Log10_mean_speed >= -0.15 & explor_Log10_mean_speed < -0.175 ~ -0.175,
    explor_Log10_mean_speed >= -0.175 & explor_Log10_mean_speed < -0.1 ~ -0.1,
    explor_Log10_mean_speed >= -0.1 & explor_Log10_mean_speed < 0 ~ 0,
    explor_Log10_mean_speed >= 0 & explor_Log10_mean_speed < 0.025 ~ 0.025,
    explor_Log10_mean_speed >= 0.025 & explor_Log10_mean_speed < 0.05 ~ 0.05,
    explor_Log10_mean_speed >= 0.05 & explor_Log10_mean_speed < 0.075 ~ 0.075,
    explor_Log10_mean_speed >= 0.075 & explor_Log10_mean_speed < 0.1 ~ 0.1,
    explor_Log10_mean_speed >= 0.1 & explor_Log10_mean_speed < 0.125 ~ 0.125,
    explor_Log10_mean_speed >= 0.125 & explor_Log10_mean_speed < 0.15 ~ 0.15,
    explor_Log10_mean_speed >= 0.15 & explor_Log10_mean_speed < 0.175 ~ 0.175,
    explor_Log10_mean_speed >= 0.175 & explor_Log10_mean_speed < 0.2 ~ 0.2,
    explor_Log10_mean_speed >= 0.2 & explor_Log10_mean_speed < 0.225 ~ 0.225,
    explor_Log10_mean_speed >= 0.225 & explor_Log10_mean_speed < 0.25 ~ 0.25,
    explor_Log10_mean_speed >= 0.25 & explor_Log10_mean_speed < 0.275 ~ 0.275,
    explor_Log10_mean_speed >= 0.275 & explor_Log10_mean_speed < 0.3 ~ 0.3,
    explor_Log10_mean_speed >= 0.3 & explor_Log10_mean_speed < 0.325 ~ 0.325,
    explor_Log10_mean_speed >= 0.325 & explor_Log10_mean_speed < 0.35 ~ 0.35,
    explor_Log10_mean_speed >= 0.35 & explor_Log10_mean_speed < 0.375 ~ 0.375,
    explor_Log10_mean_speed >= 0.375 & explor_Log10_mean_speed < 0.4 ~ 0.4,
    explor_Log10_mean_speed >= 0.4 & explor_Log10_mean_speed < 0.425 ~ 0.425,
    explor_Log10_mean_speed >= 0.425 & explor_Log10_mean_speed < 0.45 ~ 0.45,
    explor_Log10_mean_speed >= 0.45 & explor_Log10_mean_speed < 0.475 ~ 0.475,
    explor_Log10_mean_speed >= 0.475 & explor_Log10_mean_speed < 0.5 ~ 0.5,
    explor_Log10_mean_speed >= 0.5 & explor_Log10_mean_speed < 0.525 ~ 0.525,
    explor_Log10_mean_speed >= 0.525 & explor_Log10_mean_speed < 0.55 ~ 0.55,
    explor_Log10_mean_speed >= 0.55 & explor_Log10_mean_speed < 0.575 ~ 0.575,
    explor_Log10_mean_speed >= 0.575 & explor_Log10_mean_speed < 0.6 ~ 0.6,
    explor_Log10_mean_speed >= 0.6 & explor_Log10_mean_speed < 0.625 ~ 0.625,
    explor_Log10_mean_speed >= 0.625 & explor_Log10_mean_speed < 0.65 ~ 0.65,
    explor_Log10_mean_speed >= 0.65 & explor_Log10_mean_speed < 0.675 ~ 0.675,
    explor_Log10_mean_speed >= 0.675 & explor_Log10_mean_speed < 0.7 ~ 0.7,
    explor_Log10_mean_speed >= 0.7 & explor_Log10_mean_speed < 0.725 ~ 0.725,
    explor_Log10_mean_speed >= 0.725 & explor_Log10_mean_speed < 0.75 ~ 0.75,
    explor_Log10_mean_speed >= 0.75 & explor_Log10_mean_speed < 0.775 ~ 0.775,
    explor_Log10_mean_speed >= 0.775 & explor_Log10_mean_speed < 0.8 ~ 0.8,
    explor_Log10_mean_speed >= 0.8 & explor_Log10_mean_speed < 0.825 ~ 0.825,
    explor_Log10_mean_speed >= 0.825 & explor_Log10_mean_speed < 0.85 ~ 0.85,
    explor_Log10_mean_speed >= 0.85 & explor_Log10_mean_speed < 0.875 ~ 0.875,
    explor_Log10_mean_speed >= 0.875 & explor_Log10_mean_speed < 0.9 ~ 0.9,
    explor_Log10_mean_speed >= 0.9 & explor_Log10_mean_speed < 0.925 ~ 0.925,
    explor_Log10_mean_speed >= 0.925 & explor_Log10_mean_speed < 0.95 ~ 0.95,
    explor_Log10_mean_speed >= 0.95 & explor_Log10_mean_speed < 0.975 ~ 0.975,
    explor_Log10_mean_speed >= 0.975 & explor_Log10_mean_speed < 1 ~ 1,
    explor_Log10_mean_speed >= 1 & explor_Log10_mean_speed < 1.025 ~ 1.025,
    explor_Log10_mean_speed >= 1.025 & explor_Log10_mean_speed < 1.05 ~ 1.05,
    explor_Log10_mean_speed >= 1.05 & explor_Log10_mean_speed < 1.075 ~ 1.075,
    explor_Log10_mean_speed >= 1.075 & explor_Log10_mean_speed < 1.1 ~ 1.1,
    explor_Log10_mean_speed >= 1.1 & explor_Log10_mean_speed < 1.125 ~ 1.125,
    explor_Log10_mean_speed >= 1.125 & explor_Log10_mean_speed < 1.15 ~ 1.15,
    explor_Log10_mean_speed >= 1.15 & explor_Log10_mean_speed < 1.175 ~ 1.175,
    explor_Log10_mean_speed >= 1.175 & explor_Log10_mean_speed < 1.2 ~ 1.2,
    explor_Log10_mean_speed <= 1.25 ~ 1.25
  ))

tidal <- tidal2 %>% 
  group_by(expl_round, tod_start) %>% 
  summarise(se = sd(disp_in_patch_m_log)/sqrt(length(disp_in_patch_m_log)),
            disp_in_patch_m_log = mean(disp_in_patch_m_log)) #%>% 
#expl_round = format(round(explor_Log10_mean_speed, 2), nsmall = 2),
#explor_Log10_mean_speed=unique(expl_round)) %>% 
#mutate(explor_Log10_mean_speed = as.numeric(explor_Log10_mean_speed)) 

#tidal <-tidal[!duplicated(tidal), ]

colnames(tidal)[1] <- "explor_Log10_mean_speed"
tidal$explor_Log10_mean_speed <- as.numeric(tidal$explor_Log10_mean_speed)

##### MODEL 2 ########## ------------------------------------------------------------------------------------------------

## Do exploratory birds move more between patches? ####
df3$dist_bw_patch2_log <- log10(df3$dist_bw_patch)
hist(df3$dist_bw_patch)

brm_3 <- bf(dist_bw_patch2_log ~
              explor_Log10_mean_speed  * tod_start  + 
              inout + 
              night_day_ratio +
              scale_space +
              (1|tag) + (1|tideID) + (1|year))

# brm3 <- brm(brm_3, data = df3, 
#             family = "gaussian",
#             control = list(adapt_delta=0.99),
#             warmup = 500,
#             iter = 1000, thin=2,
#             chains = 1, #inits = "random",
#             seed = 12345,
#             cores = my.cores) 

#save(brm3, file= "C:/Users/sersoy/OneDrive - NIOZ/PhD/PhD_chapters/Chapter 3- Personality and patch movement/Data and R/Final_2022-06-01/brm_DISPBTWPATCH_2022-07-07-1chain.Rdata")
load(file= "brm_DISPBTWPATCH_2022-07-07-1chain.Rdata")
summary(brm3)

p4 <- pp_check(brm3, resp = "dist_bw_patch2_log")

library(emmeans)
library(rstanarm)
emmeans(brm3, pairwise ~ explor_Log10_mean_speed * tod_start )

pp_check(brm3, resp = "dist_bw_patch2_log")
mcmc_plot(brm3)
colnames(posterior_samples(brm3))[1:18]

conditional_effects(brm3, 
                    effects = "explor_Log10_mean_speed:tod_start", 
                    resp = "dist_bw_patch2_log",
                    spaghetti = T, nsamples = 150) %>% 
  plot(points = T,
       point_args = c(alpha = 1/3, size = 1), 
       mean = F, 
       rug = T, 
       theme= theme_classic())


me <- conditional_effects(brm3, effects = "explor_Log10_mean_speed:tod_start", spaghetti = T, nsamples = 150)

plot(me, plot = FALSE, points = T, mean= F, point_args = c(alpha = 1/3, size = 1))[[1]] +
  scale_color_grey() +
  #scale_color_manual(values = c("orange", "black")) +
  #scale_fill_manual(values = c("orange", "black")) +
  scale_fill_grey() +
  theme_classic()
  

mex <- as.data.frame(me$`explor_Log10_mean_speed`)  


conditional_effects(brm3, "explor_Log10_mean_speed:tod_start")
plot(me, plot = FALSE)[[1]] +
  scale_color_grey() +
  scale_fill_grey()





## Repeatability ####
df3$speed_patch <- df3$dist_in_patch_m/df3$duration_s

df_tide <- df3 %>% 
  group_by(tag,tideID) %>% 
  dplyr::summarise(speed_patch = mean(speed_patch),
            #duration_min = mean(duration_min),
            #dist_in_patch_m = mean(dist_in_patch_m),
            waterlevel_mean_cm = mean(waterlevel_mean_cm),
            #start_dist_to_tideline2 = mean(start_dist_to_tideline2),
            #mean_tod_patch = mean(mean_tod_patch)
            ) %>% 
  mutate(year = factor(substr(tideID, 1, 4)),
         tide_cont = as.numeric(substr(tideID, 5, 7)))

brm_4 <- bf(log10(speed_patch) ~
              waterlevel_mean_cm + tide_cont +
              (1|tag))

# brm4 <- brm(brm_4, data = df_tide,
#             family = "gaussian",
#             control = list(adapt_delta=0.99),
#             warmup = 500,
#             iter = 1000, thin=2,
#             chains = 1, #inits = "random",
#             seed = 12345,
#             cores = my.cores)

summary(brm4)

#save(brm4, file= "brm_REP_2022-07-07-1chain.Rdata")
load(file= "brm_REP_2022-07-07-1chain.Rdata")

#repeatability 

colnames(posterior_samples(brm4))[1:20]

var.animal.dist <- posterior_samples(brm)[,4]^2
var.res.dist <- posterior_samples(brm)[,5]^2

var.tot.dist <- var.animal.dist + var.res.dist

RDist <- var.animal.dist / var.tot.dist
mean(RDist);HPDinterval(as.mcmc(RDist),0.95)


tab_model(m)
print(VarCorr(m),comp = c("Variance","Std.Dev."))

set.seed(1)
simulated <- arm::sim(m, 1000)
describe_posterior(simulated)

posterior_tag <- apply(simulated@ranef$"tag"[ , , 1],1,var)
#posterior_tideID <- apply(simulated@ranef$"tideID"[ , , 1],1,var)
posterior_residual <- simulated@sigma^2

#among-indv
median(posterior_tag)
quantile(posterior_tag, prob=c(0.025, 0.5, 0.975))
#repeatability
quantile(posterior_tag /
           (posterior_tag +  posterior_residual),
         prob=c(0.025, 0.5, 0.975))

#within-indv
median(posterior_residual)
quantile(posterior_residual,prob=c(0.025, 0.5, 0.975))
#residual repeatability
quantile(posterior_residual /
           (posterior_tag + posterior_residual),
         prob=c(0.025, 0.5, 0.975))

CVi <- sqrt(posterior_tag) / summary(m)$coefficients[1]
quantile(CVi,prob=c(0.025, 0.5, 0.975))


# R comparison between day-night ####
#### Fit brms for group comparison
brm.formula=bf(speed_patch ~ tod_start + (0+tod_start||tag), sigma ~ 0+tod_start)

#(0+AGE||FB) estimate separate variance components by age for the among-individual variance 
#sigma ~ 0+AGE does the same for the within-individual variance. 

#run brm
fit_brm <- brm(brm.formula, data = df3,
               warmup = 500, iter = 1000, thin=4, chains=1,
               cores = 3)
summary(fit_brm)
#save(fit_brm, file= "C:/Users/sersoy/OneDrive - NIOZ/PhD/PhD_chapters/Chapter 3- Personality and patch movement/Data and R/Final_2022-06-01/brm_GroupREP_2022-07-12-1chain.Rdata")
colnames(posterior_samples(fit_brm))[1:8]

var.brms <- posterior_samples(fit_brm, 
                              c("b_Intercept", #Intercept
                                "b_tod_startnight", # Day-night coeff
                                "sd_tag__tod_startday", # among individual sd day
                                "sd_tag__tod_startnight", # among individual sd night
                                "b_sigma_tod_startday", # within individual sd day, log scale
                                "b_sigma_tod_startnight")) # within individual sd night, log scale

colnames(var.brms)=c("Intc","beta_Juv","Vi.Adult", "Vi.Juv", "Vw.Adult", "Vw.Juv")
var.brms=as.data.frame(var.brms)
var.brms %>% head() %>% knitr::kable(., digits = 2)

var.brms$beta.Adult=var.brms$Intc
var.brms$beta.Juv=var.brms$Intc+var.brms$beta_Juv

var.brms$Vi.Adult=(var.brms$Vi.Adult)^2 #need to square these values to express theme as variances
var.brms$Vi.Juv=(var.brms$Vi.Juv)^2 #need to square these values to express theme as variances
var.brms$Vw.Adult=exp(var.brms$Vw.Adult)^2 # and take the exponent for the residual standard deviation
var.brms$Vw.Juv=exp(var.brms$Vw.Juv)^2 # and take the exponent for the residual standard deviation

var.brms$tau.Adult = var.brms$Vi.Adult/(var.brms$Vi.Adult + var.brms$Vw.Adult)
var.brms$tau.Juv = var.brms$Vi.Juv/(var.brms$Vi.Juv + var.brms$Vw.Juv)

var.brms$delta.Vi = var.brms$Vi.Adult - var.brms$Vi.Juv
var.brms$delta.Vw = var.brms$Vw.Adult - var.brms$Vw.Juv
var.brms$delta.tau = var.brms$tau.Adult - var.brms$tau.Juv

#make table with variance and ratio for each group
t.1 = var.brms %>%
  dplyr::select(Vi.Adult, Vi.Juv, Vw.Adult, Vw.Juv, tau.Adult, tau.Juv) %>%
  stack() %>%
  group_by(ind) %>%
  summarise_if(is.numeric, describe_posterior)

t.1 = t.1$values
t.1$level = c("Vi.Adult", "Vi.Juv", "Vw.Adult", "Vw.Juv", "tau.Adult", "tau.Juv")
t.1 = t.1 %>%
  dplyr::select(level, Median, CI_low, CI_high, pd)
colnames(t.1) = c("Variance and ratio", "Median", "lower_CI", "Upper_CI", "Pmcmc")

t.1 %>%
  knitr::kable(digits = 2)

#make table with variance differences
t.2 = var.brms %>%
  dplyr::select(delta.Vi, delta.Vw, delta.tau) %>%
  stack() %>%
  group_by(ind) %>%
  summarise_if(is.numeric, describe_posterior)

t.2 = t.2$values
t.2$level = c("Vi", "Vw", "tau")
t.2 = t.2 %>%
  dplyr::select(level, Median, CI_low, CI_high, pd)
colnames(t.2) = c(expression(Delta), "Median", "lower_CI", "Upper_CI", "Pmcmc")

t.2 %>%
  knitr::kable(digits = 2)

#day-night repeatability do not differ

## Exploratory birds's start position to the tidal line ####
m_tideline <- lmer(data=df3, 
                   start_dist_to_tideline ~ 
                     explor_Log10_mean_speed  * inout * tod_start +
                     scale_space + 
                     (1|tag) + (1|tideID) + (1|year))

#m_tideline_sim <- arm::sim(m_tideline, 1000)
#describe_posterior(m_tideline_sim)
#summary(m_tideline)
#hist(residuals(m_tideline))
tab_model(m_tideline)
plot_model(m_tideline, type = "pred", terms= c("tod_start", "explor_Log10_mean_speed [0.1,0.9]","inout")) 

newdf <- expand.grid(explor_Log10_mean_speed=seq(-0.15,1.25,by=0.2),
                     tod_start = levels(factor(df2$tod_start)),
                     inout = levels(factor(df2$inout)),
                     scale_space=mean(df2$scale_space),
                     tag=levels(factor(df2$tag)),
                     tideID=levels(factor(df2$tideID)),
                     year=levels(factor(df2$year)))

newdf$start_dist_to_tideline2_m <- predict(m_tideline,newdata=newdf,type="response")
#newdf$dist_in_patch_m_log_unlog <- 10^(newdf$dist_in_patch_m_log)
newdf <- newdf %>% mutate(expl= case_when(explor_Log10_mean_speed < 0 ~ "Slower-explorers",
                                          explor_Log10_mean_speed > 0.9 ~ "Faster-explorers"))
newdf1 <- na.omit(newdf)

p4 <- ggplot(data=newdf1, aes(y= start_dist_to_tideline2_m, x=tod_start, col=expl)) + 
  geom_boxplot(notch =T, outlier.alpha=0.006) +
  theme_classic() + 
  facet_wrap(~inout) +
  scale_color_manual(values = c("seagreen4", "orange2"),
                     name="") +
  ylab("Est. distance from start of a patch to tidal line (m)") + 
  theme(axis.title.y = element_text(size=13),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border =element_rect(colour = "black", fill=NA, size=0.5),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=12),
        strip.background =element_blank()) +
  scale_x_discrete(labels= labels)
p4


## Exploratory birds's patch's end position to the tidal line ####
df2$end_dist_to_tideline2_m <- df2$end_dist_to_tideline2/100
m_end_tideline <- lmer(data=df2, 
                   end_dist_to_tideline2_m ~ 
                     explor_Log10_mean_speed  * inout * tod_start +
                     scale_space + 
                     (1|tag) + (1|tideID) + (1|year))

#m_tideline_sim <- arm::sim(m_tideline, 1000)
#describe_posterior(m_tideline_sim)
#summary(m_tideline)
#hist(residuals(m_tideline))
tab_model(m_end_tideline)
plot_model(m_end_tideline, type = "pred", terms= c("tod_start", "explor_Log10_mean_speed [0.1,0.9]","inout")) 

newdf <- expand.grid(explor_Log10_mean_speed=seq(-0.15,1.25,by=0.2),
                     tod_start = levels(factor(df2$tod_start)),
                     inout = levels(factor(df2$inout)),
                     scale_space=mean(df2$scale_space),
                     tag=levels(factor(df2$tag)),
                     tideID=levels(factor(df2$tideID)),
                     year=levels(factor(df2$year)))

newdf$end_dist_to_tideline2_m <- predict(m_end_tideline,newdata=newdf,type="response")
#newdf$dist_in_patch_m_log_unlog <- 10^(newdf$dist_in_patch_m_log)
newdf <- newdf %>% mutate(expl= case_when(explor_Log10_mean_speed < 0 ~ "Slower-explorers",
                                          explor_Log10_mean_speed > 0.9 ~ "Faster-explorers"))
newdf1 <- na.omit(newdf)

p5 <- ggplot(data=newdf1, aes(y= end_dist_to_tideline2_m, x=tod_start, col=expl)) + 
  geom_boxplot(notch =T, outlier.alpha=0.006) +
  theme_classic() + 
  facet_wrap(~inout) +
  scale_color_manual(values = c("seagreen4", "orange2"),
                     name="") +
  ylab("Est. distance from end of a patch to tidal line (m)") + 
  theme(axis.title.y = element_text(size=13),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border =element_rect(colour = "black", fill=NA, size=0.5),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=12),
        strip.background =element_blank()) +
  scale_x_discrete(labels= labels)
p5






####

str(df2)
summary(df2$time_bw_patch2)
m_timebtw <- lmer(data=df2, 
                  time_bw_patch2  ~ 
                         explor_Log10_mean_speed  *  tod_start + inout +
                         scale_space + 
                         (1|tag) + (1|tideID) + (1|year))

tab_model(m_timebtw)
plot_model(m_timebtw, type = "pred", terms= c("explor_Log10_mean_speed", "tod_start")) 


##put model outputs in table together ####
#Frequentistic
tab_model(m_duration, m_distance, m_disp, m_dist, m_npatches, file="modeloutputs_2022-06-07.html")

#Figures ####

newdf <- expand.grid(explor_Log10_mean_speed=seq(-0.15,1.25,by=0.2),
                     tod_start = levels(factor(df2$tod_start)),
                     inout=levels(factor(df2$inout)),
                     scale_space=mean(df2$scale_space),
                     tag=levels(factor(df2$tag)),
                     tideID=levels(factor(df2$tideID)),
                     year=levels(factor(df2$year)))

newdf$start_dist_to_tideline2 <- predict(m_tideline,newdata=newdf,type="response")

tidal <- df2 %>% 
  group_by(tag, tod_start, inout) %>% 
  summarise(se = sd(start_dist_to_tideline2_m)/sqrt(length(start_dist_to_tideline2_m)),
            start_dist_to_tideline2_m = mean(start_dist_to_tideline2_m),
            explor_Log10_mean_speed=unique(expl_round)) %>% 
  mutate(explor_Log10_mean_speed = as.numeric(explor_Log10_mean_speed))


p <- ggplot(df2, aes(x=explor_Log10_mean_speed, y=start_dist_to_tideline2_m, 
                     color=as.factor(inout))) +
  #  geom_hline(yintercept=0, linetype="dotted") +
  geom_errorbar(aes(ymin=start_dist_to_tideline2_m-se, ymax=start_dist_to_tideline2_m+se),
                data=tidal, size=0.1, width=0, alpha=0.5) +
  geom_point(data=tidal, aes(x=explor_Log10_mean_speed, y=start_dist_to_tideline2_m),
             size=1, alpha= 0.5) +
  scale_colour_manual(values=c("#999999","#0072B2")) +
  facet_wrap(~tod_start) +
  ylab("Distance to tide line (m)")+ 
  geom_smooth(method="lm", data= newdf, aes(x=explor_Log10_mean_speed,y=start_dist_to_tideline2),
              size=1) +
  theme_base() + ylim(-4,8) +
  xlab("Exploration speed (log10 cm/s)") +
  theme(legend.position= "none")
p

library(grid)
grid.arrange(p1,p2,p3,p4)

#lenght of day/night ####
library(maptools)

# these functions need the lat/lon in an unusual format
portsmouth <- matrix(c(53.252197, 5.254093), nrow=1)
for_date <- as.POSIXct("2018-08-01", tz="Amsterdam")
sunriset(portsmouth, for_date, direction="sunrise", POSIXct.out=TRUE)

ephemeris <- function(lat, lon, date, span=1, tz="UTC") {
  
  # convert to the format we need
  lon.lat <- matrix(c(lon, lat), nrow=1)
  
  # make our sequence - using noon gets us around daylight saving time issues
  day <- as.POSIXct(date, tz=tz)
  sequence <- seq(from=day, length.out=span , by="days")
  
  # get our data
  sunrise <- sunriset(lon.lat, sequence, direction="sunrise", POSIXct.out=TRUE)
  sunset <- sunriset(lon.lat, sequence, direction="sunset", POSIXct.out=TRUE)
  solar_noon <- solarnoon(lon.lat, sequence, POSIXct.out=TRUE)
  
  # build a data frame from the vectors
  data.frame(date=as.Date(sunrise$time),
             sunrise=as.numeric(format(sunrise$time, "%H%M")),
             solarnoon=as.numeric(format(solar_noon$time, "%H%M")),
             sunset=as.numeric(format(sunset$time, "%H%M")),
             day_length=as.numeric(sunset$time-sunrise$time))
  
}

d <- ephemeris(53.252197, 5.254093, "2018-08-10", 900, tz="Amsterdam")
d$night_lenght <- 24 - (d$day_length)
d$night_day_ratio <- d$night_lenght/d$day_length

df3$date <- as.Date(df3$datetime_start_UTC, "%Y-%m-%d")

df3 <- left_join(df3,d,by="date")


