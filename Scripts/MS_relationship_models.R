##########################################################
#### Fishery closures MS - Metric relationship models ####
##########################################################


# load packages ####
pacman::p_load(dplyr, ggplot2, lme4, MuMIn, DHARMa, lubridate)



# directories ####
figure_dir <- "~/PhD/Papers/Closures/"
templist <- readRDS("E:/Chapter 3 Fishery closures/surveytimes.RDS")
setwd("E:/Chapter 3 Fishery closures/")



#load files and format data ####
# data for the schools metric models
fish_n  <- read.csv("Calibrated_STEP 4 CLUSTERED SCHOOLS.csv") %>%
  filter(Exclude_below_line_depth_mean > -999) %>%
  dplyr::select(-X, -X.1) %>%
  filter(Depth_mean < 100) %>%# filter fish below 100 m, at end of pipeline we plotted them and they look erroneous
  filter(IslandYear != "Robben_2016") %>%
  mutate(Time = hms(fish$Time_S),
         Year = as.factor(Year),
         Month = as.factor(Month)) %>%
  dplyr::group_by(Year, Island, Closure, surveyID, Month) %>%
  dplyr::summarise(nschools = length(schoolID))# %>%



# include survey effort data for the number of schools models
surveytimes <- data.table::rbindlist(templist) %>%
  mutate(evi = grepl('.evi', id, fixed=TRUE)) %>%
  filter(evi == FALSE) %>%
  dplyr::select(-evi) %>%
  mutate(duration = end - start,
         surveyID = substr(id, 1,20)) %>%
  group_by(surveyID) %>%
  summarise(surveydur = sum(duration),
            surveydur_h = as.numeric(surveydur, "hours")) %>%
  ungroup() %>%
  mutate(Island = substr(surveyID, 1,6),
         Year = substr(surveyID, 8,11))

surveyeffort <- merge(fish_n, surveytimes) %>%
  mutate(schoolperhour = nschools/surveydur_h) 



# data for the trip metrics
trip <- read.csv("Trip_stats_filtered.csv") %>%
  mutate(Year_f = as.factor(Year),
         Month_f = as.factor(Month))



# data for the dive metrics including only the individuals which could be sexed
sexdives <- read.csv("dives_U.csv", header=T) %>% filter(Sex != "U") %>% 
  mutate(Sex = case_when( 
    Sex == "F" ~ "Female", 
    Sex == "M" ~ "Male"),
    Year_f = as.factor(Year),
    Month_f = as.factor(Month),
    deployID = as.factor(deployID)) %>%
  filter(Year < 2015)



# surface intervals have some very long outliers. 
# we think these are due to travel being detected as a surface interval
# use quantiles to crop outliers
sexpdsi <- subset(sexdives, postdive.dur > 0)
l_u <- quantile(sexpdsi$postdive.dur, probs = 0.975)  #2 SD from mean #bouts of travelling are included in surface intervals and so 95% quantile was used to exclude trave;
sexpdsi<- subset(sexpdsi, postdive.dur < l_u)



#turn into binary data for the dive shape model
month <- read.csv("month_tdr.csv", header=T) 
data_UV <- read.csv("Divestats.csv", header=T) %>%
  mutate(Year = as.factor(Year),
         Closure = as.factor(Closure),
         Sex = replace(Sex, is.na(Sex), "U"),
         Sex = as.factor(Sex)) %>%  
  inner_join(., month, by = "deployID") %>%
  mutate(Shape_4 = ifelse(botttim < 4, "V", "U")) %>% #define shape of dive on 4 seconds
  na.omit(.) %>%
  mutate(UorV = ifelse(Shape_4 == "U", "1", "0"),
         UorV = as.numeric(UorV),
         Month_f = as.factor(Month)) %>%
  filter(Sex != "U")



####################################
#### METRIC RELATIONSHIP MODELS ####
####################################

### school models ####

dep_hei_full <- lmerTest::lmer(Depth_mean ~ Height_mean + Closure + 
                       Height_mean*Closure + (1|Year) + (1|Month:Island), 
                     data=fish, lmerControl(optimizer = "bobyqa"), REML = FALSE)

options(na.action = "na.fail")
dredge(dep_hei_full)               

dep_hei_best <- lmerTest::lmer(Depth_mean ~ Height_mean + Closure + 
                       Height_mean*Closure + 
                       + (1|Year) + (1|Month:Island), 
                     data=fish, lmerControl(optimizer = "bobyqa"), REML = TRUE)
summary(dep_hei_best)
confint(dep_hei_best, "Height_mean:ClosureOpen", method = "Wald")
simulationOutputall.dep_hei_best <- simulateResiduals(dep_hei_best, plot = F) 

#plot residuals
plot(simulationOutputall.dep_hei_best)

testDispersion(dep_hei_best)



len_dep_full <- lmerTest::lmer(Corrected_length ~ Depth_mean + Closure + 
                       Depth_mean*Closure + (1|Year) + (1|Month:Island), 
                     data=fish, lmerControl(optimizer = "bobyqa"), REML = FALSE)


dredge(len_dep_full)               

len_dep_best <- lmerTest::lmer(Corrected_length ~ Depth_mean + Closure + 
                       Depth_mean*Closure + 
                       + (1|Year) + (1|Month:Island), 
                     data=fish, lmerControl(optimizer = "bobyqa"), REML = TRUE)
summary(len_dep_best)
confint(len_dep_best, "Depth_mean:ClosureOpen", method = "Wald")
simulationOutputall.len_dep_best <- simulateResiduals(len_dep_best, plot = F) 


#plot residuals
plot(simulationOutputall.len_dep_best)

testDispersion(len_dep_best)



len_hei_full <- lmerTest::lmer(Corrected_length ~ Height_mean + Closure + 
                       Height_mean*Closure + (1|Year) + (1|Month:Island), 
                     data=fish, lmerControl(optimizer = "bobyqa"), REML = FALSE)


dredge(len_hei_full)               

len_hei_best <- lmerTest::lmer(Corrected_length ~ Height_mean + Closure + 
                       Height_mean*Closure + (1|Year) + (1|Month:Island), 
                     data=fish, lmerControl(optimizer = "bobyqa"), REML = TRUE)
summary(len_hei_best)
confint(len_hei_best, "Height_mean:ClosureOpen", method = "Wald")

simulationOutputall.len_hei_best<- simulateResiduals(len_hei_best, plot = F) 


#plot residuals
plot(simulationOutputall.len_hei_best)

testDispersion(len_hei_best)



### penguin trip models ####

max_dur_full <- lmerTest::lmer(max_dist_km ~ duration_min + Closure + 
                       duration_min*Closure + (1|Year_f) + (1|Month_f:Colony), 
                     data=sex,lmerControl(optimizer = "bobyqa"),  REML = FALSE)


dredge(max_dur_full)               

max_dur_best <- lmerTest::lmer(max_dist_km~ duration_min + Closure + 
                       duration_min*Closure + (1|Year_f) + (1|Month_f:Colony), 
                     data=sex, lmerControl(optimizer = "bobyqa"), REML = TRUE)
summary(max_dur_best)
confint(max_dur_best, "duration_min:ClosureOpen", method = "Wald")

simulationOutputbest.max_dur_best <- simulateResiduals(max_dur_best, plot = F) 
#plot residuals
plot(simulationOutputbest.max_dur_best)



pat_dur_full <- lmerTest::lmer(path_length_km ~ duration_min + Closure + 
                       duration_min*Closure + (1|Year_f) + (1|Month_f:Colony), 
                     data=sex, REML = FALSE)

dredge(pat_dur_full)               


pat_dur_best <- lmerTest::lmer(path_length_km ~ duration_min + Closure + 
                       duration_min*Closure + (1|Year_f), 
                     data=sex, REML = TRUE)

summary(pat_dur_best)
confint(pat_dur_best, "duration_min:ClosureOpen", method = "Wald")



max_pat_full <- lmerTest::lmer(max_dist_km ~ path_length_km + Closure + 
                       path_length_km*Closure + (1|Year_f) + (1|Month_f:Colony), 
                     data=sex, REML = FALSE)

dredge(max_pat_full)               


max_pat_best <- lmerTest::lmer(max_dist_km ~ path_length_km + Closure + 
                       path_length_km*Closure + (1|Year_f) + (1|Month_f:Colony), 
                     data=sex, REML = TRUE)
summary(max_pat_best)
confint(max_pat_best, "path_length_km:ClosureOpen", method = "Wald")



### penguin dive models ####

tim_bot_full <- lmerTest::lmer(divetim ~ botttim + Closure + 
                       botttim*Closure + (1|Year_f) + (1|Month_f:Colony) + (1|deployID), 
                     data=sexpdsi, REML = FALSE)

dredge(tim_bot_full)               


tim_bot_best <- lmerTest::lmer(divetim ~ botttim + Closure + 
                       botttim*Closure + (1|Year_f) + (1|Month_f:Colony) + (1|deployID), 
                     data=sexpdsi)

summary(tim_bot_best)
confint(tim_bot_best, "botttim:ClosureOpen", method = "Wald")



bot_dep_full <- lmerTest::lmer(botttim ~ maxdep + Closure + 
                       maxdep*Closure + (1|Year_f) + (1|Month_f:Colony) + (1|deployID), 
                     data=sexpdsi, REML = FALSE, lmerControl(optimizer = "bobyqa"))

dredge(bot_dep_full)               


bot_dep_best <- lmerTest::lmer(botttim ~ maxdep + Closure + 
                       maxdep*Closure + (1|Year_f)  + (1|deployID), 
                     data=sexpdsi, REML = TRUE, lmerControl(optimizer = "bobyqa"))
summary(bot_dep_best)
confint(bot_dep_best, "maxdep:ClosureOpen", method = "Wald")


tim_sur_full <- lmerTest::lmer(divetim ~ postdive.dur + Closure + 
                       postdive.dur*Closure + (1|Year_f) + (1|Month_f:Colony) + (1|deployID), 
                     lmerControl(optimizer = "bobyqa"), 
                     data=sexpdsi, REML = FALSE)


dredge(tim_sur_full)               


tim_sur_best <- lmerTest::lmer(divetim ~ postdive.dur + Closure + 
                       postdive.dur*Closure + (1|Year_f) + (1|Month_f:Colony) + (1|deployID), 
                     lmerControl(optimizer = "bobyqa"), 
                     data=sexpdsi, REML = TRUE)
summary(tim_sur_best)
confint(tim_sur_best, "postdive.dur:ClosureOpen", method = "Wald")



bot_sur_full <- lmerTest::lmer(botttim ~ postdive.dur + Closure + 
                       postdive.dur*Closure + (1|Year_f) + (1|Month_f:Colony) + (1|deployID), 
                     data=sexpdsi, REML = FALSE)


dredge(bot_sur_full)   # saturated model not chosen as best model            



dep_tim_full <- lmerTest::lmer(maxdep ~ divetim + Closure + 
                       divetim*Closure + (1|Year_f) + (1|Month_f:Colony) + (1|deployID), 
                     data=sexpdsi, REML = FALSE, lmerControl(optimizer = "bobyqa"))


dredge(dep_tim_full)               


dep_tim_best <- lmerTest::lmer(maxdep ~ divetim + Closure + 
                       divetim*Closure + 
                       (1|Year_f) + (1|Month_f:Colony) + (1|deployID), 
                     data=sexpdsi, REML = TRUE, lmerControl(optimizer = "bobyqa"))
summary(dep_tim_best)

confint(dep_tim_best, "divetim:ClosureOpen", method = "Wald")


dep_sur_full <- lmerTest::lmer(maxdep ~ postdive.dur + Closure + 
                       postdive.dur*Closure + (1|Year_f) + (1|Month_f:Colony) + (1|deployID), 
                     data=sexpdsi, REML = FALSE, lmerControl(optimizer = "bobyqa"))


dredge(dep_sur_full)               


dep_sur_best <- lmerTest::lmer(maxdep ~ postdive.dur + Closure + 
                       postdive.dur*Closure + (1|Year_f) + (1|Month_f:Colony) + (1|deployID), 
                     data=sexpdsi, REML = TRUE, lmerControl(optimizer = "bobyqa"))
summary(dep_sur_best)
confint(dep_sur_best, "postdive.dur:ClosureOpen", method = "Wald")



#############################
## PLOTS AND PREDICTIONS ####
#############################
pal2 <- c("darkorange3","#68228B")

# SCHOOLS ####

dep_hei_plot <- 
  
  plot(ggeffects::predict_response(dep_hei_best, type = "random", 
                                   terms = c("Height_mean","Closure"), 
                                   interval = "confidence"), limit_range= T, show_data = F) +
  labs(title = "",x = "School height (m)", y = "School depth (m)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=20, y=65, size = 6, label= "a)")  + 
  annotate("text", x=50, y=65, size = 6, label= "***")+
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2) + theme_classic()  



ggeffects::predict_response(dep_hei_best, type = "random", 
                            terms = c("Height_mean","Closure"), 
                            interval = "confidence")



len_dep_plot <- 
  
  plot(ggeffects::predict_response(len_dep_best, type = "random", 
                                   terms = c("Depth_mean","Closure"), 
                                   interval = "confidence")) +
  labs(title = "",x = "School depth (m)", y = "School length (m)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=2, y=160, size = 6, label= "b)")  +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2)  + theme_classic()

ggeffects::predict_response(m2, type = "random", 
                            terms = c("Depth_mean","Closure"), 
                            interval = "confidence")



len_hei_plot <- 
  
  plot(ggeffects::predict_response(len_hei_best, type = "random", 
                                   terms = c("Height_mean","Closure"), 
                                   interval = "confidence"), show_data = F) +
  labs(title = "",x = "School height (m)", y = "School length (m)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=2, y=730, size = 6, label= "c)") +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2) + theme_classic()


figfishparameters <- ggpubr::ggarrange(dep_hei_plot, len_dep_plot, len_hei_plot, ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom")
figfishparameters

ggsave(figfishparameters, file = (paste0(figure_dir,"schoolparameters.png")), width = 8, height = 4)


# PENGUIN TRIPS ####


max_dur_plot <- 
  plot(ggeffects::predict_response(max_dur_best, type = "random", 
                                   terms = c("duration_min","Closure"), 
                                   interval = "confidence")) + labs(title = "", x = "Trip duration (min)", y = "Maximum distance from nest (km)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=1900, y=60, size = 6, label= "***") +
  annotate("text", x=50, y=60, size = 6, label= "a)")  +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2) + theme_classic()


pat_dur_plot <- 
  plot(ggeffects::predict_response(pat_dur_best, type = "random", 
                                   terms = c("duration_min","Closure"), 
                                   interval = "confidence")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=50, y=210, size = 6, label= "b)") + 
  labs(x = "Trip duration (min)", y = "Path length (km)", title = "") + theme_classic()  +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2) 


max_pat_plot <- 
  plot(ggeffects::predict_response(max_pat_best, type = "random", 
                                   terms = c("path_length_km","Closure"), 
                                   interval = "confidence")) + 
  labs(title = "", x = "Path length (km)", y = "Maximum distance from nest (km)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=90, y=55, size = 6, label= "***") +
  annotate("text", x=2, y=55, size = 6, label= "c)") +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2) + theme_classic()


figtripparameters <- ggpubr::ggarrange(max_dur_plot, pat_dur_plot, max_pat_plot, ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom")
figtripparameters
ggsave(figtripparameters, file = (paste0(figure_dir,"figtripparameters.png")), width = 8, height = 4)



# PENGUIN DIVES ####

tim_bot_plot <- 
  plot(ggeffects::predict_response(tim_bot_best, type = "random", 
                                   terms = c("botttim","Closure"), 
                                   interval = "confidence")) +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Bottom time (s)", y = "Dive duration (s)", title = "") +
  annotate("text", x=8, y=155, size = 5, label= "a)") +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2) + theme_classic()




bot_dep_plot <- 
  plot(ggeffects::predict_response(bot_dep_best, type = "random", 
                                   terms = c("maxdep","Closure"), 
                                   interval = "confidence")) +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=30, y=35, size = 5, label= "***")  +
  annotate("text", x=8, y=35, size = 5, label= "b)")  + 
  labs(x = "Maximum depth (m)", y = "Bottom time (s)", title = "")  +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2) + theme_classic()



tim_sur_plot <-   
  plot(ggeffects::predict_response(tim_sur_best, type = "random", 
                                   terms = c("postdive.dur","Closure"), 
                                   interval = "confidence")) +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=175, y=77, size = 5, label= "*")  + 
  annotate("text", x=0.8, y=77, size = 5, label= "c)")  + 
  labs(x = "Dive duration (s)", y = "Surface interval (s)", title = "")  +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2) + theme_classic()



dep_tim_plot <- 
  plot(ggeffects::predict_response(dep_tim_best, type = "random", 
                                   terms = c("divetim","Closure"), 
                                   interval = "confidence")) +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=100, y=85, size = 5, label= "***")  + 
  annotate("text", x=10, y=85, size = 5, label= "d)")  + 
  labs(x = "Dive duration (s)", y = "Maximum depth (m)", title = "")  +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2) + theme_classic()

dep_sur_plot <- 
  plot(ggeffects::predict_response(dep_sur_best, type = "random", 
                                   terms = c("postdive.dur","Closure"), 
                                   interval = "confidence")) +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=175, y=36, size = 5, label= "***")  + 
  annotate("text", x=10, y=36, size = 5, label= "e)")  + 
  labs(x = "Surface interval (s)", y = "Maximum depth (m)", title = "")  +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2) + theme_classic()


figdiveparameters <- ggpubr::ggarrange(tim_bot_plot, bot_dep_plot, tim_sur_plot, dep_tim_plot, dep_sur_plot, ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom")
figdiveparameters

ggsave(figdiveparameters, file = (paste0(figure_dir,"diveparameters.png")), width = 12, height = 8)



#### end #####