######################################
#### Fishery closures MS - Models ####
######################################


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



################
#### MODELS ####
################

### FISH  SCHOOLS ####
# number of schools model #
n_full <- glmer(nschools ~ Closure + Island + Closure:Island + 
                  (1|Year) + (1|Month:Island),
                data = surveyeffort,
                family = poisson(link = "log"), offset = log(surveydur_h))
options(na.action = "na.fail")
dredge(n_full) # closure only model selected

n_best <- glmer(nschools ~ Closure +
                  (1|Year) + (1|Month:Island),
                data = surveyeffort,
                family = poisson(link = "log"), offset = log(surveydur_h))


#plot residuals
simulation.mod<- simulateResiduals(n_best, plot = F) 
plot(simulation.mod) # qq plot looks better with a sqrt link but predictions look weird 
# and don't know what to do with the offset (shouldn't log it?)

sresid_m1 <- resid(n_best, type = "pearson")

# normal residuals
hist(n_best, freq = F)
lines(density(sresid_m1, adjust = 1))
shapiro.test(sresid_m1) # residuals are normal

# heteroscedacity
plot(residuals(simulationOutputall.m1)) # no pattern in residuals

# dispersion
testDispersion(n_best) # not significantly overdispersed



# school depth model #
depth.all <- glmer(Depth_mean ~ Closure + Island + Closure:Island +
                     (1|Year) + (1|Month:Island), 
                   data = fish, family = gaussian(link = "log"))

dredge(depth.all) # saturated model selected

depth.best <- glmer(Depth_mean ~ Closure + Island + Closure:Island +
                      (1|Year) + (1|Month:Island), 
                    data = fish, family = gaussian(link = "log"))

summary(depth.best)
confint(depth.best, "ClosureOpen", method = "Wald")


simulationOutputall.depth <- simulateResiduals(depth.best, plot = F) 

#plot residuals
plot(simulationOutputall.depth)


#model checking for depth.all
hist(residuals(depth.best))

#dispersion
testDispersion(depth.best)



# school length model #

length.all <- glmer(Corrected_length ~ Closure + Island + Closure:Island +
                      (1|Year) + (1|Month:Island), 
                    data = fish, family = gaussian(link = "log"))

dredge(length.all) # closure only

length.best <-  glmer(Corrected_length ~ Closure + 
                        (1|Year) + (1|Month:Island), 
                      data = fish, family = gaussian(link = "log"))
# remove year only random effect to get rid of singular fit
length.best <-  glmer(Corrected_length ~ Closure + 
                        (1|Month:Island), 
                      data = fish, family = gaussian(link = "log"))

simulationOutputall.length <- simulateResiduals(length.best, plot = F) 

summary(length.best)
confint(length.best, "ClosureOpen", method = "Wald")

#plot residuals
plot(simulationOutputall.length)

#dispersion
testDispersion(length.best)



# school height model #

height.all <- glmer(Height_mean ~ Closure + Island + Closure:Island +
                      (1|Year) + (1|Month:Island), 
                    data = fish, family = gaussian(link = "log"))

dredge(height.all) # closure only

height.best <- glmer(Height_mean ~ Closure + 
                       (1|Year) + (1|Month:Island), 
                     data = fish, family = gaussian(link = "log"))
#remove Year/Month random to remove singular fit
height.best <- glmer(Height_mean ~ Closure + 
                       (1|Year), 
                     data = fish, family = gaussian(link = "log"))

summary(height.best)
confint(height.best, "ClosureOpen", method = "Wald")


simulationOutputall.height <- simulateResiduals(height.best, plot = F) 

#plot residuals
plot(simulationOutputall.height)


#dispersion
testDispersion(height.best)

### PENGUIN TRIPS ####
# path length model #
PL.ALL.guassian<-  lmer(log(path_length_km) ~ Closure + Colony + Closure*Colony + Sex + 
                          (1|Year_f) + (1|Month_f:Colony),
                        data = trip, REML = FALSE, lmerControl(optimizer = "bobyqa"))

dredge(PL.ALL.guassian) #closure only

PL.best <- lmerTest::lmer(log(path_length_km) ~ Closure + 
                            (1|Year_f) + (1|Month_f:Colony),
                          data = trip, REML = TRUE, lmerControl(optimizer = "bobyqa"))

#calculate scaled residuals
simulationOutputPLguassian <- simulateResiduals(PL.best, plot = F) #they look uniformly distributed

#plot residuals
plot(simulationOutputPLguassian)


#dispersion
testDispersion(PL.best)


# trip duration model #

TD.ALL.guassian<-  lmer(log(duration_min) ~ Closure + Colony + Closure*Colony + Sex + 
                          (1|Year_f) + (1|Month_f:Colony),
                        data = trip, REML = FALSE, lmerControl(optimizer = "bobyqa"))


dredge(TD.ALL.guassian) #closure only

TD.best<-  lmerTest::lmer(log(duration_min) ~ Closure + 
                            (1|Year_f) + (1|Month_f:Colony),
                          data = trip, REML = TRUE, lmerControl(optimizer = "bobyqa"))
summary(TD.best)
confint(TD.best, "ClosureOpen")

#calculate scaled residuals
simulationOutputTDgaussian <- simulateResiduals(TD.best, plot = F) 
plot(residuals(simulationOutputTDgaussian)) #they look uniformly distributed


#plot residuals
plot(simulationOutputTDgaussian)#not great

#dispersion
testDispersion(TD.best) #fine



# maximum distance from nest model #

MD.ALL.gaussian<- lmer(log(max_dist_km) ~ Closure + Colony  + Colony*Closure + Sex + 
                         (1|Year_f) + (1|Month_f:Colony), 
                       data=trip, lmerControl(optimizer = "bobyqa"), REML = FALSE)


dredge(MD.ALL.gaussian)  # no interaction model selected             


MDbest<- lmerTest::lmer(log(max_dist_km) ~ Closure + Colony +
                          (1|Year_f) + (1|Month_f:Colony), 
                        data=trip, lmerControl(optimizer = "bobyqa"), REML = TRUE)
summary(MDbest)
confint(MDbest, "ClosureOpen")

#calculate scaled residuals
simulationOutputbest.gaussian <- simulateResiduals(MDbest, plot = F) 
plot(residuals(simulationOutputbest.gaussian))
#they look uniformly distributed

#plot residuals
plot(simulationOutputbest.gaussian)
# QQ plot looks good

#dispersion
testDispersion(MD4best) #fine


### PENGUIN DIVES ####
# bottom time model #
BT.ALL <- glmer(botttim ~ Closure + Colony + Closure*Colony + Sex + 
                  (1|Year_f) + (1|Month_f:Colony) + (1|deployID),
                data = sexdives, family = gaussian(link = "log"))

dredge(BT.ALL) # closure and sex only model chosen

BT.best <- glmer(botttim ~ Closure + Sex + 
                   (1|Year_f) + (1|deployID),
                 data = sexdives, family = gaussian(link = "log"))
summary(BT.best)
confint(BT.best, "ClosureOpen", method = "Wald")

#calculate scaled residuals
simulationOutput.BT <- simulateResiduals(BT.best, plot = F) 
plot(residuals(simulationOutput.BT))
#they look uniformly distributed

#plot residuals
plot(simulationOutput.BT)  # much better


#dispersion
testDispersion(BT.best) #looks good



# maximum depth model #

MD.ALL <- glmer(maxdep ~ Closure + Colony + Closure*Colony + Sex + 
                  (1|Year_f) + (1|Month_f:Colony) + (1|deployID),
                data = sexdives, family = gaussian(link = "log"))


dredge(MD.ALL) # closure and sex only model chosen

MD.best <-  glmer(maxdep ~ Closure + Sex + 
                    (1|Year_f) + (1|Month_f:Colony) + (1|deployID),
                  data = sexdives, family = gaussian(link = "log"))
summary(MD.best)
confint(MD.best, "ClosureOpen", method = "Wald")

#calculate scaled residuals
simulationOutput.MD <- simulateResiduals(MD.best, plot = F) 
plot(residuals(simulationOutput.MD))
#they look uniformly distributed

#plot residuals
plot(simulationOutput.MD) #much better fit

#dispersion
testDispersion(MD.best) #looks good



# dive time model #

DT.ALL <- glmer(divetim ~ Closure + Colony + Closure*Colony + Sex + 
                  (1|Year_f) + (1|Month_f:Colony) + (1|deployID),
                data = sexdives, family = gaussian(link = "log"))


dredge(DT.ALL) # closure and sex only model chosen

DT.best <- glmer(divetim ~ Closure + Sex + 
                   (1|Year_f) + (1|Month_f:Colony) + (1|deployID),
                 data = sexdives, family = gaussian(link = "log"))
summary(DT.best)
confint(DT.best, "ClosureOpen", method = "Wald")

#calculate scaled residuals
simulationOutput.DT <- simulateResiduals(DT.best, plot = F) 
plot(residuals(simulationOutput.DT))
#they look uniformly distributed

#plot residuals
plot(simulationOutput.DT)#much better


#dispersion
testDispersion(DT.best)#ok


# surface interval model #

pdsi.all<-   glmer(postdive.dur ~ Closure + Colony + Sex + Closure*Colony + 
                     (1|Year_f) + (1|Month_f:Colony) + (1|deployID),
                   data = sexpdsi, family = gaussian(link = "log"))


dredge(pdsi.all) # closure and colony only model

pdsi.best<-  glmer(postdive.dur ~ Closure + Colony + 
                     (1|Year_f) + (1|deployID),
                   data = sexpdsi, family = gaussian(link = "log"))

summary(pdsi.best)
confint(pdsi.best, "ClosureOpen", method = "Wald")

#calculate scaled residuals
simulationOutput.pdsigaussian <- simulateResiduals(pdsi.best, plot = F) 
plot(residuals(simulationOutput.pdsigaussian))
#they look uniformly distributed

#plot residuals
plot(simulationOutput.pdsigaussian)#much better


#dispersion
testDispersion(pdsi.best) #ok



# dive shape model #
data_UV$Month_f <- as.factor(data_UV$Month)

model.UV <- glmer(UorV ~  Colony + Closure + Closure*Colony + Sex + (1|Year) + (1|Month_f:Colony) + (1|deployID),
                  family = binomial(link = "logit"), data = data_UV)

dredge(model.UV) # closure and sex only model

UV.best <- glmer(UorV ~  Closure + Sex + 
                   (1|Year) + (1|Month_f:Colony) + (1|deployID),
                 family = binomial(link = "logit"), data = data_UV)

#calculate scaled residuals
simulationOutput.UV <- simulateResiduals(UV.best, plot = F) 
plot(residuals(simulationOutput.UV))
#they look uniformly distributed

#plot residuals
plot(simulationOutput.UV)

r.squaredGLMM(UV.best) 

sresidUV <- resid(UV.best, type = "pearson")
hist(sresidUV) 
fitsUV <- fitted(UV.best)
plot(sresidUV ~ fitsUV) 

#dispersion
testDispersion(UV.best) #ok


#####################################
#### MODEL PREDICTIONS AND PLOTS ####
#####################################

pal2 <- c("darkorange3","#68228B")

### FISH  SCHOOLS ####
# number of schools #
summary(n_best)
confint(n_best, "ClosureOpen", method = "Wald")

ggeffects::predict_response(n_best, type = "random", 
                            terms = c("Closure"), 
                            interval = "confidence")


nschools <- plot(ggeffects::predict_response(n_best, type = "random", 
                                             terms = c("Closure"), 
                                             interval = "confidence")) + 
  ylab("Number of schools") +
  labs(title = "") + theme_classic() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14),
        axis.title = element_text(size = 14)) +
  annotate("text", x=1.5, y=80, size = 8, label= "***") + 
  annotate("text", x=0.9, y=80, size = 8, label= "a)")   +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2)

# school depth #
summary(depth.best)
depth.int <- plot(ggeffects::predict_response(depth.best, type = "random", 
                                              terms = c("Island","Closure"), 
                                              interval = "confidence")) + ylab("School depth (m)") +
  labs(title = "") + theme_classic() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14),
        axis.title = element_text(size = 14)) +
  annotate("text", x=1.5, y=70, size = 8, label= "***") + 
  annotate("text", x=0.9, y=70, size = 8, label= "b)")   +
  scale_fill_manual(values = pal2)  +
  scale_colour_manual(values = pal2)



ggeffects::predict_response(depth.best, type = "random", 
                            terms = c("Closure","Island"), 
                            interval = "confidence")


ggeffects::predict_response(depth.best, type = "random", 
                            terms = c("Closure"), 
                            interval = "confidence")


# school length #
summary(length.best)
length.cl <- 
  plot(ggeffects::predict_response(length.best, type = "random", 
                                   terms = c("Closure"), 
                                   interval = "confidence")) + ylab("School length (m)") +
  labs(title = "") + theme_classic() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14),
        axis.title = element_text(size = 14)) +
  annotate("text", x=1.5, y=155, size = 8, label= "***") + 
  annotate("text", x=0.8, y=155, size = 8, label= "c)") 

ggsave(length.cl, file = paste0(figure_dir,"Length closure predictions.png"), width =  5, height = 6)

ggeffects::predict_response(length.best, type = "random", 
                            terms = c("Closure"), 
                            interval = "confidence")


# school height #
summary(height.best)
height.cl <- 
  plot(ggeffects::predict_response(height.best, type = "random", 
                                   terms = c("Closure"), 
                                   interval = "confidence")) +
  labs(title = "") + theme_classic() + 
  theme(axis.text=element_text(size=14),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=1.5, y=8.4, size = 8, label= "***") + 
  annotate("text", x=0.8, y=8.4, size = 8, label= "d)")  +
  labs(y = "School height (m)") 

height.cl
ggeffects::predict_response(height.best, type = "random", 
                            terms = c("Closure"), 
                            interval = "confidence")


# school figures
schoolfigs <- ggpubr::ggarrange(nschools, depth.int, length.cl , height.cl, 
                                ncol = 2, nrow = 2, 
                                legend = "bottom", common.legend = TRUE)
schoolfigs
ggsave(schoolfigs, file = paste0(figure_dir,"School closure effects.png"), width =  30, height = 30,units ="cm")




### PENGUIN TRIPS ####

# path length #
summary(PL.best)
PLpredict <- 
  plot(ggeffects::predict_response(PL.best, type = "random", 
                                   terms = c("Closure"), 
                                   interval = "confidence")) + labs(title = "",y = "Path length (km)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=0.8, y=60, size = 5, label= "a)")


ggeffects::predict_response(PL.best, type = "random", 
                            terms = c("Closure"), 
                            interval = "confidence")


# trip duration #

TDpredict <- 
  plot(ggeffects::predict_response(TD.best, type = "random", 
                                   terms = c("Closure"), 
                                   interval = "confidence")) + 
  labs(title = "",y = "Trip duration (min)") + theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=14), 
        axis.text=element_text(size=14)) +
  annotate("text", x=1.5, y=980, size = 8, label= "*")


ggeffects::predict_response(TD.best, type = "random", 
                            terms = c("Closure"), 
                            interval = "confidence")


ggsave(TDpredict, file = (paste0(figure_dir,"trip.model.closures.png")), width = 4, height = 4)

MDpredict <- 
  plot(ggeffects::predict_response(MDbest, type = "random", 
                                   terms = c("Closure"), 
                                   interval = "confidence")) + labs(title = "",y = "Maximum distance from nest (km)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=0.8, y=15, size = 5, label= "c)")



### PENGUIN DIVES ####

# bottom time #
summary(BT.best)
BTpredict <- 
  ggeffects::ggpredict(BT.best, type = "fixed", terms = c("Sex"))  %>% 
  plot(add.data = FALSE, colors = "black", show.title = FALSE, log.y = FALSE) + 
  labs(y = "Bottom time (s)") + theme_classic()  + 
  annotate("text", x=1.5, y=50, size = 10, label= "***") + 
  annotate("text", x=0.8, y=50, size = 8, label= "a)")

ggeffects::predict_response(BT.best, type = "random", 
                            terms = c("Sex"), 
                            interval = "confidence")

BTpredict2 <- 
  plot(ggeffects::predict_response(BT.best, type = "random", 
                                   terms = c("Closure"), 
                                   interval = "confidence")) + 
  labs(title = "", y = "Bottom time (s)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=1, y=17.8, size = 8, label= "a)")

ggeffects::predict_response(BT.best, type = "random", 
                            terms = c("Closure"), 
                            interval = "confidence")

# maximum depth #
summary(MD.best)
MDpredict <- 
  plot(ggeffects::predict_response(MD.best, type = "random", 
                                   terms = c("Closure"), 
                                   interval = "confidence")) + 
  labs(title = "", y = "Maximum depth (m)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=1, y=23.5, size = 8, label= "b)")

ggeffects::predict_response(MD.best, type = "random", 
                            terms = c("Closure"), 
                            interval = "confidence")

MDpredict1 <- 
  plot(ggeffects::predict_response(MD.best, type = "random", 
                                   terms = c("Sex"), 
                                   interval = "confidence")) + 
  labs(title = "", y = "Maximum depth (m)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=1.5, y=23.5, size = 6, label= "*") +
  annotate("text", x=0.8, y=23.5, size = 6, label= "b)")


ggeffects::predict_response(MD.best, type = "random", 
                            terms = c("Sex"), 
                            interval = "confidence")

# dive time #
summary(DT.best)
DTpredict <- 
  plot(ggeffects::predict_response(DT.best, type = "random", 
                                   terms = c("Closure"), 
                                   interval = "confidence")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Dive duration (s)", title = "") + 
  annotate("text", x=1, y=58, size = 8, label= "c)")

ggeffects::predict_response(DT.best, type = "random", 
                            terms = c("Closure"), 
                            interval = "confidence")


DTpredict2 <- 
  plot(ggeffects::predict_response(DT.best, type = "random", 
                                   terms = c("Sex"), 
                                   interval = "confidence")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Dive duration (s)", title = "") + 
  annotate("text", x=1.5, y=59, size = 10, label= "**") + 
  annotate("text", x=1, y=59, size = 8, label= "b)")


ggeffects::predict_response(DT.best, type = "random", 
                            terms = c("Sex"), 
                            interval = "confidence")


# surface interval #

pdsipredict <- 
  plot(ggeffects::predict_response(pdsi.best, type = "random", 
                                   terms = c("Closure"), 
                                   interval = "confidence")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Surface interval (s)", title = "") + 
  annotate("text", x=1, y=58.5, size = 8, label= "d)")

ggeffects::predict_response(pdsi.best, type = "random", 
                            terms = c("Closure"), 
                            interval = "confidence")

pdsipredict2 <- 
  plot(ggeffects::predict_response(pdsi.best, type = "random", 
                                   terms = c("Colony"), 
                                   interval = "confidence")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Surface interval (s)", title = "") +
  annotate("text", x=1.5, y=58, size = 10, label= "*") 

ggeffects::predict_response(pdsi.best, type = "random", 
                            terms = c("Colony"), 
                            interval = "confidence")


# dive figures
dive.model.closure <- ggpubr::ggarrange(BTpredict2, MDpredict, DTpredict, pdsipredict, ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")
ggsave(dive.model.closure, file = (paste0(figure_dir,"dive.model.closure.png")), width = 8, height = 8)


dive.model.sex <- ggpubr::ggarrange(BTpredict, DTpredict2, MDpredict1, ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom")
ggsave(dive.model.sex, file = (paste0(figure_dir,"dive.model.sex.png")), width = 5, height = 5)


# UV dive proportions #

UVpredict <- 
  plot(ggeffects::predict_response(UV.best, type = "random", 
                                   terms = c("Closure"), 
                                   interval = "confidence")) + 
  labs(title = "", y = "Prercentage of foraging dives") + 
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=1.5, y=0.7, size = 8, label= "**") + 
  annotate("text", x=1, y=0.7, size = 8, label= "a)") 

ggeffects::predict_response(UV.best, type = "random", 
                            terms = c("Closure"), 
                            interval = "confidence")

UVpredict2 <- 
  plot(ggeffects::predict_response(UV.best, type = "random", 
                                   terms = c("Sex"), 
                                   interval = "confidence")) + 
  labs(title = "", y = "Percentage of foraging dives") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x=1.5, y=0.72, size = 8, label= "*") + 
  annotate("text", x=1, y=0.72, size = 8, label= "b)")


# dive shape figures
UVdive.model <- ggpubr::ggarrange(UVpredict, UVpredict2,  ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")
ggsave(UVdive.model, file = (paste0(figure_dir,"UVdive.model.png")), width = 5, height = 5)


