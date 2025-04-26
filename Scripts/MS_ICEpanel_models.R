#############################################################################
#### Fishery closures MS - Models with ICE panel random effect structure ####
#############################################################################


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



#########################################
#### ICE Panel recommended structure ####
#########################################


### FISH SCHOOLS ####
# number of schools model #
n_best_ICE <- glmer(nschools ~ Closure +
                      (1|Year/Month/Island),
                    data = surveyeffort,
                    family = poisson(link = "log"), offset = log(surveydur_h))
summary(n_best_ICE)

n_best_ICE <- glmer(nschools ~ Closure +
                      (1|Year:Month:Island) + (1|Month:Island),
                    data = surveyeffort,
                    family = poisson(link = "log"), offset = log(surveydur_h))
summary(n_best_ICE)



# school depth model #
depth.best <- glmer(Depth_mean ~ Closure + Island + Closure:Island +
                      (1|Year:Month),
                    family = gaussian(link = "log"), 
                    glmerControl(optimizer = "bobyqa"),
                    data = fish)
summary(depth.best)



# school length model #
length.best <- glmer(Corrected_length ~ Closure +
                       (1|Year:Month:Island),
                     family = gaussian(link = "log"), 
                     glmerControl(optimizer = "bobyqa"),
                     data = fish)
summary(length.best)



# school height #
height.best <- glmer(Height_mean ~ Closure +
                       (1|Year),
                     family = gaussian(link = "log"), 
                     data = fish)
summary(height.best)


### PENGUIN TRIPS ####
# path length model #

PL.best <- lmerTest::lmer(log(path_length_km) ~ Closure + 
                            (1|Year:Month:Colony) + (1|Year),
                          data = trip, REML = TRUE, lmerControl(optimizer = "bobyqa"))
summary(PL.best)


# trip duration model #

TD.best<-  lmerTest::lmer(log(duration_min) ~ Closure + 
                            (1|Year:Month:Colony) + (1|Year),
                          data = trip, REML = TRUE, lmerControl(optimizer = "bobyqa"))
summary(TD.best)


# maximum distance from nest model #

MDbest<- lmerTest::lmer(log(max_dist_km) ~ Closure + Colony +
                          (1|Year:Month:Colony) + (1|Year:Month) + (1|Year), 
                        data=trip, lmerControl(optimizer = "bobyqa"), REML = TRUE)
summary(MDbest)



### PENGUIN DIVES ####
# bottom time model #

BT.best <- glmer(botttim ~ Closure + Sex + 
                   (1|Year:Month:Colony:deployID) + (1|Year:Month:Colony) + (1|Year),
                 data = sexdives, family = gaussian(link = "log"))
summary(BT.best)



# maximum depth model #

MD.best <-  glmer(maxdep ~ Closure + Sex + 
                    (1|Year:Month:Colony:deployID) + (1|Year:Month:Colony),
                  data = sexdives, family = gaussian(link = "log"))
summary(MD.best)



# dive time model #

DT.best <- glmer(divetim ~ Closure + Sex + 
                   (1|Year:Month:Colony:deployID) + (1|Year:Month:Colony) + (1|Year:Month) + (1|Year),
                 data = sexdives, family = gaussian(link = "log"))
summary(DT.best)



# surface interval model #
pdsi.best<-  glmer(postdive.dur ~ Closure + Colony + 
                     (1|Year:Month:Colony:deployID) + (1|Year:Month:Colony),
                   data = sexpdsi, family = gaussian(link = "log"))

summary(pdsi.best)



# dive shape model #

UV.best <- glmer(UorV ~  Closure + Sex + 
                   (1|Year:Month:Colony:deployID) + (1|Year:Month:Colony) + (1|Year:Month) + (1|Year),
                 family = binomial(link = "logit"), data = data_UV)
summary(UV.best)

