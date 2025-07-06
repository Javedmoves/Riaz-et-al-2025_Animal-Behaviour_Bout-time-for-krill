# Manuscript - Bout time for krill - contrasting Adelie penguin foraging behaviour during years of high and low krill availability
# Author - Javed Riaz
# Version: Final: Penguin bout script

##############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/SAERI files/Javed/Uruguay")

# Read in some relevant packages
library(ggplot2)
library(tidyverse)
library(hms)
library(viridis)
library(scales)
library(readr)
library(grid)
library(raster)
library(stringr)
library(lubridate)
library(sp)
library(diveMove) 
library(sf)
library(rnaturalearthhires)
library(data.table)
library(ggpubr)
library(rasterVis)
library(RColorBrewer)
library(devtools)
library(wesanderson)
library(glmmTMB)
library(performance)
library(gdistance)
library(mgcv)
library(gratia)
library(plotly)
library(geosphere)
library(ggOceanMaps)
library(segmented)
library(scico)
library(magick) 
library(DHARMa)
library(datawizard)
library(effects)
library(emmeans)
library(processx)
library(reticulate)


# Plotting function
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{layer(data = data, stat = StatIdentity, position = PositionIdentity,geom = ggplot2:::GeomCustomAnn,
       inherit.aes = TRUE, params = list(grob = grob,xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))}

###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

# PART 2 - PENGUIN BOUT CLASSIFICATION AND SUMMARIES

###############################################################################################################################################
##############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################



###############################################################################################################################################
# LOAD IN PENGUIN MOVEMENT AND DIVE DATA
###############################################################################################################################################


load("DiveALL_A_Hunt.RData") 

Dives <- Dive_hunt 

unique(Dives$ID) 
unique(Dives$Trip) 
unique(Dives$Year) 


Year19_20 <- Dives %>%
  filter(Year == "2019/20") 
length(unique(Year19_20$ID)) ## 19 individuals


Year20_21 <- Dives %>%
  filter(Year == "2020/21") 
length(unique(Year20_21$ID)) ## 16 individuals - not included in this study

Year21_22 <- Dives %>%
  filter(Year == "2021/22") 
length(unique(Year21_22$ID)) ## 22 individuals


DiveYearsOfInterest <-Year19_20 %>% 
  full_join(Year21_22)

unique(DiveYearsOfInterest$ID)


Dives <- DiveYearsOfInterest %>%
  mutate(randcount = 1) %>%
  group_by(ID, Trip, Year) %>%
  mutate(NoALLDives = sum(randcount)) 


## Criterion to define foraging dives
numhunts <- 4 ### 

HuntOnly <-  Dives %>% ## 
  filter(hunt > numhunts) %>% ## 
  group_by(ID, Trip) %>%
  mutate(AllHUNTCount = sum(hunt)) %>%
  ungroup()


Dives <- Dives %>% # Assigned dives as hunting on non-hunting
  mutate(HuntPresAbs = case_when(hunt > numhunts ~ 'Hunt',
                                 hunt <= numhunts ~ 'No hunt'))


Hunting <- HuntOnly %>%
  group_by(ID, Trip) %>%
  mutate(AllHuntingDives = sum(randcount)) %>%
  ungroup() %>%
  dplyr::select(ID, Trip, AllHUNTCount, AllHuntingDives) %>%
  distinct()

Dives1 <- Dives %>%
  arrange(ID, Trip, start, NoALLDives) %>%
  left_join(Hunting)


IDS <- Dives1 %>% ## Isolate the 146 trip IDs
  group_by(ID, Trip) %>%
  distinct(ID, Trip) %>%
  ungroup()


load("sum_water.RData") ## Cleaned trip durations

sum_water$Trip <- sum_water$ID_Trip

sum_water <- sum_water %>%
  left_join(IDS) %>%
  drop_na(ID) %>%
  group_by(ID, Trip) %>% 
  summarise(TripDuration = sum(trip_water)/60/60) %>%
  distinct() %>%
  ungroup()



###############################################################################################################################################
###############################################################################################################################################
# REMOVE NON-FEEDING DIVES 

## Also checks post dive diff outliers
###############################################################################################################################################
###############################################################################################################################################

AllDiveData <- Year19_20 %>% ## Merge the years of interest
  full_join(Year21_22)

Bouts <- AllDiveData %>%   
  ungroup() %>%            
  arrange(ID, start) %>%
  group_by(ID) %>%
  filter (hunt > numhunts) %>% 
  mutate(DiveEnd = start + seconds(dur)) %>%
  group_by(ID, Trip) %>%
  mutate(diff = start - lag(DiveEnd),
         diff_secs = as.numeric(diff, units = 'secs')) %>%
  drop_na(diff_secs)

outlierssurfs <- Bouts [!(Bouts$diff_secs == "0"),]
plot(outlierssurfs$diff_secs)


## Cut-off for postdive.diff

outlier_2 <- as.numeric(quantile (outlierssurfs$diff_secs, prob = 0.99))
outlier_2

###############################################################################################################################################
###############################################################################################################################################
# APPROACH 1 and 2 - CALCULATE BECs FOR EACH INDIVIDUAL AND THEN TAKE THE MEAN
# This produces BEC values ranging from 21 - 406

###############################################################################################################################################
###############################################################################################################################################


df_total = data.frame() 
df_total

for (a in unique(Bouts$ID)){
  subID <- subset(Bouts, ID == a)
  idDF <- as.data.frame(subID$ID)
  postdives.diff <- abs(diff(subID$diff_secs)) ## 
  postdives.diff <- postdives.diff[postdives.diff < 2397] # 
  lnfreq <- boutfreqs(postdives.diff, bw=1, plot=FALSE) ###
  startval <- boutinit(lnfreq, 300, plot = TRUE) ### 
  bouts2.fit <- fitNLSbouts(lnfreq, start=startval, maxiter=5000)
  bec<- bec(bouts2.fit)
  setBEC <- bec
  
  bec <- as.data.frame(print(bec))
  BoutsData <- subID %>%
    mutate(group = cumsum(diff_secs > setBEC) + 1)
  bind <- cbind(BoutsData, bec)
  df_total <- rbind(df_total,bind)
}


newdatatest <- df_total


FistSeason <- newdatatest %>%
  filter(Year == "2019/20") %>%
  distinct(`print(bec)`, .keep_all = TRUE) 
mean(FistSeason$`print(bec)`) 


SecondSeason <- newdatatest %>%
  filter(Year == "2021/22") %>%
  distinct(`print(bec)`, .keep_all = TRUE) 
mean(SecondSeason$`print(bec)`) 


PopBEC <- FistSeason %>%
  full_join(SecondSeason)

PopBEC$Approach1_ID <- PopBEC$`print(bec)`

unique(PopBEC$Approach1_ID) 
range(PopBEC$Approach1_ID) 
mean(PopBEC$Approach1_ID) 


PopBEC$Approach1_Mean <- mean(PopBEC$Approach1_ID)



###############################################################################################################################################
###############################################################################################################################################
# APPROACH 3 - CALCULATE ONE BEC FOR ALL DATA AGGREGATED -

# Produces a single BEC

###############################################################################################################################################
###############################################################################################################################################

Bouts <- Bouts %>%
  group_by(ID, Trip) %>%
  arrange(start)

postdives.diff <- abs(diff(Bouts$diff_secs)) ## 
postdives.diff <- postdives.diff[postdives.diff < 2397] # 
lnfreq <- boutfreqs(postdives.diff, bw=0.1, plot=FALSE) ###
startval <- boutinit(lnfreq, 200, plot = TRUE) ### 
bouts2.fit <- fitMLEbouts(lnfreq, start=startval, optim_opts0=NULL, optim_opts1=NULL)
bec<- bec(bouts2.fit)
setBEC <- bec

bec <- as.data.frame(print(bec))

setBEC

PopBEC <- PopBEC %>%
  mutate(Approach2 = 68)


###############################################################################################################################################
###############################################################################################################################################
# APPROACH USE IN THE MANUSCRIPT 

## This provides a BEC of 90 seconds
###############################################################################################################################################
###############################################################################################################################################


BreakPoint <- Bouts 

cumu_freq <- cbind(as.data.frame(table(BreakPoint$diff_secs)), cumsum(table(BreakPoint$diff_secs))) %>%
  mutate(CumFreq = `cumsum(table(BreakPoint$diff_secs))`) %>%
  mutate_if(is.factor, as.numeric) %>%
  dplyr::select(-3)

ggplot(cumu_freq, aes(x = Var1, y = CumFreq)) +
  geom_line(group = 1) +
  scale_y_continuous(trans = "log2")

fitted_model_lm <- lm(log(CumFreq) ~ Var1, data = cumu_freq) 
summary(fitted_model_lm)


lm_coeffs <- coef(fitted_model_lm)
lm_coeffs


ggplot(cumu_freq, aes(y = log(CumFreq), x = Var1)) +
  geom_line(group = 1) +
  geom_abline(intercept = lm_coeffs[1],
              slope = lm_coeffs[2], col = "orange")+
  geom_vline(xintercept = 40, linetype = "dashed") +
  geom_vline(xintercept = 90, linetype = "dashed")



get_segments <- segmented(fitted_model_lm, 
                      seg.Z = ~ Var1, 
                      psi = list(Var1 = c(40,90)))


summary(get_segments)

get_segments$psi  ## 90.37
slope(get_segments)


fitted_segs <- fitted(get_segments)
model_segs <- data.frame(Times = cumu_freq$Var1, Freq = fitted_segs)

bec <- 90


PopBEC <- PopBEC %>%
  mutate(Approach3 = bec)


###############################################################################################################################################
###############################################################################################################################################
# PLOT THE DATA TO MAKE COMPARISONS 
###############################################################################################################################################
###############################################################################################################################################


PopBEC_Summary <- PopBEC %>%
  ungroup() %>%
  dplyr::select(ID, Year, 24:26) %>%
  distinct()


FigS2 <- ggplot() + geom_histogram(data = PopBEC_Summary, aes(x = Approach1_ID, fill = Year), colour = "black", position = "dodge") +
  labs (x = "Bout-ending criteria (s)", y = "Number of individuals (#)") +
  theme_bw() + 
  scale_y_continuous(expand = c(0,0), breaks = c(2, 4,6,8,10,12), limits = c(0,11)) +
  geom_vline(xintercept = 77, colour = "black", linetype = "dashed", size = 1.5) +
  ggtitle("Alternate BEC calculations - approach 1 and approach 2") +
  # scale_x_continuous(breaks = c(50, 100, 150, 200, 250)) +
  scale_fill_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  # scico::scale_fill_scico_d(palette = "hawaii")+
  theme_bw(10) + theme(strip.background = element_rect(fill="gray85"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) 
FigS2   



pdf("FigS2.pdf", width = 8, height = 8, onefile = FALSE)
FigS2
dev.off()

tiff("FigS2.tiff", width = 8, height = 8, units = 'in', res = 300)
FigS2
dev.off()


# ### Some more checks
# 
# BEC_Check <- Bouts %>%
#   left_join(PopBEC_Summary)
# 
# 
# 
# for(i in unique(BEC_Check$ID)){
#   subID <- subset(BEC_Check, ID == i)
#   
#   gi <- ggplot() + 
#     geom_point(data = subID, aes(y = diff_secs, x = dur), colour = "black") +
#     geom_hline(data = subID, aes(yintercept = Approach1_ID, color = "Approach 1 ID"), size = 0.75) +
#     geom_hline(data = subID, aes(yintercept = Approach1_Mean, color = "Approach 1 Mean"), size = 0.75) +
#     # geom_hline(data = subID, aes(yintercept = Approach2, color = "Approach 2"), size = 0.75) +
#     geom_hline(data = subID, aes(yintercept = Approach3, color = "Approach 2"), size = 0.75) +
#     labs(x = "Dive duration (s)", y = "Surface interval (s)") +
#     facet_wrap(Year ~ ID) +
#     scale_y_continuous(limits = c(0, 500), breaks = c(0,50,100,150,200,250,300,350,400,450,500)) +
#     scale_color_manual(
#       name = "BEC approaches",
#       values = c("Approach 1 ID" = "#d95f02", 
#                  "Approach 1 Mean" = "red", 
#                  "Approach 2" = "purple", 
#                  "Approach 3" = "darkgreen", 
#                  "Approach 4" = "blue")) +
#     theme_bw(10) + 
#     theme(strip.background = element_rect(fill="gray85"),
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),  
#           panel.border = element_rect(colour = "black"),
#           legend.position = "right")
#   
#   gi1 <- ggplot() + 
#     geom_histogram(data = subID, aes(y = diff_secs), fill = "grey", color = "black") +
#     geom_hline(data = subID, aes(yintercept = Approach1_ID, color = "Approach 1 ID"), size = 0.75) +
#     geom_hline(data = subID, aes(yintercept = Approach1_Mean, color = "Approach 1 Mean"), size = 0.75) +
#     # geom_hline(data = subID, aes(yintercept = Approach2, color = "Approach 2"), size = 0.75) +
#     geom_hline(data = subID, aes(yintercept = Approach3, color = "Approach 2"), size = 0.75) +
#     labs(x = "Frequency", y = "Surface interval (s)") +
#     scale_y_continuous(limits = c(0, 500), breaks = c(0,50,100,150,200,250,300,350,400,450,500)) +
#     scale_color_manual(
#       name = "BEC approaches",
#       values = c("Approach 1 ID" = "#d95f02", 
#                  "Approach 1 Mean" = "red", 
#                  "Approach 2" = "purple", 
#                  "Approach 3" = "darkgreen", 
#                  "Approach 4" = "blue")) +
#     theme_bw(10) + 
#     theme(strip.background = element_rect(fill="gray85"),
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),  
#           panel.border = element_rect(colour = "black"),
#           legend.position = "right")
#   
#   # Combine plots and save
#   gi2 <- ggarrange(gi, gi1, nrow = 2, common.legend = TRUE, legend = "right")
#   ggsave(filename = sprintf('%s.png', i), plot = gi2, width = 8, height = 7)
# }



# write_rds(PopBEC_Summary, "PopBEC_Summary.rds")
# 
# write_rds(Bouts, "Bouts.rds")


###############################################################################################################################################
###############################################################################################################################################
# ASSIGN BOUTS
###############################################################################################################################################
###############################################################################################################################################

setwd("C:/Users/a46027/OneDrive - Havforskningsinstituttet/Documents/SAERI files/Javed/Uruguay")

# PopBEC_Summary <- read_rds("PopBEC_Summary.rds")
# 
# Bouts <- read_rds("Bouts.rds")

Ass_Bouts <- Bouts %>%
  left_join(PopBEC_Summary)

Ass_Bouts <- Bouts %>%
  group_by(ID, Trip) %>%
  mutate(Bout = ifelse(diff_secs < 90 | lag(diff_secs < 90),"TRUE", NA)) %>%
  drop_na() %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(group = cumsum(diff_secs >= 90) + 1) %>%
  mutate(randcount = 1)


Ass_Bouts$nestlon <- -58.933333
Ass_Bouts$nestlat <- -62.216667



###############################################################################################################################################
###############################################################################################################################################
# Bout summary calculations 
###############################################################################################################################################
###############################################################################################################################################

## Calculates some bout metrics
BoutID <- Ass_Bouts %>% 
  ungroup() %>%
  group_by(ID, Trip) %>%
  left_join(Dives1) %>%
  mutate(totaldives = NoALLDives) %>% ## Total bout-structured dives for an individual
  dplyr::select(-NoALLDives) %>%
  ungroup() %>%
  group_by(ID, Trip, group) %>%
  mutate(count=n()) %>%
  mutate(divesperbout = sum(randcount)) %>% ## Summaries count by bout group, trip and ID
  mutate(Distance=distHaversine(cbind(nestlon, nestlat),cbind(lon,lat))/1000) %>% ## Distance each dive from colony
  filter(divesperbout >= 3) %>%   # Retains only clusters of 3 or more dives
  mutate(Sumdur = sum(dur)) %>% # total bout duration in seconds
  mutate(boutduration_mins = sum(dur)/60, ## total bout duration time in minutes
         depthMean = mean(dep), ## Mean dive depth in a bout
         meanhunt = mean(hunt), # Mean hunts per dive in a bout
         meanduration = mean(dur), # Mean dive duration in a bout
         sumhunt = sum(hunt), ## Total hunts attempts per bout
         boutsuccess = sumhunt/boutduration_mins) %>% ### Hunts per minute of bouting
  ungroup() %>%
  group_by(ID, Trip) %>%
  arrange(ID, Trip, group) %>%
  mutate(TemporalBout = seq_along(group)) %>% 
  ungroup()



## Calculates summary metrics for each bout
TripID <- BoutID %>%
  group_by(ID, Trip, group) %>%
  arrange(ID, Trip, group, start) %>%
  filter(start == min(start)) %>%
  #   distinct(ID, Trip, group, .keep_all = TRUE) %>%
  left_join(sum_water) %>%
  ungroup() %>%
  group_by(ID, Trip) %>%
  mutate(allboutdives = sum(divesperbout)) %>% ## Total bout-structured dives for an individual
  mutate(boutprop = allboutdives/totaldives) %>% ### Proportion of bouts-structured dives in the foraging trip
  mutate(boutForageprop = allboutdives/AllHuntingDives) %>% ### Proportion of bouts-structured dives across foraging dives
  
  mutate(bouthuntprop = allboutdives/sum(sumhunt)) %>% 
  mutate(birdbouts = length(unique(group))) %>% ## Total number of bouts by a single individual 
  arrange(ID, Trip, start) %>%
  mutate(Boutdiff = start - lag(start),
         Bout_secs = as.numeric(Boutdiff, units = 'secs')) %>%
  mutate(Boutinterval_mins = Bout_secs - lag(Sumdur)) %>%
  mutate(Boutinterva_mins = Boutinterval_mins/60) %>%
  mutate(BoutInterDist=distHaversine(cbind(lon, lat),lag(cbind(lon,lat)))/1000)



ggplot() + geom_histogram(data = TripID, aes(x = boutprop, fill = Year), position = "dodge") 


Mean_boutForageprop <- TripID %>% 
  ungroup() %>% 
  dplyr::select(ID, boutForageprop, Year) %>% 
  distinct()

ggplot() + geom_histogram(data = Mean_boutForageprop, aes(x = boutForageprop, fill = Year), position = "dodge") 



FigS3 <-  ggplot() + geom_histogram(data = TripID, aes(x = Distance, fill = Year), colour = "black", position = "dodge") +
  # facet_wrap(~ Year) +
  scale_x_continuous(breaks = c(0,10,20,30,40,50, 60)) +
  scale_y_continuous(expand = c(0,0), breaks = c(25,50,75,100,125,150,175,200, 225), limits = c(0,225)) +
  # scico::scale_fill_scico_d(palette = "hawaii")+
  scale_fill_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  labs (y = "Frequency of bouts (#)", x = "Distance from penguin colony (km)") +
  theme_bw(10) + theme(strip.background = element_rect(fill="gray85"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) 
FigS3   


pdf("Fig S3.pdf", width = 8, height = 8, onefile = FALSE)
FigS3
dev.off()

tiff("Fig S3.tiff", width = 8, height = 8, units = "in", res = 300)
FigS3
dev.off()



Bout_Distance_Summary <- TripID %>%
  ungroup() %>%  
  mutate(Dist_cat = cut(Distance, breaks = c(0,10,20,30, 40,50, 60))) %>%
  group_by(Year, Dist_cat) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(Year) %>%
  mutate(Prop = n/sum(n))




###############################################################################################################################################
###############################################################################################################################################
# Figure 2 plot
###############################################################################################################################################
###############################################################################################################################################


AllForagingDives <- Bouts %>% 
  group_by(ID, Trip) %>%
  arrange(ID, Trip, start) %>%
  mutate(Bout = ifelse(diff_secs < 90 | lag(diff_secs < 90),"TRUE", NA)) %>%
  # drop_na() %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(group = cumsum(diff_secs >= 90) + 1) %>%
  mutate(randcount = 1)

AllForagingDives$Bout[is.na(AllForagingDives$Bout)] <- FALSE

AllForagingDives <- AllForagingDives %>%
  group_by(ID, Trip, group, Bout) %>%
  mutate(divesperbout = sum(randcount)) 

Fig2A <- AllForagingDives %>%
  ungroup() %>%
  arrange(ID, start) %>%
  filter(Trip == "A012_1920_1") %>%
  mutate(Bout = ifelse(divesperbout <= 3, FALSE, Bout)) %>%
  mutate(group = ifelse(Bout == FALSE, round(runif(n(), min = 0, max = 1), 10), group)) %>%
  group_by(group) %>%
  mutate(depthmean = mean(dep),
         sumhunt = sum(hunt)) %>%
  ungroup()


SurfaceTime <-   Fig2A %>%
  mutate(start = start + 0.1) %>%
  mutate(dep = 0)

Fig2A <- Fig2A %>% ## 
  full_join(SurfaceTime) %>%
  arrange(start) 



Fig2A_Plot <- ggplot() +
  geom_line(data = Fig2A,aes(x = start, y = dep, group = group, colour = Bout), linewidth = 0.2) +
  geom_hline(yintercept= 0) +
  # facet_wrap(~ ID, scales = "free") +
  scale_y_reverse(limits=c(136, 0), expand = c(0,0), breaks = c(0,25,50,75,100,125)) +
  scale_colour_manual(values = c("black", "orange"), labels = c("Non-bout", "Bout")) +
  scale_x_datetime(
    date_breaks = "1 hour",
    date_labels = "%H:%M") +
  labs(x = NULL, y = "Dive depth (m)", colour = "Dive\nstructure") +
  theme_bw(10) + theme(strip.background = element_rect(fill="gray85"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()) +
  # guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5, fill=NA, alpha = 1, stroke = 6))) +
  annotate("text",x=min(Fig2A$start + 2000),y=max(Fig2A$dep + 51),hjust=.2,label="BEC = 90 seconds", colour = "red", fontface = 2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(0, 0.1, -0.1, 0.1), "inches"))
# theme(legend.background = element_rect(fill = "white", colour = "black"),
#       legend.key = element_rect(fill = "white")) +
# guides(colour = guide_legend(override.aes = list(size = 5)))
Fig2A_Plot



unique(AllForagingDives$Trip)

Fig2B <- AllForagingDives %>%
  ungroup() %>%
  arrange(ID, start) %>%
  filter(Trip == "A022_2122_2") %>%
  mutate(Bout = ifelse(divesperbout <= 3, FALSE, Bout)) %>%
  mutate(group = ifelse(Bout == FALSE, round(runif(n(), min = 0, max = 10), 20), group)) %>%
  group_by(group) %>%
  mutate(depthmean = mean(dep),
         sumhunt = sum(hunt)) %>%
  ungroup()


SurfaceTimeB <-   Fig2B %>%
  mutate(start = start + 0.1) %>%
  mutate(dep = 0)

Fig2B <- Fig2B %>% ##
  full_join(SurfaceTimeB) %>%
  arrange(start) %>% 
  filter( start > as.POSIXct("2021-12-28 00:20:48") &  start < as.POSIXct("2021-12-28 11:00:00"))


Fig2B_Plot <- ggplot() +
  geom_line(data = Fig2B,aes(x = start, y = dep, group = group, colour = Bout), linewidth = 0.2) +
  geom_hline(yintercept= 0) +
  # facet_wrap(~ ID, scales = "free") +
  scale_y_reverse(limits=c(136, 0), expand = c(0,0), breaks = c(0,25,50,75,100,125)) +
  scale_colour_manual(values = c("black", "orange"), labels = c("Non-bout", "Bout")) +
  scale_x_datetime(
    date_breaks = "1 hour",
    date_labels = "%H:%M") +
  labs(x = "Time", y = "Dive depth (m)", colour = "Dive\nstructure") +
  theme_bw(10) + theme(strip.background = element_rect(fill="gray85"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()) +
  # guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5, fill=NA, alpha = 1, stroke = 6))) +
  annotate("text",x=min(Fig2B$start + 2000),y=max(Fig2B$dep -9),hjust=.2,label="BEC = 90 seconds", colour = "red", fontface = 2) +
  theme(plot.margin = unit(c(-1, 0.1, 0, 0.1), "inches")) 
# theme(legend.background = element_rect(fill = "white", colour = "black"),
#       legend.key = element_rect(fill = "white")) +
# guides(colour = guide_legend(override.aes = list(size = 5)))


Fig2B_Plot


Figure_2 <- ggarrange(Fig2A_Plot, Fig2B_Plot, nrow = 2, common.legend = TRUE, align = "h", labels = "AUTO", legend = "right")
Figure_2 <- annotate_figure(Figure_2, top = text_grob("Examples of penguin bout activity from (A) 2019/20 and (B) 2021/22",
                                                      color = "black", face = "bold", size = 11.5))
Figure_2


penguin_image <- cowplot::ggdraw() + 
  cowplot::draw_image("penguin.png", 
                      x = 1, y = 1, hjust = 1, vjust = 1, width = 0.5, height = 0.5)


Figure_2_with_image <- cowplot::ggdraw() +
  cowplot::draw_plot(Figure_2, 0, 0, 1, 1) +
  cowplot::draw_plot(penguin_image, x = 0.755, y = 0.74, width = 0.25, height = 0.25)
Figure_2_with_image


pdf("R1_Fig2.pdf", width = 9, height = 7, onefile = FALSE)
Figure_2_with_image
dev.off()


tiff("R1_Fig_2.tiff", width = 9, height = 7, units = 'in', res = 300)
Figure_2_with_image
dev.off()




###############################################################################################################################################
###############################################################################################################################################
# Figure 4 map
###############################################################################################################################################
###############################################################################################################################################


BoutDist <- BoutID %>%
  # filter(BoutOrder == "End") %>%
  ungroup() %>%
  group_by(ID, group) %>%
  slice(which.min(rowid)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(BoutDurQuant = case_when(boutduration_mins < quantile(boutduration_mins, 0.33) ~ "Low",
                                  boutduration_mins >= quantile(boutduration_mins, 0.33) & boutduration_mins <= quantile(boutduration_mins, 0.66) ~ "Medium",
                                  boutduration_mins > quantile(boutduration_mins, 0.66) ~ "High")) 

BoutDist$BoutDurQuant <- ordered(BoutDist$BoutDurQuant, levels = c("Low", "Medium", "High"))

ggplot() +geom_histogram(data = BoutDist, aes(x = BoutDurQuant), stat = "count")


coast_wgs<-read_sf("add_coastline_high_res_polygon_v7_5.shp")
st_crs(coast_wgs)  
st_geometry(coast_wgs)

AntarcticMap <- coast_wgs %>%
  st_transform(4326)


box = c(xmin = -59.35, ymin = -63, xmax = -57.6, ymax = -61.9)

cropped_sf <- st_crop(AntarcticMap, box)

cropped_sf <- cropped_sf %>%
  dplyr::select(-surface)

Antarcticamap <- ggplot() +
  geom_sf(data = cropped_sf, fill = "darkgrey")


Map2A<- Antarcticamap + 
  ggspatial:: geom_spatial_path(data = BoutDist, aes(x = lon, y = lat, group = ID), alpha = 0.5, colour = "black") +
  ggspatial:: geom_spatial_point(data = BoutDist, aes(x = lon, y = lat, colour = BoutDurQuant), alpha = 0.8, size = 0.5, shape = 1, stroke = 0.6) +
  facet_wrap(~ Year, ncol = 2) + 
  scale_colour_manual(values = c("#191900", "#F1C659", "#C34C4A")) +
  # scico::scale_fill_scico_d(palette = "lajolla", direction = 1) +
  theme_classic() + theme_bw(10) + theme(strip.background = element_rect(fill="gray85"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) +
  theme(legend.direction = "vertical") + labs(fill = "NASC\n(scaled)", colour = "Bout duration", size = "Bout duration", y = "Latitude", x = "Longitude") +
  scale_y_continuous(expand = c(0,0)) +   scale_x_continuous(expand = c(0,0)) +
  guides(colour = guide_legend(row = 1, override.aes = list(shape = 16, size = 6, fill=NA, alpha = 1, stoke = 5))) +
  guides(size = guide_legend(row = 1, override.aes = list(fill=NA, alpha = 1))) +
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.ontop = FALSE) +
  theme(legend.justification = c(0.95,0.1),
        legend.position = c(0.48,0.05),
        legend.background = element_rect(fill = "lightgrey", colour = "black"),
        legend.key = element_rect(fill = "lightgrey")) +
  coord_sf(xlim = c(-59.4, -57.6), ylim = c(-62.75, -61.9), expand = FALSE) +
  theme(plot.margin = unit(c(-0.5, 0.1, 0.1, 0.1), "inches")) 

Map2A


addSmallLegend <- function(myPlot, pointSize = 3, textSize = 7, spaceLegend = 0.5) {
  myPlot +
    guides(
      color = guide_legend(override.aes = list(size = pointSize, shape = 16, fill=NA, alpha = 1, stoke = 5))) +
    theme(legend.title = element_text(size = textSize),
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}


Map2A <- addSmallLegend(Map2A)
Map2A

BlankPlot <- ggplot() + theme_void()

Fig4A <- ggarrange(Map2A, BlankPlot, nrow = 2, heights = c(2, 1))
Fig4A



tiff("Fig4A.tiff", width = 5.5, height = 5.5, res = 300, units = "in")
Fig4A
dev.off()

pdf("Fig4A.pdf", width = 5.5, height = 5.5)
Fig4A
dev.off()




pal = c("#191900", "#F1C659", "#C34C4A")



Map4B1Data <- BoutDist %>%
  filter(Year == "2019/20")
Map4B1Data <- Map4B1Data %>%
  mutate(SuccessQuantNumeric = as.numeric(BoutDurQuant)) 


Map4B2Data <- BoutDist %>%
  filter(Year == "2021/22")
Map4B2Data <- Map4B2Data %>%
  mutate(SuccessQuantNumeric = as.numeric(BoutDurQuant))



# Create 3D  plot
Map4B1 <- plot_ly(
  data = Map4B1Data,
  x = ~lon,          
  y = ~lat,           
  z = ~depthMean,     
  color = ~BoutDurQuant, 
  size = ~SuccessQuantNumeric, 
  colors = pal,       
  type = "scatter3d",       
  mode = "markers",        
  marker = list(
    sizemode = "diameter",  
    sizeref = 10,            
    sizemin = 5,               
    # sizeref = c(10),         
    symbol = "circle",    
    line = list(width = 20, color = "black"),
    opacity = 10          
  ),
  showlegend = FALSE
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Longitude", autorange = FALSE, range = c(max(Map4B2Data$lon), min(Map4B2Data$lon))),    # Set x-axis limits
      yaxis = list(title = "Latitude", autorange = FALSE, range = c(max(Map4B2Data$lat), min(Map4B2Data$lat))),                        # Set y-axis limits
      zaxis = list(title = "Depth (m)", autorange = FALSE, range = c(130, 0))       # Set z-axis limits
    ),
    showlegend = FALSE 
  )


Map4B1




# Create  3D  plot
Map4B2 <- plot_ly(
  data = Map4B2Data,
  x = ~lon,          
  y = ~lat,          
  z = ~depthMean,     
  color = ~BoutDurQuant, 
  size = ~SuccessQuantNumeric, 
  colors = pal,       
  type = "scatter3d", 
  mode = "markers",   
  marker = list(
    sizemode = "diameter", 
    sizeref = 10,             
    sizemin = 5,               
    symbol = "circle",    
    line = list(width = 1, color = "black"), 
    opacity = 10          
  )
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Longitude", autorange = FALSE, range = c(max(Map4B2Data$lon), min(Map4B2Data$lon)),       
                   tickvals = c(-59, -58.5), 
                   ticktext = c("-59", "-58.5")), 
      yaxis = list(title = "Latitude", autorange = FALSE, range = c(max(Map4B2Data$lat), min(Map4B2Data$lat))),                        # Set y-axis limits
      zaxis = list(title = "Depth (m)", autorange = FALSE, range = c(130, 0))      
    ),
    showlegend = FALSE 
  )


Map4B2



# Export as an SVG file
orca(Map4B1, "Fig4B_First.svg")

orca(Map4B2, "Fig4B_Second.svg")




###############################################################################################################################################
###############################################################################################################################################
# Supplementary Figure 1
###############################################################################################################################################
###############################################################################################################################################


All_Acoustic <- read_rds("All_Acoustic.rds") ## Need to read All_Acoustic file from the krill script

S1 <- BoutDist %>%
  filter(Year == "2019/20")

S1$mo<-month(S1$start)
S1$day<-day(S1$start)
S1$yr<-2023
S1$yr[S1$mo<=1]<-2024
S1$date<-ymd(paste(S1$yr,"-",S1$mo,"-",S1$day))


S2 <- BoutDist %>%
  filter(Year == "2021/22")

S2$mo<-month(S2$start)
S2$day<-day(S2$start)
S2$yr<-2023
S2$yr[S2$mo<=1]<-2024
S2$date<-ymd(paste(S2$yr,"-",S2$mo,"-",S2$day))


PlotYear <- S1 %>%
  full_join(S2)


PlotYear <- PlotYear %>%
  mutate(
    ShortDate = format(date, "%d-%b"), 
    WeekBlock = paste("Week", isoweek(date)), 
    Year1 = year(date)
  )

PlotYear$WeekBlock <- as.factor(PlotYear$WeekBlock)

unique(PlotYear$Year_Guard)



All_Acoustic$start <- lubridate::ymd(All_Acoustic$Date_S)

A1 <- All_Acoustic %>%
  filter(Year == "2019/20")


A1$mo<-month(A1$start)
A1$day<-day(A1$start)
A1$yr<-2023
A1$yr[A1$mo<=1]<-2024
A1$date<-ymd(paste(A1$yr,"-",A1$mo,"-",A1$day))


A2 <- All_Acoustic %>%
  filter(Year == "2021/22")

A2$mo<-month(A2$start)
A2$day<-day(A2$start)
A2$yr<-2023
A2$yr[A2$mo<=1]<-2024
A2$date<-ymd(paste(A2$yr,"-",A2$mo,"-",A2$day))


AcousticPlotYear <- A1 %>%
  full_join(A2)


AcousticPlotYear <- AcousticPlotYear %>%
  mutate(
    ShortDate = format(date, "%d-%b"), 
    WeekBlock = paste("Week", isoweek(date)), 
    Year1 = year(date)
  )

AcousticPlotYear$WeekBlock <- as.factor(AcousticPlotYear$WeekBlock)



FigS1 <- Antarcticamap +
  ggspatial::geom_spatial_point(data = PlotYear, aes(x = lon, y = lat, color = Year), alpha = 0.5, size = 0.5) +
  ggspatial::geom_spatial_point(data = AcousticPlotYear, aes(x = Lon, y = Lat, fill = Year), alpha = 0.5, size = 0.4) +
  
  facet_grid(Year~factor(WeekBlock, levels=c("Week 49", "Week 50", "Week 51","Week 52", "Week 1"))) +
  scale_colour_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  scale_fill_viridis_d() +
  # scale_colour_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  theme_bw() + theme(strip.background = element_rect(fill="gray85"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) +
  
  labs(
    subtitle = "Spatial and temporal coverage of penguin tracking and krill acoustic data",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(legend.position = "none")  +
  theme(
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8))
# theme(legend.justification = c(0.95,0.1),
#           legend.position = c(0.9,0.7),
#           legend.background = element_rect(fill = "lightgrey", colour = "black"),
#           legend.key = element_rect(fill = "lightgrey")) 



FigS1


tiff("FigS1.tiff", width = 9, height = 4, res = 300, units = "in")
FigS1
dev.off()

pdf("FigS1.pdf", width = 9, height = 4, onefile = FALSE)
FigS1
dev.off()




###############################################################################################################################################
###############################################################################################################################################
# Save data outputs
###############################################################################################################################################
###############################################################################################################################################

# write_rds(BoutID, "BoutID_POP.rds")
# write_rds(TripID, "TripID_POP.rds")


###############################################################################################################################################
##############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

# PART 3 - PENGUIN BOUT MODELLING

###############################################################################################################################################
##############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################


BoutMOD <- TripID

BoutMOD$Year <- as.character(BoutMOD$Year)


BoutMOD <- BoutMOD %>% 
  ungroup() %>%
  dplyr::arrange(ID, Trip, start) %>%
  mutate(timeLub = lubridate::decimal_date(start)) %>%
  mutate(times = glmmTMB::numFactor(timeLub)) 


annotations <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  annotateText = c("*"),
  hjustvar = c(0) ,
  vjustvar = c(1.0))


annotationsdepth <- data.frame(
  xpos = c(-Inf),
  ypos =  c(-Inf),
  annotateText = c("*"),
  hjustvar = c(0) ,
  vjustvar = c(1.0))


###############################################################################################################################################
### Model 1 - Inter-bout Distance model
###############################################################################################################################################


ggplot() + geom_histogram(data = BoutMOD, aes(x = BoutInterDist )) +
  facet_wrap(~ Year)

InterboutModelDIST_Mod <- BoutMOD %>% 
  filter(BoutInterDist < 20) 


ggplot() + geom_histogram(data = InterboutModelDIST_Mod, aes(x = , BoutInterDist, fill = Year), colour = "black", position = "dodge") +
  labs (x = "Inter-bout period (mins)", y = "Frequency #") +
  theme_bw() + scale_y_continuous(expand = c(0,0)) 


InterboutDIST_Mod <-  glmmTMB::glmmTMB(BoutInterDist ~  Year + ou(times+0|ID), data = InterboutModelDIST_Mod, family = Gamma()) 


summary(InterboutDIST_Mod)
drop1(InterboutDIST_Mod, test = "Chisq") 
check_model(InterboutDIST_Mod) 
check_autocorrelation(InterboutDIST_Mod) 
simulationOutput <- DHARMa:: simulateResiduals(fittedModel = InterboutDIST_Mod)
plot(simulationOutput)


x <- sjPlot::get_model_data(InterboutDIST_Mod, type = "pred", terms = "Year")

x <- x %>%
  rename(Year = x,
         fit = predicted) %>%  
  mutate(Year = factor(Year, levels = c(1, 2), labels = c("2019/20", "2021/22"))) %>% 
  drop_na()


ModelPlot1 <- ggplot(data = x, aes(x = Year, y = fit)) +
  geom_point(aes(colour = Year, fill = Year), size = 3, position = position_dodge(0.4)) +
  labs(y = "Interbout distance (km)", x = "Year", fill = "Year", colour = "Year") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) + theme(legend.text = element_text(face = "italic")) +  theme(legend.title = element_blank(), legend.position = "right") +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour = Year), width=0.4, size = 1, position = position_dodge(0.4)) +
  scale_fill_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +  
  scale_colour_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  # theme(legend.position = "none") +
  theme(plot.margin = unit(c(1,1,1,1), 'lines')) +
  geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
                                    vjust=vjustvar,label=annotateText), size = 14) +
  theme(
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    # panel.border = element_rect(colour = "black"),
    legend.text = element_text(face = "plain"),  
    legend.title = element_text(face = "plain"),  
    # legend.title.align = 0.5,  
    legend.position = "right",
    legend.direction = "vertical",
    # legend.background = element_rect(fill = "white", colour = "black"),
    legend.key = element_rect(fill = "white"),
    plot.margin = unit(c(1, 1, 1, 1), 'lines')
  )
ModelPlot1




###############################################################################################################################################
### Model 2 - Inter-bout duration model
###############################################################################################################################################


InterboutModel_Mod <- BoutMOD %>%
  filter(Boutinterva_mins < 120) 



ggplot() + geom_histogram(data = InterboutModel_Mod, aes(x = , Boutinterva_mins, fill = Year), colour = "black", position = "dodge") +
  labs (x = "Interbout period (mins)", y = "Frequency #") +
  theme_bw() + scale_y_continuous(expand = c(0,0)) 



Interbout_Mod <-  glmmTMB::glmmTMB(Boutinterva_mins ~  Year + (1|ID), data = InterboutModel_Mod, family = Gamma()) 

summary(Interbout_Mod)
drop1(Interbout_Mod, test = "Chisq") 
check_model(Interbout_Mod) 
check_autocorrelation(Interbout_Mod) 
simulationOutput <- DHARMa:: simulateResiduals(fittedModel = Interbout_Mod)
plot(simulationOutput)


x <- sjPlot::get_model_data(Interbout_Mod, type = "pred", terms = "Year")

x <- x %>%
  rename(Year = x,
         fit = predicted) %>%  
  mutate(Year = factor(Year, levels = c(1, 2), labels = c("2019/20", "2021/22"))) %>% 
  drop_na()


ModelPlot2 <- ggplot(data = x, aes(x = Year, y = fit)) +
  geom_point(aes(colour = Year, fill = Year), size = 3, position = position_dodge(0.4)) +
  labs(y= "Interbout duration (min)", x= "Year")  + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) + theme(legend.text = element_text(face = "italic")) +  theme(legend.title = element_blank(), legend.position = "right") +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour = Year), width=0.4, size = 1, position = position_dodge(0.4)) +
  scale_fill_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +  
  scale_colour_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  theme(legend.position = "none") +
  geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
                                    vjust=vjustvar,label=annotateText), size = 14) + 
  theme(
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    # panel.border = element_rect(colour = "black"),
    legend.text = element_text(face = "plain"),  
    legend.title = element_text(face = "plain"),  
    legend.position = "right",
    legend.direction = "vertical",
    legend.key = element_rect(fill = "white"),
    plot.margin = unit(c(1, 1, 1, 1), 'lines')
  )
ModelPlot2



###############################################################################################################################################
### Model 3 - Bout depth Model
###############################################################################################################################################

ggplot() + geom_histogram(data = BoutMOD, aes(x = depthMean)) +
  facet_wrap(~ Year)


Depth_Mod <-  glmmTMB::glmmTMB(depthMean ~  Year + ou(times+0|ID), data = BoutMOD, family = gaussian()) 

summary(Depth_Mod)
drop1(Depth_Mod, test = "Chisq") # perform likelihood ratio tests
check_model(Depth_Mod) ## Check Assumptions
check_autocorrelation(Depth_Mod) ### 
simulationOutput <- DHARMa:: simulateResiduals(fittedModel = Depth_Mod)
plot(simulationOutput)



x <- sjPlot::get_model_data(Depth_Mod, type = "pred", terms = "Year")

x <- x %>%
  rename(Year = x,
         fit = predicted) %>%  
  mutate(Year = factor(Year, levels = c(1, 2), labels = c("2019/20", "2021/22"))) %>% 
  drop_na()

ModelPlot3 <- ggplot(data = x, aes(x = Year, y = fit)) +
  geom_point(aes(colour = Year, fill = Year), size = 3, position = position_dodge(0.4)) +
  labs(y= "Mean dive depth per bout (m)", x= "Year")  + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) + theme(legend.text = element_text(face = "italic")) +  theme(legend.title = element_blank(), legend.position = "right") +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour = Year), width=0.4, size = 1, position = position_dodge(0.4)) +
  scale_fill_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +  
  scale_colour_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  theme(legend.position = "none") +
  scale_y_reverse() +
  geom_text(data = annotationsdepth, aes(x=xpos,y=ypos,hjust=hjustvar,
                                         vjust=vjustvar,label=annotateText), size = 14)
ModelPlot3



###############################################################################################################################################
### Model 4 - Bout Duration Model
###############################################################################################################################################


ggplot() + geom_histogram(data = BoutMOD, aes(x = boutduration_mins)) +
  facet_wrap(~ Year)

Duration_Mod <-  glmmTMB::glmmTMB(boutduration_mins ~  Year + (1|ID), data = BoutMOD, family = Gamma()) 


summary(Duration_Mod)
drop1(Duration_Mod, test = "Chisq") 
check_model(Duration_Mod) 
check_autocorrelation(Duration_Mod) 
simulationOutput <- DHARMa:: simulateResiduals(fittedModel = Duration_Mod)
plot(simulationOutput)


x <- sjPlot::get_model_data(Duration_Mod, type = "pred", terms = "Year")

x <- x %>%
  rename(Year = x,
         fit = predicted) %>%  
  mutate(Year = factor(Year, levels = c(1, 2), labels = c("2019/20", "2021/22"))) %>% 
  drop_na()


ModelPlot4 <- ggplot(data = x, aes(x = Year, y = fit)) +
  geom_point(aes(colour = Year, fill = Year), size = 3, position = position_dodge(0.4)) +
  labs(y= "Bout duration (min)", x= "Year")  + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) + theme(legend.text = element_text(face = "italic")) +  theme(legend.title = element_blank(), legend.position = "right") +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour = Year), width=0.4, size = 1, position = position_dodge(0.4)) +
  scale_fill_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +  
  scale_colour_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  theme(legend.position = "none") +
  geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
                                    vjust=vjustvar,label=annotateText), size = 14)

ModelPlot4




###############################################################################################################################################
### Model 5 - PCA  Rate
###############################################################################################################################################

ggplot() + geom_histogram(data = BoutMOD, aes(x = boutsuccess)) +
  facet_wrap(~ Year)


PCAFreq_Mod <-  glmmTMB::glmmTMB(boutsuccess ~  Year + ou(times + 0|ID), data = BoutMOD, family = gaussian()) 

summary(PCAFreq_Mod)
drop1(PCAFreq_Mod, test = "Chisq") 
check_model(PCAFreq_Mod) 
check_autocorrelation(PCAFreq_Mod) 
simulationOutput <- DHARMa:: simulateResiduals(fittedModel = PCAFreq_Mod)
plot(simulationOutput)


x <- sjPlot::get_model_data(PCAFreq_Mod, type = "pred", terms = "Year")

x <- x %>%
  rename(Year = x,
         fit = predicted) %>%  
  mutate(Year = factor(Year, levels = c(1, 2), labels = c("2019/20", "2021/22"))) %>% 
  drop_na()


ModelPlot5 <- ggplot(data = x, aes(x = Year, y = fit)) +
  geom_point(aes(colour = Year, fill = Year), size = 3, position = position_dodge(0.4)) +
  labs(y= expression ("PCA rate (PCAs/min)"), x= "Year")   + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) + theme(legend.text = element_text(face = "italic")) +  theme(legend.title = element_blank(), legend.position = "right") +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour = Year), width=0.4, size = 1, position = position_dodge(0.4)) +
  scale_fill_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +  
  scale_colour_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  theme(legend.position = "none") 
# geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
#                                   vjust=vjustvar,label=annotateText), size = 14)

ModelPlot5



###############################################################################################################################################
### Model 6 - Inter-bout duration as a function of distance from colony
###############################################################################################################################################


GAMData <- InterboutModel_Mod %>%
  filter(Distance < 31) ##Max distance in 2019/20

GAMData$Year <- as.factor(GAMData$Year)

ggplot(data = GAMData, aes(x = Distance, y = Boutinterva_mins, colour = Year)) +
  geom_point() +
  geom_smooth(method = "gam")


RevisedInterBout <-  mgcv::gam(Boutinterva_mins ~  s(Distance, k = 10, bs = "tp", by = Year) + s(ID, bs = "re"), data = GAMData, method = "REML", family = Gamma(), select = FALSE)
gam.check(RevisedInterBout)
draw(RevisedInterBout, residuals = TRUE)
summary.gam(RevisedInterBout)


smooth_plot_data <- smooth_estimates(RevisedInterBout) %>%
  filter(.smooth == "s(Distance):Year2019/20") %>%
  add_confint()
smooth_plot_data1 <- smooth_estimates(RevisedInterBout) %>%
  filter(.smooth == "s(Distance):Year2021/22") %>%
  add_confint()

joinedsmooths <- smooth_plot_data %>%
  full_join(smooth_plot_data1)


ModelPlot6 <- joinedsmooths |>
  ggplot(aes(x = Distance, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .smooth),
              alpha = 0.4) +
  geom_line(aes(colour = .smooth)) + 
  scale_fill_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +   # <-- change fill scale
  scale_colour_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  # scale_x_continuous(limits = c(0,50), breaks = seq(0,50, by = 5)) +
  # geom_rug(data = GAMData,                  # <-- rug
  #          mapping = aes(x = Distance, y = NULL),
  #          sides = "b", alpha = 0.4) +  
  labs(title = NULL, y = "Interbout duration (model estimate)",
       colour = "Year", x = "Distance from penguin colony (km)", fill = "Year") + 
  theme_bw(10) + theme(strip.background = element_rect(fill="gray85"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black"))+ 
  # theme_minimal(base_size = 12) +
  theme(legend.position = "right") +
  # legend.direction = "vertical",
  #     legend.justification = c(0.95,0.1),
  #     legend.position = c(0.5,0.5),
  #     legend.background = element_rect(fill = "white", colour = "black"),
  #     legend.key = element_rect(fill = "white")) +
  theme(plot.margin = unit(c(1,1,1,1), 'lines')) +
  geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
                                    vjust=vjustvar,label=annotateText), size = 14) +
  theme(legend.direction = "vertical",
        # legend.justification = c(0.95,0.1),
        # legend.position = c(0.5,0.5),
        legend.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_rect(fill = "white")) +
  theme(plot.margin = unit(c(1,1,1,1), 'lines')) 

ModelPlot6




###############################################################################################################################################
### Extract legend just for plotting purposes
###############################################################################################################################################

LEGEND <- joinedsmooths |>
  ggplot(aes(x = Distance, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .smooth),
              alpha = 0.4) +
  geom_line(aes(colour = .smooth)) + 
  scale_fill_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +   # <-- change fill scale
  scale_colour_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  # scale_x_continuous(limits = c(0,50), breaks = seq(0,50, by = 5)) +
  # geom_rug(data = GAMData,                  # <-- rug
  #          mapping = aes(x = Distance, y = NULL),
  #          sides = "b", alpha = 0.4) +  
  labs(title = NULL, y = "Inter-bout duration (model estimate)",
       colour = "Year", x = "Distance from penguin colony (km)", fill = "Year") + 
  theme_bw(10) + theme(strip.background = element_rect(fill="gray85"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black"))+ 
  # theme_minimal(base_size = 12) +
  theme(legend.direction = "vertical",
        legend.justification = c(0.95,0.1),
        legend.position = c(0.5,0.5),
        legend.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_rect(fill = "white")) +
  theme(plot.margin = unit(c(1,1,1,1), 'lines')) +
  geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
                                    vjust=vjustvar,label=annotateText), size = 14)

LEGEND

leg <- get_legend(LEGEND)

as_ggplot(leg)


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# Figure 5 - combined model plotting

BlankPlot <- ggplot() + theme_void()

Figure5 <- ggarrange(ModelPlot2,ModelPlot1,ModelPlot3,ModelPlot4,ModelPlot5,ModelPlot6, ncol = 2, nrow = 3, common.legend = TRUE, legend = "right", align = "hv")
Figure5


pdf("Fig5.pdf", width = 7, height = 7, onefile = FALSE)
Figure5
dev.off()


tiff("Fig5.tiff", width = 7, height = 7,  res = 300, units = "in")
Figure5
dev.off()



#############################################################################################################################################
###############################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## End of Script

