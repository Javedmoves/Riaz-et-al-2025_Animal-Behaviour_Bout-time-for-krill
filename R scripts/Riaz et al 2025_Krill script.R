# Manuscript - Bout time for krill - contrasting Adelie penguin foraging behaviour during years of high and low krill availability
# Author - Javed Riaz
# Version: Final: Krill acoustic script

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
##############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

# PART 1 - KRILL DATA WRANGLING AND MODELLING

###############################################################################################################################################
##############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

# KRILL ACOUSTIC DATA

## Massage and basic spatial summaries

###############################################################################################################################################
###############################################################################################################################################

Krill2019 <- read_csv("12-2019.csv")
Krill2019$Season <- "2019/20"


Check <- Krill2019 %>%
  filter(Date_S == "20191229") %>%
  filter(VL_start < 268)

Krill2019 <- Krill2019 %>%
  filter(!Date_S == "20191229") %>%
  full_join(Check) %>%
  filter(VL_start > 160)


Krill2022 <- read_csv("01-2022.csv")
Krill2022$Season <- "2021/22"

CheckK <- Krill2022 %>%
  filter(Date_S == "20220106") %>%
  filter(VL_start > 88.9) 

Krill2022 <- Krill2022 %>%
  filter(!Date_S == "20220106") %>%
  full_join(CheckK)



All_Acoustic <- Krill2019 %>%
  # full_join(Krill2020) %>%
  full_join(Krill2022) %>%
  mutate(PRC_NASC = ifelse(PRC_NASC <= 0, 0.00000000001, PRC_NASC)) %>%
  mutate(No_Swarms = 1) %>%
  filter(Depth_mean > 0.1 )


All_Acoustic$Lon <- All_Acoustic$Lon_M
All_Acoustic$Lat <- All_Acoustic$Lat_M


ggplot() + geom_point(data = All_Acoustic, aes(x = Lon_M, y = Lat_M, colour = Season)) + 
  theme(legend.position = "right")


All_Acoustic$nestlon <- -58.933333
All_Acoustic$nestlat <- -62.216667

All_Acoustic <- All_Acoustic %>%
  mutate(Distance_Col=distHaversine(cbind(nestlon, nestlat),cbind(Lon,Lat))/1000)


Acoustic_SF <- st_as_sf(All_Acoustic, coords=c("Lon","Lat"), crs=4326)


transect_data_UTM <- st_transform(Acoustic_SF, crs = 32633)  

line_by_season <- transect_data_UTM %>%
  arrange(Season, Date_M, Time_M) %>%  
  group_by(Season) %>%
  summarise(geometry = st_cast(st_combine(geometry), "LINESTRING"))

line_by_season <- line_by_season %>%
  mutate(total_length_m = as.numeric(st_length(geometry))/1000)

hexgrid <- st_make_grid(Acoustic_SF, cellsize = c(0.05, 0.05), square = FALSE) %>% 
  st_transform(crs=4326) %>%
  st_as_sf() %>%
  mutate(GridID = c(1:NROW(.)))


ggplot()+
  geom_sf(data=hexgrid)+
  geom_sf(data=Acoustic_SF, aes(colour = Season))


joined <- st_join(hexgrid, Acoustic_SF) %>%
  group_by(Season, GridID) %>%
  summarise(DistanceColony = mean(Distance_Col),
            NASC_sum = sum(PRC_NASC),
            NASC_mean=exp(mean(log(PRC_NASC))), 
            NASC_sd=exp(sd(log(PRC_NASC))),
            No_Swarms = sum(No_Swarms),
            Depth_Mean=exp(mean(log(Depth_mean))), 
            Depth_sd=exp(sd(log(Depth_mean)))) %>%
  drop_na(NASC_sum)

joined$NASC_Sum_Standard <- datawizard::normalize(joined$NASC_sum)
joined$NASC_Mean_Standard <- datawizard::normalize(joined$NASC_mean)
joined$No_Swarms_Standard <- datawizard::normalize(joined$No_Swarms)
joined$Year <- joined$Season


All_Acoustic$Year <- All_Acoustic$Season


Sum_NASC <- ggplot() + geom_sf(data=joined, aes(fill = NASC_sum)) +
  geom_point(data = All_Acoustic, aes(x = Lon_M, y = Lat_M), size = 0.5) +
  facet_wrap(~Season, nrow = 1) +
  scale_fill_viridis()
Sum_NASC


###############################################################################################################################################
###############################################################################################################################################

### Spatial summaries of the krill data

options("scipen"=-100, "digits"=4)
options("scipen"=100, "digits"=4)


KrillSummaries<- Acoustic_SF %>%
  filter(PRC_NASC > 0.000001) %>%
  ungroup() %>%
  # mutate(Distance_Col = ifelse(Distance_Col >= 50, 50, Distance_Col)) %>%
  # mutate(Colony_Distance = cut(Distance_Col, breaks = c(0,10,20,30, 40, 50))) %>%
  group_by(Season) %>%
  summarise(
    NASC_Sum = sum(PRC_NASC),
    
    NASC_Mean = exp(mean(log(PRC_NASC[PRC_NASC > 0]))), 
    log_NASC_Mean = mean(log(PRC_NASC[PRC_NASC > 0])),
    log_NASC_SE = sd(log(PRC_NASC[PRC_NASC > 0])) / sqrt(sum(PRC_NASC > 0)),
    NASC_Lower_CI = exp(log_NASC_Mean - qnorm(0.975) * log_NASC_SE),
    NASC_Upper_CI = exp(log_NASC_Mean + qnorm(0.975) * log_NASC_SE),
    
    Depth_Mean = exp(mean(log(Depth_mean[Depth_mean > 0]))),
    log_Depth_Mean = mean(log(Depth_mean[Depth_mean > 0])),
    log_Depth_SE = sd(log(Depth_mean[Depth_mean > 0])) / sqrt(sum(Depth_mean > 0)),
    Depth_Lower_CI = exp(log_Depth_Mean - qnorm(0.975) * log_Depth_SE),
    Depth_Upper_CI = exp(log_Depth_Mean + qnorm(0.975) * log_Depth_SE)
  )




KrillSummaries <- as.data.frame(KrillSummaries)

KrillSummaries <- KrillSummaries %>%
  dplyr::select(- geometry)

KrillSummaries_Distance1 <- KrillSummaries %>%
  dplyr::filter(Season == "2019/20")
KrillSummaries_Distance2 <- KrillSummaries %>%
  dplyr::filter(Season == "2021/22")




# Binning and summarising
KrillSummaries_Distance <- All_Acoustic %>%
  filter(PRC_NASC > 0.0001) %>%
  mutate(Distance_Bin = cut(Distance_Col, 
                            breaks = seq(0, max(Distance_Col, na.rm = TRUE), by = 1),
                            labels = FALSE, 
                            include.lowest = TRUE)) %>%
  group_by(Year, Distance_Bin) %>%
  summarise(Total_NASC = sum(PRC_NASC, na.rm = TRUE), .groups = "drop")

KrillSummaries_Distance <- KrillSummaries_Distance %>%
  mutate(Distance_Bin = as.numeric(Distance_Bin))

KrillSummaries_Distance <- KrillSummaries_Distance %>%
  complete(Year, Distance_Bin = seq(min(Distance_Bin, na.rm = TRUE), max(Distance_Bin, na.rm = TRUE), by = 1), 
           fill = list(Total_NASC = 0))



write_rds(All_Acoustic, "All_Acoustic.rds")

###############################################################################################################################################
###############################################################################################################################################
##Krill Plots

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


InsetMap <- basemap(limits = c(-80, -52, -79, -51), bathymetry = FALSE , rotate = TRUE, grid.col = NA) +
  theme(legend.position = "right")  +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggspatial::geom_spatial_rect(aes(xmin = -60.5, ymin = -63, xmax = -57, ymax = -61.5),colour = "red", fill = "transparent", linewidth = 0.75) +
  theme_bw() +
  theme_void() +  # Remove all extra spacing
  theme(
    panel.border = element_rect(fill = NA, colour = "black", size = 1.1),  # Thicker black border
    plot.margin = margin(0, 0, 0, 0)  # Trim excess space
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white")  # Ensure solid white background
  ) 
InsetMap

pdf("InsetMap.pdf", width = 5, height = 5.5)
InsetMap
dev.off()


Map1<- Antarcticamap + 
  geom_sf(data=joined, aes(fill = NASC_sum)) +
  facet_wrap(~ Year) +
  scale_fill_viridis() +
  # geom_image(data = All_Acoustic, aes(x = -58.93, y = -62.22, image = image_path), size = 0.5, colour = "red") + # Adjust size as needed
  # geom_point(data = All_Acoustic, aes(x = -58.931, y = -62.22), colour = "red", shape = 17) +
  geom_point(data = All_Acoustic, aes(x = Lon_M, y = Lat_M, colour = "black"), size = 0.5) +
  facet_wrap(~ Year, ncol = 2) + 
  scale_colour_manual(values = c("#191900", "#F1C659", "#C34C4A"), labels = c("")) +
  # scico::scale_fill_scico_d(palette = "lajolla", direction = 1) +
  theme_bw(10) + theme(strip.background = element_rect(fill="gray85"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) +
  theme(legend.direction = "vertical") + labs(fill = "NASC (sum)", colour = "Acoustic\ntransect", size = "Hunt activity", y = "Latitude", x = "Longitude") +
  scale_y_continuous(expand = c(0,0)) +   scale_x_continuous(expand = c(0,0)) +
  guides(colour = guide_legend(row = 1, override.aes = list(size = 5, fill=NA, alpha = 1, stoke = 3))) +
  guides(size = guide_legend(row = 1, override.aes = list(fill=NA, alpha = 1))) +
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.ontop = FALSE) +
  theme(
    legend.justification = c(0.95, 0.05),
    legend.position = c(0.96, 0.1),
    legend.background = element_rect(fill = "lightgrey", colour = "black"),
    legend.key = element_rect(fill = "lightgrey"),
    legend.text = element_text(size = rel(0.9)),  # Reduce legend text size
    legend.title = element_text(size = rel(0.9)), # Reduce legend title size
    legend.key.size = unit(0.8, "lines")          # Reduce legend key size by 25%
  ) +
  # theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
  #       axis.title.x=element_blank()) +
  coord_sf(xlim = c(-59.4, -57.6), ylim = c(-62.75, -61.9), expand = FALSE) +
  theme(plot.margin = unit(c(0.1,0,0.2,0), 'lines')) 
# theme(plot.margin = unit(c(0, 0, 0, 0), 
#                          "cm"))

Map1


InsetMap_grob <- ggplotGrob(InsetMap)


Map1_Inset <- Map1 + annotation_custom2(
  InsetMap_grob, 
  data = data.frame(Year= "2019/20"),
  
  xmin = -58.2, xmax = Inf, 
  ymin = -62.3, ymax = Inf)
Map1_Inset

BlankPlot <- ggplot() + theme_void()

Fig1A <- ggarrange(Map1_Inset, BlankPlot, nrow = 2)
Fig1A



tiff("Fig1A.tiff", width = 5, height = 5.5, res = 300, units = "in")
Fig1A
dev.off()

pdf("Fig1A.pdf", width = 5, height = 5.5)
Fig1A
dev.off()


# ###############################################################################################################################################
# ###############################################################################################################################################
# # 3D mapping of the krill prey field

All_Acoustic$Depth_round <- plyr::round_any(All_Acoustic$Depth_mean, 20, floor)
All_Acoustic$Lat_round <- plyr::round_any(All_Acoustic$Lat_M, 0.05, floor)
All_Acoustic$Lon_round <- plyr::round_any(All_Acoustic$Lon_M, 0.05, floor)



PlotData <- All_Acoustic %>%
  filter(PRC_NASC > 0.01) %>%
  group_by(Year, Lon_round, Lat_round, Depth_round) %>%
  dplyr::summarise(PRC_NASC = sum(PRC_NASC)) %>%
  ungroup() %>%
  mutate(NASC_Quant = case_when(PRC_NASC < quantile(PRC_NASC, 0.50) ~ "Low",
                                PRC_NASC >= quantile(PRC_NASC, 0.50) & PRC_NASC <= quantile(PRC_NASC, 0.85) ~ "Medium",
                                PRC_NASC > quantile(PRC_NASC, 0.85) ~ "High"))


PlotData$SuccessQuant <- ordered(PlotData$NASC_Quant, levels = c("Low", "Medium", "High"))


PlotData <- PlotData %>%
  mutate(NASC_QuantNumeric = as.numeric(SuccessQuant))

FirstYearAccoustic <- PlotData %>%
  filter(Year == "2019/20")

SecondYearAccoustic <- PlotData %>%
  filter(Year == "2021/22") 


FirstYearAccoustic_PlotJig <- FirstYearAccoustic
FirstYearAccoustic_PlotJig$NASC_Quant <- "NA"
FirstYearAccoustic_PlotJig$SuccessQuant <- "NA"
FirstYearAccoustic_PlotJig$NASC_QuantNumeric <-  as.numeric(0)
FirstYearAccoustic_PlotJig$Depth_round <- as.numeric(199.99)
SecondYearAccoustic_PlotJig <- SecondYearAccoustic
SecondYearAccoustic_PlotJig$NASC_Quant <- "NA"
SecondYearAccoustic_PlotJig$SuccessQuant <- "NA"
SecondYearAccoustic_PlotJig$NASC_QuantNumeric <- as.numeric(0)
SecondYearAccoustic_PlotJig$Depth_round <- as.numeric(199.99)

FirstYearAccousticTest <- FirstYearAccoustic %>%
  full_join(SecondYearAccoustic_PlotJig)

SecondYearAccousticTest <- SecondYearAccoustic %>%
  full_join(FirstYearAccoustic_PlotJig)


KrillPal = c("#C34C4A","#191900", "#F1C659", "transparent")


Krill3DA <- plot_ly(
  data = FirstYearAccousticTest,
  x = ~Lon_round,          
  y = ~Lat_round,           
  z = ~Depth_round,        
  color = ~SuccessQuant,    
  size = ~NASC_QuantNumeric, 
  colors = KrillPal,       
  type = "scatter3d",       
  mode = "markers",        
  marker = list(
    sizemode = "diameter",  
    sizeref = c(10),         
    symbol = "circle",      
    opacity = 10           
  ),
  showlegend = FALSE
) 



Krill3DA <- Krill3DA %>%
  add_trace(
    x = c(-58.93), 
    y = c(-62.22), 
    z = c(0),
    type = "scatter3d", 
    mode = "markers",
    marker = list(
      size = 1,            
      symbol = "triangle", 
      color = "darkblue"   
    ),
    showlegend = FALSE     
  )

Krill3DA <- Krill3DA %>%
  layout(
    scene = list(
      xaxis = list(
        title = "Longitude",
        autorange = "reversed",
        range = c(-59.2, -58.6) 
      ),
      yaxis = list(
        title = "Latitude",
        autorange = "reversed",
        range = c(-62.75, -61.9) 
      ),
      zaxis = list(
        title = "Swarm depth (m)",
        autorange = "reversed",
        range = c(200, 0)        
      )
    )
  )

Krill3DA



Krill3DB <- plot_ly(
  data = SecondYearAccousticTest,
  x = ~Lon_round,         
  y = ~Lat_round,           
  z = ~Depth_round,         
  color = ~SuccessQuant,    
  size = ~NASC_QuantNumeric,
  colors = KrillPal,       
  type = "scatter3d",       
  mode = "markers",         
  marker = list(
    sizemode = "diameter",  
    sizeref = c(10),         
    symbol = "circle",      
    opacity = 10           
  ),
  showlegend = FALSE
) 


Krill3DB <- Krill3DB %>%
  add_trace(
    x = c(-58.93), 
    y = c(-62.22), 
    z = c(0),
    type = "scatter3d", 
    mode = "markers",
    marker = list(
      size = 1,            
      symbol = "triangle", 
      color = "darkblue"   
    ),
    showlegend = FALSE     
  )

Krill3DB <- Krill3DB %>%
  layout(
    scene = list(
      xaxis = list(
        title = "Longitude",
        autorange = "reversed",
        range = c(-59.4, -57.6)
      ),
      yaxis = list(
        title = "Latitude",
        autorange = "reversed",
        range = c(-62.75, -61.9)
      ),
      zaxis = list(
        title = "Swarm depth (m)",
        autorange = "reversed",
        range = c(200, 0),
        tickvals = seq(0, 250, 50)
      )
    )
  )


Krill3DB




# Export
orca(Krill3DA, "Fig1B_First.svg")

orca(Krill3DB, "Fig1B_Second.svg")





###############################################################################################################################################
###############################################################################################################################################
# KRILL MODELLING
###############################################################################################################################################
###############################################################################################################################################


All_Acoustic$Year <- as.factor(All_Acoustic$Year)

ModelAcoustic <- All_Acoustic %>%
  filter(Distance_Col < 38) 


KrillGAM_OVERALL <-  mgcv::gam(PRC_NASC ~  s(Distance_Col, k = 15, bs = "tp", by = Year), data = ModelAcoustic, method = "REML", family = gaussian(), select = FALSE) 

draw(KrillGAM_OVERALL, residuals = TRUE)
gam.check(KrillGAM_OVERALL)
summary.gam(KrillGAM_OVERALL)


smooth_plot_data <- smooth_estimates(KrillGAM_OVERALL) %>%
  filter(.smooth == "s(Distance_Col):Year2019/20") %>%
  add_confint()

smooth_plot_data1 <- smooth_estimates(KrillGAM_OVERALL) %>%
  filter(.smooth == "s(Distance_Col):Year2021/22") %>%
  add_confint()

joinedsmooths <- smooth_plot_data %>%
  full_join(smooth_plot_data1)


KrillMod <- joinedsmooths |>
  ggplot(aes(x = Distance_Col, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .smooth),
              alpha = 0.4) +
  geom_line(aes(colour = .smooth)) + 
  scale_fill_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +   # <-- change fill scale
  scale_colour_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  scale_x_continuous(limits = c(0,50), breaks = seq(0,50, by = 5)) +
  # geom_rug(data = BoutID,                  # <-- rug
  #          mapping = aes(x = Distance, y = NULL),
  #          sides = "b", alpha = 0.4) +  
  labs(title = NULL, y = "NASC (model estimate)", x = "Distance from penguin colony (km)",
       colour = "Year", fill = "Year") + 
  theme_bw(10) + theme(strip.background = element_rect(fill="gray85"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black"))+ 
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right") +
  # legend.direction = "vertical",
  # legend.justification = c(0.95,0.1),
  # legend.position = c(0.9,0.6),
  #  legend.text = element_text(size = rel(0.9)),  
  #  legend.title = element_text(size = rel(0.9)), 
  #  legend.key.size = unit(0.85, "lines"), 
  # legend.background = element_rect(fill = "lightgrey", colour = "black"),
  # legend.key = element_rect(fill = "lightgrey")) +
  # theme(plot.margin = unit(c(0,0,0,0), 'lines')) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), 'lines')) 

KrillMod




# Replot the histogram for Fig
KrillHist <- ggplot(KrillSummaries_Distance, aes(x = Distance_Bin, y = Total_NASC, fill = Year), alpha = 0.6) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  # scale_fill_viridis_d(name = "Year") +
  scale_fill_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +   # <-- change fill scale
  scale_colour_brewer(palette = "Set1", labels=c('2019/20', '2021/22')) +
  labs(
    # title = "Distribution of NASC in relation to distance from colony",
    x = NULL,
    y = "Total NASC"
  ) +
  scale_x_continuous(limits = c(0,50), breaks = seq(0,50, by = 5)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right") +
  theme(panel.grid.minor = element_blank()) 


KrillHist


Fig3 <- ggarrange(KrillHist, KrillMod, nrow = 2, align = "v", common.legend = TRUE, legend = "right", labels = "AUTO")
Fig3

tiff("Fig3.tiff", width = 8, height = 6, res = 300, units = "in")
Fig3
dev.off()



## Some summaries

First15 <- KrillSummaries_Distance %>%
  filter(Year == "2019/20") %>%
  mutate(Prop = Total_NASC/sum(Total_NASC)) %>% 
  filter(Distance_Bin <= 5)

sum(First15$Prop)  
sum(First15$Total_NASC) 

Second15 <- KrillSummaries_Distance %>%
  filter(Year == "2021/22") %>%
  mutate(Prop = Total_NASC/sum(Total_NASC)) %>%
  filter(Distance_Bin <= 15)

sum(Second15$Prop)  
sum(Second15$Total_NASC) 


###############################################################################################################################################
##############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

# End of krill script