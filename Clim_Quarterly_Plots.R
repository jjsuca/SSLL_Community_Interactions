library(lubridate)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(viridis)
world <- ne_countries(scale=10,returnclass = "sf")#generate high res coastlines 

Species <- c(
  "Xiphias gladius",
  "Coryphaena hippurus",
  "Lepidocybium flavobrunneum",
  "Kajikia audax",
  "Isurus oxyrinchus",
  "Lampris guttatus",
  "Ruvettus pretiosus"
)

All_Clim <- list()
Abund_Files<-list.files("~/Documents/Eco_Cast_Plus_2025/Model_Predictions/Abundance")
PA_Files<-list.files("~/Documents/Eco_Cast_Plus_2025/Model_Predictions/PA")
Env_Vars_SSLL_Quarter_Deg<-readRDS("~/Documents/LL_Data/Weekly_Quart_Deg_SSLL_GLORYS_Vars_Chl.rds")
Positions<-cbind(Env_Vars_SSLL_Quarter_Deg$Lon_Bin, Env_Vars_SSLL_Quarter_Deg$Lat_Bin)
Dates<-Env_Vars_SSLL_Quarter_Deg$Date_for_Pairing

for(sp in Species){
  
  abund_file <- grep(
    paste0("^", sp, "_"),
    Abund_Files,
    value = TRUE
  )
  
  has_pa <- sp != "Xiphias gladius"
  if(has_pa){
    
    pa_file <- grep(
      
      paste0("^", sp, "_"),
      
      PA_Files,
      
      value = TRUE
      
    )
    
  }
  Species_AB <- readRDS(
    file.path(
      "~/Documents/Eco_Cast_Plus_2025/Model_Predictions/Abundance",
      abund_file
    )
  )
  
  Mean_Species_AB <- rowMeans(
    
    exp(Species_AB),
    
    na.rm = TRUE
    
  )
  
  if(has_pa){
    
    Species_PA <- readRDS(
      
      file.path(
        
        "~/Documents/Eco_Cast_Plus_2025/Model_Predictions/PA",
        
        pa_file
        
      )
      
    )
    
    Mean_Species_PA <- rowMeans(
      
      Species_PA,
      
      na.rm = TRUE
      
    )
    
    CPUE <- Mean_Species_AB * Mean_Species_PA
    
  } else {
    
    CPUE <- Mean_Species_AB
    
  }
  df <- data.frame(
    Date = as.Date(Dates),
    Longitude = Positions[,1] - 360,
    Latitude = Positions[,2],
    CPUE = CPUE
  )
  
  df$Quarter <- quarter(df$Date)
  
  Clim <- df %>%
    filter(Quarter %in% c(1,2,4)) %>%
    group_by(
      Quarter,
      Longitude,
      Latitude
    ) %>%
    dplyr::summarize(
      CPUE = mean(CPUE, na.rm = TRUE),
      .groups = "drop"
    )
  
  Clim$Species <- sp
  
  All_Clim[[sp]] <- Clim
}

All_Clim <- bind_rows(All_Clim)
All_Clim <- All_Clim %>%
  
  group_by(Species) %>%
  
  mutate(
    
    CPUE_std = CPUE / max(CPUE, na.rm = TRUE)
    
  ) %>%
  
  ungroup()
All_Clim$Species <- factor(
  
  All_Clim$Species,
  
  levels = c(
    
    "Xiphias gladius",
    "Coryphaena hippurus",
    "Lepidocybium flavobrunneum",
    "Kajikia audax",
    "Isurus oxyrinchus",
    "Lampris guttatus",
    "Ruvettus pretiosus"
    
  ),
  
  labels = c(
    
    "Swordfish",
    "Common\ndolphinfish",
    "Escolar",
    "Striped\nmarlin",
    "Shortfin\nmako",
    "Opah",
    "Oilfish"
    
  )
  
)

All_Clim$Quarter <- factor(
  All_Clim$Quarter,
  levels = c(1,2,4),
  labels = c("Q1","Q2","Q4")
)



png(
  "~/Documents/Eco_Cast_Plus_2025/Figures/SSLL_Maps/Revisions_Species_Quarterly_Climatologies.png",
  width = 12,
  height = 14,
  units = "in",
  res = 300
)

p <- ggplot() +
  geom_raster(
    data = All_Clim,
    aes(
      x = Longitude,
      y = Latitude,
      fill = CPUE_std
    )
  ) +
  geom_sf(
    data = world,
    inherit.aes = FALSE
  ) +
  coord_sf(
    xlim = c(-178,-120),
    ylim = c(18,40),
    expand = FALSE
  ) +
  facet_grid(
    Species ~ Quarter
  ) +
  scale_fill_viridis(
    limits = c(0,1),
    name = "Relative\nCPUE"
  ) +
  theme_bw() +
  scale_x_continuous(
    
    breaks = c(-175, -160, -145, -130)
    
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_bw() +
  
  theme(
    
    strip.text.y = element_text(
      
      size = 16,
      
      face = "bold"
      
    ),
    strip.text.x = element_text(
      
      size = 20,
      
      face = "bold"
      
    ),
    
    
    # Tick labels
    axis.text.x = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 16,
      face = "bold"
    ),
    
    # Axis titles
    axis.title.x = element_text(
      size = 18,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 18,
      face = "bold"
    ),
    
    # Legend
    legend.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 16
    ),
    
    # Panel spacing
    panel.grid = element_blank(),
    
    panel.spacing = unit(
      0.3,
      "lines"
    )
  )


print(p)
dev.off()

All_Clim <- All_Clim %>%
  
  group_by(Species) %>%
  
  mutate(
    
    CPUE_std = CPUE / max(CPUE, na.rm = TRUE)
    
  ) %>%
  
  ungroup()

# Create panel labels (a-r for 6 species x 3 quarters)

panel_labels <- expand.grid(
  Species = levels(All_Clim$Species),
  Quarter = levels(All_Clim$Quarter)
)

panel_labels$Panel <- LETTERS[seq_len(nrow(panel_labels))]

# Upper-left corner of each map
panel_labels$Longitude <- -177
panel_labels$Latitude  <- 39

png(
  "~/Documents/Eco_Cast_Plus_2025/Figures/SSLL_Maps/Species_Quarterly_Climatologies_Revisions_UpdateV2.png",
  width = 12,
  height = 14,
  units = "in",
  res = 300
)

p <- ggplot() +
  
  geom_raster(
    data = All_Clim,
    aes(
      x = Longitude,
      y = Latitude,
      fill = CPUE_std
    )
  ) +
  
  geom_sf(
    data = world,
    inherit.aes = FALSE,
    linewidth = 0.2
  ) +
  
  geom_text(
    data = panel_labels,
    aes(
      x = Longitude,
      y = Latitude,
      label = Panel
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 5
  ) +
  
  coord_sf(
    xlim = c(-178, -120),
    ylim = c(18, 40),
    expand = FALSE
  ) +
  
  facet_grid(
    Species ~ Quarter,
    switch = "y"
  ) +
  
  scale_fill_viridis(
    name = "Relative\nCPUE"
  ) +
  scale_x_continuous(
    
    breaks = c(-175, -160, -145, -130)
    
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_bw() +
  
  theme(
    
    strip.text.y = element_text(
      
      size = 16,
      
      face = "bold"
      
    ),
    strip.text.x = element_text(
      
      size = 20,
      
      face = "bold"
      
    ),
  
    
    # Tick labels
    axis.text.x = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 16,
      face = "bold"
    ),
    
    # Axis titles
    axis.title.x = element_text(
      size = 18,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 18,
      face = "bold"
    ),
    
    # Legend
    legend.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 16
    ),
    
    # Panel spacing
    panel.grid = element_blank(),
    
    panel.spacing = unit(
      0.3,
      "lines"
    )
  )

print(p)

dev.off()

