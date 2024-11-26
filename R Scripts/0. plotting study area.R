#####4.5 plotting study area for predators 
library(cowplot)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)



ARS <- read.csv(file = "CoordsCovs.csv", header = TRUE) %>%  
  st_as_sf(coords = c("x", "y"), crs = "epsg:7842")

trapsuniqueplot <- trapstemp[c(1,2,3)]
trapsuniqueplot <- st_as_sf(trapsuniqueplot, coords =c("x", "y"), crs = st_crs(ARS))


coordscovsplot <- coords
studyareaplot <- plotfixed
urban_label_data <- subset(ARS, Habitat == "Urban")

st_crs(studyareaplot) <- st_crs(trapsuniqueplot)

#welcome to hell
st_write(studyareaplot, dsn = "plotshape.shp")

#reread in, i've gone arcmap route
plotfixed <- st_read("studyareaclipped.shp")
st_crs(plotfixed) <- st_crs(ARS)


plot(studyareadissolve)

st_crs(studyareaplot)
st_crs(trapsuniqueplot)
st_crs(ARS)

register_stadiamaps("f9d511a5-33ca-42aa-a527-943c07b6fe0f", write = FALSE)

st_crs(plotfixed) <-"+proj=utm +zone=56 +south +datum=WGS84 +units=m +no_defs"
st_crs(ARS) <- "+proj=utm +zone=56 +south +datum=WGS84 +units=m +no_defs"
st_crs(trapsuniqueplot) <- "+proj=utm +zone=56 +south +datum=WGS84 +units=m +no_defs"

plot(trapsuniqueplot)

ARSlat <- st_transform(ARS, crs = st_crs("+proj=longlat +datum=WGS84 +no_defs"))
plotfixedlat <- st_transform(plotfixed, crs = st_crs("+proj=longlat +datum=WGS84 +no_defs"))
trapslat <- st_transform(trapsuniqueplot, crs = st_crs("+proj=longlat +datum=WGS84 +no_defs"))

###Make AUS map

australia_map <- ne_states(country = "Australia", returnclass = "sf")


australia_map <- st_transform(australia_map, crs = st_crs("+proj=longlat +datum=WGS84 +no_defs"))

australia_map_det <- ne_countries(scale = 10,country = "Australia", returnclass = "sf")

australia_map_det <- st_transform(australia_map_det, crs = st_crs("+proj=longlat +datum=WGS84 +no_defs"))

####Create the bounding box for study area
bbox <- st_bbox(plotfixedlat) %>% 
  st_as_sfc() %>% 
  st_transform(crs = st_crs("+proj=longlat +datum=WGS84 +no_defs"))

australia_inset <- ggplot() +
  geom_sf(data = australia_map, fill = "lightgrey", color = "black") +
  geom_sf(data = bbox, fill = NA, color = "red", lwd = 0.8) +
  coord_sf(xlim = c(140, 155), ylim = c(-39, -26.6)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6),         panel.spacing = margin(10, 10, 10, 10))

australia_inset






ARSlat$Name[ARSlat$Name == "Pacific Palms 1"] <- "Pacific Palms"
ARSlat$Name[ARSlat$Name == "Hawks Nest "] <- "Hawks Nest"
ARSlat$Name[ARSlat$Name == "Smith's Lake"] <- "Smiths Lake"
ARSlat$Habitat[ARSlat$Habitat == "Urban"] <- "Urban Area"

trapsjustcoords <- trapslat
trapsPS <- trapsjustcoords[c(1:22),]
trapsPS$Habitat <- "Trail Camera Station"
trapsBA <- trapsjustcoords[c(23:61),]
trapsBA$Habitat <- "Baited Camera Station"
trapsjustcoords <- rbind(trapsPS, trapsBA)


trapsjustcoords



trapsjustcoords$Name <- "Ligma"

trapsjustcoords <- trapsjustcoords[c("Habitat", "geometry", "Name")]
ARSlatjustcoords <- ARSlat[c("Habitat", "geometry", "Name")]
alljustcoords <- rbind(trapsjustcoords, ARSlatjustcoords)


alljustcoords$Name[5] <- "River Blocking Movement West"
alljustcoords$Habitat <- factor(alljustcoords$Habitat, levels = c("Baited Camera Station", "Trail Camera Station", "Urban Area"))
levels(alljustcoords$Habitat)

alljustcoords <- alljustcoords %>% 
  arrange(Habitat)

ars_colors <- ifelse(ARSlat$Habitat == "Urban Area", 'orange', 'blue')

library(ggmap)

plotfixed_bbox <- st_bbox(plotfixedlat)
plotfixed_bbox <- as.numeric(plotfixed_bbox)

basemap <- get_stadiamap(bbox = c(left = 152.12, 
                                  bottom = -32.73,
                                  right = 152.56, 
                                  top = -32.25),
                         crop = TRUE,
                         zoom = 10, maptype = "stamen_terrain_background")

basemapgg <- ggmap(basemap)
?ggmap
basemapgg
bbox = c(left = 152.12, 
         bottom = -32.73,
         right = 152.56, 
         top = -32.25)

bbox_df <- data.frame(
  xmin = bbox["left"],
  xmax = bbox["right"],
  ymin = bbox["bottom"],
  ymax = bbox["top"]
)

habmapplot <-basemapgg+ # semi-transparent
  geom_sf(data = alljustcoords, aes(fill = Habitat), inherit.aes = FALSE, size = 3.5, shape = 21)+
  scale_fill_manual(values = c("Urban Area" = "white", "Trail Camera Station" = "purple", "Baited Camera Station" ='red'), labels = c("Urban Area" = "Urban Area", "Trail Camera Station" = "Trail Camera", "Baited Camera Station" = "Baited Camera"))+
  labs(fill = "")+
  geom_sf_label(data = subset(alljustcoords, Habitat == 'Urban Area' & Name == 'Smiths Lake'), aes(label = Name), size = 4,  nudge_x = -0.07, nudge_y = 0, label.padding = unit(0.15, "lines"), inherit.aes = FALSE)+
  geom_sf_label(data = subset(alljustcoords, Habitat == 'Urban Area' & Name == 'Pacific Palms'), aes(label = Name), size = 4,  nudge_x = -0.07, nudge_y = 0, label.padding = unit(0.15, "lines"), inherit.aes = FALSE)+
  geom_sf_label(data = subset(alljustcoords, Habitat == 'Urban Area' & Name == 'Hawks Nest'), aes(label = Name), size = 4,  nudge_x = 0.07, nudge_y = 0, label.padding = unit(0.15, "lines"), inherit.aes = FALSE)+
  geom_sf_label(data = subset(alljustcoords, Habitat == 'Trail Camera Station' & Name == 'River Blocking Movement West'), aes(label = "River and Lakes\n Blocking \nMovement West"), size = 4,  nudge_x = -0.05, nudge_y = 0.1, label.padding = unit(0.2, "lines"), inherit.aes = FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16), axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13), legend.text = element_text(size = 19), legend.key.height = unit(2, "lines"),legend.margin=margin(-10, 0, 0, 0),panel.border = element_rect(color = "red", fill = NA, size = 0.7))

habmapplot

combinedinst <- ggdraw() +
  draw_plot(habmapplot, 0, 0, 1, 1) +
  draw_plot(australia_inset, 0.65, 0.06, 0.36, 0.36)

combinedinst


png("Figures/Study Area.jpg", width = 10, height =9, res= 300, units = "in")

habmapplot
combinedinst <- ggdraw() +
  draw_plot(habmapplot, 0, 0, 1, 1) +
  draw_plot(australia_inset, 0.65, 0.08, 0.36, 0.36)

combinedinst


dev.off()
