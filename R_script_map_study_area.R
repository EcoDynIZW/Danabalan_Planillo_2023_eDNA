# Map of study area 

source("source_packages.R")


## Load Data

# sampling points
sampling_points <- st_read("data/geodata/sampling_points_25833.gpkg", crs = 25833)

extra_points <- read.csv("data/geodata/OtherStudy_points.csv") %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(crs = 25833)

# maps
berlin_border <- read_sf("data/geodata/Berlin_border_25833.gpkg", crs = 25833)
water <- read_sf("data/geodata/waterbodies_Berlin_25833.gpkg", crs = 25833)
tree <- raster("data/geodata/raster_tree_cover_berlin_25833.tif")
names(tree) <- "tree_cover"


##########################
## Prepare ggplot map

# put together all the sampling points 
sampling_points2 <- sampling_points %>% 
  as.data.frame() %>% 
  dplyr::select(id_plot, geom) %>% 
  mutate(type = case_when(id_plot == "Oh_04" ~ "only mosquitoes", 
                          TRUE ~ "flies + mosquitoes")) %>% 
  st_as_sf()

extra_points2 <- extra_points %>% 
  as.data.frame() %>% 
  dplyr::rename(id_plot = Location) %>% 
  mutate(type = "Hoffman et al. (2018)") %>% 
  rename(geom = geometry) %>% 
  st_as_sf()

all_points <- rbind(sampling_points2, extra_points2)

all_points <- all_points %>% 
  mutate(type = fct_relevel(type, "flies + mosquitoes", "only mosquitoes", "Hoffman et al. (2018)"))


# create an object to plot the tree cover
tree_rc <- reclassify(tree, rcl = 
                        c(0, 20, 10, 
                          20, 40, 30, 
                          40, 60, 50, 
                          60, 80, 70,
                          80, 100, 90))
tree2_df_tmp <- rasterToPoints(tree_rc)
tree2_df <- data.frame(tree2_df_tmp)
colnames(tree2_df)[3] <- "tree_cover"

# load images
img_mosq <- readPNG("data/PNGs/mosquito-gada61ba60_1920.png")
g_mosq <- rasterGrob(img_mosq, interpolate=TRUE)

#fly image
img_fly <- readPNG("data/PNGs/fly-gd1ec52dda_1280.png")
g_fly <- rasterGrob(img_fly, interpolate=TRUE)

# sampling picture
img_sampling <- readPNG("data/PNGs/Picture1.png")



#######################
### create map

# base map
map1 <- ggplot() +
  geom_sf(data = berlin_border) +
  geom_sf(data = water, fill = "blue", colour = "transparent") +
  geom_raster(data = tree2_df, aes(x = x, y = y, fill = as.factor(tree_cover))) +
  scale_fill_brewer(palette = "Greens", na.value = NA, 
                    name = "Tree cover percentage", 
                    labels = c("0 to 20", "20 to 40","40 to 60", 
                               "60 to 80", "80 to 100")
  ) +
  geom_sf(data = all_points, aes(pch = type), 
          fill = c("orange", "orange", "orange", "orange", "red", 
                           "grey40", "grey40", "grey40", "grey40", "grey40", "grey40"),
          size = c(4,4,4,4,4,
                   3,3,3,3,3,3)) +
  scale_shape_manual(name = "Sampling plots",
                     values = c(21,21,24)) +
  guides(shape = guide_legend(override.aes = list(fill = (c("orange", "red", "grey40")), 
                                                  size = c(4,4,3)), 
                              order = 2), 
         fill = guide_legend(order = 1, 
                             override.aes = list(size = 0.1))) +
  annotation_scale(height = unit(0.2, "cm")) +
  annotation_north_arrow(
    pad_x = unit(0.3, "cm"),
    pad_y = unit(1, "cm"),
    height = unit(0.8, "cm"),
    width = unit(0.6, "cm")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NULL, color = "grey80"), 
    axis.title = element_blank(),
    legend.position = c(0.87,0.8),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8,"lines")
  ) 


# add images
bbox(as_Spatial(berlin_border))

map_sampling <- map1 +
  annotation_custom(g_mosq, xmin=410784.6 , xmax=415784.6 , ymin=5812963, ymax=5819963) +
  annotation_custom(g_fly, xmin=410784.6 , xmax=415784.6 , ymin=5817963, ymax=5821963) +
  annotation_raster(img_sampling, xmin=370001.2, xmax=380784.6 , ymin=5830963, ymax=5838000)

# add names of sampling plots
plot_names <- sampling_points %>% 
  dplyr::select(id_plot) %>% 
  mutate(xcoord = st_coordinates(sampling_points)[,1], 
         ycoord = st_coordinates(sampling_points)[,2]) %>% 
  st_drop_geometry()
plot_names$labels <- c("Grunewald-1", 
                       "Grunewald-2", 
                       "Spandau-2",
                       "Mueggelsee", 
                       "Spandau-1")

(map_sampling2 <- map_sampling +
    geom_label_repel(data = plot_names, aes(label = labels, x = xcoord, y = ycoord)) +
    theme(plot.background = element_blank())
)


ggsave(plot = map_sampling2, 
       "output/Fig1_map_studyArea.png", 
       dpi = 600, height = 7, width = 7)

