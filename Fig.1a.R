library(ggplot2)
library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(ggpointdensity)
library(sf)
library(dplyr)
library(stringi)


data_row <- read.csv("/Users/yusha/Fig1a_wordmap.csv", stringsAsFactors = FALSE)



data_clean <- data_row %>%
  mutate(Category = trimws(Category)) %>%
  mutate(Category = gsub("[[:space:]]+", "_", Category)) %>%
  mutate(Category = stri_replace_all_regex(Category, "[\\p{C}]", "")) %>%
  filter(!is.na(Category), 
         Category != "", 
         !grepl("^[[:space:]]*$", Category)) %>%
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude)) %>%
  filter(!is.na(Latitude), !is.na(Longitude))


cat_colors <- c(
  "A" = "#BC3A24",
  "B" = "#EFB882", 
  "C" = "#06798F",
  "D" = "#7A9A01",
  "E" = "#384D73"
)


data_clean <- data_clean[data_clean$Category %in% names(cat_colors), ]
data_clean$Category <- factor(data_clean$Category, levels = names(cat_colors))


world_sf <- ne_countries(scale = "medium", returnclass = "sf")



set.seed(1234)


radius_scale <- 50000


points_sf <- st_as_sf(data_clean, coords = c("Longitude", "Latitude"), crs = 4326)


robinson_crs <- "+proj=robin +datum=WGS84"
points_proj <- st_transform(points_sf, crs = robinson_crs)
world_proj <- st_transform(world_sf, crs = robinson_crs)


coords <- st_coordinates(points_proj)
data_coords <- cbind(data_clean, coords) %>%
  group_by(X, Y) %>%
  mutate(n_dup = n(),
         dup_id = row_number()) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    
    angle = ifelse(n_dup > 1, 2 * pi * (dup_id - 1) / n_dup, 0),
    radius = ifelse(n_dup > 1, radius_scale * sqrt(n_dup), 0),
    X_jit = X + radius * cos(angle),
    Y_jit = Y + radius * sin(angle)
  ) %>%
  ungroup()



p1 <- ggplot() +
  geom_sf(data = world_proj, fill = "gainsboro", color = "white", linewidth = 0.5) +
  geom_point(
    data = data_coords,
    aes(x = X_jit, y = Y_jit, color = Category),
    size = 3,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = cat_colors,
    labels = c(
      "A" = "Invertebrate predators",
      "B" = "Invertebrate parasitoids", 
      "C" = "Invertebrate natural enemies (NEs)",
      "D" = "Invertebrate herbivores",
      "E" = "Crops"
    ),
    name = "Trophic groups:"
  ) +
  coord_sf(
    crs = robinson_crs,
    xlim = st_bbox(world_proj)[c("xmin", "xmax")] * 1.02,
    ylim = st_bbox(world_proj)[c("ymin", "ymax")] * 1.02,
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#87CEFA"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    
    legend.position = c(0.01, 0.01),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "gainsboro", color = "white", linewidth = 0.6),
    legend.margin = margin(5, 5, 5, 5),
    legend.box.margin = margin(5, 5, 5, 5),
    legend.text = element_text(size = 22),
    legend.title = element_text(face = "bold", size = 22),
    legend.spacing.y = unit(0.5, "cm"),       
    legend.key.height = unit(1, "cm")    
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 4),  
    keyheight = unit(1.5, "cm") 
  ))

p1

