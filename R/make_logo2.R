library(imager)
library(isoband)
library(sf)
library(dplyr)
library(ggplot2)

png_vers <- load.image("logo/slimr_vector.png") %>%
  flatten.alpha() %>%
  grayscale()
plot(png_vers)

png_dim <- dim(png_vers)
slimr_poly <- isobands(1:png_dim[1], png_dim[2]:1, t(png_vers[ , , 1, 1]), 
                    levels_low = -0.1, levels_high = 0.5) %>%
  iso_to_sfg() %>%
  st_as_sfc() 

plot(slimr_poly, col = "grey")

coords <- slimr_poly %>%
  st_coordinates() %>%
  as.data.frame()

y_range <- range(coords$Y)
x_range <- range(coords$X)

grid_sf <- st_make_grid(slimr_poly, cellsize = c(20, 20))
plot(grid_sf)

grid_intersect <- grid_sf %>%
  st_sf() %>%
  st_intersects(slimr_poly %>%
                  st_sf())

yes_grid <- grid_sf[sapply(grid_intersect, function(x) length(x) != 0)]
no_grid <- grid_sf[sapply(grid_intersect, function(x) length(x) == 0)]

yes_centroids <- st_centroid(yes_grid) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(X = round(X, 2)) %>%
  group_by(X) %>%
  mutate(bar_id = cur_group_id())

yes_grid <- yes_grid %>%
  st_as_sf() %>%
  mutate(bar_id = yes_centroids$bar_id) %>%
  group_by(bar_id) %>%
  summarise(count = n(), do_union = TRUE)

no_centroids <- st_centroid(no_grid) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(X = round(X, 3)) %>%
  group_by(X) %>%
  mutate(bar_id = cur_group_id())

no_grid <- no_grid %>%
  st_as_sf() %>%
  mutate(bar_id = no_centroids$bar_id) %>%
  group_by(bar_id) %>%
  summarise(count = n(), do_union = TRUE)


plot(yes_grid)
plot(no_grid)

plot(yes_grid %>% select(x), col = "red")
plot(no_grid %>% select(x), add = TRUE, col = "green")

yes_grid <- yes_grid %>%
  ungroup() %>%
  st_cast("MULTIPOLYGON") %>%
  st_cast("POLYGON", do_split = TRUE) %>%
  mutate(cents = st_coordinates(st_centroid(x))) %>%
  mutate(X = round(cents[ , 1], 3), Y = cents[ , 2],
         inside = "yes") 

no_grid <- no_grid %>%
  ungroup() %>%
  st_cast("MULTIPOLYGON") %>%
  st_cast("POLYGON", do_split = TRUE) %>%
  mutate(cents = st_coordinates(st_centroid(x))) %>%
  mutate(X = round(cents[ , 1], 3), Y = cents[ , 2],
         inside = "no") 

plot(no_grid %>% select(Y))

all_bars <- rbind(yes_grid, no_grid)


plot(all_bars)

dg <- colorspace::darken("#00843D", 0.2)
lg <- colorspace::lighten("#00843D", 0.2)

col_bars <- all_bars %>%
  group_by(X) %>%
  mutate(y_rank = rank(1 / Y),
         col = case_when(y_rank == 1 & inside == "no" ~ "grey20",
                            y_rank == 1 & inside == "yes" ~ "#FFCD00",
                            y_rank == 2 & inside == "no" ~ dg,
                            y_rank == 2 & inside == "yes" ~ "#FFCD00",
                            y_rank == 3 & inside == "no" ~ dg,
                            y_rank == 3 & inside == "yes" ~ lg,
                            y_rank == 4 & inside == "no" ~ "darkred",
                            y_rank == 4 & inside == "yes" ~ lg,
                            TRUE ~ "idunno"))

plot(col_bars %>% select(y_rank))
plot(col_bars %>% select(col))

p <- ggplot(col_bars) +
  geom_sf(aes(fill = col, colour = colorspace::darken(col, 0.3)),
          size = 0.1) +
  scale_fill_identity() +
  scale_colour_identity() +
  theme_void() +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  theme(panel.background = element_rect(colour = "grey80", fill = "black",
                                        size = 1))

p

library(hexSticker)

hexit <- sticker(p, package = "",
                 s_x = 1, s_y = 1,
                 s_width = 1.8,
                 s_height = 1,
                 h_fill = "black",
                 h_color = "#cc0000",
                 filename = "slimr_prelim_hex.png")

plot(hexit)