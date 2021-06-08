library(imager)
library(ggplot2)
library(ggquiver)
library(scico)
library(isoband)
library(sf)
library(particles)
library(tidygraph)
library(conflicted)

conflict_prefer("simulate", "particles")
conflict_prefer("position", "particles")

blur_slimr <- load.image("logo/blurred.png")

blur_slimr <- crop.borders(blur_slimr, 6, 6)

blur_slimr <- grayscale(flatten.alpha(blur_slimr)) 

plot(blur_slimr)

blur_slimr <- resize(blur_slimr, 80, 20, interpolation_type = 6)

plot(blur_slimr)

blur_slimr <- 1 - blur_slimr

grad <- imgradient(blur_slimr)

plot(grad)

grad_data <- expand_grid(x = 1:80, y = 1:20) %>%
  mutate(x_grad = as.matrix(grad[[1]])[cbind(x, y)],
         y_grad = as.matrix(grad[[2]])[cbind(x, y)],
         val = as.matrix(blur_slimr)[cbind(x, y)])

iso <- isobands(1:80, 1:20, t(as.matrix(blur_slimr)), 
                c(0, 0.2, 0.4, 0.6, 0.8), c(0.2, 0.4, 0.6, 0.8, 1.0)) %>%
  iso_to_sfg() %>%
  st_as_sfc() %>%
  st_as_sf()

ggplot(grad_data) + 
  geom_raster(aes(x, y, fill = val), alpha = 0.6) +
  geom_quiver(aes(x, y, u = x_grad, v = y_grad), vecsize = 2) +
  scale_fill_scico(palette = "bamako") +
  scale_y_reverse() +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none")


x_vel <- as.matrix(grad[[1]])
y_vel <- as.matrix(grad[[2]])

vec_scale <- 0.25

np <- 200
ni <- 1000

sim <- create_empty(np) %>% 
  simulate(alpha_decay = 0, setup = aquarium_genesis(vel_max = 0,
                                                     width = ncol(x_vel),
                                                     height = nrow(x_vel))) %>% 
  wield(reset_force, xvel = 0, yvel = 0) %>% 
  wield(random_force, xmin = -0.1, xmax = 0.1, ymin = -0.1, ymax = 0.1) %>% 
  wield(field_force, x = x_vel * vec_scale, y = y_vel * vec_scale, 
        xlim = c(-ncol(x_vel) / 2, ncol(x_vel) / 2), 
        ylim = c(-nrow(x_vel) / 2, nrow(x_vel) / 2)) %>% 
  evolve(ni, record)

traces <- data.frame(do.call(rbind, lapply(sim$history, position)))
names(traces) <- c('x', 'y')
traces$particle <- rep(1:np, ni)

ggplot(traces) +
  geom_path(aes(x, y, group = particle), size = 0.1, alpha = 0.5) + 
  theme_void() + 
  theme_minimal() +
  theme(legend.position = 'none')
