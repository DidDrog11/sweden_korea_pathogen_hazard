library(terra)
library(tidyterra)
library(dplyr)
library(rgbif)
library(geodata)
library(ranger)
library(sdmTMB)
library(biscale)
library(ggplot2)
library(tidyverse)

# 1. Environmental Data Acquisition
# Retrieve the boundary and bioclimatic variables for South Korea
kor_boundary <- geodata::gadm("KOR",
                              level = 0,
                              path = tempdir())

clim_data <- geodata::worldclim_country("KOR",
                                        var = "bio",
                                        res = 2.5,
                                        path = tempdir()) |>
  terra::mask(kor_boundary)

# Clean names for modelling compatibility
names(clim_data) <- paste0("bio_", 1:19)

# 2. Host Species Distribution Model (Apodemus agrarius)
# Retrieve occurrence data and format as a SpatVector
host_occ <- rgbif::occ_data(scientificName = "Apodemus agrarius",
                            country = "KR",
                            hasCoordinate = TRUE,
                            limit = 2000)$data |>
  dplyr::select(decimalLongitude,
                decimalLatitude) |>
  terra::vect(geom = c("decimalLongitude", "decimalLatitude"),
              crs = "EPSG:4326")

# Generate pseudo-absences and extract environmental covariates
bg_points <- terra::spatSample(kor_boundary,
                               size = 5000,
                               method = "random")

occ_extract <- terra::extract(clim_data, host_occ) |>
  dplyr::mutate(presence = 1)
bg_extract <- terra::extract(clim_data, bg_points) |>
  dplyr::mutate(presence = 0)

sdm_data <- dplyr::bind_rows(occ_extract,
                             bg_extract) |>
  dplyr::select(-ID) |>
  stats::na.omit() |>
  dplyr::mutate(presence = as.factor(presence))

# Fit Random Forest SDM and predict to the raster extent
rf_model <- ranger::ranger(presence ~ .,
                           data = sdm_data,
                           probability = TRUE)

# Custom prediction function for terra::predict with ranger probability outputs
pred_fun <- function(model, data) {
  predict(model, data)$predictions[, 2]
}

host_sdm_raster <- terra::predict(clim_data,
                                  rf_model,
                                  fun = pred_fun,
                                  na.rm = TRUE)
names(host_sdm_raster) <- "host_suitability"

# 3. Pathogen Prevalence Model (Model-Based Geostatistics)
arha_data <- read_csv(here::here("korea_prevalence_data.csv"))

# Create the spatial mesh for the SPDE approximation
mesh <- sdmTMB::make_mesh(arha_data,
                          xy_cols = c("longitude", "latitude"),
                          cutoff = 0.15)

# Fit the spatial GLMM for binomial prevalence data
fit_spde <- sdmTMB::sdmTMB(cbind(n_positive, n_tested - n_positive) ~ 1,
                           data = arha_data,
                           mesh = mesh,
                           family = binomial(link = "logit"),
                           spatial = "on")

# Extract coordinate grid from the host SDM for prediction
pred_grid <- terra::as.data.frame(host_sdm_raster,
                                  xy = TRUE,
                                  na.rm = TRUE) |>
  dplyr::rename(longitude = x,
                latitude = y)

sim_draws <- predict(fit_spde, newdata = pred_grid, nsim = 100)

prev_pred <- pred_grid |>
  dplyr::mutate(est = apply(sim_draws, 1, mean),
                uncertainty = apply(sim_draws, 1, sd),
                prevalence_prob = stats::plogis(est))

# Convert back to a SpatRaster and align with the host SDM
uncertainty_raster <- terra::rast(prev_pred, type = "xyz", crs = "EPSG:4326") |>
  terra::subset("uncertainty") |>
  terra::resample(host_sdm_raster)

# 4. Bivariate Mapping for Sampling Prioritisation
# Combine rasters and convert to dataframe for biscale classification
hazard_stack <- c(host_sdm_raster, uncertainty_raster)

hazard_df <- terra::as.data.frame(hazard_stack,
                                  xy = TRUE,
                                  na.rm = TRUE)

# Create bivariate classes (3x3 grid) based on quantiles
biv_data <- biscale::bi_class(hazard_df,
                              x = host_suitability,
                              y = uncertainty,
                              style = "quantile",
                              dim = 4)

# Generate the map
map_plot <- ggplot2::ggplot() +
  ggplot2::geom_tile(data = biv_data,
                     aes(x = x, y = y, fill = bi_class),
                     show.legend = FALSE) +
  biscale::bi_scale_fill(pal = "DkCyan2",
                         dim = 4) +
  tidyterra::geom_spatvector(data = kor_boundary,
                             fill = NA,
                             colour = "black",
                             linewidth = 0.3) +
  ggplot2::coord_sf(expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Adaptive Sampling Prioritisation",
                x = "",
                y = "") +
  ggplot2::theme(text = element_text(size = 7, family = "sans"))

# Generate the bivariate legend
legend <- biscale::bi_legend(pal = "DkCyan2",
                             dim = 4,
                             xlab = "Higher Host Suitability",
                             ylab = "Higher Uncertainty",
                             size = 8)

# Combine map and legend using cowplot or patchwork
panel_b <- cowplot::ggdraw() +
  cowplot::draw_plot(map_plot, 0, 0, 1, 1) +
  cowplot::draw_plot(legend, 0.65, 0.15, 0.25, 0.25)

# 5. Extract and Smooth True 90th Percentile Target Polygons
# Calculate the exact mathematical 90th percentiles
suit_90 <- stats::quantile(hazard_df$host_suitability, probs = 0.90, na.rm = TRUE)
uncert_90 <- stats::quantile(hazard_df$uncertainty, probs = 0.90, na.rm = TRUE)

# Filter the raw data strictly by these thresholds, ignoring the visual 3-3 class
target_df <- hazard_df |>
  dplyr::filter(host_suitability >= suit_90 & uncertainty >= uncert_90) |>
  dplyr::select(x, y, host_suitability) |>
  dplyr::mutate(target_zone = 1) # Indicator for the raster

target_rast <- terra::rast(target_df, type = "xyz", crs = "EPSG:4326")

# Convert to polygons and apply morphological closing
priority_polys <- terra::as.polygons(target_rast) |>
  terra::buffer(width = 5000) |>
  terra::buffer(width = -2500) |>
  terra::mask(kor_boundary)

# 6. Simulate Missing Pipeline Layers (WP1.2 & WP3.2)
# Panel A: Multimodal AI Ingestion (3 distinct embedding layers)
fine_rast <- terra::rast(kor_boundary, res = 0.01)

# Simulate Layer 1: ViT Structural Embedding
coarse1 <- terra::rast(kor_boundary, res = 0.05)
terra::values(coarse1) <- stats::rnorm(terra::ncell(coarse1))
vit_struct <- terra::resample(coarse1, fine_rast, method = "bilinear") |> terra::mask(kor_boundary)

# Simulate Layer 2: ViT Phenological Shift
coarse2 <- terra::rast(kor_boundary, res = 0.08)
terra::values(coarse2) <- stats::runif(terra::ncell(coarse2))
vit_pheno <- terra::resample(coarse2, fine_rast, method = "bilinear") |> terra::mask(kor_boundary)

# Simulate Layer 3: NLP Agricultural Intensity
coarse3 <- terra::rast(kor_boundary, res = 0.05)
terra::values(coarse3) <- stats::rexp(terra::ncell(coarse3))
nlp_agri <- terra::resample(coarse3, fine_rast, method = "bilinear") |> terra::mask(kor_boundary)

# Stack and name the layers for ggplot faceting
ai_stack <- c(vit_struct, vit_pheno, nlp_agri)
names(ai_stack) <- c("ViT: Structural", "ViT: Phenology", "NLP: Agriculture")
# Panel C: Logistical Constraints (Roads, Static Exclusions, Weather Mask)
precip_raster <- clim_data[["bio_16"]]
precip_thresh <- stats::quantile(terra::values(precip_raster), probs = 0.85, na.rm = TRUE)
weather_mask_rast <- terra::ifel(precip_raster >= precip_thresh, 1, NA)
dynamic_weather_mask <- terra::as.polygons(weather_mask_rast) |>
  terra::buffer(width = 2000)

# Static Exclusion: Inaccessible Mountainous Terrain (Elevation > 1000m)
elev_data <- geodata::elevation_30s(country = "KOR", path = tempdir()) |>
  terra::mask(kor_boundary)
elev_mask_rast <- terra::ifel(elev_data > 1000, 1, NA)
static_exclusion <- terra::as.polygons(elev_mask_rast) |>
  terra::buffer(width = 1000)

# Road Network & Accessibility Buffer (Major highways for reproducible plotting)
kor_roads_sf <- rnaturalearth::ne_download(scale = 10, type = "roads", category = "cultural", returnclass = "sf") |>
  sf::st_intersection(sf::st_as_sf(kor_boundary))

roads <- terra::vect(kor_roads_sf)
road_buffer <- terra::buffer(roads, width = 3000)

# Panel D: Calculate Final Feasible Sites
# Intersect smoothed priority targets with road access, then erase exclusions
feasible_sites <- terra::intersect(priority_polys, road_buffer) |>
  terra::erase(static_exclusion) |>
  terra::erase(dynamic_weather_mask)

sample_pins <- terra::spatSample(feasible_sites, size = 15, method = "random")

# 7. Construct Panels A, C, D
mainland_xlim <- c(125.7, 129.7)
mainland_ylim <- c(34.2, 38.6)

base_theme <- function() {
  list(
    ggplot2::coord_sf(xlim = mainland_xlim, ylim = mainland_ylim, expand = FALSE),
    ggplot2::theme_void(),
    ggplot2::theme(plot.title = ggplot2::element_text(size = 7, face = "bold", hjust = 0.5, family = "sans"),
                   plot.subtitle = ggplot2::element_text(size = 6, hjust = 0.5, family = "sans"),
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
                   plot.margin = ggplot2::margin(5, 5, 5, 5))
  )
}

# Sub-theme for the Panel A mini-maps
sub_theme <- function() {
  list(
    ggplot2::coord_sf(xlim = mainland_xlim, ylim = mainland_ylim, expand = FALSE),
    ggplot2::theme_void(),
    ggplot2::theme(plot.title = ggplot2::element_text(size = 5, face = "italic", hjust = 0.5),
                   panel.border = ggplot2::element_rect(colour = "grey70", fill = NA, linewidth = 0.3),
                   # Increasing the margins here shrinks the map relative to its bounding box (~80% size)
                   plot.margin = ggplot2::margin(15, 12, 15, 12, "pt"))
  )
}

# Panel A: Build 3 separate plots with unique colour scales
p_a1 <- ggplot2::ggplot() +
  tidyterra::geom_spatvector(data = kor_boundary, fill = "grey95", colour = NA) +
  tidyterra::geom_spatraster(data = vit_struct, show.legend = FALSE) +
  ggplot2::scale_fill_viridis_c(option = "magma", na.value = "transparent") +
  tidyterra::geom_spatvector(data = kor_boundary, fill = NA, colour = "grey30", linewidth = 0.3) +
  ggplot2::labs(title = "ViT: Structural") +
  sub_theme()

p_a2 <- ggplot2::ggplot() +
  tidyterra::geom_spatvector(data = kor_boundary, fill = "grey95", colour = NA) +
  tidyterra::geom_spatraster(data = vit_pheno, show.legend = FALSE) +
  ggplot2::scale_fill_viridis_c(option = "mako", na.value = "transparent") +
  tidyterra::geom_spatvector(data = kor_boundary, fill = NA, colour = "grey30", linewidth = 0.3) +
  ggplot2::labs(title = "ViT: Phenology") +
  sub_theme()

p_a3 <- ggplot2::ggplot() +
  tidyterra::geom_spatvector(data = kor_boundary, fill = "grey95", colour = NA) +
  tidyterra::geom_spatraster(data = nlp_agri, show.legend = FALSE) +
  ggplot2::scale_fill_viridis_c(option = "plasma", na.value = "transparent") +
  tidyterra::geom_spatvector(data = kor_boundary, fill = NA, colour = "grey30", linewidth = 0.3) +
  ggplot2::labs(title = "NLP: Agriculture") +
  sub_theme()

# Stitch Panel A together, apply a grouping title and background, and wrap it
library(patchwork)
panel_a_raw <- (p_a1 | p_a2 | p_a3) +
  patchwork::plot_annotation(title = "A. Multimodal AI Ingestion (WP1.2)",
                             theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5, family = "sans"),
                                                    # Adds a subtle grey box around the 3 inputs to group them
                                                    plot.background = ggplot2::element_rect(colour = "black", fill = "grey98", linewidth = 0.5),
                                                    plot.margin = ggplot2::margin(5, 5, 5, 5)))

# Wrap elements ensures patchwork treats this entire block as a single cohesive panel in the final 2x2 grid
panel_a <- patchwork::wrap_elements(panel_a_raw)

# Panel B: Rebuild and INSET the legend
legend <- biscale::bi_legend(pal = "DkCyan2",
                             dim = 4,
                             xlab = "Suitability",
                             ylab = "Uncertainty",
                             size = 6) +
  ggplot2::theme(plot.margin = ggplot2::margin(4, 4, 4, 4, "pt"), # Minimal breathing room for text
                 plot.background = ggplot2::element_rect(fill = "white", 
                                                         colour = "black", 
                                                         linewidth = 0.3), # Crisp inset border
                 panel.background = ggplot2::element_blank())

panel_b_gg <- map_plot +
  ggplot2::labs(title = "B. INLA-SPDE Hazard & Entropy (WP2.1)",
                subtitle = "Host Suitability vs. Pathogen Model Uncertainty",) +
  base_theme() +
  patchwork::inset_element(legend, 
                           left = 0.72,
                           bottom = 0.02,
                           right = 0.98,
                           top = 0.26)

# Panel C & D
panel_c <- ggplot2::ggplot() +
  tidyterra::geom_spatvector(data = kor_boundary, fill = "grey95", colour = "grey50") +
  tidyterra::geom_spatvector(data = road_buffer, fill = "#73C6B6", colour = NA, alpha = 0.5) +
  tidyterra::geom_spatvector(data = static_exclusion, fill = "black", colour = NA, alpha = 0.6) +
  tidyterra::geom_spatvector(data = dynamic_weather_mask, fill = "#3498DB", colour = NA, alpha = 0.5) +
  tidyterra::geom_spatvector(data = roads, colour = "grey20", linewidth = 0.2) +
  ggplot2::labs(title = "C. KDCA Logistics & Dynamic Masking (WP3.2)") +
  base_theme()

panel_d <- ggplot2::ggplot() +
  tidyterra::geom_spatvector(data = kor_boundary, fill = "grey95", colour = "grey50") +
  tidyterra::geom_spatvector(data = dynamic_weather_mask, fill = "#3498DB", colour = NA, alpha = 0.1, linetype = "dashed") +
  tidyterra::geom_spatvector(data = feasible_sites, fill = "#F39C12", colour = "black", linewidth = 0.4) +
  tidyterra::geom_spatvector(data = sample_pins, colour = "black", fill = "red", shape = 21, size = 1.5) +
  ggplot2::labs(title = "D. sitetool Target Deployment (WP3.3)") +
  base_theme()

# ---------------------------------------------------------
# 8. Assemble Final Figure (Pure Patchwork)
# ---------------------------------------------------------
# Using the Asymmetric Layout (Left column wider than right column)
final_4panel_figure <- (panel_a + panel_b_gg) / (panel_c + panel_d) +
  patchwork::plot_layout(widths = c(1.8, 1)) +
  patchwork::plot_annotation(title = "The Active Learning Intelligence Engine",
                             theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5, family = "sans")))

print(final_4panel_figure)

ggplot2::ggsave(filename = "active_learning_engine.png",
                plot = final_4panel_figure,
                width = 17,
                height = 19,
                units = "cm",
                dpi = 300,         # High resolution for grant peer review
                bg = "white")

ggplot2::ggsave(filename = here::here("figures", "panel_1.png"), plot = panel_a, width = 17, height = 19, units = "cm", dpi = 300, bg = "white")
ggplot2::ggsave(filename = here::here("figures", "panel_2.png"), plot = panel_b_gg, width = 17, height = 19, units = "cm", dpi = 300, bg = "white")
ggplot2::ggsave(filename = here::here("figures", "panel_3.png"), plot = panel_c, width = 17, height = 19, units = "cm", dpi = 300, bg = "white")
ggplot2::ggsave(filename = here::here("figures", "panel_4.png"), plot = panel_d, width = 17, height = 19, units = "cm", dpi = 300, bg = "white")
