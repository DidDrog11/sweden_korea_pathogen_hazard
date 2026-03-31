library(terra)
library(tidyterra)
library(tidyverse)
library(rgbif)
library(gstat)
library(automap)
library(ranger)
library(sf)
library(stars)

# 1. Define Core Parameters
hosts <- c("Rattus norvegicus", "Apodemus agrarius", "Crocidura lasiura")
viruses <- c("Orthohantavirus hantanense", "Orthohantavirus seoulense")
bounds <- ext(125.8, 129.8, 33.0, 38.6)

# 2. Geodata Acquisition and Spatial Harmonisation
kor_polygon <- geodata::gadm("KOR",
                             level = 0,
                             path = tempdir())

prk_polygon <- geodata::gadm("PRK",
                             level = 0,
                             path = tempdir())

# Extract bioclimatic variables and mask to South Korea
env_layers <- geodata::worldclim_country("KOR",
                                         var = "bio",
                                         path = tempdir()) |>
  crop(kor_polygon) |>
  mask(kor_polygon)

# 3. GBIF Occurrence Ingestion
taxon_keys <- hosts |>
  lapply(function(x) name_backbone(name = x)$usageKey) |>
  unlist()

download_req <- occ_download(
  pred_in("taxonKey", taxon_keys),
  pred("country", "KR"),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV")

occ_download_wait(download_req)

occ_records <- as.character(download_req) |>
  occ_download_get(overwrite = TRUE) |>
  occ_download_import()

occ_filtered <- occ_records |>
  filter(!is.na(decimalLongitude),
         decimalLongitude > 100,
         species %in% hosts)

# 4. Species Distribution Modelling (SDM) 
fit_sdm <- function(species_data,
                    predictors,
                    n_bg = 5000) {
  
  # 1. Format presence data and extract environmental covariates
  pres_pts <- vect(species_data,
                   geom = c("decimalLongitude", "decimalLatitude"),
                   crs = crs(predictors))
  
  pres_env <- terra::extract(predictors, pres_pts, ID = FALSE) |>
    mutate(presence = 1) |>
    na.omit()
  
  # 2. Sample background (pseudo-absences) directly from valid raster cells
  bg_env <- spatSample(predictors,
                       size = n_bg,
                       na.rm = TRUE,
                       xy = FALSE,
                       values = TRUE) |>
    mutate(presence = 0)
  
  model_df <- bind_rows(pres_env, bg_env)
  
  # 3. Fit the Random Forest
  # keep.inbag = TRUE is strictly required to calculate spatial standard errors later
  rf_model <- ranger(presence ~ .,
                     data = model_df,
                     num.trees = 500,
                     keep.inbag = TRUE)
  
  # 4. Custom prediction wrapper for terra::predict
  # Extracts both the mean prediction and the standard error
  pred_fun <- function(model, data) {
    p <- predict(model, data = data, type = "se")
    cbind(p$predictions, p$se)
  }
  
  # 5. Project across the environmental layers
  predict(predictors,
          rf_model,
          fun = pred_fun,
          na.rm = TRUE) |>
    setNames(c("prediction", "uncertainty"))
}

sdms <- map(hosts,
            ~ fit_sdm(occ_filtered[occ_filtered$species == .x, ], env_layers)) |>
  setNames(hosts)

# 5. Spatial Kriging for Pathogen Prevalence
arha <- readRDS("C:/Users/ucbtds4/R_Repositories/arenavirus_hantavirus/data/database/Project_ArHa_database_2026-01-09.rds")
arha_pathogen <- arha$host |>
  left_join(arha$pathogen, by = "host_record_id") |>
  select(host_record_id, pathogen_record_id, host_species, iso3c, latitude, longitude, pathogen_species_cleaned, pathogen_species_ncbi, number_tested, number_positive) |>
  filter(iso3c %in% c("KOR", "PRK"))

hanta_agrarius_sf <- arha_pathogen |>
  filter(pathogen_species_ncbi == "Orthohantavirus hantanense",
         host_species == "Apodemus agrarius") |>
  group_by(longitude, latitude) |>
  summarise(n_tested = sum(number_tested, na.rm = TRUE),
            n_positive = sum(number_positive, na.rm = TRUE),
            .groups = "drop") |>
  filter(n_tested > 0) |>
  mutate(prevalence = n_positive / n_tested) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# 2. Project to match the environmental template (required for gstat)
# Assuming env_layers was loaded in the previous step and is in a projected CRS
hanta_agrarius_proj <- st_transform(hanta_agrarius_sf, st_crs(env_layers))

# 3. Generate the spatial prediction and variance surfaces
krige_prevalence <- function(data_pts,
                             template_raster) {
  
  # 1. Fit the variogram
  vgm_fit <- autofitVariogram(prevalence ~ 1,
                              input_data = as(data_pts, "Spatial"))
  
  # 2. Convert to stars
  grid_stars <- st_as_stars(template_raster)
  
  # 3. Perform Kriging
  krige_out <- krige(prevalence ~ 1,
                     locations = data_pts,
                     newdata = grid_stars,
                     model = vgm_fit$var_model)
  
  # 4. Convert to SpatRaster and subset by index, NOT by name
  out_rast <- rast(krige_out)[[1:2]] |>
    setNames(c("prev_predict", "prev_var"))
  
  return(out_rast)
}

hanta_surfaces <- krige_prevalence(data_pts = hanta_agrarius_proj, 
                                   template_raster = env_layers[[1]])

w_prevalence <- 0.5
w_dose <- 0.5
w_uncertainty <- 1 - w_dose

dose_priority <- init(env_layers[[1]], fun = 0)
uncertainty_priority <- init(env_layers[[1]], fun = 0)
prevalence_priority <- init(env_layers[[1]], fun = 0)

# 2. Accumulate Host SDM Priorities
for (h in hosts) {
  P <- sdms[[h]]$prediction
  U <- sdms[[h]]$uncertainty
  
  dose_priority <- dose_priority + (w_dose * P)
  uncertainty_priority <- uncertainty_priority + (w_uncertainty * U)
}

# 3. Accumulate Pathogen Prevalence Priorities
hanta_var_rescaled <- app(hanta_surfaces$prev_var, 
                          fun = function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))

prevalence_priority <- prevalence_priority + hanta_var_rescaled


total_priority <- (w_prevalence * prevalence_priority) + ((1 - w_prevalence) * uncertainty_priority)

total_priority <- app(total_priority, 
                      fun = function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))

# 7. Spatial Sampling (Active Learning Acquisition)
# Emulating Julia's Balanced Acceptance Sampling using weighted spatSample
bon_samples <- spatSample(total_priority,
                          size = 50,
                          method = "weights",
                          replace = FALSE,
                          as.points = TRUE)

# 8. Visualisation using tidyterra
priority_plot <- ggplot() +
  geom_spatraster(data = total_priority) +
  scale_fill_viridis_c(option = "magma",
                       name = "Total Priority",
                       na.value = "transparent") +
  geom_spatvector(data = prk_polygon,
                  fill = "grey85",
                  color = "grey40",
                  linewidth = 0.5) +
  geom_spatvector(data = bon_samples,
                  color = "white",
                  fill = "black",
                  shape = 21,
                  size = 3) +
  coord_sf(xlim = c(bounds[1], bounds[2]),
           ylim = c(bounds[3], bounds[4]),
           expand = FALSE) +
  theme_void() +
  theme(legend.position = "bottom")
