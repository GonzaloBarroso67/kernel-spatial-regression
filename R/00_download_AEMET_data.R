# =========================================================
# Script 00 - Data acquisition (DOCUMENTARY ONLY)
# =========================================================
# AEMET daily precipitation data for all stations, 2024.
# Output: annual accumulated precipitation per station,
# projected to UTM 30N (EPSG:25830), restricted to mainland Spain.
#
# IMPORTANT - this script does NOT need to be run to reproduce
# the analysis. It documents how the raw data were obtained from
# the AEMET OpenData API. A full download takes several hours
# (the API is rate-limited, so the code includes long waits and
# retries between blocks) and requires a personal AEMET API key.
#
# The resulting dataset is stored in the repository at
#   data/aemet_2024_annual_precip_points_sf.rds
# and every downstream script reads from that file. To reproduce
# the results, start from script 01 and do not run this one.
#
# To run it anyway, request a free API key at
#   https://opendata.aemet.es/centrodedescargas/altaUsuario
# and set it below (section 2).
# =========================================================


# =========================================================
# 1. Load libraries
# =========================================================
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(lubridate)
library(stringr)
library(readr)
library(sf)
library(tibble)


# =========================================================
# 2. API key
# =========================================================

api_key <- "Write your API key here"


# =========================================================
# 3. Helper functions
# =========================================================

# Convert AEMET precipitation strings to numeric. AEMET uses a
# comma as decimal separator, so it is replaced by a dot first.
to_num_aemet <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, ",", ".")
  as.numeric(x)
}


# Convert AEMET coordinates (degrees-minutes-seconds packed into a
# single string, with a trailing direction letter) to decimal degrees.
# Latitudes come as DDMMSS (6 digits) and longitudes as DDDMMSS (7),
# so the number of digits is used to split the fields. S and W are
# stored as negative.
coord_aemet_to_decimal <- function(coord) {
  
  coord <- as.character(coord)
  coord <- str_trim(coord)
  
  if (is.na(coord) || coord == "") {
    return(NA_real_)
  }
  
  direction <- str_sub(coord, -1)
  numbers <- str_sub(coord, 1, -2)
  
  if (nchar(numbers) == 6) {
    degrees <- as.numeric(str_sub(numbers, 1, 2))
    minutes <- as.numeric(str_sub(numbers, 3, 4))
    seconds <- as.numeric(str_sub(numbers, 5, 6))
    
  } else if (nchar(numbers) == 7) {
    degrees <- as.numeric(str_sub(numbers, 1, 3))
    minutes <- as.numeric(str_sub(numbers, 4, 5))
    seconds <- as.numeric(str_sub(numbers, 6, 7))
    
  } else {
    return(NA_real_)
  }
  
  decimal <- degrees + minutes / 60 + seconds / 3600
  
  if (direction %in% c("S", "W")) {
    decimal <- -decimal
  }
  
  decimal
}


# Split a full year into consecutive 10-day blocks. The AEMET daily
# endpoint limits the length of each request, so the year is downloaded
# block by block rather than in a single call.
make_date_blocks <- function(year) {
  
  starts <- seq.Date(
    from = as.Date(paste0(year, "-01-01")),
    to   = as.Date(paste0(year, "-12-31")),
    by   = "10 days"
  )
  
  tibble(
    block_id = row_number(starts),
    start = starts,
    end = pmin(starts + days(9), as.Date(paste0(year, "-12-31")))
  )
}


# =========================================================
# 4. Function to download one block of daily data
# =========================================================
# The AEMET API works in two steps: the first request does not return
# the data but a URL ("datos") pointing to the actual file, which is
# then fetched in a second request. Both steps are wrapped in retries
# with long waits, because the API rate-limits requests aggressively
# and returns 429 / "caudal por minuto excedido" when overloaded.

get_aemet_daily_block <- function(start_date, end_date, api_key, max_tries = 5) {
  
  start_txt <- paste0(start_date, "T00%3A00%3A00UTC")
  end_txt   <- paste0(end_date, "T23%3A59%3A00UTC")
  
  url <- paste0(
    "https://opendata.aemet.es/opendata/api/valores/climatologicos/diarios/datos/",
    "fechaini/", start_txt,
    "/fechafin/", end_txt,
    "/todasestaciones/",
    "?api_key=", api_key
  )
  
  for (attempt in seq_len(max_tries)) {
    
    message("Attempt ", attempt, "/", max_tries, ": ", start_date, " to ", end_date)
    
    Sys.sleep(30 + 10 * attempt)
    
    first_response <- tryCatch(
      GET(url, timeout(120)),
      error = function(e) e
    )
    
    if (inherits(first_response, "error")) {
      warning("Connection error in first request: ", first_response$message)
      next
    }
    
    first_status <- status_code(first_response)
    
    if (first_status == 429) {
      warning("Rate limit in first request. Waiting longer.")
      Sys.sleep(90)
      next
    }
    
    if (first_status != 200) {
      warning("First request failed with status ", first_status)
      Sys.sleep(60)
      next
    }
    
    first_text <- content(first_response, as = "text", encoding = "UTF-8")
    first_json <- fromJSON(first_text)
    
    # The first response must contain a "datos" field with the URL of
    # the actual data file; if it is missing the block is retried.
    if (!"datos" %in% names(first_json)) {
      warning("No 'datos' URL returned.")
      next
    }
    
    data_response <- tryCatch(
      GET(first_json$datos, timeout(180)),
      error = function(e) e
    )
    
    if (inherits(data_response, "error")) {
      warning("Connection error in data request: ", data_response$message)
      next
    }
    
    data_status <- status_code(data_response)
    # AEMET data files are served in Latin-9 (ISO-8859-15).
    data_text <- content(data_response, as = "text", encoding = "ISO-8859-15")
    
    # Rate limiting can also appear inside the body of a 200 response,
    # so the text is inspected for the known limit messages.
    if (
      data_status == 429 ||
      grepl("429 Too Many Requests", data_text) ||
      grepl("caudal por minuto excedido", data_text) ||
      grepl("L.mite de peticiones", data_text)
    ) {
      warning("Rate limit in data request. Waiting longer.")
      Sys.sleep(100)
      next
    }
    
    if (data_status != 200) {
      warning("Data request failed with status ", data_status)
      Sys.sleep(60)
      next
    }
    
    data_json <- fromJSON(data_text)
    
    # A response carrying a "descripcion" field is an AEMET status
    # message (e.g. no data available), not the expected data table.
    if (is.list(data_json) && "descripcion" %in% names(data_json)) {
      warning(
        "AEMET returned message: ",
        paste(data_json$descripcion, collapse = " "),
        " | Estado: ",
        paste(data_json$estado, collapse = " ")
      )
      Sys.sleep(60)
      next
    }
    
    return(as_tibble(data_json))
  }
  
  warning("Block failed after all attempts: ", start_date, " to ", end_date)
  tibble()
}


# =========================================================
# 5. Download daily data by blocks
# =========================================================
# Each block is saved to disk as soon as it is downloaded. Blocks
# already present are skipped, so an interrupted download can be
# resumed without repeating completed requests.

year <- 2024

blocks <- make_date_blocks(year)

blocks_folder <- paste0("aemet_blocks_", year)
dir.create(blocks_folder, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(blocks))) {
  
  block_id <- blocks$block_id[i]
  start_i  <- blocks$start[i]
  end_i    <- blocks$end[i]
  
  file_i <- file.path(
    blocks_folder,
    paste0("block_", sprintf("%02d", block_id), ".rds")
  )
  
  if (file.exists(file_i)) {
    message("Block ", block_id, " already saved. Skipping.")
    next
  }
  
  message(
    "Downloading block ", block_id, " of ", nrow(blocks),
    ": ", start_i, " to ", end_i
  )
  
  block_data <- get_aemet_daily_block(
    start_date = start_i,
    end_date   = end_i,
    api_key    = api_key,
    max_tries  = 5
  )
  
  saveRDS(block_data, file_i)
  
  message("Saved block ", block_id, " with ", nrow(block_data), " rows.")
  
  Sys.sleep(90)
}


# =========================================================
# 6. Combine downloaded blocks
# =========================================================

block_files <- list.files(
  blocks_folder,
  pattern = "\\.rds$",
  full.names = TRUE
)

aemet_daily_raw <- map_dfr(block_files, readRDS)

write_csv(aemet_daily_raw, paste0("aemet_", year, "_daily_raw.csv"))


# =========================================================
# 7. Download station inventory
# =========================================================
# The inventory provides station metadata and, crucially, the
# coordinates used later to place each station in space.

get_aemet_station_inventory <- function(api_key, max_tries = 6) {
  
  url <- paste0(
    "https://opendata.aemet.es/opendata/api/valores/climatologicos/",
    "inventarioestaciones/todasestaciones/",
    "?api_key=", api_key
  )
  
  for (attempt in seq_len(max_tries)) {
    
    message("Inventory attempt ", attempt, "/", max_tries)
    
    Sys.sleep(60 + 20 * attempt)
    
    first_response <- tryCatch(
      GET(url, timeout(120)),
      error = function(e) e
    )
    
    if (inherits(first_response, "error")) {
      warning("Connection error in inventory request: ", first_response$message)
      next
    }
    
    first_status <- status_code(first_response)
    
    if (first_status == 429) {
      warning("Rate limit in inventory request. Waiting longer.")
      Sys.sleep(120)
      next
    }
    
    if (first_status != 200) {
      warning("Inventory request failed with status ", first_status)
      next
    }
    
    first_text <- content(first_response, as = "text", encoding = "UTF-8")
    first_json <- fromJSON(first_text)
    
    if (!"datos" %in% names(first_json)) {
      warning("No 'datos' URL returned for station inventory.")
      next
    }
    
    data_response <- tryCatch(
      GET(first_json$datos, timeout(120)),
      error = function(e) e
    )
    
    if (inherits(data_response, "error")) {
      warning("Connection error in inventory data request: ", data_response$message)
      next
    }
    
    data_status <- status_code(data_response)
    data_text <- content(data_response, as = "text", encoding = "ISO-8859-15")
    
    if (
      data_status == 429 ||
      grepl("429 Too Many Requests", data_text) ||
      grepl("caudal por minuto excedido", data_text) ||
      grepl("L.mite de peticiones", data_text)
    ) {
      warning("Rate limit in inventory data request. Waiting longer.")
      Sys.sleep(120)
      next
    }
    
    if (data_status != 200) {
      warning("Inventory data request failed with status ", data_status)
      next
    }
    
    return(as_tibble(fromJSON(data_text)))
  }
  
  stop("Station inventory failed after all attempts.")
}

stations_raw <- get_aemet_station_inventory(api_key)

write_csv(stations_raw, "aemet_station_inventory_raw.csv")


# =========================================================
# 8. Clean station coordinates
# =========================================================
# The packed AEMET coordinates are converted to decimal degrees and
# only the fields needed downstream are kept.

stations <- stations_raw %>%
  mutate(
    lat = map_dbl(latitud, coord_aemet_to_decimal),
    lon = map_dbl(longitud, coord_aemet_to_decimal)
  ) %>%
  select(
    indicativo,
    nombre_inventory = nombre,
    provincia_inventory = provincia,
    lat,
    lon
  )


# =========================================================
# 9. Clean daily precipitation data
# =========================================================
# AEMET encodes some daily values as text codes rather than numbers:
#   "Ip"   = trace precipitation (too small to measure); treated as 0.
#   "Acum" = the value was accumulated into a later day's reading, so
#            this day has no usable figure on its own; treated as NA.
# Empty strings are also set to NA. All remaining values are numeric.

aemet_daily <- aemet_daily_raw %>%
  mutate(
    fecha = as.Date(fecha),
    prec_original = str_trim(as.character(prec)),
    prec = case_when(
      is.na(prec_original) ~ NA_character_,
      prec_original == "" ~ NA_character_,
      prec_original == "Ip" ~ "0",
      prec_original == "Acum" ~ NA_character_,
      TRUE ~ prec_original
    ),
    prec = to_num_aemet(prec)
  )


# =========================================================
# 10. Join with station coordinates
# =========================================================

aemet_daily_points <- aemet_daily %>%
  left_join(stations, by = "indicativo")


# =========================================================
# 11. Annual accumulated precipitation per station
# =========================================================
# For each station the valid daily values are summed, and the total is
# rescaled to a full 366-day year (2024 is a leap year) as
#   annual_precip = sum_valid / n_valid * 366.
# This compensates for missing days by assuming they behave like the
# observed ones, so stations with a few gaps are not penalised. The
# distribution of missing days is examined in script 01 to justify
# this assumption.

aemet_annual_precip_all <- aemet_daily_points %>%
  group_by(
    indicativo,
    nombre_inventory,
    provincia_inventory,
    lon,
    lat
  ) %>%
  summarise(
    annual_precip_raw = if_else(
      sum(!is.na(prec)) == 0,
      NA_real_,
      sum(prec, na.rm = TRUE)
    ),
    n_valid_prec = sum(!is.na(prec)),
    n_missing_prec = sum(is.na(prec)),
    n_acum = sum(prec_original == "Acum", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    annual_precip = if_else(
      !is.na(annual_precip_raw),
      annual_precip_raw / n_valid_prec * 366,
      NA_real_
    )
  )

# Restrict the study to mainland Spain: the Balearic and Canary
# Islands, Ceuta and Melilla are excluded
aemet_annual_precip_mainland <- aemet_annual_precip_all %>%
  filter(
    !provincia_inventory %in% c(
      "ILLES BALEARS",
      "BALEARES",
      "LAS PALMAS",
      "SANTA CRUZ DE TENERIFE",
      "STA. CRUZ DE TENERIFE",
      "CEUTA",
      "MELILLA"
    ),
    !is.na(lon),
    !is.na(lat)
  )

# Quality filter: keep only stations with a positive annual total and
# at least 300 valid daily values, so the rescaling to 366 days is
# based on a sufficiently complete record.
aemet_annual_precip <- aemet_annual_precip_mainland %>%
  filter(
    !is.na(annual_precip),
    annual_precip > 0,
    n_valid_prec >= 300
  ) %>%
  mutate(
    log_annual_precip = log(annual_precip)
  )


# =========================================================
# 12. Freeze monthly missing-data counts for script 01
# =========================================================
# The full daily table is not frozen (it is large and only needed for
# diagnostics), but the monthly counts of missing and "Acum" values are
# saved so that script 01 can check how gaps are distributed over the
# year. This supports the rescaling-to-366 assumption: if gaps are
# spread fairly evenly across months rather than concentrated in wet or
# dry seasons, the rescaling introduces little bias.

final_indicativos <- aemet_annual_precip$indicativo

missing_by_station_month <- aemet_daily %>%
  filter(indicativo %in% final_indicativos) %>%
  filter(is.na(prec)) %>%
  mutate(month = month(fecha)) %>%
  group_by(indicativo, month) %>%
  summarise(
    n_missing = n(),
    n_acum = sum(prec_original == "Acum", na.rm = TRUE),   # <- esta línea
    .groups = "drop"
  )

saveRDS(
  missing_by_station_month,
  paste0("aemet_", year, "_missing_by_station_month.rds")
)


# =========================================================
# 13. Project coordinates, collapse duplicate locations,
#     and summarise the filtering steps
# =========================================================
# Coordinates are projected from geographic (EPSG:4326) to UTM 30N
# (EPSG:25830) and then divided by 1000, so that x and y are expressed
# in kilometres to match the scale used by the trend model.

aemet_annual_precip_sf <- st_as_sf(
  aemet_annual_precip,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

aemet_annual_precip_proj <- st_transform(
  aemet_annual_precip_sf,
  crs = 25830
)

coords_utm <- st_coordinates(aemet_annual_precip_proj)

aemet_annual_precip_proj <- aemet_annual_precip_proj %>%
  mutate(
    x = coords_utm[, 1] / 1000,
    y = coords_utm[, 2] / 1000
  )

# A few stations share exactly the same projected (x, y) location, and
# npsp requires unique locations for local linear estimation. These are
# collapsed into a single point: annual_precip is averaged on the
# original scale and log_annual_precip is then recomputed once from that
# average, so the two variables stay consistent. lon and lat are also
# averaged (the merged stations coincide, so the mean is representative),
# and indicativo / province are concatenated to keep a record of which
# stations were merged.
model_data <- aemet_annual_precip_proj %>%
  st_drop_geometry() %>%
  group_by(x, y) %>%
  summarise(
    annual_precip = mean(annual_precip, na.rm = TRUE),
    n_valid_prec = max(n_valid_prec, na.rm = TRUE),
    n_missing_prec = max(n_missing_prec, na.rm = TRUE),
    n_acum = max(n_acum, na.rm = TRUE),
    indicativo = paste(unique(indicativo), collapse = " / "),
    provincia_inventory = paste(unique(provincia_inventory), collapse = " / "),
    n_stations_same_location = n(),
    lon = mean(lon, na.rm = TRUE),
    lat = mean(lat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    log_annual_precip = log(annual_precip)
  )

# Rebuild the sf object from the collapsed data so that the frozen
# dataset has one row per unique location. x and y are recomputed from
# the projection to stay fully consistent with the geometry.
aemet_annual_precip_sf <- st_as_sf(
  model_data,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

aemet_annual_precip_utm <- st_transform(
  aemet_annual_precip_sf,
  crs = 25830
)

coords_utm <- st_coordinates(aemet_annual_precip_utm)

aemet_annual_precip_utm <- aemet_annual_precip_utm %>%
  mutate(
    x = coords_utm[, 1] / 1000,
    y = coords_utm[, 2] / 1000
  )

# Station counts after each step, documenting how the final sample is
# obtained.
n_initial <- n_distinct(stations$indicativo)
n_after_join_aggregation <- nrow(aemet_annual_precip_all)
n_after_scope <- nrow(aemet_annual_precip_mainland)
n_after_quality <- nrow(aemet_annual_precip)
n_final <- nrow(aemet_annual_precip_utm)

filter_summary <- tibble(
  step = c(
    "Initial station inventory",
    "After aggregation",
    "After mainland Spain and coordinate filter",
    "After precipitation quality filter",
    "After collapsing duplicate locations"
  ),
  n_stations = c(
    n_initial,
    n_after_join_aggregation,
    n_after_scope,
    n_after_quality,
    n_final
  ),
  removed_from_previous_step = c(
    NA,
    n_initial - n_after_join_aggregation,
    n_after_join_aggregation - n_after_scope,
    n_after_scope - n_after_quality,
    n_after_quality - n_final
  )
)

filter_summary


# =========================================================
# 14. Save final datasets
# =========================================================
# The .rds file keeps the sf geometry and is the frozen dataset read
# by every downstream script. The .csv is a plain copy without
# geometry, for inspection only.

write_csv(
  st_drop_geometry(aemet_annual_precip_utm),
  paste0("aemet_", year, "_annual_precip_points.csv")
)

saveRDS(
  aemet_annual_precip_utm,
  paste0("aemet_", year, "_annual_precip_points_sf.rds")
)

