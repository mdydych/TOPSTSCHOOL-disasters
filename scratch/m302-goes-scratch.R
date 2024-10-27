
# Function to get GOES files with satellite selection----
get_goes_files <- function(date,
                           hour,
                           product = "ABI-L2-AODC",
                           satellite = "goes18") {  # Add satellite parameter

  # Get date components
  date_parts <- goes_date_format(date)

  # Construct prefix
  prefix <- base::sprintf(
    "%s/%s/%s/%02d/",
    product,
    date_parts$year,
    date_parts$doy,
    hour
  )

  base::message("Searching prefix: ", prefix, " in ", satellite)

  # Get file list from bucket
  files_df <- aws.s3::get_bucket_df(
    bucket = base::paste0("noaa-", satellite),
    prefix = prefix,
    region = "us-east-1",
    max = 20
  )

  if (base::nrow(files_df) > 0) {
    files_df$LastModified <- base::as.POSIXct(
      files_df$LastModified,
      format = "%Y-%m-%dT%H:%M:%S.000Z",
      tz = "UTC"
    )
    files_df$satellite <- satellite  # Add satellite info
  }

  return(files_df)
}

# Function to construct https URL from Key
construct_url <- function(key, bucket = "noaa-goes18") {  # Default to GOES-18
  base::sprintf("https://%s.s3.amazonaws.com/%s", bucket, key)
}

# Function to read GOES data with appropriate projection----
read_goes_data <- function(url, satellite = "goes18") {
  # Set projection parameters based on satellite
  sat_params <- base::list(
    goes18 = base::list(
      lon_0 = -137.2,
      name = "GOES West"
    ),
    goes16 = base::list(
      lon_0 = -75.2,
      name = "GOES East"
    )
  )

  # Create temporary file
  temp_file <- base::tempfile(fileext = ".nc")

  # Download file
  utils::download.file(
    url = url,
    destfile = temp_file,
    mode = "wb",
    quiet = TRUE
  )

  # Set projection
  goes_projection <- base::sprintf(
    "+proj=geos +h=35786023.0 +lon_0=%f +sweep=x +ellps=GRS80",
    sat_params[[satellite]]$lon_0
  )

  # Read data
  base::message("Reading ", sat_params[[satellite]]$name, " data...")
  data <- stars::read_stars(
    temp_file,
    sub = "AOD"
  )

  # Set proper projection using sf
  sf::st_crs(data) <- goes_projection

  # Clean up
  base::unlink(temp_file)

  return(data)
}

# Function to format GOES date components
goes_date_format <- function(date) {
  base::list(
    year = base::format(date, "%Y"),
    doy = base::format(date, "%j")
  )
}

# Example of how to track a plume using both satellites----
track_smoke_plume <- function(date, hour) {
  # Get data from both satellites
  west_files <- get_goes_files(date, hour, satellite = "goes18")
  east_files <- get_goes_files(date, hour, satellite = "goes16")

  # Process most recent file from each satellite
  results <- base::list()

  if (base::nrow(west_files) > 0) {
    latest_west <- west_files |>
      dplyr::arrange(dplyr::desc(LastModified)) |>
      dplyr::slice(1)

    results$west <- read_goes_data(
      construct_url(latest_west$Key, "noaa-goes18"),
      "goes18"
    )
  }

  if (base::nrow(east_files) > 0) {
    latest_east <- east_files |>
      dplyr::arrange(dplyr::desc(LastModified)) |>
      dplyr::slice(1)

    results$east <- read_goes_data(
      construct_url(latest_east$Key, "noaa-goes16"),
      "goes16"
    )
  }

  return(results)
}

# Create buffer around Yellowknife for local analysis
yk_buffer <- sf::st_buffer(yk_boundary, 50000)  # 50km buffer

# Define analysis period (3 days before to 3 days after evacuation)
analysis_dates <- base::seq(
  base::as.Date("2023-08-13"),  # 3 days before evacuation
  base::as.Date("2023-08-19"),  # 3 days after evacuation
  by = "day"
)

# Define daily hours to analyze (daylight hours in UTC)
analysis_hours <- c(16, 18, 20, 22)  # Approximately 10am-4pm local time

# Create analysis timestamp combinations
analysis_times <- base::expand.grid(
  date = analysis_dates,
  hour = analysis_hours
)

analyze_local_aod <- function(date, hour, local_area) {
  # Get GOES files for timestamp
  files_df <- get_goes_files(
    date = date,
    hour = hour,
    product = "ABI-L2-AODC",
    satellite = "goes18"
  )

  if (base::nrow(files_df) == 0) {
    base::warning("No files found for ", date, " hour ", hour)
    return(NULL)
  }

  # Get most recent file
  latest_file <- files_df |>
    dplyr::arrange(dplyr::desc(LastModified)) |>
    dplyr::slice(1)

  # Read AOD data
  url <- construct_url(latest_file$Key, bucket = "noaa-goes18")
  aod_data <- read_goes_data(url)

  # Get local area centroid in lat/lon
  centroid <- sf::st_transform(local_area, 4326) |>
    sf::st_centroid() |>
    sf::st_coordinates()

  base::message("Area centroid (lon/lat): ",
                base::paste(centroid, collapse = ", "))

  # Convert to GOES grid coordinates
  x_grid <- base::round((centroid[1] + 137.2 + 75) * 2500/150)
  y_grid <- base::round((centroid[2] + 75) * 1500/150)

  # Define a window around the centroid
  window_size <- 100  # grid cells
  x_min <- base::max(1, x_grid - window_size)
  x_max <- base::min(2500, x_grid + window_size)
  y_min <- base::max(1, y_grid - window_size)
  y_max <- base::min(1500, y_grid + window_size)

  base::message("Grid coordinates: ",
                "x: ", x_min, "-", x_max,
                ", y: ", y_min, "-", y_max)

  # Extract the subset directly from the array
  aod_values <- aod_data[[1]][x_min:x_max, y_min:y_max]

  # Calculate statistics
  stats <- base::list(
    datetime = base::as.POSIXct(
      base::paste(date, base::sprintf("%02d:00:00", hour)),
      tz = "UTC"
    ),
    mean_aod = base::mean(aod_values, na.rm = TRUE),
    max_aod = base::max(aod_values, na.rm = TRUE),
    valid_pixels = base::sum(!base::is.na(aod_values)),
    total_pixels = base::length(aod_values)
  )

  # Verify we have valid data
  if (!base::is.nan(stats$mean_aod) &&
      !base::is.infinite(stats$max_aod) &&
      stats$valid_pixels > 0) {
    return(stats)
  }

  base::message("No valid data found")
  return(NULL)
}

# Test image focusing on Yellowknife area during peak fire
test_date <- base::as.Date("2023-08-15")  # Day before evacuation
test_hour <- 18  # Mid-afternoon local time

# Get file
files_df <- get_goes_files(
  date = test_date,
  hour = test_hour,
  product = "ABI-L2-AODF",  # Full disk product
  satellite = "goes18"
)

latest_file <- files_df |>
  dplyr::arrange(dplyr::desc(LastModified)) |>
  dplyr::slice(1)

# Create URL and download to temp file
url <- construct_url(latest_file$Key, bucket = "noaa-goes18")
temp_file <- base::tempfile(fileext = ".nc")

utils::download.file(
  url = url,
  destfile = temp_file,
  mode = "wb",
  quiet = TRUE
)

# Read AOD data from downloaded file
aod_data <- stars::read_stars(
  temp_file,
  sub = "AOD"
)

# Clean up temp file
base::unlink(temp_file)

# Create buffer for visualization extent
yk_buffer_big <- sf::st_buffer(yk_boundary, 100000)  # 100km buffer for context

# Transform to GOES projection
yk_buffer_goes <- sf::st_transform(yk_buffer_big, sf::st_crs(aod_data))

# Crop to our area of interest
aod_local <- sf::st_crop(aod_data, yk_buffer_goes)

# Create visualization
ggplot2::ggplot() +
  stars::geom_stars(data = aod_local) +
  ggplot2::geom_sf(data = sf::st_transform(yk_boundary, sf::st_crs(aod_data)),
                   fill = NA,
                   color = "red",
                   linewidth = 1) +
  ggplot2::geom_sf(data = sf::st_transform(yk_buffer_big, sf::st_crs(aod_data)),
                   fill = NA,
                   color = "blue",
                   linetype = "dashed") +
  ggplot2::scale_fill_viridis_c(
    name = "AOD",
    na.value = NA,
    option = "inferno"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "GOES-18 AOD Data Around Yellowknife",
    subtitle = format(test_date, "%B %d, %Y %H:00 UTC")
  )
