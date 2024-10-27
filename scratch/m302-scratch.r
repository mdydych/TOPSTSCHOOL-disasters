# Load required packages for this section
library(earthdatalogin)
library(sf)
library(ggplot2)

# Set up NASA Earthdata Login authentication
# This only needs to be done once per R session
# Assumes you have already run: edulog_setup()

# fire progress animate----

  # Load additional packages
  library(sf)
  library(ggplot2)
  library(httr2)
  library(gganimate)
  library(tidyterra)  # For working with raster/terra objects
  library(maptiles)   # For downloading satellite imagery

  # Create search area and get bounding box
  search_area <- sf::st_buffer(yk_boundary, 30000)
  search_area_wgs84 <- sf::st_transform(search_area, 4326)
  bbox <- sf::st_bbox(search_area_wgs84) |>
    base::round(4)

  # Download satellite imagery for our area
  # Using 'esri.worldimagery' for satellite imagery
  basemap_tiles <- maptiles::get_tiles(
    search_area_wgs84,
    provider = "Esri.WorldImagery",
    zoom = 11  # Adjust zoom level as needed
  )

  # Get FIRMS data (same as before)
  firms_url <- base::sprintf(
    "https://firms.modaps.eosdis.nasa.gov/api/area/csv/%s/VIIRS_NOAA20_NRT/%s,%s,%s,%s/10/2023-08-10",
    Sys.getenv("FIRMS_API_KEY"),
    bbox["xmin"], bbox["ymin"],
    bbox["xmax"], bbox["ymax"]
  )

  resp <- httr2::request(firms_url) |>
    httr2::req_perform()

  if (httr2::resp_status(resp) == 200) {
    # Process FIRMS data
    firms_points <- httr2::resp_body_string(resp) |>
      textConnection() |>
      utils::read.csv() |>
      sf::st_as_sf(
        coords = c("longitude", "latitude"),
        crs = 4326
      ) |>
      sf::st_transform(sf::st_crs(yk_boundary))

    # Create datetime field
    firms_points$datetime <- as.POSIXct(
      paste(firms_points$acq_date,
            sprintf("%04d", firms_points$acq_time)),
      format = "%Y-%m-%d %H%M",
      tz = "UTC"
    )

    # Create custom color palette for fire intensity
    fire_colors <- c(
      "red",         # Medium
      "darkred",    # Medium-high
      "purple"       # High intensity
    )

    # Create animated plot
    anim <- ggplot2::ggplot() +
      # Add satellite basemap
      tidyterra::geom_spatraster_rgb(data = basemap_tiles) +
      # Add buffer zone
      ggplot2::geom_sf(data = search_area,
                       fill = NA,
                       color = "white",
                       linetype = "dashed",
                       alpha = 0.7) +
      # Add city boundary
      ggplot2::geom_sf(data = yk_boundary,
                       fill = NA,
                       color = "#00BFFF",  # Deep sky blue
                       linewidth = 1.2) +
      # Add fire points
      ggplot2::geom_sf(data = firms_points,
                       ggplot2::aes(color = frp),
                       alpha = 0.8,
                       size = 2) +
      # Custom color scale
      ggplot2::scale_color_gradientn(
        colors = c(
          "red",         # Medium
          "darkred",    # Medium-high
          "purple"       # High intensity
        ),
        name = "Fire Radiative\nPower (MW)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14,
                                           face = "bold",
                                           color = "grey25"),
        plot.subtitle = ggplot2::element_text(size = 12,
                                              color = "grey25"),
        plot.caption = ggplot2::element_text(color = "grey25"),
        legend.text = ggplot2::element_text(color = "grey25"),
        legend.title = ggplot2::element_text(color = "grey25"),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        title = "Yellowknife Fire Progression",
        subtitle = "Date: {format(frame_time, '%B %d, %Y %H:%M UTC')}",
        caption = "Data: VIIRS 375m from NOAA-20 satellite | Imagery: ESRI World Imagery"
      ) +
      # Animation settings
      gganimate::transition_time(datetime) +
      gganimate::shadow_wake(wake_length = 0.2, alpha = 0.2)

    # Render animation
    gganimate::animate(
      anim,
      nframes = 200,
      fps = 8,
      width = 800,
      height = 800,
      renderer = gganimate::gifski_renderer()
    )
  }



# hotspot analysis----
  # Load required packages

# GOES ----

