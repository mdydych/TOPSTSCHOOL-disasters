# Load required packages for this section
library(earthdatalogin)
library(sf)
library(ggplot2)

# Set up NASA Earthdata Login authentication
# This only needs to be done once per R session
# Assumes you have already run: edulog_setup()

# environment creds----
Sys.setenv(
  EARTHDATA_USER = "jbrinks@isciences.com",
  EARTHDATA_PASSWORD = "fU.Xe_T56que.h_",
  FIRMS_API_KEY = "9a7258f5b04a82bccefd87cdeb2faefb"
)

earthdatalogin::edl_netrc(
  username = "jbrinks@isciences.com",
  password = "fU.Xe_T56que.h_")

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
  library(sf)
  library(dplyr)
  library(concaveman)
  library(units)
  library(ggplot2)
  library(knitr)
  library(kableExtra)

  # Create perimeters with adjusted parameters
  create_fire_perimeters <- function(points, date, concavity = 1, threshold = 1) {
    daily_points <- points |>
      dplyr::filter(date == !!date)

    if (base::nrow(daily_points) < 3) {
      base::return(NULL)
    }

    coords <- sf::st_coordinates(daily_points)
    hull <- concaveman::concaveman(coords,
                                   concavity = concavity,
                                   length_threshold = threshold)

    polygon_sf <- sf::st_polygon(list(hull)) |>
      sf::st_sfc(crs = sf::st_crs(points)) |>
      sf::st_sf()

    perimeter <- polygon_sf |>
      dplyr::mutate(
        date = date,
        n_points = base::nrow(daily_points),
        mean_frp = base::mean(daily_points$frp),
        area_km2 = base::as.numeric(units::set_units(sf::st_area(polygon_sf), "km^2"))
      )

    base::return(perimeter)
  }

  # Create daily perimeters
  unique_dates <- base::unique(firms_points$date)
  fire_perimeters <-
    base::do.call(rbind,
                  base::lapply(unique_dates,
                               function(d) create_fire_perimeters(firms_points, d, 1, 0.01)))

  # Find date with maximum area and maximum points
  max_stats <- base::list(
    area_date = fire_perimeters |>
      dplyr::arrange(dplyr::desc(area_km2)) |>
      dplyr::slice(1) |>
      dplyr::pull(date),
    points_date = fire_perimeters |>
      dplyr::arrange(dplyr::desc(n_points)) |>
      dplyr::slice(1) |>
      dplyr::pull(date)
  )

  # Create dataset for faceted plot
  comparison_perimeters <- fire_perimeters |>
    dplyr::filter(date %in% c(max_stats$area_date, max_stats$points_date)) |>
    dplyr::mutate(
      type = dplyr::case_when(
        date == max_stats$area_date ~ "Maximum Hull Area",
        date == max_stats$points_date ~ "Maximum Detection Count"
      )
    )

  # Get corresponding VIIRS points for these dates
  comparison_points <- firms_points |>
    dplyr::filter(date %in% c(max_stats$area_date, max_stats$points_date)) |>
    dplyr::mutate(
      type = dplyr::case_when(
        date == max_stats$area_date ~ "Maximum Hull Area",
        date == max_stats$points_date ~ "Maximum Detection Count"
      )
    )

  # Create faceted comparison plot
  p4b <- ggplot2::ggplot() +
    # Base layers
    ggplot2::geom_sf(data = search_area,
                     fill = NA,
                     color = "gray60",
                     linetype = "dashed") +
    ggplot2::geom_sf(data = yk_boundary,
                     fill = NA,
                     color = "#00BFFF",
                     linewidth = 1.2) +
    # Add hull polygons
    ggplot2::geom_sf(data = comparison_perimeters,
                     fill = "red",
                     alpha = 0.3) +
    # Add VIIRS points
    ggplot2::geom_sf(data = comparison_points,
                     ggplot2::aes(color = frp),
                     size = 2,
                     alpha = 0.7) +
    # Facet by type
    ggplot2::facet_wrap(~type) +
    # Style
    ggplot2::scale_color_viridis_c(
      name = "Fire Radiative\nPower (MW)",
      option = "inferno"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Comparison of Maximum Area vs. Maximum Detections",
      subtitle = base::paste("Hull Area Date:", format(max_stats$area_date, "%B %d"),
                             "  |  Detection Count Date:", format(max_stats$points_date, "%B %d")),
      caption = "Note: Convex hulls may overestimate fire extent by including unburned areas between detection points"
    )

  # Display plot
  base::print(p4b)

  # Print comparison statistics
  base::cat("\nComparison Statistics:\n")
  base::cat("\nMaximum Hull Area Date:", format(max_stats$area_date, "%B %d, %Y"), "\n")
  base::cat("- Area:", round(comparison_perimeters$area_km2[1], 2), "km²\n")
  base::cat("- Number of detections:", comparison_perimeters$n_points[1], "\n")
  base::cat("- Mean FRP:", round(comparison_perimeters$mean_frp[1], 2), "MW\n")

  base::cat("\nMaximum Detection Count Date:", format(max_stats$points_date, "%B %d, %Y"), "\n")
  base::cat("- Area:", round(comparison_perimeters$area_km2[2], 2), "km²\n")
  base::cat("- Number of detections:", comparison_perimeters$n_points[2], "\n")
  base::cat("- Mean FRP:", round(comparison_perimeters$mean_frp[2], 2), "MW\n")
