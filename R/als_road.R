#' Align and measure a road from lidar data
#'
#' From a reference road (spatial line), extracts the line with a buffer from the point cloud and computes
#' the exact positioning of the road (realignment). Then, using new the accurate shape, computes
#' road metrics including its width, its drivable width, its sinuosity as well as its state in four
#' classes. The function \link{st_snap_lines} allows to post-process the output to fix minor inaccuracies
#'  and reconnect the roads that may no longer be connected because each road is processed independently.
#'
#' @param centerline a single linestring (sf format) used as reference to search and measure the road.
#' @param roads multiple lines (sf format) used as reference to search and measure the roads
#' @param ctg a non-normalized \link[lidR:LAScatalog-class]{LAScatalog} object from lidR package
#' @param dtm RasterLayer storing the DTM with a resolution of at least of 1 m. Can be computed
#' with \link[lidR:grid_terrain]{grid_terrain}. It can be missing if a conductivity layer is provided
#' @param conductivity RasterLayer storing the pre-computed conductivity. It can be NULL in this case
#' if will be computed on the fly but the layer can be pre-computed with \link{rasterize_conductivity}
#' @param param a list of many parameters. See \link{alsroads_default_parameters}.
#' @param water a set of spatial polygons (sf format) of water bodies. This is used to mask the water
#' bodies so they cannot be mistaken as a drivable surfaces. Not mandatory but can help. It also allows
#' to detect bridges above water.
#' @param ... unused

#' @return An sf object similar to the input with additional attributes and an updated geometry. If
#' the class is 3 or 4 the original geometry is preserved to prevent adding more error. The new attributes
#' are ROADWITH, DRIVABLEWIDTH, PERCABOVEROAD (percentage of points between 0.5 and 5 meter above the road)
#' SHOULDERS (average number of shoulders found), SINUOSITY, CONDUCTIVITY (conductivity per linear meters)
#' SCORE (a road state score) and CLASS (4 classes derived from the SCORE). See references
#'
#' @references Roussel, J.-R., Bourdon, J.-F., Morley, I. D., Coops, N. C., & Achim, A. (2022).
#' Correction , update , and enhancement of vectorial forestry road maps using ALS data a pathfinder
#' and seven metrics. International Journal of Applied Earth Observation and Geoinformation, 114(September),
#' 103020. https://doi.org/10.1016/j.jag.2022.103020
#' @export
#' @examples
#' library(lidR)
#' library(sf)
#' library(raster)
#'
#' dir  <- system.file("extdata", "", package="ALSroads")
#' road <- system.file("extdata", "j5gr_centerline_971487.gpkg", package="ALSroads")
#' dtm  <- system.file("extdata", "j5gr_dtm.tif", package="ALSroads")
#' ctg  <- readLAScatalog(dir)
#' road <- st_read(road, "original", quiet = TRUE)
#' dtm  <- raster(dtm)
#'
#' # Voluntarily add more error to the road
#' crs <- st_crs(road)
#' st_geometry(road) <- st_geometry(road) + st_sfc(st_point(c(-8, 0)))
#' st_crs(road) <- crs
#'
#' plot(dtm, col = gray(1:50/50))
#' plot(ctg, add = TRUE)
#' plot(st_geometry(road), add = TRUE, col = "red")
#'
#' res <- measure_road(ctg, road, dtm = dtm)
#' res
#' poly <- sf::st_buffer(res, res$ROADWIDTH/2)
#'
#' plot(dtm, col = gray(1:50/50))
#' plot(st_geometry(road), col = "red", add = TRUE) # Inaccurate road track
#' plot(st_geometry(res), col = "blue", add = TRUE) # Corrected road track
#'
#' domain <- "https://servicesmatriciels.mern.gouv.qc.ca:443"
#' path <- "/erdas-iws/ogc/wmts/Inventaire_Ecoforestier/Inventaire_Ecoforestier/default/"
#' tiles <- "GoogleMapsCompatibleExt2:epsg:3857/{z}/{y}/{x}.jpg"
#' url <- paste0(domain, path, tiles)
#' m = mapview::mapview(list(road, poly),
#'   layer.name = c("Inaccurate", "Corrected"),
#'   color = c("red", "blue"), map.type = "Esri.WorldImagery")
#' leaflet::addTiles(m@map, url)
#'
#' \dontrun{
#' conductivity <- system.file("extdata", "j5gr_conductivity.tif", package="ALSroads")
#' conductivity <- raster(conductivity)
#' plot(conductivity, col = viridis::viridis(50))
#'
#' res <- measure_road(ctg, road, conductivity = conductivity)
#'
#' plot(st_geometry(road), col = "red") # Inaccurate road track
#' plot(st_geometry(res), col = "blue", add = TRUE) # Corrected road track
#' }
#' @useDynLib ALSroads, .registration = TRUE
#' @import data.table
measure_road = function(ctg, centerline, dtm = NULL, conductivity = NULL, water = NULL, param = alsroads_default_parameters,
                        return_all = NULL, return_stack = NULL, dl_model = NULL, ...)
{
  # Plenty of checks before to run anything
  dots <- list(...)

  if (!is.null(water)) {
    water <- st_transform(water, st_crs(centerline))
  }

  lidR::opt_progress(ctg) <- getOption("ALSroads.debug.verbose")
  geometry_type_road <- sf::st_geometry_type(centerline)
  if (geometry_type_road != "LINESTRING") stop(glue::glue("Expecting LINESTRING geometry for 'centerline' but found {geometry_type_road} geometry instead."), call. = FALSE)
  if (nrow(centerline) > 1) stop("Expecting a single LINESTRING", call. = FALSE)
  if (!methods::is(ctg, "LAScatalog")) stop("Expecting a LAScatalog", call. = FALSE)
  if (!is.null(water)) { if (any(!sf::st_geometry_type(water) %in% c("MULTIPOLYGON", "POLYGON"))) stop("Expecting POLYGON geometry type for 'water'", call. = FALSE) }
  if (sf::st_is_longlat(centerline)) stop("Expecting a projected CRS for 'centerline' but found geographic CRS instead.", call. = FALSE)
  if (param[["extraction"]][["road_max_width"]] > param[["extraction"]][["road_buffer"]]) stop("'road_max_width' parameter must be smaller than 'road_buffer' parameter", call. = FALSE)
  if (is.null(dtm) & is.null(conductivity)) stop("'dtm' and 'conductivity' cannot be both NULL", call. = FALSE)
  #if (!is.null(dtm) & !is.null(conductivity)) stop("'dtm' or 'conductivity' must be NULL", call. = FALSE)
  crs <- sf::st_crs(centerline)
  # If the collection is not indexed we throw a warning and even an error if the density is high
  # Without spatial indexation the centerline extraction is terribly long
  # (No warning when called from measure_roads because warning are thrown once only by measure_roads)
  if (!isFALSE(dots$Windex)) alert_no_index(ctg)

  # Display progress
  if (getOption("ALSroads.debug.progress")) cat("Progress: ")

  # Check for weird roads. We already found a road with a 90 degrees node
  # at the end of the road and it broke the output
  warn_weird_road(centerline)

  # We need to handle loop in special way so we define a bool
  is_loop <- st_is_loop(centerline)

  # If the confidence on the location of road is 1 we do not need to relocate
  # We know that the input is already perfect and we only need to measure the road
  relocate <- param[["constraint"]][["confidence"]] < 1

  # Get the units to display informative messages
  dist_unit  <- crs$units

  # We generate a default output in case we should exit early the function.
  # Basically a road with the original geometry and NA metrics
  new_road <- road_class0(centerline)

  # Check the distance between start and end points of the road. If they are too close (e.g. closer than
  # the size of the buffer) we need to reduce the size of the buffer otherwise we are likely to do
  # something wrong.
  if (!is_loop)
  {
    bu <- param[["extraction"]][["road_buffer"]]
    p1 <- lwgeom::st_startpoint(centerline)
    p2 <- lwgeom::st_endpoint(centerline)
    d  <- as.numeric(sf::st_distance(p1, p2))
    if (d < 2 * param[["extraction"]][["road_buffer"]])
      param[["extraction"]][["road_buffer"]] <- (d / (2 * bu + 5) * 2 * bu)/2
  }

  # Return the original geometry without computing anything if the centerline is too short to compute anything
  length_min <- 4*param[["extraction"]][["section_length"]]
  len <- as.numeric(sf::st_length(centerline))
  if (len < length_min)
  {
    warning(glue::glue("Road too short (< {length_min} {dist_unit}) to compute anything. Original road returned."), call. = FALSE)
    verbose("Done\n") ; cat("\n")
    return(new_road)
  }

  # Estimate if the road must be split in subsections. This have to goals
  # - It is an optimization to avoid processing giant conductivity raster from long road
  # - It allows to handle loop roads cases
  cut <- floor(len/param[["extraction"]][["road_max_len"]])

  if (cut > 0) { message(sprintf("Long road detected. Splitting the roads in %d chunks of %d %s to process.", cut+1, round(len/(cut+1)), dist_unit)) }
  if (cut == 0 & is_loop) { message(sprintf("Loop detected. Splitting the roads in 2 chunks of %d %s to process.", round(len/2,0), dist_unit)) ; cut = 1 }
  # If a split is required do the cut
  if (cut > 0)
  {
    # If we need to cut the road, the road is spitted and we recursively send each
    # piece recursively into this function
    cuts  <- seq(0,1, length.out = cut+2)
    from  <- cuts[-length(cuts)]
    to    <- cuts[-1]
    roads <- lapply(seq_along(from), function(i) { lwgeom::st_linesubstring(centerline, from[i], to[i]) })
    roads <- do.call(rbind, roads)

    res   <- measure_roads(ctg, roads, dtm, conductivity=conductivity, water=water, param=param, return_all=return_all , return_stack=return_stack, dl_model=dl_model)
    # geom  <- st_merge_line(res)
    #
    # # Because we split the line we have multiple independent results
    # # We recombine the multiple output in a single one
    # new_road$ROADWIDTH     <- mean(res$ROADWIDTH)
    # new_road$DRIVABLEWIDTH <- mean(res$ROADWIDTH)
    # new_road$PERCABOVEROAD <- mean(res$PERCABOVEROAD)
    # new_road$SHOULDERS     <- mean(res$SHOULDERS)
    # new_road$SINUOSITY     <- sinuosity(new_road)
    # new_road$ROADWIDTH     <- mean(res$ROADWIDTH)
    # new_road$CONDUCTIVITY  <- mean(res$CONDUCTIVITY)
    # new_road$SCORE         <- road_score(new_road, param)
    # new_road$CLASS         <- get_class(new_road$SCORE)
    # if (!is.na(new_road$CLASS) && new_road$CLASS < 4) sf::st_geometry(new_road) <- geom
    # verbose("Done\n") ; cat("\n")
    #
    # # We exit the function. The next code being the regular case we no splitting
    # new_road <- rename_sf_column(new_road, centerline)
    new_road = res
    return(new_road)
  }
  # Query the roads in a collection of files
  las <- extract_road(ctg, centerline, param)

  # Exit early. This should never happen
  if (lidR::is.empty(las))
  {
    warning("No point found.", call. = FALSE)
    verbose("Done\n") ; cat("\n")
    return(new_road)
  }

  # If we can't assume that the road is correctly positioned we need to recompute its
  # location accurately
  if (relocate)
  {
    # This is maybe the most important function of the code. It takes the point cloud, the original
    # road, the dtm or conductivity layer, water bodies and param to draw a new accurate line
    dots = list(...)
    dots$return_all=return_all
    dots$return_stack=return_stack

    result <- least_cost_path(las, centerline, dtm, conductivity, water, param,dl_model, dots)
    res1     <- result[[1]]
    layers_lidar    <- result[[2]]
    if (!is.null(res1$path)) {
      res      <- res1$path
      vertices <- res1$vertices
      res$CONDUCTIVITY = mean(vertices$CONDUCTIVITY)
    } else {
      res <- res1
    }
    if (sf::st_geometry_type(res) == "MULTILINESTRING")
      stop("Internal error. A MULTILINESTRING has been returned. Please report", call. = FALSE)

    # If the conductivity of the result is 0 it means that the path finder was not able to reach the
    # point B. We assume the road does not exist.
    if (res$CONDUCTIVITY == 0)
    {
      warning("Impossible to travel to the end of the road. Road does not exist.", call. = FALSE)
      new_road <- road_class4(centerline)
      verbose("Done\n") ; cat("\n")
      return(new_road)
    }

    # The new centerline geometry is very short? It is likely a bug that may arise for curved
    # and super short roads. The case is handled to avoid failure but in practice it should
    # not happen for regular road. It is an exception that must be handled
    if (as.numeric(sf::st_length(res)) < length_min)
    {
      warning(glue::glue("The computed road is too short (< {length_min} {dist_unit}) to compute anything. Original road returned."), call. = FALSE)
      new_road <- road_class0(centerline)
      new_road$CLASS <- 4
      verbose("Done\n") ; cat("\n")
      return(new_road)
    }

    # Eventually, if the result exist and is not too short the new road is the
    # one returned by least_cost_path
    new_road <- res
  }
  else
  {
    new_road$CONDUCTIVITY <- 1
  }
  # We now have an accurate road (hopefully). We can make measurement on it.
  # This step extracts the width profiles of the road, the percentage of points
  # and relocate more accurately the centerline.
  if (!is.null(param$segment_class) &&
      as.numeric(sf::st_length(new_road)) >= as.numeric(param$segment_class))
  {
    segment_length <- param$segment_class
    total_len <- as.numeric(sf::st_length(new_road))
    seg_len   <- as.numeric(segment_length)
    ratio     <- seg_len / total_len

    starts <- seq(0, 1 - ratio, by = ratio)
    ends   <- seq(ratio, 1, by = ratio)
    ends[length(ends)] <- 1

    n      <- min(length(starts), length(ends))
    starts <- starts[1:n]
    ends   <- ends[1:n]

    segments_list <- mapply(function(f, t) {
      lwgeom::st_linesubstring(new_road, from = f, to = t)
    }, starts, ends, SIMPLIFY = FALSE)

    road_segments <- do.call(rbind, segments_list)
    results       <- vector("list", length = nrow(road_segments))

    for (i in seq_len(nrow(road_segments)))
    {
      road_i <- road_segments[i, ]

      endpoints <- sf::st_cast(road_i, "POINT")
      idx       <- sf::st_nearest_feature(endpoints, vertices)

      road_i$CONDUCTIVITY <- mean(
        vertices$CONDUCTIVITY[idx],
        na.rm = TRUE
      )

      slice_metrics <- road_measure(las, road_i, param)
      metrics       <- road_metrics(road_i, slice_metrics)
      metrics[["SCORE"]] <- road_score(metrics, param)
      metrics[["CLASS"]] <- get_class(metrics[["SCORE"]])

      # Add a segment index so you can trace which segment is which
      metrics[["SEG_ID"]] <- i

      ngeom         <- attr(road_i, "sf_column")
      new_geometry  <- sf::st_geometry(road_i)

      attribute_table <- sf::st_drop_geometry(centerline)
      if ("length_m" %in% names(attribute_table)) {
        attribute_table["length_m"] <- segment_length
      }
      attribute_table <- cbind(attribute_table, metrics)
      attribute_table[[ngeom]] <- new_geometry

      road_i <- sf::st_as_sf(attribute_table)
      sf::st_crs(road_i) <- crs

      results[[i]] <- cbind(road_i, metrics)
    }

    # --- Simply bind all segments — no grouping, no averaging ---
    new_road <- do.call(rbind, results)

    # Optional: create a unique ID per segment if OGF_ID or ROADID_unique exists
    if ("OGF_ID" %in% names(new_road)) {
      new_road$OGF_ID <- paste0(new_road$OGF_ID, "_", new_road$SEG_ID)
    }
    if ("ROADID_unique" %in% names(new_road)) {
      new_road$ROADID_unique <- paste0(new_road$ROADID_unique, "_", new_road$SEG_ID)
    }
  } else {

    slice_metrics <- road_measure(las, new_road, param)

    # Smooth the centerline using a spline adjustment
    if (relocate && nrow(slice_metrics) > 4L)
    {
      spline <- adjust_spline(slice_metrics)
      spline <- sf::st_simplify(spline, dTolerance = 1)
      spline <- sf::st_set_crs(spline, crs)
      sf::st_geometry(new_road) <- spline
    }


    metrics <- road_metrics(new_road, slice_metrics)
    class_result   <- road_class_full(new_road, metrics, slice_metrics, layers_lidar, param)
    metrics[["CLASS"]] = class_result$class
    metrics[["SCORE"]] = class_result$score
    #
    # metrics[["SCORE"]] <- road_score(metrics, param)
    # metrics[["CLASS"]] <- get_class(metrics[["SCORE"]])
  }
  is_single <- nrow(new_road) == 1L
  if (is_single)
  {
    # Merge the tables of attributes
    ngeom <- attr(new_road, "sf_column")
    new_geometry <- sf::st_geometry(new_road)

    attribute_table <- sf::st_drop_geometry(centerline)
    attribute_table <- cbind(attribute_table, metrics)
    attribute_table[[ngeom]] <- new_geometry


    new_road <- sf::st_as_sf(attribute_table)

    sf::st_crs(new_road) <- crs
    keep_fields <- c(
      "geometry",
      "ROADWIDTH",
      "DRIVABLEWIDTH",
      "PERCABOVEROAD",
      "SHOULDERS",
      "SINUOSITY",
      "CONDUCTIVITY",
      "OGF_ID",
      "SCORE",
      "CLASS"
    )

    if ("length_m" %in% names(new_road)) {
      keep_fields <- c(keep_fields, "length_m")
    }

    if ("ROADID_unique" %in% names(new_road)) {
      keep_fields <- c(keep_fields, "ROADID_unique")
    }

    new_road <- dplyr::select(new_road, dplyr::any_of(keep_fields))
  }

  # Hidden option for JF Bourdons
  keep_class = dots$keep_class
  if (is.null(keep_class))
    keep_class = 4L
  else
    stopifnot(is.numeric(keep_class), length(keep_class) == 1)

  mean_class <- mean(new_road$CLASS, na.rm = TRUE)
  if (mean_class > keep_class)
  {
    new_road <- centerline
    sf::st_crs(new_road) <- crs

    verbose("Done (reverted to centerline)\n")
    new_road <- rename_sf_column(new_road, centerline)
    return(new_road)
  }

  verbose("Done\n") ; cat("\n")
  new_road <- rename_sf_column(new_road, centerline)
  return(new_road)

}


#' @export
#' @rdname measure_road
measure_roads = function(ctg, roads, dtm, conductivity = NULL, water = NULL, param = alsroads_default_parameters,  return_all = NULL, return_stack = NULL, dl_model = NULL, ...)
{

  alert_no_index(ctg)
  i <- 1:nrow(roads)
  res <- lapply(i, function(j)
  {
    if (getOption("ALSroads.debug.verbose") | getOption("ALSroads.debug.progress")) cat("Road", j, "of", nrow(roads), " ")
    print(j)
    tryCatch(
    {

      withCallingHandlers(

        measure_road(ctg, roads[j,], dtm, conductivity, water, param,
                     return_all = return_all, return_stack = return_stack, dl_model = dl_model,Windex = FALSE, ...),
        warning = function(w) {
          warning(paste0("Warning in road ", j, ": ", conditionMessage(w)))
          invokeRestart("muffleWarning") # Suppress warning after handling
        }
      )
    },
    error = function(e)
    {
      warning(paste0("Error in road ", i, ": NULL returned.\nThe error was: ", e))
      return(NULL)
    })
  })
  # do.call(rbind, res)
  bind_rows(res)
}

measure_roads_par <- function(ctg, roads, dtm, conductivity = NULL, water = NULL, param = alsroads_default_parameters) {
  alert_no_index(ctg)
  library(parallel)
  # Prepare inputs for parallel processing
  inputs <- lapply(1:nrow(roads), function(i) {
    list(
      ctg = ctg,
      road = roads[i, ],
      dtm = dtm,
      conductivity = conductivity,
      water = water,
      param = param
    )
  })

  cl <- parallel::makeCluster(parallel::detectCores()-20)
  on.exit(parallel::stopCluster(cl))  # ensure the cluster shuts down

  clusterEvalQ(cl, {
    library(lidR)
    library(sf)
    library(raster)
    library(mapview)
    library(leaflet)
    library(ggplot2)
    library(leaflet.esri)
    library(ALSroads)
  })

  # Export the function you are using in parLapply
  parallel::clusterExport(cl, varlist = c("measure_road_par"))

  # Run in parallel
  results <- parallel::parLapply(cl, inputs, function(args) do.call(measure_road_par, args))

  return(results)
}

alert_no_index <- function(ctg)
{
  is_copc = substr(ctg$filename, nchar(ctg$filename)-8, nchar(ctg$filename))
  if (all(is_copc == ".copc.laz"))
  {
    if (utils::packageVersion("rlas") >= "1.7.0")
      return(invisible())
    else
      message(paste0("copc files are supported using package rlas >= 1.7.0. Currently installed: ", utils::packageVersion("rlas")))
  }

  if (!lidR::is.indexed(ctg))
  {
    d <- lidR::density(ctg)
    if (d < 5)
    {
      message("No spatial index for LAS/LAZ files in this collection.")
      return(invisible())
    }
    else if (d < 10)
    {
      warning("No spatial index for LAS/LAZ files in this collection.", call. = FALSE)
      return(invisible())
    }
    else
    {
      stop("No spatial index for LAS/LAZ files in this collection.")
    }
  }
}

road_class0 <- function(centerline)
{
  new_road <- centerline
  new_road$ROADWIDTH     <- NA
  new_road$DRIVABLEWIDTH <- NA
  new_road$PERCABOVEROAD <- NA
  new_road$SHOULDERS     <- NA
  new_road$SINUOSITY     <- NA
  new_road$CONDUCTIVITY  <- NA
  new_road$SCORE         <- NA
  new_road$CLASS         <- 0

  # reorder the columns so outputs are consistent even if exiting early
  ngeom <- attr(new_road, "sf_column")
  names <- names(new_road)
  names <- names[names != ngeom]
  names <- append(names, ngeom)
  data.table::setcolorder(new_road, names)

  new_road <- rename_sf_column(new_road, centerline)

  return(new_road)
}

road_class4 <- function(centerline)
{
  new_road <- centerline
  new_road$ROADWIDTH     <- 0
  new_road$DRIVABLEWIDTH <- 0
  new_road$PERCABOVEROAD <- 100
  new_road$SHOULDERS     <- 0
  new_road$SINUOSITY     <- NA
  new_road$CONDUCTIVITY  <- 0
  new_road$SCORE         <- 0
  new_road$CLASS         <- 4

  # reorder the columns so outputs are consistent even if exiting early
  ngeom <- attr(new_road, "sf_column")
  names <- names(new_road)
  names <- names[names != ngeom]
  names <- append(names, ngeom)
  data.table::setcolorder(new_road, names)

  new_road <- rename_sf_column(new_road, centerline)

  return(new_road)
}

warn_weird_road <- function(centerline)
{
  angles <- st_angles(centerline)
  if (any(angles > 90))
  {
    if (any(angles[c(1, length(angles))] > 90))
      warning("Sharp turn (< 90 degrees) at one or both ends of the input road. This is weird and may lead to invalid outputs.", call. = FALSE)
    else
      warning("Sharp turn (< 90 degrees) between two consecutive segments of the input road. This is weird and may lead to invalid outputs.", call. = FALSE)
  }
}

rename_sf_column <- function(x,as)
{
  # Ensure the sf_colum is the same than the input
  current <- attr(x, "sf_column")
  name    <- attr(as, "sf_column")
  names(x)[names(x) == current] = name
  sf::st_geometry(x) = name
  x
}

