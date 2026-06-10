#' Extracts sequentially slices perpendicularly to the road
#'
#' Extracts sequentially slices perpendicularly to the road and makes operation on each slice
#' to return profiles along the road
#'
#' @noRd
road_measure = function(las, road, param)
{
  if (is.null(las@index[["quadtree"]]))  stop("The point cloud is not spatially indexed", call. = TRUE)
  # Retrieve start and end points of the road
  coords_road <- sf::st_coordinates(road)
  start <- coords_road[1,1:2]
  end <- coords_road[nrow(coords_road),1:2]

  # Split the path in n sections of ~ same length
  path_lenght <- as.numeric(sf::st_length(road))
  each <- param[["extraction"]][["section_length"]]/path_lenght
  at <- round(seq(0,1, each),4)
  at[length(at)] = 1
  points <- sf::st_line_sample(road, sample = at)
  points <- sf::st_cast(points, "POINT")
  points <- sf::st_as_sf(points)

  # Retrieve the "temporal" position along the line so we can order the points later
  CC  <- sf::st_coordinates(points)
  CCX <- CC[,1]
  CCY <- CC[,2]
  dx  <- c(CCX, end[1]) - c(start[1], CCX)
  dy  <- c(CCY, end[2]) - c(start[2], CCY)
  dd  <- sqrt(dx*dx+dy*dy)
  dist <- cumsum(dd[-length(dd)])
  pdist <- dist/dist[length(dist)]
  points$DISTANCE <- path_lenght*pdist

  # We can now loop through each segment between two consecutive points
  n <- nrow(CC)
  ids <- 1:(n-1)
  .segment_metrics <- vector("list", n)
  print_at = round(seq(1, n, length.out = 20)) # prints every ~5%
  for (i in ids)
  {
    if (i %in% print_at)
    {
      verbose("Computing road metrics... ", round(i/length(ids)*100,0), "%\r", sep = "")
      utils::flush.console()
    }

    # Extraction of the segment
    p1 <- CC[i,1:2]
    p2 <- CC[i+1,1:2]
    dx <- p1[1] - p2[1]
    dy <- p1[2] - p2[2]
    angle  <- atan(dy/dx)
    if (angle < 0) angle <- angle + pi
    center <- (p1+p2)/2
    height <- sqrt(sum((p1-p2)^2))
    width <- param[["extraction"]][["road_max_width"]]
    xmin <- center[1] - width/2
    xmax <- center[1] + width/2
    ymin <- center[2] - height/2
    ymax <- center[2] + height/2
    las_slice <- clip_orectangle_with_index(las, xmin, ymin, xmax, ymax, angle-pi/2)

    # Normalize the segment. In practice both normalized and raw point cloud will be used
    # but we can unnormalized on-the-fly and thus later we will pass only the normalized version
    # /!\ Normalization can fail because we never know what point cloud we get for a specific
    # segment. No ground point, no point, not enough ground points and so on. Few error
    # along the road are allowed with the try-catch block
    nlas_slice <- tryCatch(
    {
      lidR::normalize_height(las_slice, lidR::tin(), Wdegenerated = FALSE, na.rm = TRUE)
    },
    error = function(e)
    {
      if (isTRUE(getOption("ALSroads.debug")))
      {
        f <- tempfile(fileext = ".las")
        lidR::writeLAS(las_slice, f)
        message(glue::glue("Normalization impossible in segment {i}. Segment skiped with error : {e} "))
        message(glue::glue("The LAS objects that caused the failure has been saved in {f}."))
      }
      return(NULL)
    })

    if (is.null(nlas_slice)) next

    # We need to rotate the point cloud along a single axis maintaining a valid LAS object.
    # This allow to work in 2D (see figure 7)
    nlas_slice <- las_rotate(nlas_slice, angle, center)

    # We now have a 2D-like points cloud of a slice perpendicular to the road. We can perform measurements
    # such as road width, percentage of points and so on...
    # /!\ Here again there are many rare cases were it could fail such as water only + bridge,
    # missing points. Hard to find each specific case. Few error along the road are allowed
    # with the try-catch block
    m <- tryCatch(
    {
      slice_metrics(nlas_slice, param)
    },
    error = function(e)
    {
      if (isTRUE(getOption("ALSroads.debug")))
      {
        f <- tempfile(fileext = ".las")
        nlas_slice <- lidR::add_lasattribute(nlas_slice, name = "Zref", desc = "Absolute Elvation")
        lidR::writeLAS(nlas_slice, f)
        cat("\n")
        cat(glue::glue("Computation impossible in segment {i}. segment_road_metrics() failed with error : {e}\n"))
        cat(glue::glue("The LAS objects that caused the failure has been saved in {f}.\n"))
        cat("\n")
      }
      return(NULL)
    },
    warning = function(e)
    {
      if (isTRUE(getOption("ALSroads.debug")))
      {
        f <- tempfile(fileext = ".las")
        nlas_slice <- lidR::add_lasattribute(nlas_slice, name = "Zref", desc = "Absolute Elvation")
        lidR::writeLAS(nlas_slice, f)
        cat("\n")
        cat(glue::glue("Computation warning in segment {i}. segment_road_metrics() failed with error : {e}\n"))
        cat(glue::glue("The LAS objects that caused the failure has been saved in {f}.\n"))
        cat("\n")
      }
    })

    if (is.null(m)) next

    # We add some additional metrics to know the "temporal" location of this segment
    # and record the original position of the road + we updated previous locations
    m$distance_to_start <- points[["DISTANCE"]][i]
    .segment_metrics[[i]] <- m
  }

  verbose("Computing road metrics... 100%\n")

  segment_metrics <- data.table::rbindlist(.segment_metrics)
  segment_metrics <- sf::st_as_sf(segment_metrics, coords = c("xroad", "yroad"), crs = sf::st_crs(las))

  # Insert back end points because road_measure discards the ending of the road
  # This piece of code should go in road_measure I guess
  end <- segment_metrics[nrow(segment_metrics),]
  end$distance_to_start <- sf::st_length(road)
  end <- sf::st_set_geometry(end, lwgeom::st_endpoint(road))
  end <- sf::st_set_crs(end, sf::NA_crs_)
  end <- sf::st_set_crs(end, sf::st_crs(las))
  segment_metrics <- rbind(segment_metrics, end)
  return(segment_metrics)
}

#' Compute the "average" metrics for the road from all the metrics of each segment
#'
#' @param road a single spatial line in sf format
#' @param segment_metrics the metric of each segment as outputted by road_generic_loop(..., mode = "measure")
#'
#' @return An updated road with the metrics including the average state, width, drivable width, ...
#' @noRd
road_metrics =  function(road, segment_metrics)
{
  # Average percentage of point above x
  avg_percentage_above_05 <- round(mean(segment_metrics[["pzabove05"]], na.rm = TRUE), 1)
  avg_percentage_above_2 <- round(mean(segment_metrics[["pzabove2"]], na.rm = TRUE), 1)
  avg_shoulders <- round(mean(segment_metrics[["number_accotements"]], na.rm = TRUE)/2*100, 1)

  # Average width
  road_width <- round(mean(segment_metrics[["road_width"]], na.rm = TRUE), 1)

  # Average drivable width.
  drivable_width <- round(mean(segment_metrics[["drivable_width"]], na.rm = TRUE), 1)

  # Measure of the right of way not available yet
  right_of_way <- NA_real_

  # Sinuosity on XY and on Z
  z <- segment_metrics[["zroad"]]
  x <- segment_metrics[["distance_to_start"]]
  S <- round(sinuosity(road), 2)

  road_metrics = data.frame(
    ROADWIDTH = road_width,
    DRIVABLEWIDTH = drivable_width,
    #RIGHTOFWAY = right_of_way,
    PERCABOVEROAD = avg_percentage_above_05,
    SHOULDERS = avg_shoulders,
    SINUOSITY = S,
    CONDUCTIVITY = road$CONDUCTIVITY)

  if (getOption("ALSroads.debug.metrics")) plot_road_metrics(road, road_metrics, segment_metrics)

  return(road_metrics)
}

road_score = function(metrics, param)
{
  P     <- metrics$PERCABOVEROAD
  W     <- metrics$DRIVABLEWIDTH
  S     <- metrics$SHOULDERS
  Sigma <- metrics$CONDUCTIVITY

  verbose("P = ", P, "\n")
  verbose("W = ", W, "\n")
  verbose("S = ", S, "\n")
  verbose("Sigma = ", Sigma, "\n")

  # The overall state is the median state
  # Estimate the existence of a road in the segment based on percentage of veg points and width
  p <- param[["state"]][["percentage_veg_thresholds"]]
  road_exist <- activation(P, p, "piecewise-linear",  asc = FALSE)

  p <- param[["state"]][["drivable_width_thresholds"]]
  drivable_exist <- activation(W, p, "piecewise-linear",  asc = TRUE)

  p <- param[["state"]][["conductivity_thresholds"]]
  score_exist <- activation(Sigma, p, "piecewise-linear",  asc = TRUE)

  p <- param[["state"]][["shoulder_thresholds"]]
  shoulder_exist <- activation(S, p, "piecewise-linear",  asc = TRUE)

  pexist <- (road_exist + drivable_exist + score_exist + shoulder_exist) / 4*100

  verbose("Estimating the state of the road...\n")
  verbose("   - Estimated probability based on vegetation (P):", round(road_exist,1), "\n")
  verbose("   - Estimated probability based on road size (W):", round(drivable_exist,1), "\n")
  verbose("   - Estimated probability based on conductivity (sigma):", round(score_exist,1), "\n")
  verbose("   - Estimated probability based on road shoulders (S):", round(shoulder_exist,1), "\n")
  verbose("   - Estimated probability:", round(pexist,1), "\n")

  return(round(pexist,1))
}

get_class =  function(score)
{
  5 - as.integer(cut(score, breaks = c(-1,25,50,75,101), right = FALSE))
}

las_rotate <- function(las, angle, center)
{
  xfactor <- las@header@PHB[["X scale factor"]]
  yfactor <- las@header@PHB[["Y scale factor"]]
  rot <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), ncol = 2)
  coords <- as.matrix(lidR:::coordinates(las))
  coords[,1] <- coords[,1] - center[1]
  coords[,2] <- coords[,2] - center[2]
  coords <- coords %*% rot
  X <- coords[,1]
  Y <- coords[,2]
  lidR:::fast_quantization(X, xfactor, 0)
  lidR:::fast_quantization(Y, yfactor, 0)
  las@data[["X"]] <- X
  las@data[["Y"]] <- Y
  las@header@PHB[["X offset"]] <- 0
  las@header@PHB[["Y offset"]] <- 0
  las <- lidR::las_update(las)
  data.table::setattr(las, "rotation", rot)
  data.table::setattr(las, "offset", center)
  return(las)
}

road_class_full <- function(road, metrics, segment_metrics, layers_lidar, param) {

  # ── Step 1: extract raster values within actual road width corridor ──────
  corridor_vals <- extract_corridor_values(
    road         = road,
    road_width_m = metrics$ROADWIDTH,    # from road_metrics()
    layers_lidar = layers_lidar
  )

  # ── Step 2: Class 5 gate — must happen before anything else ─────────────
  ghost <- detect_ghost_road(corridor_vals, param)

  if (ghost$is_ghost) {
    return(data.frame(
      class           = 5L,
      class_name      = "Ghost / reclaimed",
      score           = 0,
      D1_surface      = NA, D2_width  = NA, D3_canopy = NA,
      D4_detect       = NA, D5_drain  = NA,
      ghost_reason    = ghost$reason,
      chm_median      = ghost$evidence["chm_median"],
      density_mean    = ghost$evidence["density_mean"],
      canopy_fraction = ghost$evidence["canopy_fraction"]
    ))
  }

  # ── Step 3: compute dimensions for Class 1-4 ────────────────────────────
  # Augment segment_metrics with per-segment raster values now that we know
  # the road exists and is detectable
  segment_metrics <- augment_with_layers(segment_metrics, layers_lidar)

  dims   <- compute_dimensions(metrics, segment_metrics,corridor_vals, param)

  result <- classify_road(dims, param)

  data.frame(
    class           = result$class,
    score           = round(result$score, 1),
    D1_surface      = round(result$dims["D1"], 1),
    D2_width        = round(result$dims["D2"], 1),
    D3_canopy       = round(result$dims["D3"], 1),
    D4_detect       = round(result$dims["D4"], 1),
    D5_drain        = round(result$dims["D5"], 1),
    ghost_reason    = NA_character_,
    chm_median      = NA_real_,
    density_mean    = NA_real_,
    canopy_fraction = NA_real_
  )

}
default_classification_param <- function() {
  list(
    # Class 5 hard-gate thresholds
    ghost_chm_threshold    = 3.0,    # metres; CHM median above this → ghost
    ghost_max_density      = 0.20,   # normalised density; below → ghost
    ghost_min_canopy_cover = 0.85,   # fraction; above → ghost

    # Class 1-4 thresholds (unchanged)
    min_half_width_m       = 1.5,
    lateral_buffer_m       = 3.0,
    class1_min_all_dims    = 70,
    class2_min_detect      = 60,
    class2_min_width       = 50,
    class2_min_any_dim     = 30,
    class3_min_detect      = 40,
    class4_min_detect      = 20
  )
}

detect_ghost_road <- function(corridor_vals, param) {

  p <- param[["classification"]] %||% default_classification_param()

  # Safe fallback if extraction failed
  if (is.null(corridor_vals) || nrow(corridor_vals) == 0) {
    return(list(
      is_ghost = TRUE,
      reason   = "No raster pixels found inside road corridor",
      evidence = c(chm_median = NA, density_mean = NA, canopy_fraction = NA)
    ))
  }

  chm     <- corridor_vals[["chm"]]
  density <- corridor_vals[["density"]]

  chm     <- chm[!is.na(chm)]
  density <- density[!is.na(density)]

  # ── Trigger 1: tall vegetation on road surface
  chm_median     <- if (length(chm) > 0) median(chm) else NA_real_
  ghost_chm      <- !is.na(chm_median) && chm_median > (p$ghost_chm_threshold %||% 3.0)

  # ── Trigger 2: ground point density too low to trust measurements
  density_mean   <- if (length(density) > 0) mean(density) else NA_real_
  ghost_density  <- !is.na(density_mean) && density_mean < (p$ghost_max_density %||% 0.20)

  # ── Trigger 3: canopy closure fraction > 85% of corridor pixels
  canopy_fraction <- if (length(chm) > 0) mean(chm > 2.0) else NA_real_
  ghost_canopy    <- !is.na(canopy_fraction) &&
                     canopy_fraction > (p$ghost_min_canopy_cover %||% 0.85)

  is_ghost <- ghost_chm | ghost_density | ghost_canopy

  reason <- dplyr::case_when(
    ghost_chm     ~ sprintf("CHM median %.1f m > threshold inside road corridor", chm_median),
    ghost_density ~ sprintf("Ground density %.3f below minimum detection threshold", density_mean),
    ghost_canopy  ~ sprintf("Canopy closure %.0f%% > 85%% of road corridor", canopy_fraction * 100),
    TRUE          ~ "Not a ghost road"
  )

  list(
    is_ghost = is_ghost,
    reason   = reason,
    evidence = c(
      chm_median      = round(chm_median,    3),
      density_mean    = round(density_mean,  3),
      canopy_fraction = round(canopy_fraction, 3)
    )
  )
}
augment_with_layers <- function(slice_metrics, layers_lidar) {

  if (is.null(layers_lidar)) return(slice_metrics)

  # Extract coordinates of each segment centroid
  coords <- sf::st_coordinates(slice_metrics)

  # Extract raster values at those points
  vals <- raster::extract(layers_lidar, coords[, c("X", "Y")])
  vals <- as.data.frame(vals)

  # Map layer names to the column names dim_* functions look for
  if ("density"   %in% names(vals)) slice_metrics$point_density <- vals$density
  if ("roughness" %in% names(vals)) slice_metrics$roughness     <- vals$roughness
  if ("chm"       %in% names(vals)) slice_metrics$chm_val       <- vals$chm
  if ("slope"     %in% names(vals)) slice_metrics$slope_val     <- vals$slope

  return(slice_metrics)
}
detect_ghost_road <- function(corridor_vals, param) {

  p <- param[["classification"]] %||% default_classification_param()

  # Safe fallback if extraction failed
  if (is.null(corridor_vals) || nrow(corridor_vals) == 0) {
    return(list(
      is_ghost = TRUE,
      reason   = "No raster pixels found inside road corridor",
      evidence = c(chm_median = NA, density_mean = NA, canopy_fraction = NA)
    ))
  }

  chm     <- corridor_vals[["chm"]]
  density <- corridor_vals[["density"]]

  chm     <- chm[!is.na(chm)]
  density <- density[!is.na(density)]

  # ── Trigger 1: tall vegetation on road surface
  chm_median     <- if (length(chm) > 0) median(chm) else NA_real_
  ghost_chm      <- !is.na(chm_median) && chm_median > (p$ghost_chm_threshold %||% 3.0)

  # ── Trigger 2: ground point density too low to trust measurements
  density_mean   <- if (length(density) > 0) mean(density) else NA_real_
  ghost_density  <- !is.na(density_mean) && density_mean < (p$ghost_max_density %||% 0.20)

  # ── Trigger 3: canopy closure fraction > 85% of corridor pixels
  canopy_fraction <- if (length(chm) > 0) mean(chm > 2.0) else NA_real_
  ghost_canopy    <- !is.na(canopy_fraction) &&
    canopy_fraction > (p$ghost_min_canopy_cover %||% 0.85)

  is_ghost <- ghost_chm | ghost_density | ghost_canopy

  reason <- dplyr::case_when(
    ghost_chm     ~ sprintf("CHM median %.1f m > threshold inside road corridor", chm_median),
    ghost_density ~ sprintf("Ground density %.3f below minimum detection threshold", density_mean),
    ghost_canopy  ~ sprintf("Canopy closure %.0f%% > 85%% of road corridor", canopy_fraction * 100),
    TRUE          ~ "Not a ghost road"
  )

  list(
    is_ghost = is_ghost,
    reason   = reason,
    evidence = c(
      chm_median      = round(chm_median,    3),
      density_mean    = round(density_mean,  3),
      canopy_fraction = round(canopy_fraction, 3)
    )
  )
}
extract_corridor_values <- function(road, road_width_m, layers_lidar) {

  if (is.null(layers_lidar) || is.na(road_width_m) || road_width_m <= 0)
    return(NULL)

  # Buffer the centerline by half the road width to get the corridor polygon
  half_width <- max(road_width_m / 2, 1.0)   # minimum 1 m to avoid empty polygon
  corridor   <- sf::st_buffer(road, dist = half_width)

  # Convert to Spatial for raster::extract
  corridor_sp <- sf::as_Spatial(corridor)

  # Extract all pixel values within the corridor (not just centroids)
  vals <- raster::extract(
    layers_lidar,
    corridor_sp,
    df       = TRUE,
    cellnumbers = TRUE
  )

  if (is.null(vals) || nrow(vals) == 0) return(NULL)

  # Also compute each pixel's distance to centerline centre
  # (useful for distinguishing road surface vs shoulder pixels later)
  cell_coords <- raster::xyFromCell(layers_lidar[[1]], vals[["cell"]])
  road_center <- sf::st_centroid(road) |> sf::st_coordinates()
  vals$distance_to_center <- sqrt(
    (cell_coords[, 1] - road_center[1, 1])^2 +
      (cell_coords[, 2] - road_center[1, 2])^2
  )

  return(vals)
}

compute_dimensions <- function(metrics, segment_metrics, corridor_vals, param) {

  # Ensure classification sub-list exists
  if (is.null(param[["classification"]])) {
    param[["classification"]] <- default_classification_param()
  }


  D1 = dim_surface(metrics, segment_metrics, corridor_vals, param)
  D2 = dim_width(metrics, segment_metrics, param)
  D3 = dim_canopy(metrics, segment_metrics, param)
  D4 = dim_detectability(segment_metrics, param)
  D5 = dim_drainage(metrics, segment_metrics, param)

  list(
    D1,
    D2,
    D3,
    D4,
    D5
  )
}


# ── Default classification parameters ────────────────────────────────────────

default_classification_param <- function() {
  list(
    # Class 5 hard-gate thresholds (applied before any scoring)
    ghost_max_density      = 0.20,   # normalised point density; below → ghost
    ghost_max_skel_cont    = 0.30,   # skeleton continuity 0-1; below → ghost
    ghost_min_canopy_cover = 0.85,   # canopy fraction; above → ghost

    # Width thresholds (metres, actual road widths)
    min_half_width_m       = 1.5,    # below this in p10 → pinch-point penalty
    lateral_buffer_m       = 3.0,    # buffer around road for lateral canopy check

    # Class boundary thresholds on the composite score (0-100)
    # A road reaches Class N only if ALL dimension scores >= threshold_N
    # AND the composite score is in the band for class N.
    class1_min_all_dims    = 70,     # every dimension must exceed this for Class 1
    class2_min_detect      = 60,     # detectability floor for Class 2
    class2_min_width       = 50,     # width floor for Class 2
    class2_min_any_dim     = 30,     # no dimension may fall below this for Class 2
    class3_min_detect      = 40,     # detectability floor for Class 3
    class4_min_detect      = 20      # detectability floor for Class 4
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  HELPER UTILITIES
# ══════════════════════════════════════════════════════════════════════════════

#' Clamp a value to [0, 100], replacing NA/NaN/Inf with a neutral 50
safe_score <- function(x, lo = 0, hi = 100) {
  if (!is.finite(x)) return(50)
  as.numeric(max(lo, min(hi, x)))
}

#' Piecewise-linear activation (same API as existing alsroads activation())
#' Ascending: low x → 0, high x → 100
#' Descending: low x → 100, high x → 0
activation_score <- function(x, thresholds, asc = TRUE) {
  lo <- thresholds[1]
  hi <- thresholds[2]
  if (asc) {
    val <- (x - lo) / (hi - lo + 1e-9) * 100
  } else {
    val <- (hi - x) / (hi - lo + 1e-9) * 100
  }
  safe_score(val)
}

#' Weighted mean of a vector, ignoring NA
wmean <- function(x, w = NULL) {
  if (all(is.na(x))) return(NA_real_)
  if (is.null(w)) return(mean(x, na.rm = TRUE))
  sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
}


# ══════════════════════════════════════════════════════════════════════════════
#  DIMENSION 1 — SURFACE QUALITY
#  Sources: segment_metrics columns produced by slice_metrics()
#           + road-level conductivity from road_metrics()
# ══════════════════════════════════════════════════════════════════════════════

#' D1: Surface quality score (0-100)
#'
#' Three sub-signals averaged:
#'   planarity   — low roughness on road pixels  (1 - normalised roughness)
#'   slope_ok    — road surface is near-flat
#'   hard_surface — intensity/conductivity proxy for paved/gravel vs soil
#'
#' @param metrics       output of road_metrics()  (one row data.frame)
#' @param segment_metrics output of road_measure() (one row per segment)
#' @param param         classification param list
#' @return named list: score (0-100), planarity, slope_ok, hard_surface
dim_surface <- function(metrics, segment_metrics, corridor_vals, param) {

  # ── Planarity: roughness from segment_metrics or corridor raster
  roughness_col <- if ("roughness" %in% names(segment_metrics)) "roughness"
  else if ("zsd"  %in% names(segment_metrics)) "zsd"
  else NULL

  if (!is.null(roughness_col)) {
    r         <- mean(segment_metrics[[roughness_col]], na.rm = TRUE)
    planarity <- activation_score(r, c(0.15, 0.02), asc = FALSE)  # c(bad, good)
  } else if (!is.null(corridor_vals) && "roughness" %in% names(corridor_vals)) {
    r         <- mean(corridor_vals[["roughness"]], na.rm = TRUE)
    planarity <- activation_score(r, c(0.15, 0.02), asc = FALSE)
  } else {
    planarity <- 50
  }

  # ── Slope: from slope raster pixels inside the road corridor
  # 0° = flat = good (100),  15°+ = too steep = bad (0)
  if (!is.null(corridor_vals) && "slope" %in% names(corridor_vals)) {
    s        <- mean(corridor_vals[["slope"]], na.rm = TRUE)
    slope_ok <- activation_score(s, c(15, 0), asc = FALSE)   # c(bad, good)
  } else {
    slope_ok <- 50   # unknown → neutral
  }

  # ── Hard surface: conductivity from road_metrics()
  if ("CONDUCTIVITY" %in% names(metrics) && !is.na(metrics$CONDUCTIVITY)) {
    p            <- c(0.2, 0.8)
    hard_surface <- activation_score(metrics$CONDUCTIVITY, p, asc = TRUE)  # c(bad, good)
  } else {
    hard_surface <- 50
  }

  score <- safe_score(mean(c(planarity, slope_ok, hard_surface)))

  list(
    score        = score,
    planarity    = planarity,
    slope_ok     = slope_ok,
    hard_surface = hard_surface
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  DIMENSION 2 — WIDTH CONTINUITY
#  Key improvement over existing road_score():
#  We track the 10th-percentile width (pinch points) and the fraction of
#  segments where width is below the minimum drivable threshold.
# ══════════════════════════════════════════════════════════════════════════════

#' D2: Width continuity score (0-100)
#'
#' Three sub-signals:
#'   width_median    — median drivable width along road
#'   width_p10       — 10th-pctile drivable width (worst pinch point)
#'   width_cont      — fraction of segments above minimum threshold
#'
#' @param metrics         road_metrics() output
#' @param segment_metrics road_measure() output
#' @param param           classification param list
#' @return named list: score, width_median, width_p10, width_continuity
dim_width <- function(metrics, segment_metrics, param) {

  dw <- segment_metrics[["drivable_width"]]
  dw <- dw[!is.na(dw)]

  if (length(dw) == 0) {
    return(list(score = 0, width_median = 0, width_p10 = 0, width_continuity = 0))
  }

  min_w <- param[["classification"]][["min_half_width_m"]] * 2   # full width
  if (is.null(min_w)) min_w <- 2.0

  # Thresholds: 0.2 m → score 0, 4 m+ → score 100
  width_median <- activation_score(median(dw), c(0.2, 4), asc = TRUE)
  width_p10    <- activation_score(quantile(dw, 0.10), c(0.2, min_w * 2), asc = TRUE)

  # Continuity: fraction of segments wide enough to drive
  width_cont   <- safe_score(mean(dw >= min_w, na.rm = TRUE) * 100)

  score <- safe_score(mean(c(width_median, width_p10, width_cont)))

  list(
    score           = score,
    width_median    = width_median,
    width_p10       = width_p10,
    width_continuity = width_cont
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  DIMENSION 3 — CANOPY ENCROACHMENT  (stratified into 3 levels)
#
#  The existing PERCABOVEROAD collapses all vegetation into one number.
#  Here we split it into:
#    overhead  — canopy directly above road (blocks sky, damages surface)
#    lateral   — encroachment into the road corridor from the sides
#    ground    — vegetation growing on the road surface itself (reclamation)
# ══════════════════════════════════════════════════════════════════════════════

#' D3: Canopy encroachment score (0-100)
#'
#' @param metrics         road_metrics() output
#' @param segment_metrics road_measure() output  — must have:
#'   pzabove05   fraction of points > 0.5 m above road surface (overhead veg)
#'   pzabove2    fraction of points > 2 m above road (canopy)
#'   number_accotements  shoulder count (used to infer lateral openness)
#' @param param           classification param list
#' @return named list: score, overhead, lateral, ground_level
dim_canopy <- function(metrics, segment_metrics, param) {

  # ── Overhead canopy: pzabove2 is the cleanest signal for canopy closure
  pzabove2 <- mean(segment_metrics[["pzabove2"]], na.rm = TRUE)
  overhead  <- activation_score(pzabove2,
                                c(60, 5),
                                asc = TRUE)

  # ── Lateral encroachment: use pzabove05 minus pzabove2
  # Points between 0.5 m and 2 m are lateral understorey / shrubs at vehicle height
  pzabove05 <- mean(segment_metrics[["pzabove05"]], na.rm = TRUE)
  lateral_veg <- max(0, pzabove05 - pzabove2)    # shrub-height band
  lateral <- activation_score(lateral_veg, c(40, 5), asc = TRUE)

  # ── Ground-level vegetation: pzabove05 with a low threshold
  # If even low vegetation (0.5 m) covers > 80% of the road, the surface is reclaimed
  ground_level <- activation_score(pzabove05, c(80, 10), asc = TRUE)

  score <- safe_score(mean(c(overhead, lateral, ground_level)))

  list(
    score        = score,
    overhead     = overhead,
    lateral      = lateral,
    ground_level = ground_level
  )
}

# Null-coalescing operator (base R does not have one)
`%||%` <- function(a, b) if (!is.null(a)) a else b


# ══════════════════════════════════════════════════════════════════════════════
#  DIMENSION 4 — DETECTABILITY
#
#  This gates everything else. A road that LiDAR cannot reliably detect
#  cannot be reliably measured on any other dimension.
#
#  Three sub-signals:
#    point_density  — returns per m² inside the road corridor
#    skel_cont      — skeleton continuity (longest connected run / total length)
#    return_ratio   — single-return fraction (open sky → clean ground signal)
# ══════════════════════════════════════════════════════════════════════════════

#' D4: Detectability score (0-100)
#'
#' @param segment_metrics road_measure() output — should have:
#'   point_density  (if available)
#'   n              number of points per segment (fallback density proxy)
#'   road_width     for area normalisation
#' @param param           classification param list
#' @return named list: score, point_density, skel_continuity, return_ratio
dim_detectability <- function(segment_metrics, param) {

  section_length <- param[["extraction"]][["section_length"]] %||% 10   # metres

  # ── Point density: points / m² in road corridor
  if ("point_density" %in% names(segment_metrics)) {
    dens <- mean(segment_metrics[["point_density"]], na.rm = TRUE)
  } else if ("n" %in% names(segment_metrics) &&
             "road_width" %in% names(segment_metrics)) {
    # Approximate density from count and corridor area
    area <- segment_metrics[["road_width"]] * section_length
    dens <- mean(segment_metrics[["n"]] / (area + 1e-6), na.rm = TRUE)
  } else {
    dens <- NA_real_
  }
  # Typical ALS density: 2-4 returns/m² = adequate, 8+ = dense
  point_density <- if (!is.na(dens)) activation_score(dens, c(.2, .8), asc = TRUE) else 0.5

  # ── Skeleton continuity: longest run of non-NA segments / total segments
  # A gap in the road detection → NA in segment_metrics
  # Count consecutive non-NA drivable_width values
  dw <- segment_metrics[["drivable_width"]]
  if (length(dw) > 0) {
    runs     <- rle(!is.na(dw))
    valid_runs <- runs$lengths[runs$values]
    longest  <- if (length(valid_runs) > 0) max(valid_runs) else 0
    skel_cont_raw <- longest / length(dw)
    skel_cont <- safe_score(skel_cont_raw * 100)
  } else {
    skel_cont <- 0
  }

  # ── Return ratio: first-return fraction proxy
  # If segment_metrics has return info use it; otherwise fall back to pzabove2:
  # very high canopy cover means many multi-return pulses → low single-return fraction
  if ("return_ratio" %in% names(segment_metrics)) {
    rr <- mean(segment_metrics[["return_ratio"]], na.rm = TRUE)
    return_ratio <- activation_score(rr, c(0.3, 0.9), asc = TRUE)
  } else {
    # Use 1 - pzabove2/100 as a rough proxy: open sky → low canopy fraction → high score
    pz2 <- mean(segment_metrics[["pzabove2"]], na.rm = TRUE)
    return_ratio <- if (!is.na(pz2)) safe_score((1 - pz2/100) * 100) else 50
  }

  score <- safe_score(mean(c(point_density, skel_cont, return_ratio)))

  list(
    score         = score,
    point_density = point_density,
    skel_cont     = skel_cont,
    return_ratio  = return_ratio
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  DIMENSION 5 — DRAINAGE & SHOULDERS
#
#  A maintained road sheds water.  Three sub-signals:
#    shoulder      — presence and continuity of road shoulders
#    conductivity  — alsroads sigma (existing metric, reused here)
#    camber        — cross-slope symmetry (water drains off both sides)
# ══════════════════════════════════════════════════════════════════════════════

#' D5: Drainage and shoulder score (0-100)
#'
#' @param metrics         road_metrics() output
#' @param segment_metrics road_measure() output — uses number_accotements
#' @param param           classification param list
#' @return named list: score, shoulder, conductivity, camber
dim_drainage <- function(metrics, segment_metrics, param) {

  # ── Shoulders: number_accotements is 0, 1, or 2 per segment
  # Average × 50 gives a 0-100 scale (2 shoulders = 100)
  acc <- segment_metrics[["number_accotements"]]
  acc <- acc[!is.na(acc)]
  shoulder <- if (length(acc) > 0) safe_score(mean(acc) / 2 * 100) else 0

  # ── Continuity of shoulders: fraction of segments with ≥ 1 shoulder
  shoulder_cont <- if (length(acc) > 0) safe_score(mean(acc >= 1) * 100) else 0
  shoulder_score <- safe_score((shoulder + shoulder_cont) / 2)

  # ── Conductivity (directly from existing metrics)
  p <- param[["state"]][["conductivity_thresholds"]] %||% c(0.2, 0.8)
  conductivity <- if (!is.na(metrics$CONDUCTIVITY)) {
    activation_score(metrics$CONDUCTIVITY, p, asc = TRUE)
  } else 50

  # ── Camber: cross-slope symmetry
  # If segment_metrics has left_slope and right_slope columns use them directly.
  # Otherwise proxy: road_width vs drivable_width ratio
  # (if much of the width is not drivable, the cross-profile is irregular)
  if (all(c("left_slope", "right_slope") %in% names(segment_metrics))) {
    ls <- segment_metrics[["left_slope"]]
    rs <- segment_metrics[["right_slope"]]
    diff <- abs(ls - rs)
    camber <- activation_score(mean(diff, na.rm = TRUE), c(0.10, 0.01), asc = FALSE)
  } else if (all(c("road_width", "drivable_width") %in% names(segment_metrics))) {
    rw <- segment_metrics[["road_width"]]
    dw <- segment_metrics[["drivable_width"]]
    ratio <- mean(dw / (rw + 1e-6), na.rm = TRUE)   # 1 = fully drivable → good camber
    camber <- safe_score(ratio * 100)
  } else {
    camber <- 50
  }

  score <- safe_score(mean(c(shoulder_score, conductivity, camber)))

  list(
    score        = score,
    shoulder     = shoulder_score,
    conductivity = conductivity,
    camber       = camber
  )
}

classify_road <- function(dims, param) {

  p   <- param[["classification"]]
  if (is.null(p)) p <- default_classification_param()

  # Unpack dimension scores
  dims$D1 = dims[[1]]
  dims$D2 = dims[[2]]
  dims$D3 = dims[[3]]
  dims$D4 = dims[[4]]
  dims$D5 = dims[[5]]

  D1 <- dims$D1$score
  D2 <- dims$D2$score
  D3 <- dims$D3$score
  D4 <- dims$D4$score
  D5 <- dims$D5$score

  all_d  <- c(D1 = D1, D2 = D2, D3 = D3, D4 = D4, D5 = D5)
  composite <- safe_score(mean(all_d))

  # Flat feature vector for logging / downstream ML
  features <- c(
    # D1
    surf_planarity    = dims$D1$planarity,
    surf_slope_ok     = dims$D1$slope_ok,
    surf_hard_surface = dims$D1$hard_surface,
    # D2
    width_median      = dims$D2$width_median,
    width_p10         = dims$D2$width_p10,
    width_continuity  = dims$D2$width_continuity,
    # D3
    canopy_overhead   = dims$D3$overhead,
    canopy_lateral    = dims$D3$lateral,
    canopy_ground     = dims$D3$ground_level,
    # D4
    detect_density    = dims$D4$point_density,
    detect_skel_cont  = dims$D4$skel_cont,
    detect_return_ratio = dims$D4$return_ratio,
    # D5
    drain_shoulder    = dims$D5$shoulder,
    drain_conductivity = dims$D5$conductivity,
    drain_camber      = dims$D5$camber
  )

  # ── Class 5 hard gate  ────────────────────────────────────────────────────
  if (is_ghost_road(dims$D4, dims$D3, param)) {
    return(list(class = 5L, score = composite, dims = all_d, features = features))
  }

  # ── Rules (ordered, first match wins)  ───────────────────────────────────
  cls <- if (all(all_d >= (p$class1_min_all_dims %||% 70))) {
    1L
  } else if (D4 >= (p$class2_min_detect  %||% 60) &&
             D2 >= (p$class2_min_width   %||% 50) &&
             all(all_d >= (p$class2_min_any_dim %||% 30))) {
    2L
  } else if (D4 >= (p$class3_min_detect  %||% 40) &&
             (D1 < 50 | D2 < 50)) {
    3L
  } else if (D4 >= (p$class4_min_detect  %||% 20)) {
    4L
  } else {
    5L
  }

  list(class = cls, score = composite, dims = all_d, features = features)
}

is_ghost_road <- function(d4, d3, param) {
  p <- param[["classification"]]

  ghost_density  <- (d4$point_density / 100) < (p$ghost_max_density      %||% 0.20)
  ghost_skel     <- (d4$skel_cont     / 100) < (p$ghost_max_skel_cont    %||% 0.30)
  ghost_canopy   <- (d3$overhead      / 100) < (1 - (p$ghost_min_canopy_cover %||% 0.85))

  ghost_density | ghost_skel | ghost_canopy
}

