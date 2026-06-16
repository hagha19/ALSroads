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

safe_score <- function(x, lo = 0, hi = 100) {
  if (!is.finite(x)) return(50)
  as.numeric(max(lo, min(hi, x)))
}

# thresholds always c(bad_value, good_value)
# direction encoded in threshold order, asc=TRUE means low→high is good
activation_score <- function(x, thresholds, asc = TRUE) {
  lo  <- thresholds[1]
  hi  <- thresholds[2]
  raw <- (x - lo) / (hi - lo + 1e-9) * 100
  val <- if (asc) raw else 100 - raw
  safe_score(val)
}

class_name <- function(cls) {
  c("1" = "Maintained",
    "2" = "Passable",
    "3" = "Degraded",
    "4" = "Barely passable",
    "5" = "Ghost / reclaimed")[as.character(cls)]
}

default_classification_param <- function() {
  list(

    # ── Ghost gate (2/3 signals must agree) ──────────────────────
    ghost_max_shrub_layer        = 50,    # % points 0.5-2m on road surface
    ghost_min_drivable_fraction  = 0.20,  # fraction of segments >= min width
    ghost_max_density            = 0.10,  # ground return fraction

    # ── Five core features — thresholds c(bad, good) ─────────────
    # entropy_norm: 0 = all points near ground, 1 = spread across all heights
    entropy_bad   = 0.80,   # fully vegetated distribution
    entropy_good  = 0.30,   # clean road surface distribution

    # width_median: median drivable width in metres
    width_bad     = 1.0,    # barely a track
    width_good    = 6.0,    # wide logging road

    # lateral vegetation: pzabove05_drive - pzabove2_drive (shrubs at vehicle height)
    lateral_bad   = 40,     # % — heavily encroached
    lateral_good  = 5,      # % — essentially clear

    # pzabove2_drive: canopy closing above driveway
    canopy_bad    = 70,     # % — nearly fully closed
    canopy_good   = 10,     # % — mostly open

    # density: ground return fraction
    density_bad   = 0.05,   # almost no ground signal
    density_good  = 0.40,   # strong ground signal

    # ── Class boundaries on composite score (0-100) ──────────────
    class1_min    = 75,
    class2_min    = 55,
    class3_min    = 35,
    class4_min    = 15,

    # ── Individual feature floors per class ──────────────────────
    # A road cannot be Class 1 if any single feature is catastrophic
    class1_min_each  = 60,   # every feature must exceed this for Class 1
    class2_min_each  = 30,   # every feature must exceed this for Class 2

    # ── Width ─────────────────────────────────────────────────────
    min_half_width_m = 1.5,

    # ── Drainage assurance upgrade ────────────────────────────────
    drainage_upgrade_threshold = 75
  )
}


# ══════════════════════════════════════════════════════════════════
#  SCORE THE FIVE FEATURES
# ══════════════════════════════════════════════════════════════════

score_features <- function(segment_metrics, corridor_vals, param) {

  p <- param[["classification"]] %||% default_classification_param()

  dw <- segment_metrics[["drivable_width"]]
  dw <- dw[!is.na(dw)]

  # ── 1. Entropy (vertical point distribution)
  # Low = points near ground = clear road
  # High = spread across heights = vegetated / reclaimed
  entropy_raw <- if ("entropy_norm" %in% names(segment_metrics)) {
    mean(segment_metrics[["entropy_norm"]], na.rm = TRUE)
  } else NA_real_

  f_entropy <- if (!is.na(entropy_raw)) {
    activation_score(entropy_raw,
                     c(p$entropy_good %||% 0.3
                       ,p$entropy_bad  %||% 0.80),
                     asc = FALSE)   # low entropy = good
  } else 50

  # ── 2. Width median (drivable width)
  # Wide = good, narrow = bad
  width_raw <- if (length(dw) > 0) median(dw) else NA_real_

  f_width <- if (!is.na(width_raw)) {
    activation_score(width_raw,
                     c(p$width_bad  %||% 1.0,
                       p$width_good %||% 6.0),
                     asc = TRUE)   # wide = good
  } else 50

  # ── 3. Lateral vegetation above driveway (pzabove05_drive - pzabove2_drive)
  # Shrubs at vehicle height inside driveable corridor
  # Low = clear = good, high = encroached = bad
  pz05_drive <- mean(segment_metrics[["pzabove05_drive"]], na.rm = TRUE)
  pz2_drive  <- mean(segment_metrics[["pzabove2_drive"]],  na.rm = TRUE)
  lateral_raw <- max(0, pz05_drive - pz2_drive)

  f_lateral <- activation_score(lateral_raw,
                                c(p$lateral_good %||% 5,
                                  p$lateral_bad  %||% 40),
                                asc = FALSE)   # low lateral veg = good

  # ── 4. Canopy above driveway (pzabove2_drive)
  # Canopy closing over the drivable corridor
  # Low = open sky = good, high = closed canopy = bad
  f_canopy <- activation_score(pz2_drive,
                               c(p$canopy_good %||% 10,
                                 p$canopy_bad  %||% 70),
                                 asc = FALSE)   # low canopy = good

  # ── 5. Ground return density
  # High fraction of ground returns = LiDAR sees the surface clearly = good
  density_raw <- if (!is.null(corridor_vals) &&
                     "density" %in% names(corridor_vals)) {
    mean(corridor_vals[["density"]], na.rm = TRUE)
  } else NA_real_

  f_density <- if (!is.na(density_raw)) {
    activation_score(density_raw,
                     c(p$density_bad  %||% 0.05,
                       p$density_good %||% 0.40),
                     asc = TRUE)   # high density = good
  } else 50

  list(
    entropy  = f_entropy,
    width    = f_width,
    lateral  = f_lateral,
    canopy   = f_canopy,
    density  = f_density,
    # raw values for logging
    raw = list(
      entropy  = entropy_raw,
      width    = width_raw,
      lateral  = lateral_raw,
      canopy   = pz2_drive,
      density  = density_raw
    )
  )
}


# ══════════════════════════════════════════════════════════════════
#  CLASSIFIER
# ══════════════════════════════════════════════════════════════════

classify_road <- function(features,  param) {

  p <- param[["classification"]] %||% default_classification_param()

  # Unpack feature scores
  f <- c(
    entropy = features$entropy,
    width   = features$width,
    lateral = features$lateral,
    canopy  = features$canopy,
    density = features$density
  )

  # Composite score: simple mean of all five
  composite <- safe_score(mean(f, na.rm = TRUE))

  # ── Class rules ─────────────────────────────────────────────────
  # Primary rule: composite score band
  # Secondary rule: no individual feature below per-class floor
  # This prevents a road scoring Class 1 on four features but
  # being catastrophically bad on one

  cls <- if (
    composite >= (p$class1_min %||% 75) &&
    all(f >= (p$class1_min_each %||% 60))
  ) {
    1L

  } else if (
    composite >= (p$class2_min %||% 55) &&
    all(f >= (p$class2_min_each %||% 30))
  ) {
    2L

  } else if (
    composite >= (p$class3_min %||% 35)
  ) {
    3L

  } else if (
    composite >= (p$class4_min %||% 15)
  ) {
    4L

  } else {
    5L
  }

  # ── Drainage assurance upgrade ───────────────────────────────────
  # upgraded <- FALSE
  # if (cls %in% c(2L, 3L, 4L) &&
  #     drainage$assurance >= (p$drainage_upgrade_threshold %||% 75)) {
  #   cls      <- cls - 1L
  #   upgraded <- TRUE
  # }

  list(
    class     = cls,
    score     = composite,
    features  = f
    )
}


# ══════════════════════════════════════════════════════════════════
#  MAIN ENTRY POINT
# ══════════════════════════════════════════════════════════════════

road_class_full <- function(road, metrics, segment_metrics,
                            layers_lidar, param) {

  # ── Step 1: extract raster values within actual road width ───────
  corridor_vals <- extract_corridor_values(
    road         = road,
    road_width_m = metrics$ROADWIDTH,
    layers_lidar = layers_lidar
  )

  # ── Step 2: ghost gate ───────────────────────────────────────────
  ghost <- detect_ghost_road(segment_metrics, corridor_vals, param)

  if (ghost$is_ghost) {
    return(data.frame(
      class             = 5L,
      class_name        = "Ghost / reclaimed",
      score             = 0,
      f_entropy         = NA_real_,
      f_width           = NA_real_,
      f_lateral         = NA_real_,
      f_canopy          = NA_real_,
      f_density         = NA_real_,
      ghost_reason      = ghost$reason,
      shrub_layer       = ghost$evidence["shrub_layer"],
      drivable_fraction = ghost$evidence["drivable_fraction"],
      density           = ghost$evidence["density"]
    ))
  }

  # ── Step 3: score features and classify ─────────────────────────
  features <- score_features(segment_metrics, corridor_vals, param)
  # drainage <- compute_drainage_assurance(segment_metrics, param)
  result   <- classify_road(features, param)

  return(data.frame(
    class              = result$class,
    class_name         = class_name(result$class),
    score              = round(result$score,            1),
    f_entropy          = round(result$features["entropy"], 1),
    f_width            = round(result$features["width"],   1),
    f_lateral          = round(result$features["lateral"], 1),
    f_canopy           = round(result$features["canopy"],  1),
    f_density          = round(result$features["density"], 1),
    ghost_reason       = NA_character_,
    shrub_layer        = NA_real_,
    drivable_fraction  = NA_real_,
    density            = NA_real_
  ))
}
# # ══════════════════════════════════════════════════════════════════
# #  PARAMETERS
# # ══════════════════════════════════════════════════════════════════
#
# default_classification_param <- function() {
#   list(
#
#     # ── Ghost gate (2/3 signals must agree) ──────────────────────
#     ghost_max_shrub_layer        = 50,    # % points 0.5-2m on road surface
#     ghost_min_drivable_fraction  = 0.20,  # fraction of segments drivable
#     ghost_max_density            = 0.10,  # ground return fraction
#
#     # ── D1 Surface thresholds ─────────────────────────────────────
#     roughness_bad  = 0.2,   # m — roughness above this → score 0
#     roughness_good = 0.05,   # m — roughness below this → score 100
#     slope_bad      = 15,     # degrees — slope above this → score 0
#     slope_good     = 3,      # degrees — slope below this → score 100
#
#     # ── Class 1 ──────────────────────────────────────────────────
#     class1_min_surface  = 65,
#     class1_min_width    = 65,
#     class1_min_detect   = 65,
#     class1_max_canopy   = 60,   # looser — forest roads always have canopy
#
#     # ── Class 2 ──────────────────────────────────────────────────
#     class2_min_detect   = 45,
#     class2_min_width    = 40,
#     class2_min_any_dim  = 20,
#
#     # ── Class 3 ──────────────────────────────────────────────────
#     class3_min_detect   = 25,
#
#     # ── Class 4 ──────────────────────────────────────────────────
#     class4_min_detect   = 10,
#
#     # ── D5 Drainage assurance ─────────────────────────────────────
#     drainage_upgrade_threshold = 75,  # upgrades class by 1 if above this
#
#     # ── Width ─────────────────────────────────────────────────────
#     min_half_width_m    = 1.5,
#     lateral_buffer_m    = 3.0
#   )
# }
#
#
# # ══════════════════════════════════════════════════════════════════
# #  GHOST GATE
# # ══════════════════════════════════════════════════════════════════
#
detect_ghost_road <- function(segment_metrics, corridor_vals, param) {

  p <- param[["classification"]] %||% default_classification_param()

  if (is.null(segment_metrics) || nrow(segment_metrics) == 0) {
    return(list(
      is_ghost = TRUE,
      reason   = "No segment metrics available",
      evidence = c(shrub_layer = NA, drivable_fraction = NA, density = NA)
    ))
  }

  # ── Signal 1: points between 0.5m and 2m
  # shrubs/saplings growing ON the road surface
  shrub_layer <- mean(
    segment_metrics[["pzabove05"]] - segment_metrics[["pzabove2"]],
    na.rm = TRUE
  )
  ghost_shrub <- !is.na(shrub_layer) &&
    shrub_layer > (p$ghost_max_shrub_layer %||% 50)

  # ── Signal 2: drivable width collapsed
  min_w <- (p$min_half_width_m %||% 1.5) * 2

  drivable_fraction <- mean(
    segment_metrics[["drivable_width"]] >= min_w, na.rm = TRUE
  )
  ghost_width <- !is.na(drivable_fraction) &&
    drivable_fraction < (p$ghost_min_drivable_fraction %||% 0.20)

  # ── Signal 3: ground return density (% ground / total points)
  density <- if (!is.null(corridor_vals) && "density" %in% names(corridor_vals)) {
    mean(corridor_vals[["density"]], na.rm = TRUE)
  } else NA_real_

  ghost_density <- !is.na(density) &&
    density < (p$ghost_max_density %||% 0.10)

  # ── 2 out of 3 must agree
  n_triggered <- sum(c(ghost_shrub, ghost_width, ghost_density))
  is_ghost    <- n_triggered >= 2

  reason <- dplyr::case_when(
    is_ghost & ghost_shrub & ghost_width & ghost_density ~
      sprintf("All 3: shrub=%.1f%%, drivable=%.0f%%, density=%.3f",
              shrub_layer, drivable_fraction * 100, density),
    is_ghost & ghost_shrub & ghost_density ~
      sprintf("Shrub=%.1f%% + density=%.3f",
              shrub_layer, density),
    is_ghost & ghost_width & ghost_density ~
      sprintf("Drivable=%.0f%% + density=%.3f",
              drivable_fraction * 100, density),
    is_ghost & ghost_shrub & ghost_width ~
      sprintf("Shrub=%.1f%% + drivable=%.0f%%",
              shrub_layer, drivable_fraction * 100),
    TRUE ~ "Not a ghost road"
  )

  list(
    is_ghost = is_ghost,
    reason   = reason,
    evidence = c(
      shrub_layer       = round(shrub_layer,       2),
      drivable_fraction = round(drivable_fraction, 2),
      density           = round(density,           3)
    )
  )
}
#
#
# # ══════════════════════════════════════════════════════════════════
# #  D1 — SURFACE QUALITY
# # ══════════════════════════════════════════════════════════════════
#
# dim_surface <- function(metrics, segment_metrics, corridor_vals, param) {
#
#   p <- param[["classification"]] %||% default_classification_param()
#
#   # ── Planarity: roughness from segment_metrics, fallback to corridor raster
#   roughness_col <- if ("roughness" %in% names(segment_metrics)) "roughness"
#   else if ("zsd"  %in% names(segment_metrics)) "zsd"
#   else NULL
#
#   if (!is.null(roughness_col)) {
#     r <- mean(segment_metrics[[roughness_col]], na.rm = TRUE)
#   } else if (!is.null(corridor_vals) && "roughness" %in% names(corridor_vals)) {
#     r <- mean(corridor_vals[["roughness"]], na.rm = TRUE)
#   } else {
#     r <- NA_real_
#   }
#
#   planarity <- if (!is.na(r)) {
#     activation_score(r,
#                      c(p$roughness_good  %||% 0.06,
#                        p$roughness_bad %||% 0.2),
#                      asc = TRUE)   # c(bad, good) — lower roughness is better
#   } else 50
#
#   # ── Slope: from slope raster pixels inside road corridor (degrees)
#   # 0° = flat = good (100),  15°+ = too steep = bad (0)
#   if (!is.null(corridor_vals) && "slope" %in% names(corridor_vals)) {
#     s        <- mean(corridor_vals[["slope"]], na.rm = TRUE)
#     slope_ok <- activation_score(s,
#                                  c(p$slope_bad  %||% 15,
#                                    p$slope_good %||% 5),
#                                  asc = FALSE)   # c(bad, good)
#   } else {
#     slope_ok <- 50
#   }
#
#   score <- safe_score(mean(c(planarity, slope_ok), na.rm = TRUE))
#
#   list(
#     score     = score,
#     planarity = planarity,
#     slope_ok  = slope_ok
#   )
# }
#
#
# # ══════════════════════════════════════════════════════════════════
# #  D2 — WIDTH CONTINUITY
# # ══════════════════════════════════════════════════════════════════
#
# dim_width <- function(metrics, segment_metrics, param) {
#
#   p      <- param[["classification"]] %||% default_classification_param()
#   min_w  <- (p$min_half_width_m %||% 1.5) * 2   # full width
#
#   dw <- segment_metrics[["drivable_width"]]
#   dw <- dw[!is.na(dw)]
#
#   if (length(dw) == 0) {
#     return(list(score = 0, width_median = 0, width_p10 = 0, width_continuity = 0))
#   }
#
#   # Median width: 0m = bad, 8m+ = good
#   width_median <- activation_score(median(dw), c(1, 6), asc = TRUE)
#
#   # 10th percentile: captures worst pinch points
#   width_p10 <- activation_score(
#     quantile(dw, 0.10, na.rm = TRUE),
#     c(0.5, min_w),
#     asc = TRUE
#   )
#
#   # Continuity: fraction of segments wide enough to drive
#   # e.g. 0.80 = road drivable in 80% of cross-sections → score 80
#   width_cont <- safe_score(mean(dw >= min_w, na.rm = TRUE) * 100)
#
#   score <- safe_score(mean(c(width_median, width_p10, width_cont)))
#
#   list(
#     score            = score,
#     width_median     = width_median,
#     width_p10        = width_p10,
#     width_continuity = width_cont
#   )
# }
#
#
# # ══════════════════════════════════════════════════════════════════
# #  D3 — CANOPY ENCROACHMENT  (stratified by height)
# # ══════════════════════════════════════════════════════════════════
#
# dim_canopy <- function(metrics, segment_metrics, param) {
#
#   p <- param[["classification"]] %||% default_classification_param()
#
#   pzabove05 <- mean(segment_metrics[["pzabove05_drive"]], na.rm = TRUE)
#   pzabove2  <- mean(segment_metrics[["pzabove2_drive"]],  na.rm = TRUE)
#
#   # ── Overhead: canopy closing above road (> 2m)
#   # Uses existing param thresholds if available
#   overhead_thresh <- param[["state"]][["percentage_veg_thresholds"]] %||% c(60, 5)
#   overhead <- activation_score(pzabove2, overhead_thresh, asc = FALSE)
#
#   # ── Lateral: shrubs at vehicle height (0.5m - 2m) encroaching from sides
#   lateral_veg <- max(0, pzabove05 - pzabove2)
#   lateral     <- activation_score(lateral_veg, c(40, 10), asc = FALSE)
#
#   # ── Ground level: any vegetation > 0.5m on road surface
#   ground_level <- activation_score(pzabove05, c(80, 20), asc = FALSE)
#
#
#   score <- safe_score(mean(c(overhead, lateral, ground_level)))
#
#   list(
#     score        = score,
#     overhead     = overhead,
#     lateral      = lateral,
#     ground_level = ground_level
#   )
# }
#
#
# # ══════════════════════════════════════════════════════════════════
# #  D4 — DETECTABILITY
# # ══════════════════════════════════════════════════════════════════
#
# dim_detectability <- function(segment_metrics, corridor_vals, param) {
#
#   # ── Signal 1: ground return density (% ground points / total)
#   density <- if (!is.null(corridor_vals) && "density" %in% names(corridor_vals)) {
#     mean(corridor_vals[["density"]], na.rm = TRUE)
#   } else NA_real_
#
#   # 0% ground returns = bad, 40%+ = open road with good ground signal
#   density_score <- if (!is.na(density)) {
#     activation_score(density, c(0.1, 0.40), asc = TRUE)
#   } else 50
#
#   # ── Signal 2: hillshade texture (coefficient of variation)
#   # Multidirectional hillshade (8 azimuths) from rasterize_conductivity4()
#   # stored as "dtm" layer in your raster stack
#   # Clear road surface → uniform hillshade → low CV → good
#   # Vegetated/degraded → noisy hillshade → high CV → bad
#   hillshade_score <- if (!is.null(corridor_vals) && "dtm" %in% names(corridor_vals)) {
#     hvals <- corridor_vals[["dtm"]]
#     hvals <- hvals[!is.na(hvals)]
#     if (length(hvals) > 5) {
#       cv <- sd(hvals) / (mean(hvals) + 1e-6)
#       # CV: 0.05 = very smooth (good), 0.40 = very rough (bad)
#       activation_score(cv, c(0.40, 0.05), asc = FALSE)
#     } else 50
#   } else 50
#
#   if ("entropy_norm" %in% names(segment_metrics)) {
#     ent <- mean(segment_metrics[["entropy_norm"]], na.rm = TRUE)
#     # 0.30 = low entropy = good, 0.90 = high entropy = bad
#     entropy_score <- activation_score(ent, c(0.80, 0.30), asc = FALSE)
#   } else {
#     entropy_score <- 50
#   }
#
#   score <- safe_score(mean(c(density_score, hillshade_score, entropy_score), na.rm = TRUE))
#
#   list(
#     score     = score,
#     density   = density_score,
#     hillshade = hillshade_score
#   )
# }
#
#
# # ══════════════════════════════════════════════════════════════════
# #  D5 — DRAINAGE ASSURANCE  (modifier, not a standalone dimension)
# # ══════════════════════════════════════════════════════════════════
#
# compute_drainage_assurance <- function(segment_metrics, param) {
#
#   # ── Shoulders: presence and continuity along road
#   acc           <- segment_metrics[["number_accotements"]]
#   shoulder_mean <- safe_score(mean(acc, na.rm = TRUE) / 2 * 100)
#   shoulder_cont <- safe_score(mean(acc >= 1, na.rm = TRUE) * 100)
#   shoulder      <- mean(c(shoulder_mean, shoulder_cont))
#
#   # ── Camber: drivable width / road width ratio
#   # High ratio → crown intact → water drains both sides → maintained road
#   ratio  <- mean(
#     segment_metrics[["drivable_width"]] /
#       (segment_metrics[["road_width"]] + 1e-6),
#     na.rm = TRUE
#   )
#   camber <- safe_score(ratio * 100)
#
#   assurance <- safe_score(mean(c(shoulder, camber)))
#
#   list(
#     assurance = assurance,
#     shoulder  = shoulder,
#     camber    = camber
#   )
# }
#
#
# # ══════════════════════════════════════════════════════════════════
# #  COMPUTE ALL DIMENSIONS
# # ══════════════════════════════════════════════════════════════════
#
# compute_dimensions <- function(metrics, segment_metrics, corridor_vals, param) {
#
#   if (is.null(param[["classification"]])) {
#     param[["classification"]] <- default_classification_param()
#   }
#
#   list(
#     D1 = dim_surface(metrics, segment_metrics, corridor_vals, param),
#     D2 = dim_width(metrics, segment_metrics, param),
#     D3 = dim_canopy(metrics, segment_metrics, param),
#     D4 = dim_detectability(segment_metrics, corridor_vals, param)
#   )
# }
#
#
# # ══════════════════════════════════════════════════════════════════
# #  CLASSIFIER
# # ══════════════════════════════════════════════════════════════════
#
# classify_road <- function(dims, drainage, param) {
#
#   p  <- param[["classification"]] %||% default_classification_param()
#
#   D1 <- dims$D1$score
#   D2 <- dims$D2$score
#   D3 <- dims$D3$score
#   D4 <- dims$D4$score
#   DA <- drainage$assurance
#
#   all_d     <- c(D1 = D1, D2 = D2, D3 = D3, D4 = D4)
#   composite <- safe_score(mean(all_d))
#
#   # ── Base classification on D1-D4 ────────────────────────────────
#   cls <- if (
#     D1 >= (p$class1_min_surface %||% 65) &&
#     D2 >= (p$class1_min_width   %||% 65) &&
#     D4 >= (p$class1_min_detect  %||% 65) &&
#     D3 >= (p$class1_max_canopy  %||% 60)
#   ) {
#     1L
#   } else if (
#     D4 >= (p$class2_min_detect  %||% 45) &&
#     D2 >= (p$class2_min_width   %||% 40) &&
#     all(all_d >= (p$class2_min_any_dim %||% 20))
#   ) {
#     2L
#   } else if (
#     # ── Class 3: degraded but still a functioning road ─────────────
#     # D4: detectable enough to trust measurements (floor)
#     D4 >= (p$class3_min_detect  %||% 25) &&
#     # D1: surface is degraded but road surface still exists
#     D1 >= (p$class3_min_surface %||% 20) &&
#     D1 <  (p$class3_max_surface %||% 50) &&
#     # D2: road is narrow but still passable in some segments
#     D2 >= (p$class3_min_width   %||% 20) &&
#     D2 <  (p$class3_max_width   %||% 50)
#   ) {
#     3L
#   } else if (
#     D4 >= (p$class4_min_detect  %||% 10)
#   ) {
#     4L
#   } else {
#     5L
#   }
#
#   # ── Drainage assurance upgrade ───────────────────────────────────
#   # Strong shoulder + camber evidence upgrades one class up
#   # Never downgrades, never upgrades Class 1, never upgrades Class 5
#   upgraded <- FALSE
#   if (cls %in% c(2L, 3L, 4L) &&
#       DA  >= (p$drainage_upgrade_threshold %||% 75)) {
#     cls      <- cls - 1L
#     upgraded <- TRUE
#   }
#
#   list(
#     class     = cls,
#     score     = composite,
#     upgraded  = upgraded,
#     dims      = all_d,
#     drainage  = DA
#   )
# }
#
#
# # ══════════════════════════════════════════════════════════════════
# #  CORRIDOR EXTRACTION
# # ══════════════════════════════════════════════════════════════════
#
extract_corridor_values <- function(road, road_width_m, layers_lidar) {

  if (is.null(layers_lidar) || is.na(road_width_m) || road_width_m <= 0)
    return(NULL)

  half_width  <- max(road_width_m / 2, 1.0)
  corridor    <- sf::st_buffer(road, dist = half_width)
  corridor_sp <- sf::as_Spatial(corridor)

  vals <- raster::extract(
    layers_lidar,
    corridor_sp,
    df          = TRUE,
    cellnumbers = TRUE
  )

  if (is.null(vals) || nrow(vals) == 0) return(NULL)

  return(vals)
}
#
#
# # ══════════════════════════════════════════════════════════════════
# #  MAIN ENTRY POINT
# # ══════════════════════════════════════════════════════════════════
#
# road_class_full <- function(road, metrics, segment_metrics, layers_lidar, param) {
#
#   # ── Step 1: extract raster values within actual road width ───────
#   corridor_vals <- extract_corridor_values(
#     road         = road,
#     road_width_m = metrics$ROADWIDTH,
#     layers_lidar = layers_lidar
#   )
#
#   # ── Step 2: ghost gate ───────────────────────────────────────────
#   ghost <- detect_ghost_road(segment_metrics, corridor_vals, param)
#
#   if (ghost$is_ghost) {
#     return(data.frame(
#       class             = 5L,
#       class_name        = "Ghost / reclaimed",
#       score             = 0,
#       D1_surface        = NA_real_,
#       D2_width          = NA_real_,
#       D3_canopy         = NA_real_,
#       D4_detect         = NA_real_,
#       drainage_assurance = NA_real_,
#       upgraded          = FALSE,
#       ghost_reason      = ghost$reason,
#       shrub_layer       = ghost$evidence["shrub_layer"],
#       drivable_fraction = ghost$evidence["drivable_fraction"],
#       density           = ghost$evidence["density"]
#     ))
#   }
#
#   # ── Step 3: compute dimensions ───────────────────────────────────
#   dims     <- compute_dimensions(metrics, segment_metrics, corridor_vals, param)
#   drainage <- compute_drainage_assurance(segment_metrics, param)
#   result   <- classify_road(dims, drainage, param)
#
#   return(data.frame(
#     class             = result$class,
#     score             = round(result$score,        1),
#     D1_surface        = round(result$dims["D1"],   1),
#     D2_width          = round(result$dims["D2"],   1),
#     D3_canopy         = round(result$dims["D3"],   1),
#     D4_detect         = round(result$dims["D4"],   1),
#     drainage_assurance = round(result$drainage,    1),
#     upgraded          = result$upgraded,
#     ghost_reason      = NA_character_,
#     shrub_layer       = NA_real_,
#     drivable_fraction = NA_real_,
#     density           = NA_real_
#   ))
# }
