#' @export
#' @rdname rasterize_conductivity
rasterize_conductivity2 <- function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  UseMethod("rasterize_conductivity2", las)
}

#' @export
rasterize_conductivity2.LAS <- function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  use_intensity <- "Intensity" %in% names(las)
  display <- getOption("ALSroads.debug.finding")
  pkg <- if (is.null(dtm)) getOption("lidR.raster.default") else lidR:::raster_pkg(dtm)

  if (is.null(dtm))
  {
    dtm <- lidR::rasterize_terrain(las, 1, lidR::tin(), pkg = "raster")
  }
  else if (lidR:::is_raster(dtm))
  {
    res <- round(raster::res(dtm)[1], 2)
    if (res > 1) stop("The DTM must have a resolution of 1 m or less.")

    bb = lidR::st_bbox(las)
    dtm <- raster::crop(dtm, raster::extent(bb))

    if (res < 1)
      dtm <- raster::aggregate(dtm, fact = 1/res, fun = mean)
  }
  else
  {
    stop("dtm must be a RasterLayer or must be NULL")
  }

  # plot(dtm, col = gray(1:30/30))

  mask = NULL
  if (!is.null(water) && length(sf::st_geometry(water)) > 0)
  {
    id <- NULL
    water <- sf::st_geometry(water)
    bbox <- suppressWarnings(sf::st_bbox(las))
    bbox <- sf::st_set_crs(bbox, sf::st_crs(water))
    mask <- sf::st_crop(water, bbox)
    if (length(mask) > 0)
    {
      las = lidR::classify_poi(las, lidR::LASWATER, roi = water)
      mask <- sf::as_Spatial(mask)
    }
    else
    {
      mask = NULL
    }
  }

  # Force to use raster
  if (lidR:::raster_pkg(dtm) == "terra")
    dtm <- raster::raster(dtm)

  nlas <- lidR::normalize_height(las, dtm) |> suppressMessages() |> suppressWarnings()
  nlas@data[Classification == LASWATER & Z > 2, Classification := LASBRIGDE]
  bridge = lidR::filter_poi(nlas, Classification == LASBRIGDE)
  bridge = sf::st_coordinates(bridge, z = FALSE)

  # Terrain metrics using the raster package (slope, roughness)
  slope <- terra::terrain(dtm, opt = c("slope"), unit = "degrees")
  if (!is.null(mask)) slope <- raster::mask(slope, mask, inverse = T)
  if(display) raster::plot(slope, col = gray(1:30/30), main = "slope")

  smoothdtm = raster::focal(dtm, matrix(1,5,5), mean)
  roughdtm = dtm - smoothdtm
  roughness <- terra::terrain(roughdtm, opt = c("roughness"), unit = "degrees")

  if(display) {
    raster::plot(roughdtm, col = gray(1:30/30), main = "Residual roughness")
    raster::plot(roughness, col = gray(1:30/30), main = "Roughness")
  }

  # Slope-based conductivity
  s <- param$conductivity$s
  sigma_s <- activation(slope, s, "piecewise-linear", asc = FALSE)
  sigma_s <- raster::aggregate(sigma_s, fact = 2, fun = mean)

  if (display) raster::plot(sigma_s, col = viridis::viridis(25), main = "Conductivity slope")
  verbose("   - Slope conductivity map \n")

  # Roughness-based conductivity
  r <- param$conductivity$r
  sigma_r <- activation(roughness, r, "piecewise-linear", asc = FALSE)
  sigma_r <- raster::aggregate(sigma_r, fact = 2, fun = mean)

  if (display) raster::plot(sigma_r, col = viridis::viridis(25), main = "Conductivity roughness")
  verbose("   - Roughness conductivity map\n")

  # Edge-based conductivity
  e    <- param$conductivity$e
  sobl <- sobel.RasterLayer(slope)
  #plot(sigma_e, col = gray(1:30/30))
  sigma_e <- activation(sobl, e, "thresholds", asc = FALSE)
  sigma_e <- raster::aggregate(sigma_e, fact = 2, fun = mean)

  if (display) raster::plot(sigma_e, col = viridis::viridis(25), main = "Conductivity Sobel edges")
  verbose("   - Sobel conductivity map\n")

  # Intensity-based conductivity
  sigma_i   <- dtm
  sigma_i[] <- 0
  if (use_intensity)
  {
    template2m <- raster::aggregate(dtm, fact = 2)
    template2m[] <- 0
    irange = intensity_range_by_flightline(las, template2m)
    #irange = terra::focal(irange, matrix(1,3,3), mean, na.rm = T)
    if (!is.null(mask)) irange <- raster::mask(irange, mask, inverse = T)

    if (display) raster::plot(irange, col = heat.colors(20), main = "Intensity range")

    #th <- stats::quantile(irange[], probs = q, na.rm = TRUE)
    th <- c(0.25,0.35)
    #sigma_i <- template2m
    sigma_i <- activation(irange, th, "piecewise-linear", asc = FALSE)
    #sigma_i = activation2(irange)

    if (display) raster::plot(sigma_i, col = viridis::viridis(20), main = "Conductivity intensity")
    verbose("   - Intensity conductivity map\n")
  }


  # CHM-based conductivity
  h <- param$conductivity$h
  chm <- lidR::grid_canopy(nlas, dtm, lidR::p2r())
  if (display) raster::plot(chm, col = height.colors(25),  main = "CHM")
  sigma_h <- dtm
  sigma_h <- activation(chm, h, "piecewise-linear", asc = FALSE)
  sigma_h <- raster::aggregate(sigma_h, fact = 2, fun = mean)

  if (display) raster::plot(sigma_h, col = viridis::inferno(25), main = "Conductivity CHM")
  verbose("   - CHM conductivity map\n")

  # Lowpoints-based conductivity
  # Check the presence/absence of lowpoints
  z1 <- 0.5
  z2 <- 3
  th <- 0.01
  lp <- density_lp_by_flightline(nlas, dtm)
  lp[is.na(lp)] = 0
  sigma_lp <- activation(lp, th, "thresholds", asc = FALSE)
  sigma_lp <- raster::aggregate(sigma_lp, fact = 2, fun = min)

  if (display) raster::plot(lp, col = viridis::inferno(20), main = "Number of low point")
  if (display) raster::plot(sigma_lp, col = viridis::inferno(2), main = "Bottom layer")
  verbose("   - Bottom layer conductivity map\n")

  # Density-based conductivity
  # Notice that the paper makes no mention of smoothing
  q <- param$conductivity$d

  d <- density_gnd_by_flightline1(las, sigma_lp, drop_angles = 0)
  d[is.na(d)] = 0
  M   <- matrix(1,3,3)
  d = terra::focal(d, M, mean, na.rm = T)
  th <- mean(d[d>0], na.rm = T)
  d[is.na(d)] = 0

  if (!is.null(mask)) d <- raster::mask(d, mask, inverse = T, updatevalue = 0)
  if (display) raster::plot(d, col = viridis::inferno(15), main = "Density of ground points")

  sigma_d <- activation(d, c(0.33, 0.66), "piecewise-linear")

  if (display)  raster::plot(sigma_d, col = viridis::inferno(25), main = "Conductivity density")
  verbose("   - Density conductivity map\n")

  hard_slope = slope < 25
  hard_slope <- raster::aggregate(hard_slope, fact = 2, fun = min)


  # Final conductivity sigma
  alpha = param$conductivity$alpha
  alpha$i = alpha$i * as.numeric(use_intensity)
  max_coductivity <- sum(unlist(alpha))+1
  sigma <- (alpha$s* sigma_s + alpha$e * sigma_e + alpha$d * sigma_d + alpha$h * sigma_h + alpha$r * sigma_r + alpha$i * sigma_i)
  sigma <- sigma/max_coductivity

  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Conductivity")

  sigma <- edge_enhancement(sigma, interation = 50, lambda = 0.05, k = 20, pass = 2)

  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")

  sigma[sigma < 0.1] = 0.1
  #raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")
  sigma[hard_slope == 0] = 0.1
  #raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")

  sigma[is.na(sigma)] <- 0.1 # lakes
  cells = raster::cellFromXY(sigma, bridge)
  sigma[cells] = 0.75
  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity with bridge")

  if (pkg == "terra") sigma <- terra::rast(sigma)

  dots = list(...)

  if (!is.null(dots$rawlayer))
  {
    dtm2 = raster::aggregate(dtm, fact = 2, fun = mean)
    slope2 = raster::aggregate(slope, fact = 2, fun = mean)
    roughness2 =  raster::aggregate(roughness, fact = 2, fun = mean)
    chm2 =  raster::aggregate(chm, fact = 2, fun = mean)
    lp2 =  raster::aggregate(lp, fact = 2, fun = mean)
    sobl2 =  raster::aggregate(sobl, fact = 2, fun = mean)
    u = raster::stack(slope2, roughness2, chm2, sigma_i, lp2, d, sobl2)
    names(u) = c("slope", 'rought', "chm", "intensity", "low", "density", "sobel")
    return(u)
  }


  return(sigma)
}

#' @export
rasterize_conductivity2.LAScluster = function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  x <- lidR::readLAS(las)
  if (lidR::is.empty(x)) return(NULL)

  sigma <- rasterize_conductivity2(x, dtm, water, param, ...)
  sigma <- lidR:::raster_crop(sigma, lidR::st_bbox(las))
  return(sigma)
}

#' @export
rasterize_conductivity2.LAScatalog = function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  # Enforce some options
  if (lidR::opt_select(las) == "*") lidR::opt_select(las) <- "xyzciap"

  # Compute the alignment options including the case when res is a raster/stars/terra
  alignment <- lidR:::raster_alignment(1)

  if (lidR::opt_chunk_size(las) > 0 && lidR::opt_chunk_size(las) < 2*alignment$res)
    stop("The chunk size is too small. Process aborted.", call. = FALSE)

  # Processing
  options <- list(need_buffer = TRUE, drop_null = TRUE, raster_alignment = alignment, automerge = TRUE)
  output  <- lidR::catalog_apply(las, rasterize_conductivity2, dtm = dtm, water = water, param = param, ..., .options = options)
  return(output)
}

rasterize_conductivity3 <- function(las, centerline, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  UseMethod("rasterize_conductivity3", las)
}

#' @export
rasterize_conductivity3.LAS <- function(las, centerline, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  use_intensity <- "Intensity" %in% names(las)
  display <- getOption("ALSroads.debug.finding")
  pkg <- if (is.null(dtm)) getOption("lidR.raster.default") else lidR:::raster_pkg(dtm)

  if (is.null(dtm))
  {
    dtm <- lidR::rasterize_terrain(las, 1, lidR::tin(), pkg = "raster")
  }
  else if (lidR:::is_raster(dtm))
  {
    res <- round(raster::res(dtm)[1], 2)
    if (res > 1) stop("The DTM must have a resolution of 1 m or less.")

    bb = lidR::st_bbox(las)
    dtm <- raster::crop(dtm, raster::extent(bb))

    if (res < 1)
      dtm <- raster::aggregate(dtm, fact = 1/res, fun = mean)
  }
  else
  {
    stop("dtm must be a RasterLayer or must be NULL")
  }

  # plot(dtm, col = gray(1:30/30))

  mask = NULL
  if (!is.null(water) && length(sf::st_geometry(water)) > 0)
  {
    id <- NULL
    water <- sf::st_geometry(water)
    bbox <- suppressWarnings(sf::st_bbox(las))
    bbox <- sf::st_set_crs(bbox, sf::st_crs(water))
    mask <- sf::st_crop(water, bbox)
    if (length(mask) > 0)
    {
      las = lidR::classify_poi(las, lidR::LASWATER, roi = water)
      mask <- sf::as_Spatial(mask)
    }
    else
    {
      mask = NULL
    }
  }
  # Force to use raster
  if (lidR:::raster_pkg(dtm) == "terra")
    dtm <- raster::raster(dtm)

  nlas <- lidR::normalize_height(las, dtm) |> suppressMessages() |> suppressWarnings()
  nlas@data[Classification == LASWATER & Z > 2, Classification := LASBRIGDE]
  bridge = lidR::filter_poi(nlas, Classification == LASBRIGDE)
  bridge = sf::st_coordinates(bridge, z = FALSE)

  # Terrain metrics using the raster package (slope, roughness)
  slope <- terra::terrain(dtm, opt = c("slope"), unit = "degrees")
  if (!is.null(mask)) slope <- raster::mask(slope, mask, inverse = T)
  if(display) raster::plot(slope, col = gray(1:30/30), main = "slope")
  area = lidR::area(las)
  n_points <- nrow(las@data)
  densitypoints = n_points/area
  w <- matrix(1, nrow = 5, ncol = 5)
  dsm <- lidR::grid_canopy(nlas, dtm, lidR::p2r())
  smoothdsm = raster::focal(dsm, matrix(1,5,5), mean)
  diff_dsm = dsm - smoothdsm
  roughness_dsm <- terra::terrain(diff_dsm, v = "roughness", unit = "degrees")
  if (densitypoints<10) {roughness_dsm<- terra::focal(roughness_dsm, w = w, fun = function(x) {
        if (is.na(x[13])) {  # Center of 5x5 window is element 13 in vectorized form
               return(mean(x, na.rm = TRUE))
           } else {
                 return(x[13])
             }
     }, na.policy = "omit", fillvalue = NA)}
  sigma_dsm <- activation(roughness_dsm, c(0.01,20), "piecewise-linear", asc = FALSE)
  sigma_dsm <- raster::aggregate(sigma_dsm, fact = 2, fun = mean)
  if(display) {
    raster::plot(roughness_dsm, col = gray(1:30/30), main = "dsm roughness")
    raster::plot(sigma_dsm, col = gray(1:30/30), main = "conductivity dsm Roughness")
  }

  smoothdtm = raster::focal(dtm, matrix(1,5,5), mean)
  roughdtm = dtm - smoothdtm
  roughness <- terra::terrain(roughdtm, opt = c("roughness"), unit = "degrees")

  if(display) {
    raster::plot(roughdtm, col = gray(1:30/30), main = "Residual roughness")
    raster::plot(roughness, col = gray(1:30/30), main = "Roughness")
  }

  # Slope-based conductivity
  # s <- param$conductivity$s
  slope_segment <- extract_raster_along_centerline(centerline, slope, pixel_buffer = 5)
  s <- compute_segment_percentiles(slope_segment)

  sigma_s <- activation(slope, s, "piecewise-linear", asc = FALSE)
  sigma_s <- raster::aggregate(sigma_s, fact = 2, fun = mean)

  if (display) raster::plot(sigma_s, col = viridis::viridis(25), main = "Conductivity slope")
  verbose("   - Slope conductivity map \n")

  # Roughness-based conductivity
  # r <- param$conductivity$r
  roughness_segment <- extract_raster_along_centerline(centerline, roughness, pixel_buffer = 5)
  r <- compute_segment_percentiles(roughness_segment)
  sigma_r <- activation(roughness, r, "piecewise-linear", asc = FALSE)
  # sigma_r[sigma_r >= 0.98] <- 0
  sigma_r <- raster::aggregate(sigma_r, fact = 2, fun = max)

  if (display) raster::plot(sigma_r, col = viridis::viridis(25), main = "Conductivity roughness")
  verbose("   - Roughness conductivity map\n")

  # Edge-based conductivity
  e    <- param$conductivity$e
  sobl <- sobel.RasterLayer(slope)
  #plot(sigma_e, col = gray(1:30/30))
  sigma_e <- activation(sobl, e, "thresholds", asc = FALSE)
  sigma_e <- raster::aggregate(sigma_e, fact = 2, fun = mean)

  if (display) raster::plot(sigma_e, col = viridis::viridis(25), main = "Conductivity Sobel edges")
  verbose("   - Sobel conductivity map\n")

  # Intensity-based conductivity
  sigma_i   <- dtm
  sigma_i[] <- 0

  if (use_intensity)
  {
    template2m <- raster::aggregate(dtm, fact = 2)
    template2m[] <- 0
    irange = intensity_range_by_flightline(las, template2m)
    if (densitypoints<10) {irange<- terra::focal(irange, w = w, fun = function(x) {
           if (is.na(x[13])) {  # Center of 5x5 window is element 13 in vectorized form
                 return(mean(x, na.rm = TRUE))
             } else {
                   return(x[13])
               }
       }, na.policy = "omit", fillvalue = NA)}
    #irange = terra::focal(irange, matrix(1,3,3), mean, na.rm = T)
    if (!is.null(mask)) irange <- raster::mask(irange, mask, inverse = T)

    if (display) raster::plot(irange, col = heat.colors(20), main = "Intensity range")

    #th <- stats::quantile(irange[], probs = q, na.rm = TRUE)
    # th <- c(0.15,0.25)
    # irange_segment <- extract_raster_along_centerline(centerline, irange, pixel_buffer = 8)
    # max_th = compute_segment_percentiles(irange_segment)[1]
    th = c(0.05 , 0.35)
    #sigma_i <- template2m
    sigma_i <- activation(irange, th, "piecewise-linear", asc = FALSE)
    #sigma_i = activation2(irange)

    if (display) raster::plot(sigma_i, col = viridis::viridis(20), main = "Conductivity intensity")
    verbose("   - Intensity conductivity map\n")
  }


  # CHM-based conductivity
  # h <- param$conductivity$h
  chm <- lidR::grid_canopy(nlas, dtm, lidR::p2r())
  if (densitypoints<10) {chm<- terra::focal(chm, w = w, fun = function(x) {
         if (is.na(x[13])) {  # Center of 5x5 window is element 13 in vectorized form
              return(mean(x, na.rm = TRUE))
          } else {
                 return(x[13])
             }
     }, na.policy = "omit", fillvalue = NA)}
  if (display) raster::plot(chm, col = height.colors(25),  main = "CHM")
  # chm_segment <- extract_raster_along_centerline(centerline, chm, pixel_buffer = 8)
  # max_h = find_trough(chm_segment)
  mask_chm = activation(chm, c(5 ,10), "piecewise-linear", asc = FALSE) + 0.1
  if (densitypoints<25) { h = c(0 ,3)} else {h = c(0 ,1)}
  sigma_h <- dtm
  sigma_h <- activation(chm, h, "piecewise-linear", asc = FALSE)
  sigma_h <- raster::aggregate(sigma_h, fact = 2, fun = max)
  mask_chm <- raster::aggregate(mask_chm, fact = 2, fun = max)

  if (display) raster::plot(sigma_h, col = viridis::inferno(25), main = "Conductivity CHM")
  if (display) raster::plot(mask_chm, col = viridis::inferno(25), main = "mask_chm CHM")
  verbose("   - CHM conductivity map\n")

  # Lowpoints-based conductivity
  # Check the presence/absence of lowpoints
  z1 <- 0.5
  z2 <- 3
  th <- 0.2
  lp <- density_lp_by_flightline(nlas, dtm)
  if (densitypoints<10) {lp<- terra::focal(lp, w = w, fun = function(x) {
         if (is.na(x[13])) {  # Center of 5x5 window is element 13 in vectorized form
               return(mean(x, na.rm = TRUE))
           } else {
                 return(x[13])
             }
     }, na.policy = "omit", fillvalue = NA)}
  lp[is.na(lp)] = 0
  sigma_lp <- activation(lp, th, "thresholds", asc = FALSE)
  sigma_lp <- raster::aggregate(sigma_lp, fact = 2, fun = min)

  if (display) raster::plot(lp, col = viridis::inferno(20), main = "Number of low point")
  if (display) raster::plot(sigma_lp, col = viridis::inferno(2), main = "Bottom layer")
  verbose("   - Bottom layer conductivity map\n")

  # Density-based conductivity
  # Notice that the paper makes no mention of smoothing
  q <- param$conductivity$d

  d <- density_gnd_by_flightline(las, sigma_lp, drop_angles = 0)
  if (densitypoints<10) {d<- terra::focal(d, w = w, fun = function(x) {
         if (is.na(x[13])) {  # Center of 5x5 window is element 13 in vectorized form
               return(mean(x, na.rm = TRUE))
           } else {
                 return(x[13])
             }
     }, na.policy = "omit", fillvalue = NA)}
  d[is.na(d)] = 0
  M   <- matrix(1,3,3)
  d = terra::focal(d, M, mean, na.rm = T)
  th <- mean(d[d>0], na.rm = T)
  d[is.na(d)] = 0
  d_segment <- extract_raster_along_centerline(centerline, d, pixel_buffer = 8)
  th <- compute_segment_percentiles(d_segment)
  if (!is.null(mask)) d <- raster::mask(d, mask, inverse = T, updatevalue = 0)
  if (display) raster::plot(d, col = viridis::inferno(15), main = "Density of ground points")

  sigma_d <- activation(d, th, "piecewise-linear")

  if (display)  raster::plot(sigma_d, col = viridis::inferno(25), main = "Conductivity density")
  verbose("   - Density conductivity map\n")

  hard_slope = slope < 25
  hard_slope <- raster::aggregate(hard_slope, fact = 2, fun = min)

  # Final conductivity sigma
  alpha = param$conductivity$alpha
  alpha$i = alpha$i * as.numeric(use_intensity)
  max_coductivity <- sum(unlist(alpha))+1
  if (densitypoints<10) {

  sigma <- (mask_chm+0.1) * (alpha$s*sigma_s + alpha$e*sigma_e + 0.5*alpha$h * sigma_h + alpha$r * sigma_r + 0.5*alpha$i * sigma_i)
  } else if (densitypoints<25) {
    sigma <- (mask_chm+0.1) * (sigma_dsm + alpha$s*sigma_s + alpha$e*sigma_e +alpha$d * sigma_d + alpha$h * sigma_h + alpha$r * sigma_r + alpha$i * sigma_i)
  } else {
    sigma <- (sigma_d* mask_chm+0.1) * (sigma_dsm + alpha$s*sigma_s + alpha$e*sigma_e +alpha$d * sigma_d + alpha$h * sigma_h + alpha$r * sigma_r + alpha$i * sigma_i)
  }
  sigma <- sigma/max_coductivity

  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Conductivity")

  # sigma <- edge_enhancement(sigma, interation = 50, lambda = 0.05, k = 20, pass = 2)

  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")

  sigma[sigma < 0.1] = 0.1
  #raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")
  sigma[hard_slope == 0] = 0.1
  #raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")

  sigma[is.na(sigma)] <- 0.1 # lakes
  cells = raster::cellFromXY(sigma, bridge)
  sigma[cells] = 0.75
  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity with bridge")

  if (pkg == "terra") sigma <- terra::rast(sigma)

  dots = list(...)


  if (!is.null(dots$rawlayer))
  {
    dtm2 = raster::aggregate(dtm, fact = 1, fun = mean)
    dtm2 <- terra::rast(dtm2)  # convert if needed
    # 8 azimuth directions
    azimuths <- c(0, 45, 90, 135, 180, 225, 270, 315)
    altitude <- 45   # sun elevation
    slope1  <- terrain(dtm2, v = "slope",  unit = "radians")
    aspect1 <- terrain(dtm2, v = "aspect", unit = "radians")
    # Compute hillshades
    shades <- lapply(azimuths, function(a) {
      terra::shade(slope = slope1, aspect = aspect1, angle = altitude, direction = a)
    })
    all_shades <- terra::rast(shades)
    dtm2 <- raster::raster(terra::app(all_shades, fun = mean, na.rm = TRUE))
    slope2 = raster::aggregate(slope, fact = 1, fun = mean)
    dsm2 =  raster::aggregate(dsm, fact = 1, fun = mean)
    roughness_dsm2 =  raster::aggregate(roughness_dsm, fact = 1, fun = mean)
    roughness2 =  raster::aggregate(roughness, fact = 1, fun = mean)
    sobl2 =  raster::aggregate(sobl, fact = 1, fun = mean)
    irange2 =  raster::aggregate(irange, fact = 1, fun = mean)
    chm2 =  raster::aggregate(chm, fact = 1, fun = mean)
    lp2 =  raster::aggregate(lp, fact = 1, fun = mean)
    d2 =  raster::aggregate(d, fact = 1, fun = mean)
    ref <- dtm2
    slope2 <- raster::resample(slope2, ref)
    dsm2 <- raster::resample(dsm2, ref)
    roughness_dsm2 <- raster::resample(roughness_dsm2, ref)
    roughness2 <- raster::resample(roughness2, ref)
    sobl2 <- raster::resample(sobl2, ref)
    irange2 <- raster::resample(irange2, ref)
    chm2 <- raster::resample(chm2, ref)
    lp2 <- raster::resample(lp2, ref)
    d2 <- raster::resample(d2, ref)

    hull <- sf::st_buffer(centerline, param$extraction$road_buffer)
    p <- sf::st_buffer(centerline, 1)
    f <- fasterize::fasterize(p, ref)
    f <- raster::distance(f)
    f <- raster::mask(f, hull)
    fmin <- min(f[], na.rm = T)
    fmax <- max(f[], na.rm = T)
    target_min <- 1-param[["constraint"]][["confidence"]]
    f <- (1-(((f - fmin) * (1 - target_min)) / (fmax - fmin)))
    f <- raster::resample(f, ref)

    u = raster::stack(ref, slope2, dsm2, roughness_dsm2, roughness2, sobl2, irange2, chm2,
                      lp2, d2,f)
    names(u) = c("dtm", 'slope', "dsm", "roughness_dsm", "roughness", "sobel", "intensity",
                 "chm", "low", "density", "distance")
    return(u)
  }

  return(sigma)
}

#' @export
rasterize_conductivity3.LAScluster = function(las, centerline, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  x <- lidR::readLAS(las)
  if (lidR::is.empty(x)) return(NULL)

  sigma <- rasterize_conductivity3(x, centerline, dtm, water, param, ...)
  sigma <- lidR:::raster_crop(sigma, lidR::st_bbox(las))
  return(sigma)
}

#' @export
rasterize_conductivity3.LAScatalog = function(las, centerline, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  # Enforce some options
  if (lidR::opt_select(las) == "*") lidR::opt_select(las) <- "xyzciap"

  # Compute the alignment options including the case when res is a raster/stars/terra
  alignment <- lidR:::raster_alignment(1)

  if (lidR::opt_chunk_size(las) > 0 && lidR::opt_chunk_size(las) < 2*alignment$res)
    stop("The chunk size is too small. Process aborted.", call. = FALSE)

  # Processing
  options <- list(need_buffer = TRUE, drop_null = TRUE, raster_alignment = alignment, automerge = TRUE)
  output  <- lidR::catalog_apply(las, rasterize_conductivity3, centerline= centerline, dtm = dtm, water = water, param = param, ..., .options = options)
  return(output)
}


rasterize_conductivity4 <- function(las, centerline, dtm = NULL, water = NULL, param = alsroads_default_parameters2, dl_model= NULL)
{
  UseMethod("rasterize_conductivity4", las)
}

#' @export
rasterize_conductivity4 <- function(las, centerline, dtm = NULL, water = NULL, param = alsroads_default_parameters2, dl_model= NULL)
{
  use_intensity <- "Intensity" %in% names(las)
  display <- getOption("ALSroads.debug.finding")
  pkg <- if (is.null(dtm)) getOption("lidR.raster.default") else lidR:::raster_pkg(dtm)

  if (is.null(dtm))
  {
    dtm <- lidR::rasterize_terrain(las, 1, lidR::tin(), pkg = "raster")
  }
  else if (lidR:::is_raster(dtm))
  {
    res <- round(raster::res(dtm)[1], 2)
    if (res > 1) stop("The DTM must have a resolution of 1 m or less.")

    bb = lidR::st_bbox(las)
    dtm <- raster::crop(dtm, raster::extent(bb))

    if (res < 1)
      dtm <- raster::aggregate(dtm, fact = 1/res, fun = mean)
  }
  else
  {
    stop("dtm must be a RasterLayer or must be NULL")
  }

  # plot(dtm, col = gray(1:30/30))

  mask = NULL
  if (!is.null(water) && length(sf::st_geometry(water)) > 0)
  {
    id <- NULL
    water <- sf::st_geometry(water)
    bbox <- suppressWarnings(sf::st_bbox(las))
    bbox <- sf::st_set_crs(bbox, sf::st_crs(water))
    mask <- sf::st_crop(water, bbox)
    if (length(mask) > 0)
    {
      las = lidR::classify_poi(las, lidR::LASWATER, roi = water)
      mask <- sf::as_Spatial(mask)
    }
    else
    {
      mask = NULL
    }
  }
  # Force to use raster
  if (lidR:::raster_pkg(dtm) == "terra")
    dtm <- raster::raster(dtm)

  nlas <- lidR::normalize_height(las, dtm) |> suppressMessages() |> suppressWarnings()
  nlas@data[Classification == LASWATER & Z > 2, Classification := LASBRIGDE]
  bridge = lidR::filter_poi(nlas, Classification == LASBRIGDE)
  bridge = sf::st_coordinates(bridge, z = FALSE)

  # Terrain metrics using the raster package (slope, roughness)
  slope <- terra::terrain(dtm, opt = c("slope"), unit = "degrees")
  if (!is.null(mask)) slope <- raster::mask(slope, mask, inverse = T)
  if(display) raster::plot(slope, col = gray(1:30/30), main = "slope")
  area = lidR::area(las)
  n_points <- nrow(las@data)
  densitypoints = n_points/area
  w <- matrix(1, nrow = 5, ncol = 5)
  dsm <- lidR::grid_canopy(nlas, dtm, lidR::p2r())
  smoothdsm = raster::focal(dsm, matrix(1,5,5), mean)
  diff_dsm = dsm - smoothdsm
  roughness_dsm <- terra::terrain(diff_dsm, v = "roughness", unit = "degrees")
  if (densitypoints<10) {roughness_dsm<- terra::focal(roughness_dsm, w = w, fun = function(x) {
    if (is.na(x[13])) {  # Center of 5x5 window is element 13 in vectorized form
      return(mean(x, na.rm = TRUE))
    } else {
      return(x[13])
    }
  }, na.policy = "omit", fillvalue = NA)}
  sigma_dsm <- activation(roughness_dsm, c(0.01,20), "piecewise-linear", asc = FALSE)
  sigma_dsm <- raster::aggregate(sigma_dsm, fact = 1, fun = mean)
  if(display) {
    raster::plot(roughness_dsm, col = gray(1:30/30), main = "dsm roughness")
    raster::plot(sigma_dsm, col = gray(1:30/30), main = "conductivity dsm Roughness")
  }

  smoothdtm = raster::focal(dtm, matrix(1,5,5), mean)
  roughdtm = dtm - smoothdtm
  roughness <- terra::terrain(roughdtm, opt = c("roughness"), unit = "degrees")

  if(display) {
    raster::plot(roughdtm, col = gray(1:30/30), main = "Residual roughness")
    raster::plot(roughness, col = gray(1:30/30), main = "Roughness")
  }

  # Slope-based conductivity
  # s <- param$conductivity$s
  slope_segment <- extract_raster_along_centerline(centerline, slope, pixel_buffer = 5)
  s <- compute_segment_percentiles(slope_segment)

  sigma_s <- activation(slope, s, "piecewise-linear", asc = FALSE)
  sigma_s <- raster::aggregate(sigma_s, fact = 1, fun = mean)

  if (display) raster::plot(sigma_s, col = viridis::viridis(25), main = "Conductivity slope")
  verbose("   - Slope conductivity map \n")

  # Roughness-based conductivity
  # r <- param$conductivity$r
  roughness_segment <- extract_raster_along_centerline(centerline, roughness, pixel_buffer = 5)
  r <- compute_segment_percentiles(roughness_segment)
  sigma_r <- activation(roughness, r, "piecewise-linear", asc = FALSE)
  # sigma_r[sigma_r >= 0.98] <- 0
  sigma_r <- raster::aggregate(sigma_r, fact = 1, fun = max)

  if (display) raster::plot(sigma_r, col = viridis::viridis(25), main = "Conductivity roughness")
  verbose("   - Roughness conductivity map\n")

  # Edge-based conductivity
  e    <- param$conductivity$e
  sobl <- sobel.RasterLayer(slope)
  #plot(sigma_e, col = gray(1:30/30))
  sigma_e <- activation(sobl, e, "thresholds", asc = FALSE)
  sigma_e <- raster::aggregate(sigma_e, fact = 1, fun = mean)

  if (display) raster::plot(sigma_e, col = viridis::viridis(25), main = "Conductivity Sobel edges")
  verbose("   - Sobel conductivity map\n")

  # Intensity-based conductivity
  sigma_i   <- dtm
  sigma_i[] <- 0

  if (use_intensity)
  {
    template2m <- raster::aggregate(dtm, fact = 1)
    template2m[] <- 0
    irange = intensity_range_by_flightline(las, template2m)
    if (densitypoints<10) {irange<- terra::focal(irange, w = w, fun = function(x) {
      if (is.na(x[13])) {  # Center of 5x5 window is element 13 in vectorized form
        return(mean(x, na.rm = TRUE))
      } else {
        return(x[13])
      }
    }, na.policy = "omit", fillvalue = NA)}
    #irange = terra::focal(irange, matrix(1,3,3), mean, na.rm = T)
    if (!is.null(mask)) irange <- raster::mask(irange, mask, inverse = T)

    if (display) raster::plot(irange, col = heat.colors(20), main = "Intensity range")

    #th <- stats::quantile(irange[], probs = q, na.rm = TRUE)
    # th <- c(0.15,0.25)
    # irange_segment <- extract_raster_along_centerline(centerline, irange, pixel_buffer = 8)
    # max_th = compute_segment_percentiles(irange_segment)[1]
    th = c(0.05 , 0.35)
    #sigma_i <- template2m
    sigma_i <- activation(irange, th, "piecewise-linear", asc = FALSE)
    #sigma_i = activation2(irange)

    if (display) raster::plot(sigma_i, col = viridis::viridis(20), main = "Conductivity intensity")
    verbose("   - Intensity conductivity map\n")
  }


  # CHM-based conductivity
  # h <- param$conductivity$h
  chm <- lidR::grid_canopy(nlas, dtm, lidR::p2r())
  if (densitypoints<10) {chm<- terra::focal(chm, w = w, fun = function(x) {
    if (is.na(x[13])) {  # Center of 5x5 window is element 13 in vectorized form
      return(mean(x, na.rm = TRUE))
    } else {
      return(x[13])
    }
  }, na.policy = "omit", fillvalue = NA)}
  if (display) raster::plot(chm, col = height.colors(25),  main = "CHM")
  # chm_segment <- extract_raster_along_centerline(centerline, chm, pixel_buffer = 8)
  # max_h = find_trough(chm_segment)
  mask_chm = activation(chm, c(5 ,10), "piecewise-linear", asc = FALSE) + 0.1
  if (densitypoints<25) { h = c(0 ,3)} else {h = c(0 ,1)}
  sigma_h <- dtm
  sigma_h <- activation(chm, h, "piecewise-linear", asc = FALSE)
  sigma_h <- raster::aggregate(sigma_h, fact = 1, fun = max)
  mask_chm <- raster::aggregate(mask_chm, fact = 1, fun = max)

  if (display) raster::plot(sigma_h, col = viridis::inferno(25), main = "Conductivity CHM")
  if (display) raster::plot(mask_chm, col = viridis::inferno(25), main = "mask_chm CHM")
  verbose("   - CHM conductivity map\n")

  # Lowpoints-based conductivity
  # Check the presence/absence of lowpoints
  z1 <- 0.5
  z2 <- 3
  th <- 0.2
  lp <- density_lp_by_flightline(nlas, dtm)
  if (densitypoints<10) {lp<- terra::focal(lp, w = w, fun = function(x) {
    if (is.na(x[13])) {  # Center of 5x5 window is element 13 in vectorized form
      return(mean(x, na.rm = TRUE))
    } else {
      return(x[13])
    }
  }, na.policy = "omit", fillvalue = NA)}
  lp[is.na(lp)] = 0
  sigma_lp <- activation(lp, th, "thresholds", asc = FALSE)
  sigma_lp <- raster::aggregate(sigma_lp, fact = 1, fun = min)

  if (display) raster::plot(lp, col = viridis::inferno(20), main = "Number of low point")
  if (display) raster::plot(sigma_lp, col = viridis::inferno(2), main = "Bottom layer")
  verbose("   - Bottom layer conductivity map\n")

  # Density-based conductivity
  # Notice that the paper makes no mention of smoothing
  q <- param$conductivity$d

  d <- density_gnd_by_flightline(las, sigma_lp, drop_angles = 0)
  if (densitypoints<10) {d<- terra::focal(d, w = w, fun = function(x) {
    if (is.na(x[13])) {  # Center of 5x5 window is element 13 in vectorized form
      return(mean(x, na.rm = TRUE))
    } else {
      return(x[13])
    }
  }, na.policy = "omit", fillvalue = NA)}
  d[is.na(d)] = 0
  M   <- matrix(1,3,3)
  d = terra::focal(d, M, mean, na.rm = T)
  th <- mean(d[d>0], na.rm = T)
  d[is.na(d)] = 0
  d_segment <- extract_raster_along_centerline(centerline, d, pixel_buffer = 8)
  th <- compute_segment_percentiles(d_segment)
  if (!is.null(mask)) d <- raster::mask(d, mask, inverse = T, updatevalue = 0)
  if (display) raster::plot(d, col = viridis::inferno(15), main = "Density of ground points")

  sigma_d <- activation(d, th, "piecewise-linear")

  if (display)  raster::plot(sigma_d, col = viridis::inferno(25), main = "Conductivity density")
  verbose("   - Density conductivity map\n")

  hard_slope = slope < 25
  hard_slope <- raster::aggregate(hard_slope, fact = 1, fun = min)

  # Final conductivity sigma
  alpha = param$conductivity$alpha
  alpha$i = alpha$i * as.numeric(use_intensity)
  max_coductivity <- sum(unlist(alpha))+1
  if (densitypoints<10) {

    sigma <- (mask_chm+0.1) * (alpha$s*sigma_s + alpha$e*sigma_e + 0.5*alpha$h * sigma_h + alpha$r * sigma_r + 0.5*alpha$i * sigma_i)
  } else if (densitypoints<25) {
    sigma <- (mask_chm+0.1) * (sigma_dsm + alpha$s*sigma_s + alpha$e*sigma_e +alpha$d * sigma_d + alpha$h * sigma_h + alpha$r * sigma_r + alpha$i * sigma_i)
  } else {
    sigma <- (sigma_d* mask_chm+0.1) * (sigma_dsm + alpha$s*sigma_s + alpha$e*sigma_e +alpha$d * sigma_d + alpha$h * sigma_h + alpha$r * sigma_r + alpha$i * sigma_i)
  }
  sigma <- sigma/max_coductivity

  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Conductivity")

  # sigma <- edge_enhancement(sigma, interation = 50, lambda = 0.05, k = 20, pass = 2)

  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")

  sigma[sigma < 0.1] = 0.1
  #raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")
  sigma[hard_slope == 0] = 0.1
  #raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")

  sigma[is.na(sigma)] <- 0.1 # lakes
  cells = raster::cellFromXY(sigma, bridge)
  sigma[cells] = 0.75
  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity with bridge")

  if (pkg == "terra") sigma <- terra::rast(sigma)
  dtm2 = raster::aggregate(dtm, fact = 1, fun = mean)
  dtm2 <- terra::rast(dtm2)  # convert if needed
  # 8 azimuth directions
  azimuths <- c(0, 45, 90, 135, 180, 225, 270, 315)
  altitude <- 45   # sun elevation
  slope1  <- terrain(dtm2, v = "slope",  unit = "radians")
  aspect1 <- terrain(dtm2, v = "aspect", unit = "radians")
  # Compute hillshades
  shades <- lapply(azimuths, function(a) {
    terra::shade(slope = slope1, aspect = aspect1, angle = altitude, direction = a)
  })
  all_shades <- terra::rast(shades)
  dtm2 <- raster::raster(terra::app(all_shades, fun = mean, na.rm = TRUE))
  slope2 = raster::aggregate(slope, fact = 1, fun = mean)
  dsm2 =  raster::aggregate(dsm, fact = 1, fun = mean)
  roughness_dsm2 =  raster::aggregate(roughness_dsm, fact = 1, fun = mean)
  roughness2 =  raster::aggregate(roughness, fact = 1, fun = mean)
  sobl2 =  raster::aggregate(sobl, fact = 1, fun = mean)
  irange2 =  raster::aggregate(irange, fact = 1, fun = mean)
  chm2 =  raster::aggregate(chm, fact = 1, fun = mean)
  lp2 =  raster::aggregate(lp, fact = 1, fun = mean)
  d2 =  raster::aggregate(d, fact = 1, fun = mean)
  ref <- dtm2
  slope2 <- raster::resample(slope2, ref)
  dsm2 <- raster::resample(dsm2, ref)
  roughness_dsm2 <- raster::resample(roughness_dsm2, ref)
  roughness2 <- raster::resample(roughness2, ref)
  sobl2 <- raster::resample(sobl2, ref)
  irange2 <- raster::resample(irange2, ref)
  chm2 <- raster::resample(chm2, ref)
  lp2 <- raster::resample(lp2, ref)
  d2 <- raster::resample(d2, ref)

  hull <- sf::st_buffer(centerline, param$extraction$road_buffer)
  p <- sf::st_buffer(centerline, 1)
  f <- fasterize::fasterize(p, ref)
  f <- raster::distance(f)
  f <- raster::mask(f, hull)
  fmin <- min(f[], na.rm = T)
  fmax <- max(f[], na.rm = T)
  target_min <- 1-param[["constraint"]][["confidence"]]
  f <- (1-(((f - fmin) * (1 - target_min)) / (fmax - fmin)))
  f <- raster::resample(f, ref)


  # dtm2 = raster::aggregate(dtm, fact = 2, fun = mean)
  # slope2 = raster::aggregate(slope, fact = 2, fun = mean)
  # dsm2 =  raster::aggregate(dsm, fact = 2, fun = mean)
  # roughness_dsm2 =  raster::aggregate(roughness_dsm, fact = 2, fun = mean)
  # roughness2 =  raster::aggregate(roughness, fact = 2, fun = mean)
  # sobl2 =  raster::aggregate(sobl, fact = 2, fun = mean)
  # irange2 =  raster::aggregate(irange, fact = 2, fun = mean)
  # chm2 =  raster::aggregate(chm, fact = 2, fun = mean)
  # lp2 =  raster::aggregate(lp, fact = 2, fun = mean)
  # d2 =  raster::aggregate(d, fact = 2, fun = mean)
  # ref <- dtm2
  # slope2 <- raster::resample(slope2, ref)
  # dsm2 <- raster::resample(dsm2, ref)
  # roughness_dsm2 <- raster::resample(roughness_dsm2, ref)
  # roughness2 <- raster::resample(roughness2, ref)
  # sobl2 <- raster::resample(sobl2, ref)
  # irange2 <- raster::resample(irange2, ref)
  # chm2 <- raster::resample(chm2, ref)
  # lp2 <- raster::resample(lp2, ref)
  # d2 <- raster::resample(d2, ref)
  #
  # hull <- sf::st_buffer(centerline, param$extraction$road_buffer)
  # p <- sf::st_buffer(centerline, 1)
  # f <- fasterize::fasterize(p, ref)
  # f <- raster::distance(f)
  # f <- raster::mask(f, hull)
  # fmin <- min(f[], na.rm = T)
  # fmax <- max(f[], na.rm = T)
  # target_min <- 1-param[["constraint"]][["confidence"]]
  # f <- (1-(((f - fmin) * (1 - target_min)) / (fmax - fmin)))
  # f <- raster::resample(f, ref)

  # u = raster::stack(ref, slope2, dsm2, roughness_dsm2, roughness2, sobl2, irange2, chm2,
  #                   lp2, d2,f)
  # names(u) = c("dtm", 'slope', "dsm", "roughness_dsm", "roughness", "sobel", "intensity",
  #              "chm", "low", "density", "distance")
  u = raster::stack(ref, slope2, roughness_dsm2, roughness2, sobl2, irange2, chm2,
                    lp2, d2)
  names(u) = c("dtm", 'slope',  "roughness_dsm", "roughness", "sobel", "intensity",
               "chm", "low", "density")
  u <- raster::mask(u, f)
  # u <- reclassify(u, cbind(NA, 0))
  # u[[11]] = 0
  tiles = tile_raster(terra::rast(u) ,  target_size = 512)
  prediction_tiles <- list()

  for (i in 1:length(tiles$tiles)) {
    cat("Processing tile", i, "of", length(tiles$tiles), "\n")

    prediction_tiles[[i]] <- process_tile(tiles$tiles[[i]], dl_model)
  }
  combined <- terra::mosaic(terra::sprc(prediction_tiles))
  aligned <- terra::resample(combined, terra::rast(u), method = "near")
  result <- terra::mask(terra::crop(aligned, terra::rast(u[[3]])), terra::rast(u[[3]])) + 0.05
  r = raster::raster(result) * sigma
  r = r/raster::cellStats(r, stat = 'max')
  layers <- raster::stack(dtm2, slope2, roughness2, chm2, d2)
  names(layers) <- c("hillshade", "slope", "roughness", "chm", "density")

  return(list(r, layers))

}

#' @export
rasterize_conductivity4.LAScluster = function(las, centerline, dtm = NULL, water = NULL, param = alsroads_default_parameters2, dl_model= NULL)
{
  x <- lidR::readLAS(las)
  if (lidR::is.empty(x)) return(NULL)

  sigma <- rasterize_conductivity4(x, centerline, dtm, water, param, dl_model)
  sigma <- lidR:::raster_crop(sigma, lidR::st_bbox(las))
  return(sigma)
}

#' @export
rasterize_conductivity4.LAScatalog = function(las, centerline, dtm = NULL, water = NULL, param = alsroads_default_parameters2, dl_model= NULL)
{
  # Enforce some options
  if (lidR::opt_select(las) == "*") lidR::opt_select(las) <- "xyzciap"

  # Compute the alignment options including the case when res is a raster/stars/terra
  alignment <- lidR:::raster_alignment(1)

  if (lidR::opt_chunk_size(las) > 0 && lidR::opt_chunk_size(las) < 2*alignment$res)
    stop("The chunk size is too small. Process aborted.", call. = FALSE)

  # Processing
  options <- list(need_buffer = TRUE, drop_null = TRUE, raster_alignment = alignment, automerge = TRUE)
  output  <- lidR::catalog_apply(las, rasterize_conductivity4, centerline= centerline, dtm = dtm, water = water, param = param, dl_model, .options = options)
  return(output)
}


anisotropic_diffusion_filter = function(x, interation = 50, lambda = 0.2, k = 10)
{
  M = raster::as.matrix(x)
  M2 = anisotropic_diffusion(M*255, interation, lambda, k)
  M2 = M2/255
  y = x
  y[] = M2
  y
}

edge_enhancement = function(x, interation = 50, lambda = 0.05, k = 20, pass = 2)
{
  #smx = quantile(x[], na.rm = T, probs = 0.99)
  for (i in 1:pass)
    x = anisotropic_diffusion_filter(x, interation, lambda, k)

  x = raster::stretch(x, minv = 0, maxv = 1, maxq = 0.995)
  x
}

intensity_range_by_flightline <- function(las, res)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(las$PointSourceID)

  if (length(ids) == 1)
  {
    ans = rasterize_intensityrange(las, res)
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(las, PointSourceID == i)
    q = quantile(psi$Intensity, probs = 0.95)
    psi@data[Intensity>q, Intensity := q]
    ans[[as.character(i)]] <- rasterize_intensityrange(psi, res)
  }

  ans =terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans = terra::stretch(ans, maxv = 1)
  #plot(ans, col = heat.colors(50))
  ans <- terra::app(ans, fun = max, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

extract_raster_along_centerline <- function(centerline, data_raster, pixel_buffer = 8) {
  # Convert pixel buffer to actual distance using raster resolution
  res_val <- mean(res(data_raster))  # Assuming square or near-square cells
  dist_buffer <- pixel_buffer * res_val

  # Buffer the centerline
  buffered_line <- st_buffer(centerline, dist = dist_buffer)

  # Mask the CHM raster to the buffered area
  raster_crop <- crop(data_raster, buffered_line)
  raster_masked <- mask(raster_crop, buffered_line)

  return(raster_masked)
}

compute_segment_percentiles <- function(segment) {
  # Extract raster values and remove NA
  values <- raster::values(segment)
  values <- values[!is.na(values)]

  # Check if values are sufficient
  if (length(values) < 10) {
    warning("Not enough data to compute percentiles.")
    return(c(p5 = NA, p95 = NA))
  }

  # Return vector
  return(as.numeric(quantile(values, probs = c(0.05, 0.95), na.rm = TRUE)))
}

find_trough <- function(data, breaks = 50, degree = 2) {
  # Extract non-NA CHM values

  values <- data[]
  values <- values[!is.na(values)]
  hist(values, breaks = 50, probability = TRUE,
       main = " Histogram with Density", xlab = "Height (m)", col = "skyblue")
  lines(density(values, na.rm = TRUE), col = "red", lwd = 2)
  hist_data <- hist(values, breaks = breaks, plot = FALSE)
  counts <- hist_data$counts
  mids <- hist_data$mids

  # Find index of histogram peak
  max_idx <- which.max(counts[5:length(counts)]) + 4

  # Fit polynomial to counts before peak
  fit <- lm(counts[1:max_idx] ~ poly(mids[1:max_idx], degree = degree))

  # Predict on same mids
  y_pred <- predict(fit, newdata = data.frame(mids[1:max_idx]))

  # Find index of minimum
  min_idx <- which.min(y_pred)+1
  trough_x <- mids[min_idx]

  return(trough_x)
}

density_gnd_by_flightline1 <- function(las, res, drop_angles = 3)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(las$PointSourceID)
  gnd = filter_ground(las)
  angle_max = max(abs(las$ScanAngle))
  gnd = filter_poi(gnd, abs(ScanAngle) < angle_max-drop_angles)

  if (length(ids) == 1)
  {
    ans = rasterize_density(gnd, res)
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(gnd, PointSourceID == i)
    d <- rasterize_density(psi, res)
    d[d == 0] = NA
    q = quantile(d[], probs = 0.95, na.rm = TRUE)
    d[d>q] = q
    ans[[as.character(i)]] = d
  }

  ans = terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans = terra::stretch(ans, maxv = 1)
  #plot(ans, col = heat.colors(50))
  ans <- terra::app(ans, fun = max, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

density_gnd_by_flightline <- function(las, res, drop_angles = 0)
{
  .N <- PointSourceID <- NULL
  las_filtered <- filter_poi(las, abs(ScanAngle) < (max(abs(las$ScanAngle)) - drop_angles))

  # Separate ground and non-ground
  gnd <- filter_ground(las_filtered)
  veg <- filter_poi(las_filtered, Classification != 2)  # Adjust based on your classification

  ids <- unique(las_filtered$PointSourceID)
  res <- terra::rast(res)

  ans <- vector("list", 0)

  for (i in ids) {
    # Get points for this flightline
    psi_gnd <- filter_poi(gnd, PointSourceID == i)
    psi_veg <- filter_poi(veg, PointSourceID == i)

    # Calculate densities
    d_gnd <- rasterize_density(psi_gnd, res)
    d_veg <- rasterize_density(psi_veg, res)
    d_total <- d_gnd + d_veg

    # Calculate ground ratio
    d_ratio <- d_gnd / d_total

    # Handle division by zero
    d_ratio[d_total == 0] <- NA

    # Clip outliers (optional)
    d_ratio[d_ratio > 1] <- 1  # Shouldn't happen but just in case
    d_ratio[d_ratio < 0] <- 0

    ans[[as.character(i)]] <- d_ratio
  }

  ans <- terra::rast(ans)
  # Take the maximum ratio across flightlines
  ans <- terra::app(ans, fun = max, na.rm = TRUE)

  return(raster::raster(ans))
}

density_lp_by_flightline <- function(las, res)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(las$PointSourceID)
  z1  <- 0.5
  z2  <- 3
  tmp <- lidR::filter_poi(las, Z > z1, Z < z2)

  if (length(ids) == 1)
  {
    ans = rasterize_density(tmp, res)
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids)
  {
    psi <- lidR::filter_poi(tmp, PointSourceID == i)
    d <- rasterize_density(psi, res)
    d[d == 0] = NA
    q = quantile(d[], probs = 0.95, na.rm = TRUE)
    d[d>q] = q
    ans[[as.character(i)]] = d
  }

  ans = terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans = terra::stretch(ans, maxv = 1)
  #plot(ans, col = heat.colors(50))
  ans <- terra::app(ans, fun = max, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

low_points_by_flightline <- function(nlas, res, th = 2)
{
  .N <- PointSourceID <- NULL

  z1 <- 0.5
  z2 <- 3
  res <- terra::rast(res)
  ids <- unique(nlas$PointSourceID)

  if (length(ids) == 1)
  {
    tmp <- lidR::filter_poi(nlas, Z > z1, Z < z2)
    ans = lidR:::rasterize_fast(tmp, res, 0, "count", pkg = "terra")
    return(raster::raster(ans > th))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(nlas, PointSourceID == i, Z > z1, Z < z2)
    d <- lidR:::rasterize_fast(psi, res, 0, "count", pkg = "terra")
    d[is.na(d)] = 0
    ans[[as.character(i)]] = d <= 3
  }

  ans = terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans <- terra::app(ans, fun = min, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

chm_by_flightline <- function(nlas, res)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(nlas$PointSourceID)

  if (length(ids) == 1)
  {
    ans = rasterize_canopy(nlas, res, p2r())
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(nlas, PointSourceID == i)
    d <- rasterize_canopy(psi, res, lidR::p2r())
    ans[[as.character(i)]] = d
  }

  ans = terra::rast(ans)
  ans = terra::mask(ans, terra::vect(water), inverse = T)
  #plot(ans, col = height.colors(50))
  ans = terra::stretch(ans, minv = 0)
  #plot(ans, col = height.colors(50))
  ans <- terra::app(ans, fun = min, na.rm = TRUE)
  #plot(ans, col = height.colors(50))
  ans =  ans - min(ans[], na.rm = TRUE)
  #plot(ans, col = height.colors(50))

  return(raster::raster(ans))
}


rasterize_intensityrange <- function(las, res)
{
  # Detect outliers of intensity and change their value. This is not perfect but avoid many troubles
  Intensity <- NULL
  outliers <- as.integer(stats::quantile(las$Intensity, probs = 0.98))
  las@data[Intensity > outliers, Intensity := outliers]

  # Switch Z and Intensity trick to use fast lidR internal function
  Z <- las[["Z"]]
  las@data[["Z"]] <-  las@data[["Intensity"]]
  imax <- lidR:::rasterize_fast(las, res, 0, "max")
  imin <- lidR:::rasterize_fast(las, res, 0, "min")
  irange <- imax - imin
  return(irange)
}
tile_raster <- function(r, target_size = 256) {
  nr <- nrow(r)
  nc <- ncol(r)

  # number of tiles needed
  n_tiles_rows <- ceiling((nr / target_size)-0.1)
  n_tiles_cols <- ceiling((nc / target_size)-0.1)

  tiles <- list()
  tile_coords <- list()

  # loop over rows and cols of tiles
  for (i in 1:n_tiles_rows) {
    for (j in 1:n_tiles_cols) {

      # compute pixel indices
      row_start <- (i - 1) * target_size + 1
      col_start <- (j - 1) * target_size + 1
      row_end <- min(row_start + target_size , nr)
      col_end <- min(col_start + target_size , nc)

      # compute spatial extent
      # row/col -> xy
      xy_top_left     <- xyFromCell(r, cellFromRowCol(r, row_start, col_start))
      xy_bottom_right <- xyFromCell(r, cellFromRowCol(r, row_end, col_end))

      tile_extent <- ext(
        xy_top_left[1],     # xmin
        xy_bottom_right[1], # xmax
        xy_bottom_right[2], # ymin
        xy_top_left[2]      # ymax
      )

      # crop to extent
      tile <- terra::crop(r, tile_extent)

      # store
      tiles[[length(tiles) + 1]] <- tile
      tile_coords[[length(tile_coords) + 1]] <- list(
        row_start = row_start,
        row_end = row_end,
        col_start = col_start,
        col_end = col_end,
        extent = tile_extent
      )
    }
  }

  return(list(tiles = tiles, coords = tile_coords, original_raster = r))
}
resize_conductivity <- function(r, target_size = 256) {
  # Replace NA with 0
  r[is.na(r)] <- 0

  nr <- nrow(r)
  nc <- ncol(r)
  nlyrs <- terra::nlyr(r)  # Get number of layers - corrected function name
  # r <- terra::rast(r)
  # Case 1: Larger → crop center
  if (nr >= target_size && nc >= target_size) {
    if (nr == target_size && nc == target_size) {
      # Return the original raster if already the target size
      result <- r
    } else {
      row_start <- floor((nr - target_size) / 2) + 1
      col_start <- floor((nc - target_size) / 2) + 1
      row_end <- row_start + target_size - 1
      col_end <- col_start + target_size - 1

      # Get spatial coordinates of crop window
      xy1 <- terra::xyFromCell(r, cellFromRowCol(r, row_start, col_start))
      xy2 <- terra::xyFromCell(r, cellFromRowCol(r, row_end, col_end))

      # Build extent (xmin, xmax, ymin, ymax)
      e <- terra::ext(xy1[1], xy2[1], xy2[2], xy1[2])

      result <- terra::crop(r, e)
    }
  } else {
    # Case 2: Smaller → pad with zeros
    new_r <- terra::rast(nrow = target_size, ncol = target_size, nlyr = nlyrs,
                  crs = st_crs(r))

    # Align extent to be same as original, centered
    ext_orig <- ext(r)
    resx <- xres(r)
    resy <- yres(r)
    width <- target_size * resx
    height <- target_size * resy
    cx <- (xmin(ext_orig) + xmax(ext_orig)) / 2
    cy <- (ymin(ext_orig) + ymax(ext_orig)) / 2
    ext_new <- ext(cx - width/2, cx + width/2, cy - height/2, cy + height/2)
    terra::ext(new_r) <- ext_new

    # Reproject original into new canvas (no resampling, nearest neighbor)
    # new_r_raster <- raster(new_r)
    result <- terra::resample(r, new_r, method = "near")

    # Replace NA with 0 again (from padding)
    result[is.na(result)] <- 0
  }

  return(result)
}
prepare_raster_for_prediction <- function(raster_obj) {
  # Resize to 256x256


  # Convert raster to array
  raster_array <- terra::as.array(raster_obj)

  # Check dimensions and reshape if needed
  dims <- dim(raster_array)

  # Expected shape: (256, 256, 11)
  if (length(dims) == 2) {
    # Single layer raster - add third dimension
    raster_array <- array(raster_array, dim = c(dims[1], dims[2], 1))
  } else if (length(dims) == 3 && dims[3] != 9) {
    # Wrong number of layers - take first 11 or pad
    if (dims[3] > 9) {
      raster_array <- raster_array[, , 1:9]
    } else {
      # Pad with zeros
      padded_array <- array(0, dim = c(dims[1], dims[2], 9))
      padded_array[, , 1:dims[3]] <- raster_array
      raster_array <- padded_array
    }
  }

  # Add batch dimension: (1, 256, 256, 11)
  raster_array <- array(raster_array, dim = c(1, dim(raster_array)))

  return(raster_array)
}
convert_predictions_to_raster <- function(predictions_array, original_raster) {
  # Remove batch dimension and get the prediction data
  # Assuming predictions have shape: (1, 256, 256, 1) or similar
  prediction_data <- predictions_array[1,1 , , ]  # Remove batch dim and channel dim

  # Create a template raster with the same spatial properties as original
  # but with 256x256 dimensions (from model output)
  template <- terra::rast(
    nrows = 512,
    ncols = 512,
    extent = terra::ext(original_raster),  # Keep original extent
    crs = st_crs(original_raster)      # Keep original CRS
  )

  # Fill the template with prediction data
  values(template) <- prediction_data

  return(template)
}
process_tile <- function(tile, model_combined1) {

  # Global normalisation constants — must match Python training script exactly
  SLOPE_GLOBAL_MAX <- 25.0   # degrees — forestry road operational limit
  CHM_GLOBAL_MAX   <- 40.0   # metres  — max canopy height
  SLOPE_BAND       <- 2      # 2nd band in R (1-indexed) = slope
  CHM_BAND         <- 7      # 7th band in R (1-indexed) = CHM

  u_normalized <- tile

  for (i in 1:terra::nlyr(u_normalized)) {
    layer_vals <- terra::values(u_normalized[[i]])

    if (i == SLOPE_BAND) {
      # Global fixed normalisation — preserves absolute slope meaning
      u_normalized[[i]] <- terra::clamp(u_normalized[[i]] / SLOPE_GLOBAL_MAX,
                                        lower = 0, upper = 1)

    } else if (i == CHM_BAND) {
      # Global fixed normalisation — preserves absolute canopy height meaning
      u_normalized[[i]] <- terra::clamp(u_normalized[[i]] / CHM_GLOBAL_MAX,
                                        lower = 0, upper = 1)

    } else {
      # Per-tile min-max for all other bands
      layer_min <- min(layer_vals, na.rm = TRUE)
      layer_max <- max(layer_vals, na.rm = TRUE)
      u_normalized[[i]] <- (u_normalized[[i]] - layer_min) / (layer_max - layer_min + 1e-5)
    }

    u_normalized[[i]] <- terra::app(u_normalized[[i]], fun = as.numeric)
  }

  # Zero out NAs after normalisation — matches training script behaviour
  u_normalized[is.na(u_normalized)] <- 0

  resized    <- resize_conductivity(u_normalized, 512)
  input      <- prepare_raster_for_prediction(resized)
  predictions <- model_combined1$predict(input)
  prediction_raster <- convert_predictions_to_raster(predictions, resized)
  return(prediction_raster)
}
# process_tile <- function(tile, model_combined1) {
#   # Resize to exactly 256x256 (pad if smaller)
#
#   u_normalized <- tile
#
#   for (i in 1:terra::nlyr(u_normalized)) {
#     layer_vals <- values(u_normalized[[i]])
#     layer_max <- max(layer_vals, na.rm = TRUE)
#     layer_min <- min(layer_vals, na.rm = TRUE)
#
#     u_normalized[[i]] <- (u_normalized[[i]]-layer_min) / (layer_max-layer_min + 1e-5)
#
#     # if (layer_max > 1) {
#     #   u_normalized[[i]] <- (u_terra[[i]]-layer_min) / (layer_max-layer_min + 1e-5)
#     # }
#
#     # Convert to float
#     u_normalized[[i]] <- terra::app(u_normalized[[i]], fun = as.numeric)
#   }
#
#   tile = u_normalized
#   tile[is.na(tile)] <- 0
#   resized <- resize_conductivity(tile,  512)
#   input <- prepare_raster_for_prediction(resized)
#   predictions <- model_combined1$predict(input)
#   prediction_raster <- convert_predictions_to_raster(predictions, resized)
#   return(prediction_raster)
# }
reassemble_predictions <- function(prediction_tiles, tile_coords, original_raster) {
  # Create empty result raster with same dimensions as original
  result <- rast(
    nrows = nrow(original_raster),
    ncols = ncol(original_raster),
    nlyrs = 1,
    extent = ext(original_raster),
    crs = st_crs(original_raster)
  )

  # Initialize with NA
  result[] <- NA

  # Place each prediction tile in the correct position
  for (i in 1:length(prediction_tiles)) {
    coords <- tile_coords[[i]]
    pred_tile <- prediction_tiles[[i]]

    # If tile was padded, crop back to original tile size
    tile_rows <- coords$row_end - coords$row_start
    tile_cols <- coords$col_end - coords$col_start

    if (tile_rows < 512 || tile_cols < 512) {
      pred_tile <- pred_tile[1:tile_rows, 1:tile_cols]
    }

    # Insert into result raster
    result[coords$row_start:coords$row_end-1,
           coords$col_start:coords$col_end-1] <- values(pred_tile)
  }

  return(result)
}

