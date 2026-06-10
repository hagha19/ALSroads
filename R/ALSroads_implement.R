ibrary(lidR)
library(sf)
library(raster)
library(mapview)
library(leaflet)
library(ggplot2)
library(reticulate)
library(parallel)
library(dplyr)

Sys.setlocale("LC_ALL", "English_United States.UTF-8")
devtools::load_all()

alsroads_default_parameters2 <- alsroads_default_parameters

# --- Data loading (unchanged) ---
ctg   <- readLAScatalog("E:/Dryden/las/drydenFMU/",
                        filter = "-max_depth 3 -drop_class 7 18", progress = TRUE)
roads <- read_sf("H:/Sarah/dryden_test_roads.shp")
roads <- st_transform(roads, st_crs(ctg))
dtm   <- raster("E:/Dryden/dtm/dtm.vrt")
water <- st_read("E:/Dryden/maps/water/Dryden_Waterbodies.gpkg")
water <- st_transform(water, st_crs(roads))

# --- Parameters (unchanged) ---
options(ALSroads.debug.finding = FALSE, ALSroads.debug.metrics = FALSE,
        ALSroads.debug.measuring = FALSE)
alsroads_default_parameters2$state$percentage_veg_thresholds  <- c(55, 70)
alsroads_default_parameters2$state$conductivity_thresholds     <- c(0.15, 0.25)
alsroads_default_parameters2$state$shoulder_thresholds         <- c(20, 40)
alsroads_default_parameters2$state$drivable_width_thresholds   <- c(0.5, 2)
alsroads_default_parameters2$extraction$road_buffer            <- 80
alsroads_default_parameters2$conductivity$s                    <- c(3, 15)
alsroads_default_parameters2$conductivity$e                    <- 20
alsroads_default_parameters2$conductivity$alpha$d              <- 2
alsroads_default_parameters2$conductivity$alpha$h              <- 1.5
alsroads_default_parameters2$conductivity$alpha$r              <- 1.5
alsroads_default_parameters2$conductivity$alpha$e              <- 1
alsroads_default_parameters2$conductivity$alpha$s              <- 1
alsroads_default_parameters2$conductivity$alpha$i              <- 0.25
alsroads_default_parameters2$constraint$confidence             <- 0.4

conductivity  <- "v4"
param         <- alsroads_default_parameters2
return_all    <- TRUE
return_stack  <- TRUE
roads         <- st_cast(roads, "LINESTRING")

# =============================================================================
# STEP 1 — Load DL model in main session (Python only needed here)
# =============================================================================
Sys.setenv(RETICULATE_PYTHON = "C:/Users/HAGHA19/AppData/Local/anaconda3/envs/tf_env/python.exe")
use_condaenv("tf_env", required = TRUE)

correct_path    <- normalizePath("C:/Users/HAGHA19/OneDrive - Université Laval/01_Projects/01_Ontario_roads_project/outputs/ALSroads_DL/code/", winslash = "/")
checkpoint_path <- normalizePath("C:/Users/HAGHA19/OneDrive - Université Laval/01_Projects/01_Ontario_roads_project/outputs/ALSroadsdl run output1/saved_models/best_road_model.pth", winslash = "/")

py_run_string(paste0(
  "import sys\nsys.path.append('", correct_path, "')\n",
  "import torch\nimport numpy as np\n",
  "sys.modules.pop('pytorch_model_loader_ssl', None)\n",
  "from pytorch_model_loader_ssl_new import AttentionResUNet\n"
))
py_run_string(paste0(
  "checkpoint_filepath = '", checkpoint_path, "'\n",
  "model_base = AttentionResUNet(in_ch=9, nc=1, drop=0.0, bn=True)\n",
  "model_base.load_state_dict(torch.load(checkpoint_filepath, map_location='cpu', weights_only=True))\n",
  "class ModelWrapper:\n",
  "    def __init__(self, model, device='cpu'):\n",
  "        self.model = model.to(device); self.device = torch.device(device)\n",
  "    def predict(self, x):\n",
  "        self.model.eval()\n",
  "        with torch.no_grad():\n",
  "            if isinstance(x, np.ndarray): x = torch.from_numpy(x).float()\n",
  "            if x.ndim == 4 and x.shape[3] == 9: x = x.permute(0,3,1,2)\n",
  "            seg, *_ = self.model(x)\n",
  "            return seg.cpu().numpy()\n",
  "model_combined = ModelWrapper(model_base, device='cpu')\n"
))
dl_model <- py$model_combined

# =============================================================================
# STEP 2 — Precompute DL prediction over ALL road buffers (run ONCE)
# =============================================================================

# Set paths — change these to suit your project
tmp_dir     <- "H:/Sarah/dl_chunks_dryden"     # intermediate per-tile tifs
dl_pred_tif <- "H:/Sarah/dl_prediction_dryden.tif"  # final merged output

if (!file.exists(dl_pred_tif)) {
  rasterize_conductivity4_precompute(
    ctg      = ctg,
    roads    = roads,
    dtm      = dtm,
    dl_model = dl_model,
    water    = water,
    param    = alsroads_default_parameters2,
    out_tif  = dl_pred_tif,
    tmp_dir  = tmp_dir,
    chunk_buffer = 50
  )
} else {
  message("DL prediction raster already exists — skipping precomputation.")
}

# =============================================================================
# STEP 3 — Parallel road processing (NO Python needed in workers)
# =============================================================================

cl <- parallel::makeCluster(9)
on.exit(parallel::stopCluster(cl))

clusterEvalQ(cl, {
  library(lidR); library(sf); library(raster); library(ALSroads)
})

clusterExport(cl, varlist = c(
  "return_all", "ctg", "roads", "dtm", "water",
  "param", "conductivity", "return_stack",
  "dl_pred_tif"     # just a file path string — tiny to serialize
))

measure_roads1 <- function(i) {
  Sys.setlocale("LC_ALL", "French_France.UTF-8")
  tryCatch({
    withCallingHandlers({
      ALSroads::measure_road(
        ctg, roads[i, , drop = FALSE],
        dtm,
        conductivity  = conductivity,   # still "v4" — triggers the v4 code path
        water         = water,
        param         = param,
        return_all    = return_all,
        return_stack  = return_stack,
        dl_model      = NULL,           # not needed — precomputed raster is used
        dl_pred_raster = dl_pred_tif    # NEW: path to precomputed tif
      )
    }, warning = function(w) {
      warning(sprintf("Warning in road %d: %s", i, conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    warning(sprintf("Error in road %d: %s", i, conditionMessage(e)))
    NULL
  })
}

clusterExport(cl, varlist = "measure_roads1")

system.time(res <- parLapply(cl, 1:nrow(roads), measure_roads1))

res_clean   <- res[!sapply(res, is.null)]
common_cols <- Reduce(intersect, lapply(res_clean, names))
res1        <- do.call(rbind, lapply(res_clean, function(x) x[, common_cols]))
st_write(res1, "H:/Sarah/Dryden_DL_output.shp", delete_dsn = TRUE)
