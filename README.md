![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)

# ALSroads

ALSroads provides tools to relocate, measure, and estimate the condition of forestry roads from an inaccurate map using airborne LiDAR point clouds.

An early version of the algorithm is described in:

> Roussel, J., Bourdon, J., Morley, I. D., Coops, N. C., & Achim, A. (2022). Correction, update, and enhancement of vectorial forestry road maps using ALS data, a pathfinder, and seven metrics. *International Journal of Applied Earth Observation and Geoinformation*, 114, 103020. https://doi.org/10.1016/j.jag.2022.103020

---

## Installation

```r
remotes::install_github("hagha19/ALSroads")
```

### R dependencies

```r
install.packages(c("lidR", "sf", "raster", "terra", "lwgeom",
                   "dplyr", "data.table", "glue", "reticulate", "httr"))
```

---

## Deep Learning Extension (optional)

The package supports an Attention ResUNet deep learning model (`conductivity = "v4"`) for improved road detection. This requires additional setup.

### 1 — Download the model files

Two files must be downloaded from HuggingFace and saved locally:

**Model checkpoint** (~140 MB):
```
https://huggingface.co/hagha-19/alsroads-AttResUNet-model/resolve/main/best_road_model.pth
```

**Python model loader script**:
```
https://huggingface.co/hagha-19/alsroads-AttResUNet-model/resolve/main/pytorch_model_loader_ssl_new.py
```

Or download both automatically in Python:

```python
from huggingface_hub import hf_hub_download

hf_hub_download(
    repo_id   = "hagha-19/alsroads-AttResUNet-model",
    filename  = "best_road_model.pth",
    repo_type = "model",
    local_dir = "./models"
)

hf_hub_download(
    repo_id   = "hagha-19/alsroads-AttResUNet-model",
    filename  = "pytorch_model_loader_ssl_new.py",
    repo_type = "model",
    local_dir = "./scripts"
)
```

### 2 — Set up a Python environment

```bash
# Create and activate a virtual environment
python -m venv env_alsroads
source env_alsroads/bin/activate        # Linux / Mac
env_alsroads\Scripts\activate           # Windows

# Install PyTorch and numpy
pip install torch torchvision numpy
```

### 3 — Set paths in your R script

```r
PYTHON_EXE      <- "path/to/your/python.exe"
PYTHON_CODE_DIR <- "path/to/folder/containing/pytorch_model_loader_ssl_new.py"
CHECKPOINT_PATH <- "path/to/best_road_model.pth"
```

---

## Usage

```r
library(ALSroads)
library(lidR)
library(sf)
library(raster)

# ── Load data ────────────────────────────────────────────────────
ctg   <- readLAScatalog("path/to/las/folder/",
                         filter = "-max_depth 3 -drop_class 7 18")
roads <- st_read("path/to/roads.shp") |>
           st_cast("LINESTRING") |>
           st_transform(st_crs(ctg))
dtm   <- raster("path/to/dtm.vrt")

# ── Parameters ───────────────────────────────────────────────────
param <- alsroads_default_parameters
param$constraint$confidence <- 0.4
param$extraction$road_buffer <- 80

# ── Classical conductivity (no deep learning) ─────────────────────
res <- measure_roads(ctg, roads, dtm,
                     conductivity = "v2",
                     param        = param)

st_write(res, "roads_corrected.shp")
```

### With deep learning (`v4`)

```r
library(reticulate)

# Paths to downloaded files
PYTHON_EXE      <- "path/to/python.exe"
PYTHON_CODE_DIR <- "path/to/scripts/"
CHECKPOINT_PATH <- "path/to/models/best_road_model.pth"

use_python(PYTHON_EXE, required = TRUE)

py_run_string(paste0(
  "import sys, torch, numpy as np\n",
  "sys.path.append('", PYTHON_CODE_DIR, "')\n",
  "from pytorch_model_loader_ssl_new import AttentionResUNet\n",
  "model = AttentionResUNet(in_ch=9, nc=1, drop=0.0, bn=True)\n",
  "model.load_state_dict(torch.load('", CHECKPOINT_PATH, "',\n",
  "  map_location='cpu', weights_only=True))\n",
  "class ModelWrapper:\n",
  "    def __init__(self, model):\n",
  "        self.model = model.eval()\n",
  "    def predict(self, x):\n",
  "        with torch.no_grad():\n",
  "            if isinstance(x, np.ndarray): x = torch.from_numpy(x).float()\n",
  "            if x.ndim == 4 and x.shape[3] == 9: x = x.permute(0,3,1,2)\n",
  "            seg, *_ = self.model(x)\n",
  "            return seg.cpu().numpy()\n",
  "dl_model = ModelWrapper(model)\n"
))

res <- measure_roads(ctg, roads, dtm,
                     conductivity = "v4",
                     param        = param,
                     dl_model     = py$dl_model)

st_write(res, "roads_corrected_dl.shp")
```

---

## Road condition classes

| Class | Label | Description |
|---|---|---|
| 1 | Maintained | Year-round access |
| 2 | Passable | Passable with care |
| 3 | Degraded | Seasonal access only |
| 4 | Barely passable | Specialist vehicle required |
| 5 | Ghost / reclaimed | Road no longer exists |

---

## Output fields

| Field | Description |
|---|---|
| `ROADWIDTH` | Total road width (m) |
| `DRIVABLEWIDTH` | Drivable width (m) |
| `PERCABOVEROAD` | % LiDAR points above road surface |
| `SHOULDERS` | Shoulder presence (%) |
| `SINUOSITY` | Road sinuosity index |
| `CONDUCTIVITY` | Road surface conductivity score |
| `SCORE` | Overall road condition score (0–100) |
| `CLASS` | Road condition class (1–5) |
