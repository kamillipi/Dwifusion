# DWIFusion – MATLAB App for IVIM/IVCM Analysis

**DWIFusion** is a MATLAB App Designer application for fitting and visualizing **IntraVoxel Incoherent Motion (IVIM)**  model from diffusion-weighted MRI (DWI) data.

The app integrates fitting routines with an interactive GUI for exploring raw DWI data, applying masks, running parameter estimation, and visualizing fitted IVIM maps.

---

## Features

- **Data loading**
  - Import NIfTI (`.nii`, `.nii.gz`) diffusion datasets
  - Automatically or manually load b-values (`.bval` or `.txt`)

- **Fitting methods**
  - `1step`: Nonlinear simultaneous fit of IVIM parameters
  - `segmented`: Two-step fit (tissue diffusion first, then perfusion)
  - `grid`: Two-step Grid search fitting with parameter constraints
  - (extra - available only in command line) `IVCM_grid`/`IVCM_seg` IntraVoxelCoherentMotion velocity parameter `V` fitting.
        reuqires  gradient timing (`large_delta` - separation, `small_delta` - duration)

- **Interactive viewer**
  - Orthoslice visualization of raw data and parameter maps (`S0`, `f`, `D*`, `D`, `f·D*`)
  - Adjustable window/level (WL/WW) controls
  - Crosshair-linked parameter readout

- **Analysis tools**
  - On-click **signal plot** with measured data vs. fitted IVIM model
  - Export fitted results as **NIfTI** maps for external use
  - Supports optional ROI masking

- **Parameter outputs**
  - `S0`   – baseline signal  
  - `f`    – perfusion fraction  
  - `D*`   – pseudo-diffusion coefficient  
  - `D`    – diffusion coefficient  
  - `V`    – vascular velocity (IVCM only)  

---

## Installation

1. Clone or download this repository:
   ```bash
   git clone https://github.com/kamillipi/DWIFusion.git
   ```

2. Open MATLAB and add the repository folder (and subfolders) to your MATLAB path:
   ```matlab
   addpath(genpath('path/to/DWIFusion'))
   ```

3. Launch the app:
   ```matlab
   app = DWIfusion;
   ```

---

## Usage

**Load Data**
   - From the *File → Load raw* menu, select a diffusion-weighted NIfTI file and corresponding b-values file.
   - Optionally load a binary mask (*Process → Load mask*).

**Choose Fitting Method**
   - Select `1step`, `segmented`, or `grid` from the **Method** dropdown.
   - Adjust the **b-split value** (default: 250 s/mm²) if using segmented/grid methods.

**Run Fitting**
   - Click **Calculate** to estimate parameter maps from the raw DWI data.
   - Results will appear in the **Analysis** tab and can be explored interactively.

**OR**

**Import Precomputed Results**
   - If you already have IVIM/IVCM parameter maps saved as a NIfTI file, load them via  
     *File → Load results*.  
   - This allows you to:
     - Skip fitting and directly explore parameter maps (`S0`, `f`, `D*`, `D`, `f·D*`).
     - Compare different fitting strategies or datasets.
     - Reuse previously computed results without recalculating.

**After these steps you will be able to**

4. **Export Results**
   - Save fitted or imported parameter maps as NIfTI files via  
     *File → Export results → To NIfTI*.

5. **Visualize Results**
   - Switch between parameter maps (`S0`, `f`, `D*`, `D`, `f·D*`, raw data).
   - Adjust WL/WW interactively.
   - Click **Show signal plot** to view model fit vs. measured data at the crosshair (only if raw data is available).

---

## Project Structure

- `DWIfusion.m`  
  MATLAB App Designer class for the GUI.

- `get_ivim_v5.m`  
  Main entry point for IVIM/IVCM fitting (handles parsing, masking, and dispatching to fitting routines).

- `fit_IVIM_1step.m`  
  Nonlinear least-squares IVIM fit.

- `fit_IVIM_segmented.m`  
  Segmented IVIM fitting (D/S0 first, then D*, f).

- `fit_IVIM_grid.m`  
  Grid search IVIM fitting.

- `fit_IVCM_segmented.m`  
  Segmented IVCM fitting.

- `fit_IVCM_grid.m`  
  Grid search IVCM fitting.

- `process_mask.m`  
  Utility for applying or computing masks.

- `save_IVIM.m`  
  Utility for saving fitted maps as NIfTI.

---

## Notes

- **Data requirements**
  - Ensure the number of b-values matches the number of DWI volumes.
  - Units:
    - b-values in s/mm²
    - D, D* in mm²/s
    - V in mm/s
    - Δ, δ in seconds

- **Performance**
  - Grid search methods can be computationally expensive; masking is recommended.
  - Parallelization (`parfor`) is used to speed up voxel-wise fitting.

---

## Requirements

- MATLAB R2021a or newer (App Designer, Parallel Computing Toolbox recommended)
- Image Processing Toolbox (for NIfTI I/O and visualization)


---

## Acknowledgements

This app builds upon standard IVIM and novel IVCM models in diffusion MRI literature and was developed to support both research and teaching applications.
