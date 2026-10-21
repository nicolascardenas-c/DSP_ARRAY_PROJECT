# Direction of Arrival Estimation — Bartlett, MVDR and MUSIC

**TU Delft | EE4795 — Statistical Digital Signal Processing**  
**Project:** Array Processing for Communication Systems

## Overview

This project compares three Direction of Arrival (DoA) estimators for a uniform linear antenna array (ULA):

- **Bartlett** — conventional beamforming.
- **MVDR (Capon)** — adaptive minimum-variance beamforming.
- **MUSIC** — noise-subspace-based high-resolution estimation.

The goal is to compare angular resolution and robustness to additive noise for two source separations: **10°** and **50°**.

## Data and experimental setup

The supplied datasets are `spcom_10sep.mat` and `spcom_50sep.mat` (included in `spcom.zip`, if redistribution is permitted). Each contains received complex-valued array observations with:

- **5 antennas**, **20 snapshots**, **2 sources** and **half-wavelength antenna spacing**.
- Reference angles inferred from the data: **(20°, 30°)** and **(-20°, 30°)**, respectively.

The script performs three experiments:

1. **Original data:** spatial spectra, estimated DoAs and covariance eigenvalues.
2. **Added noise:** the same analyses at nominal SNRs of **30, 20, 10 and 0 dB** using spatially white complex Gaussian noise.
3. **Monte Carlo:** **100 independent noise realizations** per SNR, reporting DoA success rate and angular RMSE.

A trial is successful if both estimated DoAs are within **2°** of their respective reference angles. RMSE is calculated only for trials where both peaks are detected; unsuccessful two-peak estimates are still included. A **3 dB minimum peak prominence** is used for peak selection.

## Requirements

- MATLAB (with support for local functions in scripts).
- No additional toolboxes are required.
- `spcom.zip` **or** both extracted `.mat` files, in the same folder as the MATLAB script.

## How to run

1. Download or clone this repository.
2. Ensure `DSP_DoA_no_snapshots.m` and `spcom.zip` (or both `.mat` files) are in the same folder.
3. Open **`DSP_DoA_no_snapshots.m`** in MATLAB and click **Run**, or execute:

   ```matlab
   run('DSP_DoA_no_snapshots.m')
   ```

The script automatically extracts `spcom.zip` if necessary. It generates MATLAB figures of the spatial spectra and covariance eigenvalues, and prints the DoA estimates and Monte Carlo metrics in the Command Window. Fixed random seeds are used for reproducibility.

## Implementation notes

- The sample spatial covariance is computed as `R = X*X'/N`.
- MVDR uses `R\A` (a linear-system solve), **without diagonal loading**.
- MUSIC estimates the noise subspace from the eigenvectors associated with the smallest covariance eigenvalues.
- Angular scan step: **0.01°** for plots and **0.1°** for Monte Carlo, to reduce execution time.
- Peaks are detected without `findpeaks` or the Signal Processing Toolbox.

The report provides the mathematical formulation, results and interpretation of the experiments.
