# utils/visualization.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

def plot_capillary_traces(
        heatmap: np.ndarray,
        extent: tuple,
        capillary_amount: int,
        smoothing_factor: float
    ):
    """
    Plot summed capillary traces from the 2D‐gel heatmap.

    - First shows all capillaries overlaid.
    - Then shows each in its own figure.

    Parameters
    ----------
    heatmap : 2D numpy array
    extent  : (pI_min, pI_max, lnMW_min, lnMW_max)
    capillary_amount : number of vertical segments to split into
    smoothing_factor : sigma for Gaussian smoothing
    """
    pi_min, pi_max, lnmw_min, lnmw_max = extent
    n_rows, n_cols = heatmap.shape

    # reconstruct MW axis in kDa
    lnmw = np.linspace(lnmw_min, lnmw_max, n_rows)
    mw_vals = np.exp((lnmw + 6.4014) / 5.3779)

    # define column ranges
    per = n_cols // capillary_amount
    segments = [
        range(i*per, (i+1)*per if i<capillary_amount-1 else n_cols)
        for i in range(capillary_amount)
    ]

    # combined overlay
    plt.figure(figsize=(8, 5))
    for idx, cols in enumerate(segments, 1):
        prof = heatmap[:, list(cols)].sum(axis=1)
        prof_s = gaussian_filter1d(prof, sigma=smoothing_factor)
        plt.plot(mw_vals, prof_s, label=f"Capillary {idx}")
    plt.xlabel("Molecular Weight (kDa)")
    plt.ylabel("Aggregated Abundance")
    plt.title("Combined Capillary Traces")
    plt.legend(loc='best', fontsize='small', ncol=2)
    plt.tight_layout()
    plt.show()

    # individual
    for idx, cols in enumerate(segments, 1):
        prof = heatmap[:, list(cols)].sum(axis=1)
        prof_s = gaussian_filter1d(prof, sigma=smoothing_factor)
        plt.figure(figsize=(6, 4))
        plt.plot(mw_vals, prof_s)
        plt.xlabel("Molecular Weight (kDa)")
        plt.ylabel("Aggregated Abundance")
        plt.title(f"Capillary {idx}")
        plt.tight_layout()
        plt.show()
