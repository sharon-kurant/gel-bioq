# utils/visualization.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

def plot_capillary_traces(
        heatmap: np.ndarray,
        extent: tuple,
        capillary_amount: int,
        smoothing_factor: float
    ) -> tuple[list, list, np.ndarray]:
    """
    Compute and plot capillary traces, returning figures and data.

    Parameters
    ----------
    heatmap : 2D numpy array
        The gel intensity grid.
    extent : (pI_min, pI_max, lnMW_min, lnMW_max)
    capillary_amount : int
        Number of vertical pI segments.
    smoothing_factor : float
        Sigma for Gaussian smoothing.

    Returns
    -------
    figs : list of matplotlib.figure.Figure
        [combined_overlay_fig, cap1_fig, cap2_fig, …]
    capillary_data : list of 1D numpy arrays
        Smoothed abundance profiles for each capillary.
    x_values : 1D numpy array
        Molecular weight axis (kDa) corresponding to profiles.
    """
    pi_min, pi_max, lnmw_min, lnmw_max = extent
    n_rows, n_cols = heatmap.shape

    # rebuild MW axis in kDa
    lnmw = np.linspace(lnmw_min, lnmw_max, n_rows)
    x_values = np.exp((lnmw + 6.4014) / 5.3779)

    # split into capillary segments
    per = n_cols // capillary_amount
    segments = [
        range(i*per, (i+1)*per if i<capillary_amount-1 else n_cols)
        for i in range(capillary_amount)
    ]

    figs = []
    capillary_data = []

    # 1) Combined overlay
    fig_all, ax_all = plt.subplots(figsize=(8, 5))
    for idx, cols in enumerate(segments, 1):
        raw = heatmap[:, list(cols)].sum(axis=1)
        smooth = gaussian_filter1d(raw, sigma=smoothing_factor)
        capillary_data.append(smooth)
        ax_all.plot(x_values, smooth, label=f"Capillary {idx}")
    ax_all.set(
        xlabel="Molecular Weight (kDa)",
        ylabel="Aggregated Abundance",
        title="Combined Capillary Traces"
    )
    ax_all.legend(loc='best', fontsize='small', ncol=2)
    figs.append(fig_all)

    # 2) Individual plots
    for idx, smooth in enumerate(capillary_data, 1):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(x_values, smooth)
        ax.set(
            xlabel="Molecular Weight (kDa)",
            ylabel="Aggregated Abundance",
            title=f"Capillary {idx}"
        )
        figs.append(fig)

    return figs, capillary_data, x_values
