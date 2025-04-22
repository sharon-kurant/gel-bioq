import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

def plot_capillary_traces(
        heatmap: np.ndarray,
        extent: tuple,
        capillary_amount: int,
        smoothing_factor: float
    ) -> list:
    """
    Returns a list of Figures: one overlay + one per capillary.
    """
    pi_min, pi_max, lnmw_min, lnmw_max = extent
    n_rows, n_cols = heatmap.shape

    # recover kDa axis
    lnmw = np.linspace(lnmw_min, lnmw_max, n_rows)
    mw_vals = np.exp((lnmw + 6.4014) / 5.3779)

    # split columns
    per = n_cols // capillary_amount
    segments = [
        range(i*per, (i+1)*per if i<capillary_amount-1 else n_cols)
        for i in range(capillary_amount)
    ]

    figs = []

    # combined overlay
    fig_all, ax_all = plt.subplots(figsize=(8,5))
    for idx, cols in enumerate(segments, 1):
        prof = heatmap[:, list(cols)].sum(axis=1)
        prof_s = gaussian_filter1d(prof, sigma=smoothing_factor)
        ax_all.plot(mw_vals, prof_s, label=f"Capillary {idx}")
    ax_all.set(
        xlabel="Molecular Weight (kDa)",
        ylabel="Aggregated Abundance",
        title="Combined Capillary Traces"
    )
    ax_all.legend(loc='best', fontsize='small', ncol=2)
    figs.append(fig_all)

    # individual
    for idx, cols in enumerate(segments, 1):
        fig, ax = plt.subplots(figsize=(6,4))
        prof = heatmap[:, list(cols)].sum(axis=1)
        prof_s = gaussian_filter1d(prof, sigma=smoothing_factor)
        ax.plot(mw_vals, prof_s)
        ax.set(
            xlabel="Molecular Weight (kDa)",
            ylabel="Aggregated Abundance",
            title=f"Capillary {idx}"
        )
        figs.append(fig)

    return figs
