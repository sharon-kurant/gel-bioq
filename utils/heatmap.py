# utils/heatmap.py

import numpy as np
from scipy.stats import multivariate_normal
from scipy.ndimage import gaussian_filter, map_coordinates
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

def create_custom_blue_cmap():
    """
    Build and return a LinearSegmentedColormap mimicking Coomassie‑blue gels.
    """
    # RGBA for background, start, end
    background = np.array([187/255, 222/255, 251/255, 0.3])
    start_color = np.array([10/255, 100/255, 170/255, 1.0])
    end_color   = np.array([0/255,  41/255, 107/255, 1.0])

    # non‑linear positions for sensitivity at low values
    alpha = 0.3
    positions = (np.linspace(0, 1, 256) ** alpha)
    positions /= positions.max()

    # allocate arrays
    r = np.zeros(256); g = np.zeros(256)
    b = np.zeros(256); a = np.zeros(256)

    # first entry = background
    r[0], g[0], b[0], a[0] = background

    # interpolate rest
    for i in range(1, 256):
        pos = (i - 1) / 254
        r[i] = start_color[0] * (1-pos) + end_color[0] * pos
        g[i] = start_color[1] * (1-pos) + end_color[1] * pos
        b[i] = start_color[2] * (1-pos) + end_color[2] * pos
        a[i] = start_color[3] * (1-pos) + end_color[3] * pos

    colors = np.vstack((r, g, b, a)).T
    return LinearSegmentedColormap.from_list("gel_blues", colors)


def generate_combined_heatmap(
        data_list,
        grid_size=(500, 500),
        *,
        # σ‑multipliers
        sigma_x_factor=0.005,
        sigma_y_factor=0.005,
        # spot morphology
        spot_irregularity=True,
        spot_spread_jitter=(0.75, 1.35),
        spot_angle_jitter=True,
        # streaks
        apply_streaks=True,
        streak_orient='both',
        streak_prob=0.12,
        # spot trains
        apply_spot_trains=True,
        train_N=(10, 50),
        train_decay=0.55,
        train_offset=0.9,
        # smile distortion
        apply_smile=True,
        smile_rel_amp=0.02,
        smile_curve='S',
        smile_pow=2.0,
        smile_s_coef=0.3,
        # edge fading
        apply_edges=True,
        edge_prob=0.6,
        edge_strength=(1.35, 1.0),
        # stain & haze
        background_intensity=0.02,
        texture_strength=0.18,
        texture_scale=0.85,
        haze_strength=0.35,
        haze_sigma=6.0,
        dynamic_range_gamma=0.65,
        # dropout & abundance noise
        apply_dropout=False,
        drop_frac_range=(0.0, 0.5),
        apply_abundance_variation=False,
        abundance_var_range=(0.5, 2.0),
        abundance_var_sd=0.3,
        # reproducibility
        random_seed=None
    ):
    """
    Return (heatmap, extent) for a simulated 2D‑gel with optional artefacts.
    """
    rng = np.random.default_rng(random_seed)

    # 1) Dropout
    data = list(data_list)
    if apply_dropout:
        frac = rng.uniform(*drop_frac_range)
        keep_n = int(len(data) * (1 - frac))
        idxs = rng.choice(len(data), keep_n, replace=False)
        data = [data[i] for i in idxs]

    # 2) Abundance variation
    if apply_abundance_variation:
        varied = []
@@ -96,60 +107,82 @@ def generate_combined_heatmap(
                factor = rng.normal(1.0, abundance_var_sd)
            varied.append((pid, mw, pi, abund * factor))
        data = varied

    # 3) Grid setup
    pi_vals = np.array([d[2] for d in data])
    mw_vals = np.array([d[1] / 1000 for d in data])  # kDa
    pi_min, pi_max = pi_vals.min(), pi_vals.max()
    pad = 0.05 * (pi_max - pi_min)
    pi_min, pi_max = pi_min - pad, pi_max + pad

    def ln_transform(x): return 5.3779 * np.log(x) - 6.4014
    y_min, y_max = ln_transform(mw_vals.min() * 0.9), ln_transform(mw_vals.max() * 1.1)

    X = np.linspace(pi_min, pi_max, grid_size[1])
    Y = np.linspace(y_min, y_max, grid_size[0])
    X, Y = np.meshgrid(X, Y)
    heat = np.zeros(grid_size)

    sig_x = sigma_x_factor * (pi_max - pi_min)
    sig_y = sigma_y_factor * (y_max - y_min)

    # 4) Spots + streaks
    for pid, mw, pi, abund in data:
        cy = ln_transform(mw / 1000)

        # spot‑specific width & angle
        spot_sig_x = sig_x
        spot_sig_y = sig_y
        if spot_irregularity:
            jitter = rng.uniform(*spot_spread_jitter, size=2)
            spot_sig_x *= jitter[0]
            spot_sig_y *= jitter[1]

        stretch_x = 1.0
        stretch_y = 1.0
        if apply_streaks and rng.random() < streak_prob:
            orient = streak_orient
            if orient == 'both':
                orient = 'horizontal' if rng.random() < 0.5 else 'vertical'
            if orient == 'horizontal':
                stretch_x = 6.0
            else:
                stretch_y = 6.0

        cov = np.array([
            [(spot_sig_x * stretch_x) ** 2, 0],
            [0, (spot_sig_y * stretch_y) ** 2]
        ])

        if spot_irregularity and spot_angle_jitter:
            theta = rng.uniform(0, np.pi)
            rot = np.array([
                [np.cos(theta), -np.sin(theta)],
                [np.sin(theta),  np.cos(theta)]
            ])
            cov = rot @ cov @ rot.T

        Z = multivariate_normal(mean=[pi, cy], cov=cov).pdf(np.dstack((X, Y)))
        heat += Z * (abund / Z.sum())

    # 5) Spot trains
    if apply_spot_trains:
        topN = train_N if isinstance(train_N, int) else rng.integers(train_N[0], train_N[1]+1)
        idxs = np.argsort([-d[3] for d in data])[:topN]
        for i in idxs:
            pid, mw, pi, abund = data[i]
            repeats = rng.integers(1, 8)
            for k in range(1, repeats+1):
                shift = k * sig_x * train_offset
                new_pi = pi + shift
                sat   = abund * (train_decay**k)
                Z = multivariate_normal(
                        mean=[new_pi, ln_transform(mw/1000)],
                        cov=[[sig_x**2,0],[0,sig_y**2]]
                    ).pdf(np.dstack((X, Y)))
                heat += Z * (sat / Z.sum())

    # 6) Gel distortion
    h, w = grid_size
    # 6a) smile
    if apply_smile:
@@ -161,58 +194,96 @@ def generate_combined_heatmap(
            smile = amp * (np.abs(x_norm) ** smile_pow) * (1 - smile_s_coef * x_norm)

        src_y = np.arange(h)[:, None] + smile            # shape (h, w)
        src_x = np.tile(np.arange(w), (h, 1))            # shape (h, w)
        heat = map_coordinates(heat, [src_y, src_x], order=1, mode='nearest')


    # 6b) edge fading
    if apply_edges:
        start_f, end_f = edge_strength
        for edge in ['left','right','top','bottom']:
            if rng.random() < edge_prob:
                if edge in ('left','right'):
                    width = int(0.08 * w)
                    for i in range(width):
                        f   = start_f - (start_f - end_f)*i/width
                        col = i if edge=='left' else -i-1
                        heat[:, col] *= f
                else:
                    height = int(0.08 * h)
                    for j in range(height):
                        f   = start_f - (start_f - end_f)*j/height
                        row = j if edge=='top' else -j-1
                        heat[row, :] *= f

    # 7) Background staining & texture
    if background_intensity > 0 or texture_strength > 0:
        grad_x = np.linspace(-1, 1, w)
        grad_y = np.linspace(-1, 1, h)[:, None]
        gradient = 0.35 * grad_x + 0.35 * grad_y

        mottled = gaussian_filter(
            rng.normal(0, 1, size=grid_size),
            sigma=max(1, texture_scale * min(grid_size) / 40)
        )
        mottled = (mottled - mottled.min()) / (np.ptp(mottled) + 1e-9) - 0.5

        heat += background_intensity * (1 + gradient)
        heat += texture_strength * mottled

    # 8) Haze/halos for stronger blotting
    if haze_strength > 0:
        blurred = gaussian_filter(heat, sigma=haze_sigma)
        heat = (1 - haze_strength) * heat + haze_strength * blurred

    # 9) Dynamic range shaping (γ < 1 brightens low intensities)
    heat = np.maximum(heat, 0)
    if dynamic_range_gamma != 1.0 and heat.max() > 0:
        norm = heat / heat.max()
        heat = (norm ** dynamic_range_gamma) * heat.max()

    extent = (pi_min, pi_max, y_min, y_max)
    return heat, extent


def plot_heatmap(heatmap, extent, *, cmap=None, gel_like=True):
    """
    Display the heatmap (pI × MW) and return the Matplotlib Figure.
    """
    if cmap is None:
        cmap = create_custom_blue_cmap()

    fig_face = '#f7f4f7' if gel_like else 'white'
    fig, ax = plt.subplots(figsize=(10, 10), facecolor=fig_face)
    im = ax.imshow(
        heatmap,
        origin='lower',
        extent=extent,
        aspect='auto',
        cmap=cmap
    )

    ax.set_xlabel("Isoelectric Point (pI)")
    ax.set_ylabel("Molecular Weight (kDa)")
    ax.set_title("2D Gel Electrophoresis Simulation", fontsize=13, pad=10)

    def inv_ln(v): return np.exp((v + 6.4014)/5.3779)
    ticks = [t for t in ax.get_yticks() if extent[2] <= t <= extent[3]]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{inv_ln(t):.0f}" for t in ticks])

    if gel_like:
        ax.grid(False)
        ax.set_facecolor('#ece9f1')
        ax.tick_params(colors='#3f3d56', labelsize=9)
        for spine in ax.spines.values():
            spine.set_alpha(0.25)
        cbar = fig.colorbar(im, ax=ax, label="Protein Abundance")
        cbar.outline.set_visible(False)
        cbar.ax.tick_params(labelsize=9)
    else:
        ax.grid(True, which='major', linestyle='-', linewidth=0.5, alpha=0.4)
        ax.set_facecolor('#bbdefb')
        fig.colorbar(im, ax=ax, label="Protein Abundance")

    return fig
