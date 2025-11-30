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
    grid_size=(700, 700),
    *,
    sigma_x_factor=0.006,
    sigma_y_factor=0.004,
    bleed_sigma=3.5,
    smear_strength=0.45,
    paper_texture=0.28,
    scan_noise=0.2,
    tone_gamma=0.8,
    haze_strength=0.25,
    random_seed=None,
):
    """Return (heatmap, extent) for a simulated 2D gel.

    This version focuses on a photographic gel feel by layering spot bleed,
    smear, paper texture, scanner noise, and tone‑curve shaping onto the
    analytical spot layout.
    """

    rng = np.random.default_rng(random_seed)

    # 1) Grid + coordinate transforms
    data = list(data_list)
    pi_vals = np.array([d[2] for d in data])
    mw_vals = np.array([d[1] / 1000 for d in data])  # kDa
    pi_min, pi_max = pi_vals.min(), pi_vals.max()
    pad = 0.06 * (pi_max - pi_min)
    pi_min, pi_max = pi_min - pad, pi_max + pad

    def ln_transform(x):
        return 5.3779 * np.log(x) - 6.4014

    y_min, y_max = ln_transform(mw_vals.min() * 0.85), ln_transform(mw_vals.max() * 1.1)

    X = np.linspace(pi_min, pi_max, grid_size[1])
    Y = np.linspace(y_min, y_max, grid_size[0])
    X, Y = np.meshgrid(X, Y)
    heat = np.zeros(grid_size)

    sig_x = sigma_x_factor * (pi_max - pi_min)
    sig_y = sigma_y_factor * (y_max - y_min)

    # 2) Core spots with irregularity
    blot_noise = gaussian_filter(rng.normal(0, 1, size=grid_size), sigma=1.4)
    blot_noise = (blot_noise - blot_noise.min()) / (np.ptp(blot_noise) + 1e-9)
    for pid, mw, pi, abund in data:
        cy = ln_transform(mw / 1000)
        jitter = rng.uniform(0.65, 1.5, size=2)
        theta = rng.uniform(0, np.pi)
        rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        cov = rot @ np.array([[ (sig_x * jitter[0]) ** 2, 0], [0, (sig_y * jitter[1]) ** 2]]) @ rot.T

        Z = multivariate_normal(mean=[pi, cy], cov=cov).pdf(np.dstack((X, Y)))
        # spot‑local blotchiness
        spot_mask = 0.7 + 0.6 * blot_noise
        blob = Z * spot_mask
        blob_sum = blob.sum() + 1e-9
        heat += blob * (abund / blob_sum)

    # 3) Bleed / bloom and smear
    if bleed_sigma > 0:
        bleed = gaussian_filter(heat, sigma=bleed_sigma)
        heat = 0.6 * heat + 0.4 * bleed

    if smear_strength > 0:
        smear_y = gaussian_filter(heat, sigma=(3 * smear_strength + 0.5, 0.4))
        smear_x = gaussian_filter(heat, sigma=(0.6, 3 * smear_strength + 0.5))
        heat = (1 - smear_strength) * heat + 0.35 * smear_strength * smear_y + 0.35 * smear_strength * smear_x

    # 4) Gentle smile warp for gel curve
    h, w = grid_size
    x_norm = np.linspace(-1, 1, w)
    smile = (0.012 * h) * (x_norm ** 2)
    src_y = np.arange(h)[:, None] + smile
    src_x = np.tile(np.arange(w), (h, 1))
    heat = map_coordinates(heat, [src_y, src_x], order=1, mode="nearest")

    # 5) Background stain + paper texture
    grad_x = np.linspace(-1, 1, w)
    grad_y = np.linspace(-1, 1, h)[:, None]
    gradient = 0.25 * grad_y + 0.15 * grad_x
    paper = gaussian_filter(rng.normal(0, 1, size=grid_size), sigma=28)
    paper = (paper - paper.min()) / (np.ptp(paper) + 1e-9) - 0.5
    heat += 0.08 + 0.05 * gradient + paper_texture * paper

    # 6) Scanner noise overlay
    if scan_noise > 0:
        scan = gaussian_filter(rng.normal(0, 1, size=grid_size), sigma=1.2)
        banding = np.sin(np.linspace(0, np.pi * 4, h))[:, None]
        heat += scan_noise * (0.6 * scan + 0.4 * banding)

    # 7) Bloom / haze
    if haze_strength > 0:
        blurred = gaussian_filter(heat, sigma=3.0)
        heat = (1 - haze_strength) * heat + haze_strength * blurred

    # 8) Tone curve + normalization
    heat = np.maximum(heat, 0)
    if heat.max() > 0:
        norm = heat / heat.max()
        norm = norm ** tone_gamma
        heat = norm / (norm.max() + 1e-9)

    extent = (pi_min, pi_max, y_min, y_max)
    return heat, extent


def plot_heatmap(heatmap, extent, *, cmap=None):
    """Render the simulated gel with photographic styling."""

    if cmap is None:
        cmap = create_custom_blue_cmap()

    fig, ax = plt.subplots(figsize=(8, 10), facecolor="#f8f5f1")

    # subtle paper vignette
    h, w = heatmap.shape
    xv = np.linspace(-1, 1, w)
    yv = np.linspace(-1, 1, h)[:, None]
    vignette = 1 - 0.12 * (xv ** 2 + yv ** 2)

    rgba = cmap(np.clip(heatmap * vignette, 0, 1))
    ax.imshow(rgba, origin="lower", extent=extent, aspect="auto")

    ax.set_xlabel("Isoelectric Point (pI)", color="#1f1c2e")
    ax.set_ylabel("Molecular Weight (kDa)", color="#1f1c2e")
    ax.set_title("2D Gel Electrophoresis Simulation", fontsize=13, pad=10, color="#1f1c2e")

    def inv_ln(v):
        return np.exp((v + 6.4014) / 5.3779)

    ticks = [t for t in ax.get_yticks() if extent[2] <= t <= extent[3]]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{inv_ln(t):.0f}" for t in ticks])

    ax.tick_params(colors="#1f1c2e", labelsize=9)
    for spine in ax.spines.values():
        spine.set_visible(False)

    return fig
