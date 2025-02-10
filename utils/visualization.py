import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d

def create_gel_plot(
    properties1,
    properties2,
    normalized_abundance1,
    normalized_abundance2,
    organism1,
    organism2,
    ratio1,
    ratio2,
    augmentation_type=None
):
    """Create 2D gel plot with optional augmentations."""
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111)
    
    # Store data for export
    plot_data = {
        'protein_ids': [],
        'pI_values': [],
        'mw_values': [],
        'abundance1': [],
        'abundance2': []
    }
    
    def plot_organism(properties, abundances, color, label, ratio):
        if properties:
            pi, mw = zip(*[(prop[2], prop[1]) for prop in properties])
            log_mw = [logarithmic_transform(mw_val) for mw_val in mw]
            sizes = [abundances.get(prop[0], 1) for prop in properties]
            
            # Store data for export
            plot_data['protein_ids'].extend([prop[0] for prop in properties])
            plot_data['pI_values'].extend(pi)
            plot_data['mw_values'].extend(mw)
            plot_data['abundance1' if label == organism1 else 'abundance2'].extend(sizes)
            
            ax.scatter(pi, log_mw, alpha=0.5, color=color, s=sizes, 
                      label=f'{label} ({ratio}%)')
    
    plot_organism(properties1, normalized_abundance1, 'blue', organism1, ratio1)
    plot_organism(properties2, normalized_abundance2, 'red', organism2, ratio2)
    
    ax.set_title('2-D Gel Plot')
    ax.set_xlabel('Isoelectric Point (pI)')
    ax.set_ylabel('Molecular Weight (kDa)')
    ax.legend()
    ax.grid(True, linestyle='--', linewidth=0.5)
    
    # Set y-axis labels
    yticks = ax.get_yticks()
    ytick_labels = [f'{int(np.exp((tick + 6.4014) / 5.3779)):.0f}' for tick in yticks]
    ax.set_yticklabels(ytick_labels)
    
    plt.tight_layout()
    return fig, plot_data

def create_capillary_plot(
    filtered_props1,
    filtered_props2,
    normalized_abundance1,
    normalized_abundance2,
    organism1,
    organism2,
    smoothing_sigma,
    gaussian_std,
    show_organism1,
    show_organism2,
    show_sum
):
    """Create capillary plot with updated parameters."""
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    
    x_values = np.linspace(0, 20, 2000)
    y1 = np.zeros_like(x_values)
    y2 = np.zeros_like(x_values)
    
    # Calculate distributions with custom gaussian std
    for props, abundances, y_values in [
        (filtered_props1, normalized_abundance1, y1),
        (filtered_props2, normalized_abundance2, y2)
    ]:
        for prop in props:
            protein_id, mw, pi = prop
            abundance = abundances.get(protein_id, 0)
            if abundance > 0 and mw > 0:
                gaussian = norm.pdf(x_values, loc=mw/1000, 
                                 scale=gaussian_std)
                y_values += gaussian * abundance
    
    # Apply smoothing
    y1_smooth = gaussian_filter1d(y1, sigma=smoothing_sigma)
    y2_smooth = gaussian_filter1d(y2, sigma=smoothing_sigma)
    y_sum = gaussian_filter1d(y1 + y2, sigma=smoothing_sigma)
    
    if show_organism1:
        ax.plot(x_values, y1_smooth, color='blue', 
                label=organism1, linewidth=1.5)
    if show_organism2:
        ax.plot(x_values, y2_smooth, color='red', 
                label=organism2, linewidth=1.5)
    if show_sum:
        ax.plot(x_values, y_sum, color='green', 
                label='Sum', linewidth=2, alpha=0.7)
    
    ax.set_title('Capillary View')
    ax.set_xlabel('Molecular Weight (kDa)')
    ax.set_ylabel('Volume')
    if show_organism1 or show_organism2 or show_sum:
        ax.legend()
    ax.grid(True, linestyle='--', linewidth=0.5)
    
    # Set reasonable y-axis limits
    max_y = max(
        max(y1_smooth) if show_organism1 else 0,
        max(y2_smooth) if show_organism2 else 0,
        max(y_sum) if show_sum else 0
    )
    ax.set_ylim(0, max_y * 1.1)
    
    plt.tight_layout()
    return fig, {'y1': y1_smooth, 'y2': y2_smooth, 'y_sum': y_sum}, x_values

def logarithmic_transform(mw):
    """Transform molecular weight for visualization."""
    return 5.3779 * np.log(mw) - 6.4014