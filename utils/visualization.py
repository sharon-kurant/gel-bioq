import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d

def create_gel_plot(
    protein_properties1,
    protein_properties2,
    normalized_abundance1,
    normalized_abundance2,
    organism1,
    organism2,
    ratio1,
    ratio2,
    augmentation_type=None
):
    """Create 2D gel plot."""
    fig, ax = plt.subplots(figsize=(12, 8))
    
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
            if label == organism1:
                plot_data['protein_ids'].extend([prop[0] for prop in properties])
                plot_data['pI_values'].extend(pi)
                plot_data['mw_values'].extend([mw_val/1000 for mw_val in mw])  # Convert to kDa
                plot_data['abundance1'].extend(sizes)
                plot_data['abundance2'].extend([0] * len(properties))  # Fill with zeros for alignment
            else:
                plot_data['protein_ids'].extend([prop[0] for prop in properties])
                plot_data['pI_values'].extend(pi)
                plot_data['mw_values'].extend([mw_val/1000 for mw_val in mw])  # Convert to kDa
                plot_data['abundance1'].extend([0] * len(properties))  # Fill with zeros for alignment
                plot_data['abundance2'].extend(sizes)
            
            # Store data for export
            plot_data['protein_ids'].extend([prop[0] for prop in properties])
            plot_data['pI_values'].extend(pi)
            plot_data['mw_values'].extend([mw_val/1000 for mw_val in mw])  # Convert to kDa for export
            plot_data['abundance1' if label == organism1 else 'abundance2'].extend(sizes)
            
            plt.scatter(pi, log_mw, alpha=0.5, color=color, s=sizes, 
                      label=f'{label} ({ratio}%)')
    
    plot_organism(protein_properties1, normalized_abundance1, 'blue', organism1, ratio1)
    plot_organism(protein_properties2, normalized_abundance2, 'red', organism2, ratio2)
    
    plt.title('2-D Gel Plot')
    plt.xlabel('Isoelectric Point (pI)')
    plt.ylabel('Molecular Weight (kDa)')
    plt.legend()
    plt.grid(True, linestyle='--', linewidth=0.5)
    
    # Set y-axis labels
    yticks = plt.gca().get_yticks()
    ytick_labels = [f'{int(np.exp((tick + 6.4014) / 5.3779) / 1000):.0f}' for tick in yticks]
    plt.gca().set_yticklabels(ytick_labels)
    
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
    """Create capillary plot."""
    fig, ax = plt.subplots(figsize=(6, 4))
    x_values = np.linspace(0, 20, 2000)
    y1 = np.zeros_like(x_values)
    y2 = np.zeros_like(x_values)
    
    # Calculate distributions
    for props, abundances, y_values in [
        (filtered_props1, normalized_abundance1, y1),
        (filtered_props2, normalized_abundance2, y2)
    ]:
        for prop in props:
            protein_id, mw, pi = prop
            abundance = abundances.get(protein_id, 0)
            if abundance > 0 and mw > 0:
                gaussian = norm.pdf(x_values, loc=mw/1000, 
                                 scale=max(abundance/100, 0.01))  # Revert to original scale calculation
                y_values += gaussian * abundance
    
    # Apply smoothing
    y1_smooth = gaussian_filter1d(y1, sigma=smoothing_sigma)
    y2_smooth = gaussian_filter1d(y2, sigma=smoothing_sigma)
    y_sum = gaussian_filter1d(y1 + y2, sigma=smoothing_sigma)
    
    # Plot lines based on visibility settings
    if show_organism1:
        plt.plot(x_values, y1_smooth, color='blue', 
                label=organism1, linewidth=1.5)
    if show_organism2:
        plt.plot(x_values, y2_smooth, color='red', 
                label=organism2, linewidth=1.5)
    if show_sum:
        plt.plot(x_values, y_sum, color='green', 
                label='Sum', linewidth=2, alpha=0.7)
    
    plt.xlabel('Molecular Weight (kDa)')
    plt.ylabel('Volume')
    if show_organism1 or show_organism2 or show_sum:
        plt.legend()
    plt.grid(True, linestyle='--', linewidth=0.5)
    
    # Set y-axis limits
    max_y = max(
        max(y1_smooth) if show_organism1 else 0,
        max(y2_smooth) if show_organism2 else 0,
        max(y_sum) if show_sum else 0
    )
    plt.ylim(0, max_y * 1.1)
    
    return fig, {'y1': y1_smooth, 'y2': y2_smooth, 'y_sum': y_sum}, x_values

def logarithmic_transform(mw):
    """Transform molecular weight for visualization."""
    return 5.3779 * np.log(mw) - 6.4014