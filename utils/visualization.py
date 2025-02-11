import numpy as np
import matplotlib
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
    augmentation_type=None,
    mw_scale=1.0
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
            pi = [prop[2] for prop in properties]
            mw = [prop[1] for prop in properties]
            log_mw = [logarithmic_transform(mw_val) for mw_val in mw]
            sizes = [abundances.get(prop[0], 1) for prop in properties]
            
            # Create marker that will be scaled into an ellipse
            if mw_scale != 1.0:
                # Create elliptical markers by scaling circle in y direction
                marker = matplotlib.markers.MarkerStyle(
                    marker='o',
                    transform=matplotlib.transforms.Affine2D().scale(1, mw_scale)
                )
            else:
                marker = 'o'
            
            # Store data for export
            if label == organism1:
                plot_data['protein_ids'].extend([prop[0] for prop in properties])
                plot_data['pI_values'].extend(pi)
                plot_data['mw_values'].extend([mw_val/1000 for mw_val in mw])
                plot_data['abundance1'].extend(sizes)
                plot_data['abundance2'].extend([0] * len(properties))
            else:
                plot_data['protein_ids'].extend([prop[0] for prop in properties])
                plot_data['pI_values'].extend(pi)
                plot_data['mw_values'].extend([mw_val/1000 for mw_val in mw])
                plot_data['abundance1'].extend([0] * len(properties))
                plot_data['abundance2'].extend(sizes)
            
            plt.scatter(pi, log_mw, alpha=0.5, color=color, s=sizes, 
                      marker=marker, label=f'{label} ({ratio}%)')
    
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
    show_sum,
    cap_start,
    cap_end,
    mw_scale=1.0,
    enable_spillage=False,
    spillage_width=0.1
):
    """Create capillary plot with spillage handling."""
    fig, ax = plt.subplots(figsize=(6, 4))
    x_values = np.linspace(0, 20, 2000)
    y1 = np.zeros_like(x_values)
    y2 = np.zeros_like(x_values)
    
    def calculate_spillage_ratio(pI, std_dev, cap_start, cap_end):
        """
        Calculate what portion of the protein belongs in this capillary.
        Using error function (erf) to calculate the area under the normal distribution.
        """
        from scipy.special import erf
        
        # Calculate the proportion of the protein that falls within this capillary
        proportion = 0.5 * (erf((cap_end - pI)/(std_dev * np.sqrt(2))) - 
                          erf((cap_start - pI)/(std_dev * np.sqrt(2))))
        return max(0, min(1, proportion))  # Ensure value is between 0 and 1
    
    # Calculate distributions
    for props, abundances, y_values in [
        (filtered_props1, normalized_abundance1, y1),
        (filtered_props2, normalized_abundance2, y2)
    ]:
        for prop in props:
            protein_id, mw, pI_tuple = prop
            original_pI, current_pI = pI_tuple if isinstance(pI_tuple, tuple) else (pI_tuple, pI_tuple)
            
            if enable_spillage:
                # Calculate spillage contribution
                spillage_ratio = calculate_spillage_ratio(current_pI, spillage_width, cap_start, cap_end)
                # Only proceed if there's any contribution to this capillary
                if spillage_ratio > 0:
                    abundance = abundances.get(protein_id, 0)
                    if abundance > 0 and mw > 0:
                        base_std = max(abundance/100, 0.01)
                        scaled_std = base_std * mw_scale if mw_scale != 1.0 else base_std
                        gaussian = norm.pdf(x_values, loc=mw/1000, scale=scaled_std)
                        y_values += gaussian * abundance * spillage_ratio
            else:
                # Original behavior without spillage
                if cap_start <= original_pI < cap_end:
                    abundance = abundances.get(protein_id, 0)
                    if abundance > 0 and mw > 0:
                        base_std = max(abundance/100, 0.01)
                        scaled_std = base_std * mw_scale if mw_scale != 1.0 else base_std
                        gaussian = norm.pdf(x_values, loc=mw/1000, scale=scaled_std)
                        y_values += gaussian * abundance
    
    # Apply smoothing - adjust based on MW scale
    adjusted_sigma = smoothing_sigma * mw_scale
    y1_smooth = gaussian_filter1d(y1, sigma=adjusted_sigma)
    y2_smooth = gaussian_filter1d(y2, sigma=adjusted_sigma)
    y_sum = gaussian_filter1d(y1 + y2, sigma=adjusted_sigma)
    
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
    
    plt.title(f'Capillary: pI range ({cap_start:.2f} - {cap_end:.2f})')
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