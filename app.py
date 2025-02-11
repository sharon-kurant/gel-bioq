import streamlit as st
from datetime import datetime

# Import utility functions
from utils.data_processing import *
from utils.visualization import create_gel_plot, create_capillary_plot
from utils.export import create_download_package

# Page configuration
st.set_page_config(
    page_title="2D Gel Electrophoresis Analysis",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Get available organisms
ORGANISMS = scan_available_organisms()

if not ORGANISMS:
    st.error("No organism data found. Please check the data/abundance directory.")
    st.stop()

# Title and description
st.title("2D Gel Electrophoresis Analysis")
st.markdown("Select organisms and adjust parameters to visualize the results.")

# Sidebar controls
with st.sidebar:
    st.header("Settings")
    
    # Organism selection
    st.subheader("Organism Selection")
    available_organisms = list(ORGANISMS.keys())
    
    if len(available_organisms) < 2:
        st.error("At least two organisms are required for comparison.")
        st.stop()
    
    organism1 = st.selectbox(
        "First Organism",
        options=available_organisms,
        index=0,
        key="organism1"
    )
    
    # Filter out the first selected organism from the second dropdown
    remaining_organisms = [org for org in available_organisms if org != organism1]
    organism2 = st.selectbox(
        "Second Organism",
        options=remaining_organisms,
        index=0,
        key="organism2"
    )
    
    # Analysis parameters
    st.subheader("Analysis Parameters")
    min_mw = st.number_input("Minimum Molecular Weight", value=0.000001, format="%.6f", key="min_mw")
    max_mw = st.number_input("Maximum Molecular Weight", value=500000.0, key="max_mw")
    min_normalized_abundance = st.number_input("Minimum Normalized Abundance", value=2.0, key="min_norm_abundance")
    
    # Gaussian parameter
    gaussian_std = st.slider("Gaussian Standard Deviation", 0.01, 1.0, 0.01, 0.01, key="gaussian_std")
    
    # Sample ratios
    st.subheader("Sample Ratios")
    ratio1 = st.slider(f"{organism1} Ratio (%)", 0, 100, 95, key="ratio1")
    ratio2 = st.slider(f"{organism2} Ratio (%)", 0, 100, 5, key="ratio2")
    
    # Augmentation settings
    st.subheader("Augmentation Settings")
    augmentation_type = st.radio(
        "Select Augmentation Type",
        options=["None", "pI Shift", "MW Scale"],
        key="augmentation_type"
    )
    
    if augmentation_type == "pI Shift":
        pI_shift = st.slider(
            "pI Shift Amount",
            min_value=-1.0,
            max_value=1.0,
            value=0.0,
            step=0.1,
            help="Shift pI values left (negative) or right (positive)",
            key="pi_shift"
        )
    elif augmentation_type == "MW Scale":
        mw_scale = st.slider(
            "MW Scale Factor",
            min_value=0.5,
            max_value=2.0,
            value=1.0,
            step=0.1,
            help="Scale molecular weights (0.5 = squeeze, 2.0 = stretch)",
            key="mw_scale"
        )
    
    st.markdown("---")
    
    # Noise Settings
    st.subheader("Noise Settings")
    enable_mw_noise = st.checkbox(
        "Add Random Noise to Molecular Weights (±50%)", 
        value=False,
        help="Randomly increase or decrease each protein's molecular weight by up to 50%",
        key="enable_noise"
    )
    
    
    # Capillary settings
    st.markdown("---")
    st.subheader("Capillary Analysis")    
    num_capillaries = st.slider("Number of Capillaries", 1, 100, 8, key="num_capillaries")
    
    # Add spillage controls
    enable_spillage = st.sidebar.checkbox(
        "Enable Capillary Spillage", 
        value=False, 
        help="Allow proteins to contribute to adjacent capillaries based on their circle size"
    )
    
    # Visualization toggles
    st.subheader("Show/Hide Plots")
    show_organism1 = st.sidebar.checkbox(f"Show {organism1}", value=True, key="show_org1")
    show_organism2 = st.sidebar.checkbox(f"Show {organism2}", value=True, key="show_org2")
    show_sum = st.sidebar.checkbox("Show Sum", value=True, key="show_sum")
    
    # Smoothing control
    st.subheader("Visualization Settings")
    smoothing_sigma = st.slider("Smoothing Factor", 1, 20, 5, key="smoothing")
    
    st.markdown("---")
    
    # Download section
    st.subheader("Download Results")
    download_button = st.button("Download All Results", key="download")
    
# Main app logic
try:
    with st.spinner("Processing data..."):
        # Process data for both organisms
        properties1, abundance1 = process_organism_data(
            ORGANISMS[organism1], min_mw, max_mw
        )
        properties2, abundance2 = process_organism_data(
            ORGANISMS[organism2], min_mw, max_mw
        )
        
        # Apply noise if enabled
        if enable_mw_noise:
            properties1 = add_mw_noise(properties1)
            properties2 = add_mw_noise(properties2)
        
        # Store original properties for capillary calculations
        original_properties1 = properties1.copy()
        original_properties2 = properties2.copy()
        
        # Apply pI shift if selected
        if augmentation_type == "pI Shift":
            properties1, actual_shift = shift_pI(properties1, pI_shift)
            properties2, _ = shift_pI(properties2, pI_shift)
            st.info(f"Applied pI shift of {actual_shift:.2f} units")
    
        
        # Normalize abundances
        normalized_abundance1, normalized_abundance2 = normalize_abundance(
            abundance1, abundance2, ratio1, ratio2
        )
        
        # Create visualization tabs
        tab1, tab2 = st.tabs(["2D Gel Plot", "Capillary Analysis"])
        
        with tab1:
            st.subheader("2D Gel Plot")
            gel_fig, gel_data = create_gel_plot(
                properties1, properties2,
                normalized_abundance1, normalized_abundance2,
                organism1, organism2,
                ratio1, ratio2,
                augmentation_type,
                mw_scale=mw_scale if augmentation_type == "MW Scale" else 1.0
            )
            st.pyplot(gel_fig)
        
        with tab2:
            st.subheader("Capillary Analysis")
            capillary_figs = []
            capillary_data = []
            x_values = None
            
            # Calculate capillary ranges
            min_pI, max_pI = get_pI_range(original_properties1, original_properties2)
            capillary_ranges = calculate_capillary_ranges(min_pI, max_pI, num_capillaries)
            
            # Create capillary plots
            cols = st.columns(2)
            for i, (cap_start, cap_end) in enumerate(capillary_ranges):
                with cols[i % 2]:
                    st.write(f"Capillary {i + 1}: pI {cap_start:.2f} - {cap_end:.2f}")
                    
                    # Filter proteins based on their original positions
                    filtered_props1 = filter_by_pI_range(properties1, cap_start, cap_end)
                    filtered_props2 = filter_by_pI_range(properties2, cap_start, cap_end)
                    
                    # create_capillary_plot
                    fig, data, x_vals = create_capillary_plot(
                        filtered_props1, filtered_props2,
                        normalized_abundance1, normalized_abundance2,
                        organism1, organism2,
                        smoothing_sigma,
                        gaussian_std,
                        show_organism1, show_organism2, show_sum,
                        cap_start, cap_end,
                        mw_scale=mw_scale if augmentation_type == "MW Scale" else 1.0,
                        enable_spillage=enable_spillage  # Remove spillage_width parameter
                    )
                    
                    capillary_figs.append(fig)
                    capillary_data.append(data)
                    if x_values is None:
                        x_values = x_vals
                    
                    st.pyplot(fig)
                    
        # Handle download if requested
        if download_button:
            # Update the parameters dictionary in the download section
            parameters = {
                'Analysis_Parameters': {
                    'min_molecular_weight': min_mw,
                    'max_molecular_weight': max_mw,
                    'min_normalized_abundance': min_normalized_abundance,
                    'gaussian_std': gaussian_std
                },
                'Sample_Ratios': {
                    organism1: ratio1,
                    organism2: ratio2
                },
                'Augmentations': {
                    'type': augmentation_type,
                    'pI_shift': pI_shift if augmentation_type == "pI Shift" else 0,
                    'mw_scale': mw_scale if augmentation_type == "MW Scale" else 1.0,
                    'mw_noise_enabled': enable_mw_noise
                },
                'Capillary_Settings': {
                    'number_of_capillaries': num_capillaries,
                    'spillage_enabled': enable_spillage,
                    'ranges': [
                        {
                            f'capillary_{i+1}': {
                                'pI_range': {
                                    'start': ranges[0],
                                    'end': ranges[1]
                                }
                            }
                        } for i, ranges in enumerate(capillary_ranges)
                    ]
                },
                'Visualization_Settings': {
                    'smoothing_factor': smoothing_sigma,
                    'show_organism1': show_organism1,
                    'show_organism2': show_organism2,
                    'show_sum': show_sum
                }
            }
            
            zip_buffer = create_download_package(
                parameters=parameters,
                gel_plot_fig=gel_fig,
                gel_plot_data=gel_data,
                capillary_figs=capillary_figs,
                capillary_data=capillary_data,
                x_values=x_values,
                organism1=organism1,
                organism2=organism2
            )
            
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            st.sidebar.download_button(
                label="Click to Download",
                data=zip_buffer,
                file_name=f'gel_analysis_{timestamp}.zip',
                mime="application/zip"
            )


except Exception as e:
    st.error(f"An error occurred: {str(e)}")
    st.write("Please check that all data files are in the correct locations:")
    st.code("""
    gel/
    ├── data/
    │   ├── abundance/
    │   │   ├── 4932-WHOLE_ORGANISM-integrated.txt
    │   │   └── 511145-WHOLE_ORGANISM-integrated.txt
    │   └── fasta/
    │       ├── fasta.v11.5.4932.fa
    │       └── fasta.v11.5.511145.fa
    """)

# Instructions at the bottom
st.markdown("""
---
### How to Use
1. Select organisms from the dropdown menus in the sidebar

2. Adjust Analysis Parameters:
   * Set molecular weight range (kDa)
   * Adjust minimum normalized abundance
   * Configure Gaussian standard deviation for peak shapes

3. Set Sample Ratios:
   * Adjust relative proportions of each organism (0-100%)

4. Choose Augmentations (optional):
   * pI Shift: Shift proteins left/right on pI axis (-1 to +1 units)
   * MW Scale: Stretch/squeeze protein spots vertically (0.5x to 2x)
   * MW Noise: Add random variation to molecular weights (±50%)

5. Configure Capillary Analysis:
   * Set number of capillaries (1-100)
   * Enable/disable capillary spillage
   * Toggle visibility of individual organisms and sum
   * Adjust smoothing factor for line plots

6. View Visualizations:
   * 2D Gel Plot tab: Shows complete protein distribution
   * Capillary Analysis tab: Shows protein distribution in pI ranges
   * Use toggles to show/hide individual organisms and sum

7. Export Results:
   * Click "Download Results" to get:
     - All visualization parameters
     - 2D gel plot image and data
     - Capillary analysis plots and data
     - Complete parameter settings

Note: Augmentations and settings can be combined to explore different aspects of the protein distribution patterns.
""")