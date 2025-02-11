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
    min_mw = st.number_input("Minimum Molecular Weight", value=0.000001, format="%.6f")
    max_mw = st.number_input("Maximum Molecular Weight", value=500000.0)
    min_normalized_abundance = st.number_input("Minimum Normalized Abundance", value=2.0)
    
    # Gaussian parameter
    gaussian_std = st.slider("Gaussian Standard Deviation", 0.01, 1.0, 0.01, 0.01)
    
    # Sample ratios
    st.subheader("Sample Ratios")
    ratio1 = st.slider(f"{organism1} Ratio (%)", 0, 100, 95)
    ratio2 = st.slider(f"{organism2} Ratio (%)", 0, 100, 5)
    
    # Capillary settings
    st.subheader("Capillary Analysis")
    num_capillaries = st.slider("Number of Capillaries", 1, 100, 8)
    
    # Augmentation settings
    st.subheader("Augmentation Settings")
    augmentation_type = st.radio(
        "Select Augmentation Type",
        options=["None", "pI Shift", "MW Stretching"]
    )
    
    if augmentation_type == "pI Shift":
        pI_shift = st.slider(
            "pI Shift Amount",
            min_value=-1.0,
            max_value=1.0,
            value=0.0,
            step=0.1,
            help="Shift pI values left (negative) or right (positive)"
        )
        if abs(pI_shift) > 0:
            st.info("Points will be shifted while ensuring all remain visible in the plot")
    elif augmentation_type == "MW Stretching":
        st.info("Molecular weights will be stretched to fill the visualization space")
        
    st.sidebar.markdown("---")
    st.sidebar.subheader("Noise Settings")
    enable_mw_noise  = st.sidebar.checkbox("Add Random Noise to Molecular Weights (±50%)", 
                                   value=False,
                                   help="Randomly increase or decrease each protein's molecular weight by up to 50%")
    
    # Visualization settings
    st.subheader("Visualization Settings")
    smoothing_sigma = st.slider("Smoothing Factor", 1, 20, 5)
    
    # Line visibility controls
    st.subheader("Capillary View Settings")
    show_organism1 = st.checkbox(f"Show {organism1}", value=True)
    show_organism2 = st.checkbox(f"Show {organism2}", value=True)
    show_sum = st.checkbox("Show Sum", value=True)
    
    # Download section
    st.markdown("---")
    st.subheader("Download Results")
    download_button = st.button("Download All Results")

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
        
        # Store original properties for capillary ranges
        original_properties1 = properties1
        original_properties2 = properties2
        
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
                augmentation_type
            )
            st.pyplot(gel_fig)
        
        with tab2:
            with tab2:
                st.subheader("Capillary Analysis")
                capillary_figs = []
                capillary_data = []
                x_values = None
                
                # Calculate capillary ranges
                min_pI, max_pI = get_pI_range(properties1, properties2)
                capillary_ranges = calculate_capillary_ranges(min_pI, max_pI, num_capillaries)
                
                # Create capillary plots
                cols = st.columns(2)
                for i, (cap_start, cap_end) in enumerate(capillary_ranges):
                    with cols[i % 2]:
                        st.write(f"Capillary {i + 1}: pI {cap_start:.2f} - {cap_end:.2f}")
                        
                        # Filter proteins for this capillary
                        filtered_props1 = filter_by_pI_range(properties1, cap_start, cap_end)
                        filtered_props2 = filter_by_pI_range(properties2, cap_start, cap_end)
                        
                        fig, data, x_vals = create_capillary_plot(
                            filtered_props1, filtered_props2,
                            normalized_abundance1, normalized_abundance2,
                            organism1, organism2,
                            smoothing_sigma,
                            gaussian_std,
                            show_organism1, show_organism2, show_sum,
                            cap_start, cap_end  # Add these parameters
                        )
                        
                        capillary_figs.append(fig)
                        capillary_data.append(data)
                        if x_values is None:
                            x_values = x_vals
                        
                        st.pyplot(fig)
                    
        # Handle download if requested
        if download_button:
            parameters = {
                'Analysis_Parameters': {
                    'min_molecular_weight': min_mw,
                    'max_molecular_weight': max_mw,
                    'min_normalized_abundance': min_normalized_abundance,
                    'gaussian_std': gaussian_std,
                    'mw_noise_enabled': enable_mw_noise,
                    'pI_shift': pI_shift if augmentation_type == "pI Shift" else 0
                },
                'Sample_Ratios': {
                    organism1: ratio1,
                    organism2: ratio2
                },
                'Capillaries': {
                    f'capillary_{i+1}': {
                        'pI_range': {
                            'start': ranges[0],
                            'end': ranges[1]
                        }
                    } for i, ranges in enumerate(capillary_ranges)
                },
                'Visualization_Settings': {
                    'smoothing_factor': smoothing_sigma
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

# Add instructions at the bottom
st.markdown("""
---
### How to Use
1. Select organisms from the dropdown menus in the sidebar
2. Adjust the analysis parameters as needed:
   - Molecular weight range
   - Sample ratios
   - Abundance threshold
   - Number of capillaries
   - Smoothing factor for capillary views
3. View the results in the tabs above:
   - 2D Gel Plot tab shows the main visualization
   - Capillary Analysis tab shows individual capillary distributions
4. Toggle visibility of different lines in capillary views
5. Download all results using the download button
""")