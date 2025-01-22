import streamlit as st
from datetime import datetime

# Import utility functions
from utils.data_processing import (
    process_organism_data,
    normalize_abundance,
    calculate_capillary_ranges,
    filter_by_pI_range,
    get_pI_range
)
from utils.visualization import create_gel_plot, create_capillary_plot
from utils.export import create_download_package

# Page configuration
st.set_page_config(
    page_title="2D Gel Electrophoresis Analysis",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Define organism mapping
ORGANISMS = {
    "S. cerevisiae": "4932",
    "E. coli": "511145"
}

# Title and description
st.title("2D Gel Electrophoresis Analysis")
st.markdown("Select organisms and adjust parameters to visualize the results.")

# Sidebar controls
with st.sidebar:
    st.header("Settings")
    
    # Organism selection
    st.subheader("Organism Selection")
    organism1 = st.selectbox(
        "First Organism",
        options=list(ORGANISMS.keys()),
        index=0,
        key="organism1"
    )
    organism2 = st.selectbox(
        "Second Organism",
        options=list(ORGANISMS.keys()),
        index=1,
        key="organism2"
    )
    
    # Analysis parameters
    st.subheader("Analysis Parameters")
    min_mw = st.number_input("Minimum Molecular Weight", value=0.000001, format="%.6f")
    max_mw = st.number_input("Maximum Molecular Weight", value=500000.0)
    min_normalized_abundance = st.number_input("Minimum Normalized Abundance", value=2.0)
    
    # Sample ratios
    st.subheader("Sample Ratios")
    ratio1 = st.slider(f"{organism1} Ratio (%)", 0, 100, 95)
    ratio2 = st.slider(f"{organism2} Ratio (%)", 0, 100, 5)
    
    # Capillary settings
    st.subheader("Capillary Analysis")
    num_capillaries = st.slider("Number of Capillaries", 1, 20, 8)
    
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
        
        # Normalize abundances
        normalized_abundance1, normalized_abundance2 = normalize_abundance(
            abundance1, abundance2, ratio1, ratio2
        )
        
        # Create visualization tabs
        tab1, tab2 = st.tabs(["2D Gel Plot", "Capillary Analysis"])
        
        with tab1:
            st.subheader("2D Gel Plot")
            gel_fig = create_gel_plot(
                properties1, properties2,
                normalized_abundance1, normalized_abundance2,
                organism1, organism2,
                ratio1, ratio2
            )
            st.pyplot(gel_fig)
        
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
                        show_organism1, show_organism2, show_sum
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
                    'min_normalized_abundance': min_normalized_abundance
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
                parameters,
                gel_fig,
                capillary_figs,
                capillary_data,
                x_values,
                organism1,
                organism2
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