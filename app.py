import streamlit as st
from utils.data_processing import (
    scan_available_organisms,
    parse_fasta,
    compute_protein_properties,
    filter_by_molecular_weight,
    filter_by_abundance_threshold,
    parse_abundance,
    normalize_abundance
)
from utils.heatmap import (
    generate_combined_heatmap,
    plot_heatmap
)
from utils.capillary import plot_capillary_traces
from utils.export import create_download_package


def main():
    st.set_page_config(page_title="2D Gel Dashboard", layout="wide")
    st.sidebar.title("2‑D Gel Settings")

    # — Organism selectors —
    orgs = scan_available_organisms()
    if len(orgs) < 2:
        st.error("Need at least two organism datasets in data/")
        return

    org1 = st.sidebar.selectbox("Organism 1", orgs, index=0)
    org2 = st.sidebar.selectbox("Organism 2", orgs, index=1)

    if org1 == org2:
        st.sidebar.error("Please pick two different organisms")

    # — Filters —
    st.sidebar.subheader("Filters")
    min_mw = st.sidebar.number_input("Min MW (Da)", min_value=1_000, value=1_000)
    max_mw = st.sidebar.number_input("Max MW (Da)", min_value=min_mw, value=500_000)
    min_abund = st.sidebar.slider("Min normalized abundance", 0.0, 1.0, 0.0)

    # — Sample ratios —
    st.sidebar.subheader("Sample ratios (%)")
    ratio1 = st.sidebar.slider("Organism 1 ratio", 0, 100, 80)
    ratio2 = 100 - ratio1
    st.sidebar.markdown(f"**Organism 2 ratio:** {ratio2}%")

    # — Artefact toggles & parameters —
    st.sidebar.subheader("Heatmap artefacts")

    # 1) Streaks
    apply_streaks = st.sidebar.checkbox("Apply streaks", value=True)
    if apply_streaks:
        streak_orient = st.sidebar.selectbox(
            "Orientation", ["both", "horizontal", "vertical"], index=0
        )
        streak_prob = st.sidebar.slider("Streak probability", 0.0, 1.0, 0.12)
    else:
        streak_orient, streak_prob = "both", 0.12

    # 2) Spot trains
    apply_spot_trains = st.sidebar.checkbox("Apply spot trains", value=True)
    if apply_spot_trains:
        train_N = st.sidebar.slider("N parents", 1, 100, (10, 50))
        train_decay = st.sidebar.slider("Decay factor", 0.0, 1.0, 0.55)
        train_offset = st.sidebar.slider("pI offset (σ)", 0.0, 2.0, 0.9)
    else:
        train_N, train_decay, train_offset = (10, 50), 0.55, 0.9

    # 3) Smile distortion
    apply_smile = st.sidebar.checkbox("Apply smile", value=True)
    if apply_smile:
        smile_rel_amp = st.sidebar.slider("Smile amplitude", 0.0, 0.1, 0.02)
        smile_curve = st.sidebar.selectbox("Curve shape", ["S", "quadratic"], index=0)
        smile_pow = st.sidebar.slider("Smile power", 1.0, 5.0, 2.0)
        smile_s_coef = st.sidebar.slider("S‑coef", 0.0, 1.0, 0.3)
    else:
        smile_rel_amp, smile_curve, smile_pow, smile_s_coef = 0.02, "S", 2.0, 0.3

    # 4) Edge fading
    apply_edges = st.sidebar.checkbox("Apply edges", value=True)
    if apply_edges:
        edge_prob = st.sidebar.slider("Edge probability", 0.0, 1.0, 0.60)
        edge_strength = st.sidebar.slider(
            "Edge strength", 0.5, 2.0, (1.35, 1.0)
        )
    else:
        edge_prob, edge_strength = 0.60, (1.35, 1.0)

    # 5) Dropout
    apply_dropout = st.sidebar.checkbox("Random dropout", value=False)
    if apply_dropout:
        drop_frac_range = st.sidebar.slider(
            "Drop fraction range", 0.0, 0.5, (0.0, 0.2)
        )
    else:
        drop_frac_range = (0.0, 0.0)

    # 6) Abundance variation
    apply_abundance_variation = st.sidebar.checkbox(
        "Abundance variation", value=False
    )
    if apply_abundance_variation:
        abundance_var_range = st.sidebar.slider(
            "Abund. range", 0.1, 5.0, (0.5, 2.0)
        )
        abundance_var_sd = st.sidebar.slider(
            "Abund. σ", 0.0, 1.0, 0.3
        )
    else:
        abundance_var_range, abundance_var_sd = (1.0, 1.0), 0.0

    # — Gaussian blob widths —
    st.sidebar.subheader("Blob σ‑multipliers")
    sigma_x_factor = st.sidebar.slider("σ_x factor", 0.001, 0.02, 0.005)
    sigma_y_factor = st.sidebar.slider("σ_y factor", 0.001, 0.02, 0.005)

    # — Grid size & seed —
    st.sidebar.subheader("Grid & Seed")
    grid_w = st.sidebar.number_input("Grid width", 100, 2000, 500)
    grid_h = st.sidebar.number_input("Grid height", 100, 2000, 500)
    random_seed = st.sidebar.number_input("Random seed", value=42, step=1)

    # — Capillary analysis —
    st.sidebar.subheader("Capillary analysis")
    cap_amount = st.sidebar.number_input("Number of capillaries", 1, 50, 10)
    smoothing = st.sidebar.slider("Smoothing σ", 0.0, 10.0, 1.0)

    # — Generate button —
    if st.sidebar.button("Generate"):
        # 1) Load & preprocess both organisms
        rec1 = parse_fasta(f"data/fasta/fasta.v11.5.{org1}.fa")
        props1 = compute_protein_properties(rec1)
        props1 = filter_by_molecular_weight(props1, min_mw, max_mw)
        props1 = filter_by_abundance_threshold(props1, min_abund)
        abund1 = parse_abundance(f"data/abundance/{org1}-WHOLE_ORGANISM-integrated.txt")

        rec2 = parse_fasta(f"data/fasta/fasta.v11.5.{org2}.fa")
        props2 = compute_protein_properties(rec2)
        props2 = filter_by_molecular_weight(props2, min_mw, max_mw)
        props2 = filter_by_abundance_threshold(props2, min_abund)
        abund2 = parse_abundance(f"data/abundance/{org2}-WHOLE_ORGANISM-integrated.txt")

        norm1, norm2 = normalize_abundance(
            abund1, abund2, ratio1, ratio2
        )

        # 2) Merge for heatmap
        merged1 = [(*t, norm1[t[0]]) for t in props1 if t[0] in norm1]
        merged2 = [(*t, norm2[t[0]]) for t in props2 if t[0] in norm2]
        combined = merged1 + merged2

        # 3) Generate & show heatmap
        heat, extent = generate_combined_heatmap(
            combined,
            grid_size=(grid_h, grid_w),
            apply_streaks=apply_streaks,
            streak_orient=streak_orient,
            streak_prob=streak_prob,
            apply_spot_trains=apply_spot_trains,
            train_N=train_N,
            train_decay=train_decay,
            train_offset=train_offset,
            apply_smile=apply_smile,
            smile_rel_amp=smile_rel_amp,
            smile_curve=smile_curve,
            smile_pow=smile_pow,
            smile_s_coef=smile_s_coef,
            apply_edges=apply_edges,
            edge_prob=edge_prob,
            edge_strength=edge_strength,
            apply_dropout=apply_dropout,
            drop_frac_range=drop_frac_range,
            apply_abundance_variation=apply_abundance_variation,
            abundance_var_range=abundance_var_range,
            abundance_var_sd=abundance_var_sd,
            sigma_x_factor=sigma_x_factor,
            sigma_y_factor=sigma_y_factor,
            random_seed=random_seed
        )
        fig = plot_heatmap(heat, extent)
        st.pyplot(fig)

        # 4) Capillary tab
        st.markdown("### Capillary analysis")
        plot_capillary_traces(heat, extent, cap_amount, smoothing)

        # 5) Download package
        params = dict(
            org1=org1, org2=org2, min_mw=min_mw, max_mw=max_mw,
            min_abund=min_abund, ratio1=ratio1, ratio2=ratio2,
            streaks=(apply_streaks, streak_orient, streak_prob),
            trains=(apply_spot_trains, train_N, train_decay, train_offset),
            smile=(apply_smile, smile_rel_amp, smile_curve, smile_pow, smile_s_coef),
            edges=(apply_edges, edge_prob, edge_strength),
            dropout=(apply_dropout, drop_frac_range),
            abund_var=(apply_abundance_variation, abundance_var_range, abundance_var_sd),
            sigma_factors=(sigma_x_factor, sigma_y_factor),
            grid=(grid_h, grid_w), seed=random_seed,
            capillaries=(cap_amount, smoothing)
        )
        buf = create_download_package(
            params, fig, heat, None, None, None
        )
        st.download_button("Download Results", data=buf, file_name="results.zip")


if __name__ == "__main__":
    main()

# import streamlit as st
# from datetime import datetime

# # Import utility functions
# from utils.data_processing import *
# from utils.visualization import create_gel_plot, create_capillary_plot
# from utils.export import create_download_package

# # Page configuration
# st.set_page_config(
#     page_title="2D Gel Electrophoresis Analysis",
#     layout="wide",
#     initial_sidebar_state="expanded"
# )

# # Get available organisms
# ORGANISMS = scan_available_organisms()

# if not ORGANISMS:
#     st.error("No organism data found. Please check the data/abundance directory.")
#     st.stop()

# # Title and description
# st.title("2D Gel Electrophoresis Analysis")
# st.markdown("Select organisms and adjust parameters to visualize the results.")

# # Sidebar controls
# with st.sidebar:
#     st.header("Settings")
    
#     # Organism selection
#     st.subheader("Organism Selection")
#     available_organisms = list(ORGANISMS.keys())
    
#     if len(available_organisms) < 2:
#         st.error("At least two organisms are required for comparison.")
#         st.stop()
    
#     organism1 = st.selectbox(
#         "First Organism",
#         options=available_organisms,
#         index=0,
#         key="organism1"
#     )
    
#     # Filter out the first selected organism from the second dropdown
#     remaining_organisms = [org for org in available_organisms if org != organism1]
#     organism2 = st.selectbox(
#         "Second Organism",
#         options=remaining_organisms,
#         index=0,
#         key="organism2"
#     )
    
#     # Analysis parameters
#     st.subheader("Analysis Parameters")
#     min_mw = st.number_input("Minimum Molecular Weight", value=0.000001, format="%.6f", key="min_mw")
#     max_mw = st.number_input("Maximum Molecular Weight", value=500000.0, key="max_mw")
#     min_normalized_abundance = st.number_input("Minimum Normalized Abundance", value=2.0, key="min_norm_abundance")
    
#     # Gaussian parameter
#     gaussian_std = st.slider("Gaussian Standard Deviation", 0.01, 1.0, 0.01, 0.01, key="gaussian_std")
    
#     # Sample ratios
#     st.subheader("Sample Ratios")
#     ratio1 = st.slider(f"{organism1} Ratio (%)", 0, 100, 95, key="ratio1")
#     ratio2 = st.slider(f"{organism2} Ratio (%)", 0, 100, 5, key="ratio2")
    
#     # Augmentation settings
#     st.subheader("Augmentation Settings")
#     augmentation_type = st.radio(
#         "Select Augmentation Type",
#         options=["None", "pI Shift", "MW Scale"],
#         key="augmentation_type"
#     )
    
#     if augmentation_type == "pI Shift":
#         pI_shift = st.slider(
#             "pI Shift Amount",
#             min_value=-1.0,
#             max_value=1.0,
#             value=0.0,
#             step=0.1,
#             help="Shift pI values left (negative) or right (positive)",
#             key="pi_shift"
#         )
#     elif augmentation_type == "MW Scale":
#         mw_scale = st.slider(
#             "MW Scale Factor",
#             min_value=0.5,
#             max_value=2.0,
#             value=1.0,
#             step=0.1,
#             help="Scale molecular weights (0.5 = squeeze, 2.0 = stretch)",
#             key="mw_scale"
#         )
    
#     st.markdown("---")
    
#     # Noise Settings
#     st.subheader("Noise Settings")
#     enable_mw_noise = st.checkbox(
#         "Add Random Noise to Molecular Weights (±50%)", 
#         value=False,
#         help="Randomly increase or decrease each protein's molecular weight by up to 50%",
#         key="enable_noise"
#     )
    
    
#     # Capillary settings
#     st.markdown("---")
#     st.subheader("Capillary Analysis")    
#     num_capillaries = st.slider("Number of Capillaries", 1, 100, 8, key="num_capillaries")
    
#     # Add spillage controls
#     enable_spillage = st.sidebar.checkbox(
#         "Enable Capillary Spillage", 
#         value=False, 
#         help="Allow proteins to contribute to adjacent capillaries based on their circle size"
#     )
    
#     # Visualization toggles
#     st.subheader("Show/Hide Plots")
#     show_organism1 = st.sidebar.checkbox(f"Show {organism1}", value=True, key="show_org1")
#     show_organism2 = st.sidebar.checkbox(f"Show {organism2}", value=True, key="show_org2")
#     show_sum = st.sidebar.checkbox("Show Sum", value=True, key="show_sum")
    
#     # Smoothing control
#     st.subheader("Visualization Settings")
#     smoothing_sigma = st.slider("Smoothing Factor", 1, 20, 5, key="smoothing")
    
#     st.markdown("---")
    
#     # Download section
#     st.subheader("Download Results")
#     download_button = st.button("Download All Results", key="download")
    
# # Main app logic
# try:
#     with st.spinner("Processing data..."):
#         # Process data for both organisms
#         properties1, abundance1 = process_organism_data(
#             ORGANISMS[organism1], min_mw, max_mw
#         )
#         properties2, abundance2 = process_organism_data(
#             ORGANISMS[organism2], min_mw, max_mw
#         )
        
#         # Apply noise if enabled
#         if enable_mw_noise:
#             properties1 = add_mw_noise(properties1)
#             properties2 = add_mw_noise(properties2)
        
#         # Store original properties for capillary calculations
#         original_properties1 = properties1.copy()
#         original_properties2 = properties2.copy()
        
#         # Apply pI shift if selected
#         if augmentation_type == "pI Shift":
#             properties1, actual_shift = shift_pI(properties1, pI_shift)
#             properties2, _ = shift_pI(properties2, pI_shift)
#             st.info(f"Applied pI shift of {actual_shift:.2f} units")
    
        
#         # Normalize abundances
#         normalized_abundance1, normalized_abundance2 = normalize_abundance(
#             abundance1, abundance2, ratio1, ratio2
#         )
        
        
#         # Create visualization tabs
#         tab1, tab2 = st.tabs(["2D Gel Plot", "Capillary Analysis"])
        
#         with tab1:
#             st.subheader("2D Gel Plot")
#             gel_fig, gel_data = create_gel_plot(
#                 properties1, properties2,
#                 normalized_abundance1, normalized_abundance2,
#                 organism1, organism2,
#                 ratio1, ratio2,
#                 augmentation_type,
#                 mw_scale=mw_scale if augmentation_type == "MW Scale" else 1.0
#             )
#             st.pyplot(gel_fig)
        
#         with tab2:
#             st.subheader("Capillary Analysis")
#             capillary_figs = []
#             capillary_data = []
#             x_values = None
            
#             # Calculate capillary ranges
#             min_pI, max_pI = get_pI_range(original_properties1, original_properties2)
#             capillary_ranges = calculate_capillary_ranges(min_pI, max_pI, num_capillaries)
            
#             # Create capillary plots
#             cols = st.columns(2)
#             for i, (cap_start, cap_end) in enumerate(capillary_ranges):
#                 with cols[i % 2]:
#                     st.write(f"Capillary {i + 1}: pI {cap_start:.2f} - {cap_end:.2f}")
                    
#                     # Filter proteins based on their original positions
#                     filtered_props1 = filter_by_pI_range(properties1, cap_start, cap_end)
#                     filtered_props2 = filter_by_pI_range(properties2, cap_start, cap_end)
                    
#                     # create_capillary_plot
#                     fig, data, x_vals = create_capillary_plot(
#                         filtered_props1, filtered_props2,
#                         normalized_abundance1, normalized_abundance2,
#                         organism1, organism2,
#                         smoothing_sigma,
#                         gaussian_std,
#                         show_organism1, show_organism2, show_sum,
#                         cap_start, cap_end,
#                         mw_scale=mw_scale if augmentation_type == "MW Scale" else 1.0,
#                         enable_spillage=enable_spillage  # Remove spillage_width parameter
#                     )
                    
#                     capillary_figs.append(fig)
#                     capillary_data.append(data)
#                     if x_values is None:
#                         x_values = x_vals
                    
#                     st.pyplot(fig)
                    
#         # Handle download if requested
#         if download_button:
#             # Update the parameters dictionary in the download section
#             parameters = {
#                 'Analysis_Parameters': {
#                     'min_molecular_weight': min_mw,
#                     'max_molecular_weight': max_mw,
#                     'min_normalized_abundance': min_normalized_abundance,
#                     'gaussian_std': gaussian_std
#                 },
#                 'Sample_Ratios': {
#                     organism1: ratio1,
#                     organism2: ratio2
#                 },
#                 'Augmentations': {
#                     'type': augmentation_type,
#                     'pI_shift': pI_shift if augmentation_type == "pI Shift" else 0,
#                     'mw_scale': mw_scale if augmentation_type == "MW Scale" else 1.0,
#                     'mw_noise_enabled': enable_mw_noise
#                 },
#                 'Capillary_Settings': {
#                     'number_of_capillaries': num_capillaries,
#                     'spillage_enabled': enable_spillage,
#                     'ranges': [
#                         {
#                             f'capillary_{i+1}': {
#                                 'pI_range': {
#                                     'start': ranges[0],
#                                     'end': ranges[1]
#                                 }
#                             }
#                         } for i, ranges in enumerate(capillary_ranges)
#                     ]
#                 },
#                 'Visualization_Settings': {
#                     'smoothing_factor': smoothing_sigma,
#                     'show_organism1': show_organism1,
#                     'show_organism2': show_organism2,
#                     'show_sum': show_sum
#                 }
#             }
            
#             zip_buffer = create_download_package(
#                 parameters=parameters,
#                 gel_plot_fig=gel_fig,
#                 gel_plot_data=gel_data,
#                 capillary_figs=capillary_figs,
#                 capillary_data=capillary_data,
#                 x_values=x_values,
#                 organism1=organism1,
#                 organism2=organism2
#             )
            
#             timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
#             st.sidebar.download_button(
#                 label="Click to Download",
#                 data=zip_buffer,
#                 file_name=f'gel_analysis_{timestamp}.zip',
#                 mime="application/zip"
#             )


# except Exception as e:
#     st.error(f"An error occurred: {str(e)}")
#     st.write("Please check that all data files are in the correct locations:")
#     st.code("""
#     gel/
#     ├── data/
#     │   ├── abundance/
#     │   │   ├── 4932-WHOLE_ORGANISM-integrated.txt
#     │   │   └── 511145-WHOLE_ORGANISM-integrated.txt
#     │   └── fasta/
#     │       ├── fasta.v11.5.4932.fa
#     │       └── fasta.v11.5.511145.fa
#     """)

# # Instructions at the bottom
# st.markdown("""
# ---
# ### How to Use
# 1. Select organisms from the dropdown menus in the sidebar

# 2. Adjust Analysis Parameters:
#    * Set molecular weight range (kDa)
#    * Adjust minimum normalized abundance
#    * Configure Gaussian standard deviation for peak shapes

# 3. Set Sample Ratios:
#    * Adjust relative proportions of each organism (0-100%)

# 4. Choose Augmentations (optional):
#    * pI Shift: Shift proteins left/right on pI axis (-1 to +1 units)
#    * MW Scale: Stretch/squeeze protein spots vertically (0.5x to 2x)
#    * MW Noise: Add random variation to molecular weights (±50%)

# 5. Configure Capillary Analysis:
#    * Set number of capillaries (1-100)
#    * Enable/disable capillary spillage
#    * Toggle visibility of individual organisms and sum
#    * Adjust smoothing factor for line plots

# 6. View Visualizations:
#    * 2D Gel Plot tab: Shows complete protein distribution
#    * Capillary Analysis tab: Shows protein distribution in pI ranges
#    * Use toggles to show/hide individual organisms and sum

# 7. Export Results:
#    * Click "Download Results" to get:
#      - All visualization parameters
#      - 2D gel plot image and data
#      - Capillary analysis plots and data
#      - Complete parameter settings

# Note: Augmentations and settings can be combined to explore different aspects of the protein distribution patterns.
# """)