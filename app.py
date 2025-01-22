import streamlit as st
import os
from datetime import datetime
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

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

# Function to get file paths
def get_file_paths(organism_id):
    fasta_path = os.path.join('data', 'fasta', f'fasta.v11.5.{organism_id}.fa')
    abundance_path = os.path.join('data', 'abundance', f'{organism_id}-WHOLE_ORGANISM-integrated.txt')
    return fasta_path, abundance_path

# Parse FASTA files
def parse_fasta(filepath):
    sequences = []
    if os.path.exists(filepath):
        for record in SeqIO.parse(filepath, "fasta"):
            sequences.append(record)
    return sequences

# Filter sequences
def filter_sequences(sequences):
    filtered_sequences = []
    for record in sequences:
        if 'X' not in str(record.seq):
            filtered_sequences.append(record)
    return filtered_sequences

# Calculate protein properties
def calculate_properties(sequences):
    properties = []
    for record in sequences:
        seq = str(record.seq)
        analysis = ProteinAnalysis(seq)
        mw = analysis.molecular_weight()
        pi = analysis.isoelectric_point()
        properties.append((record.id, mw, pi))
    return properties

# Parse abundance files
def parse_abundance(filepath):
    abundance_dict = {}
    if os.path.exists(filepath):
        with open(filepath, 'r') as file:
            for line in file:
                if not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        name, abundance = parts
                        abundance_dict[name] = float(abundance)
                    elif len(parts) == 3:
                        _, name, abundance = parts
                        abundance_dict[name] = float(abundance)
    return abundance_dict

# Filter by molecular weight
def filter_by_molecular_weight(properties, min_mw=None, max_mw=None):
    filtered_properties = []
    for prop in properties:
        protein_id, mw, pi = prop
        if (min_mw is None or mw >= min_mw) and (max_mw is None or mw <= max_mw):
            filtered_properties.append(prop)
    return filtered_properties

# Normalize abundance
def normalize_abundance(abundance1, abundance2, ratio1=50, ratio2=50, min_size=1, max_size=1000):
    all_values = [v * ratio1 / 100 for v in abundance1.values()] + [v * ratio2 / 100 for v in abundance2.values()]
    min_abundance = min(all_values) if all_values else 0
    max_abundance = max(all_values) if all_values else 1

    def normalize(value, ratio):
        if max_abundance == min_abundance:
            return max_size
        return min_size + (max_size - min_size) * (value * ratio / 100 - min_abundance) / (max_abundance - min_abundance)

    normalized_abundance1 = {k: normalize(v, ratio1) for k, v in abundance1.items()}
    normalized_abundance2 = {k: normalize(v, ratio2) for k, v in abundance2.items()}

    return normalized_abundance1, normalized_abundance2

# Logarithmic transform
def logarithmic_transform(mw):
    return 5.3779 * np.log(mw) - 6.4014

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
    
    st.sidebar.markdown("---")
    st.sidebar.subheader("Download Results")
    download_button = st.sidebar.button("Download All Results")

# Main app logic
try:
    # Get file paths
    fasta_path1, abundance_path1 = get_file_paths(ORGANISMS[organism1])
    fasta_path2, abundance_path2 = get_file_paths(ORGANISMS[organism2])

    with st.spinner("Processing data..."):
        # Process FASTA files
        sequences1 = parse_fasta(fasta_path1)
        sequences2 = parse_fasta(fasta_path2)
        
        if not sequences1 or not sequences2:
            st.error("Error reading FASTA files. Please check the file paths.")
            st.stop()
        
        # Filter and calculate properties
        filtered_sequences1 = filter_sequences(sequences1)
        filtered_sequences2 = filter_sequences(sequences2)
        
        protein_properties1 = calculate_properties(filtered_sequences1)
        protein_properties2 = calculate_properties(filtered_sequences2)
        
        # Parse abundance data
        abundance1 = parse_abundance(abundance_path1)
        abundance2 = parse_abundance(abundance_path2)
        
        if not abundance1 or not abundance2:
            st.error("Error reading abundance files. Please check the file paths.")
            st.stop()

        # Create visualization tabs
        tab1, tab2 = st.tabs(["2D Gel Plot", "Capillary Analysis"])
    
    with tab1:
        st.subheader("2D Gel Plot")
        gel_fig = create_gel_plot(
            protein_properties1, protein_properties2,
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
        min_pI = min(min(prop[2] for prop in protein_properties1),
                    min(prop[2] for prop in protein_properties2))
        max_pI = max(max(prop[2] for prop in protein_properties1),
                    max(prop[2] for prop in protein_properties2))
        capillary_width = (max_pI - min_pI) / num_capillaries
        
        # Create capillary plots
        cols = st.columns(2)
        for i in range(num_capillaries):
            cap_start = min_pI + i * capillary_width
            cap_end = min_pI + (i + 1) * capillary_width
            
            with cols[i % 2]:
                st.write(f"Capillary {i + 1}: pI {cap_start:.2f} - {cap_end:.2f}")
                
                filtered_props1 = [prop for prop in protein_properties1
                                 if cap_start <= prop[2] < cap_end]
                filtered_props2 = [prop for prop in protein_properties2
                                 if cap_start <= prop[2] < cap_end]
                
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
                        'start': min_pI + i * capillary_width,
                        'end': min_pI + (i + 1) * capillary_width
                    }
                } for i in range(num_capillaries)
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
3. View the results in the tabs above:
   - 2D Gel Plot tab shows the main visualization
   - Capillary Analysis tab shows individual capillary distributions
4. The visualizations update automatically when you change any parameters
""")