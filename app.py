import streamlit as st
import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d

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
            
            # Filter and process data
            protein_properties1 = filter_by_molecular_weight(protein_properties1, min_mw, max_mw)
            protein_properties2 = filter_by_molecular_weight(protein_properties2, min_mw, max_mw)
            
            normalized_abundance1, normalized_abundance2 = normalize_abundance(
                abundance1, abundance2, ratio1, ratio2
            )
            
            # Create scatter plot
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Extract and transform data
            if protein_properties1:
                pi1, mw1 = zip(*[(prop[2], prop[1]) for prop in protein_properties1])
                log_mw1 = [logarithmic_transform(mw) for mw in mw1]
                sizes1 = [normalized_abundance1.get(prop[0], 1) for prop in protein_properties1]
                plt.scatter(pi1, log_mw1, alpha=0.5, color='blue', s=sizes1, 
                          label=f'{organism1} ({ratio1}%)')

            if protein_properties2:
                pi2, mw2 = zip(*[(prop[2], prop[1]) for prop in protein_properties2])
                log_mw2 = [logarithmic_transform(mw) for mw in mw2]
                sizes2 = [normalized_abundance2.get(prop[0], 1) for prop in protein_properties2]
                plt.scatter(pi2, log_mw2, alpha=0.5, color='red', s=sizes2, 
                          label=f'{organism2} ({ratio2}%)')
            
            plt.title('2-D Gel Plot')
            plt.xlabel('Isoelectric Point (pI)')
            plt.ylabel('Molecular Weight (kDa)')
            plt.legend()
            plt.grid(True, linestyle='--', linewidth=0.5)
            
            # Set y-axis labels
            yticks = plt.gca().get_yticks()
            ytick_labels = [f'{int(np.exp((tick + 6.4014) / 5.3779) / 1000):.0f}' for tick in yticks]
            plt.gca().set_yticklabels(ytick_labels)
            
            st.pyplot(fig)
            plt.close()
        
        with tab2:
            st.subheader("Capillary Analysis")
            
            # Add line visibility controls in sidebar
            st.sidebar.subheader("Capillary View Settings")
            show_organism1 = st.sidebar.checkbox(f"Show {organism1}", value=True)
            show_organism2 = st.sidebar.checkbox(f"Show {organism2}", value=True)
            show_sum = st.sidebar.checkbox("Show Sum", value=True)
            
            # Add smoothing control
            smoothing_sigma = st.sidebar.slider("Smoothing Factor", 1, 20, 5)
            
            # Calculate capillary ranges
            min_pI = min(
                min(prop[2] for prop in protein_properties1),
                min(prop[2] for prop in protein_properties2)
            )
            max_pI = max(
                max(prop[2] for prop in protein_properties1),
                max(prop[2] for prop in protein_properties2)
            )
            capillary_width = (max_pI - min_pI) / num_capillaries
            
            # Create capillary plots
            cols = st.columns(2)
            for i in range(num_capillaries):
                cap_start = min_pI + i * capillary_width
                cap_end = min_pI + (i + 1) * capillary_width
                
                with cols[i % 2]:
                    st.write(f"Capillary {i + 1}: pI {cap_start:.2f} - {cap_end:.2f}")
                    
                    # Filter proteins for this capillary
                    filtered_props1 = [
                        prop for prop in protein_properties1
                        if cap_start <= prop[2] < cap_end
                    ]
                    filtered_props2 = [
                        prop for prop in protein_properties2
                        if cap_start <= prop[2] < cap_end
                    ]
                    
                    # Create capillary plot with higher resolution
                    fig, ax = plt.subplots(figsize=(6, 4))
                    x_values = np.linspace(0, 20, 2000)  # Increased resolution
                    y1 = np.zeros_like(x_values)
                    y2 = np.zeros_like(x_values)
                    
                    # Calculate distributions
                    for prop in filtered_props1:
                        protein_id, mw, pi = prop
                        abundance = normalized_abundance1.get(protein_id, 0)
                        if abundance > 0 and mw > 0:
                            gaussian = norm.pdf(x_values, loc=mw/1000, 
                                            scale=max(abundance/100, 0.01))
                            y1 += gaussian * abundance
                    
                    for prop in filtered_props2:
                        protein_id, mw, pi = prop
                        abundance = normalized_abundance2.get(protein_id, 0)
                        if abundance > 0 and mw > 0:
                            gaussian = norm.pdf(x_values, loc=mw/1000, 
                                            scale=max(abundance/100, 0.01))
                            y2 += gaussian * abundance
                    
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
                    
                    # Set plot properties
                    plt.title(f'Capillary {i + 1}')
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
                    plt.ylim(0, max_y * 1.1)  # Add 10% padding
                    
                    # Display plot
                    st.pyplot(fig)
                    plt.close()

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