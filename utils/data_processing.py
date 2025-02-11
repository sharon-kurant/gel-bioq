import os
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def get_file_paths(organism_id):
    """Get file paths for FASTA and abundance files."""
    fasta_path = os.path.join('data', 'fasta', f'fasta.v11.5.{organism_id}.fa')
    abundance_path = os.path.join('data', 'abundance', f'{organism_id}-WHOLE_ORGANISM-integrated.txt')
    return fasta_path, abundance_path

def parse_fasta(filepath):
    """Parse FASTA file and return sequences."""
    sequences = []
    if os.path.exists(filepath):
        for record in SeqIO.parse(filepath, "fasta"):
            sequences.append(record)
    return sequences

def filter_sequences(sequences):
    """Filter out sequences containing ambiguous characters."""
    filtered_sequences = []
    for record in sequences:
        if 'X' not in str(record.seq):
            filtered_sequences.append(record)
    return filtered_sequences

def calculate_properties(sequences):
    """Calculate molecular weight and isoelectric point for each sequence."""
    properties = []
    for record in sequences:
        seq = str(record.seq)
        analysis = ProteinAnalysis(seq)
        mw = analysis.molecular_weight()
        pi = analysis.isoelectric_point()
        properties.append((record.id, mw, pi))
    return properties

def parse_abundance(filepath):
    """Parse abundance file and return dictionary of protein abundances."""
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

def filter_by_molecular_weight(properties, min_mw=None, max_mw=None):
    """Filter properties by molecular weight range."""
    filtered_properties = []
    for prop in properties:
        protein_id, mw, pi = prop
        if (min_mw is None or mw >= min_mw) and (max_mw is None or mw <= max_mw):
            filtered_properties.append(prop)
    return filtered_properties

def normalize_abundance(abundance1, abundance2, ratio1=50, ratio2=50, min_size=1, max_size=1000):
    """Normalize abundance values with ratio."""
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

def process_organism_data(organism_id, min_mw=None, max_mw=None):
    """Process all data for a single organism."""
    fasta_path, abundance_path = get_file_paths(organism_id)
    
    # Process sequences
    sequences = parse_fasta(fasta_path)
    filtered_sequences = filter_sequences(sequences)
    properties = calculate_properties(filtered_sequences)
    
    # Apply molecular weight filter if specified
    if min_mw is not None or max_mw is not None:
        properties = filter_by_molecular_weight(properties, min_mw, max_mw)
    
    # Get abundance data
    abundance = parse_abundance(abundance_path)
    
    return properties, abundance

def filter_by_abundance(properties, abundance, min_abundance):
    """Filter properties by minimum abundance threshold."""
    return [prop for prop in properties if abundance.get(prop[0], 0) >= min_abundance]

# def get_pI_range(properties1, properties2):
#     """Calculate the overall pI range using original pI values."""
#     all_pIs = []
#     for props in [properties1, properties2]:
#         if props and len(props[0]) > 3:  # If properties contain both original and shifted pI
#             all_pIs.extend([prop[2] for prop in props])  # Use original pI (index 2)
#         else:
#             all_pIs.extend([prop[2] for prop in props])  # Use the only pI value
#     return min(all_pIs), max(all_pIs)

def calculate_capillary_ranges(min_pI, max_pI, num_capillaries):
    """Calculate the pI ranges for each capillary."""
    capillary_width = (max_pI - min_pI) / num_capillaries
    return [(min_pI + i * capillary_width, min_pI + (i + 1) * capillary_width)
            for i in range(num_capillaries)]

def filter_by_pI_range(properties, pI_start, pI_end, use_shifted=True):
    """
    Filter properties by pI range.
    If use_shifted is True, use the shifted pI values for filtering.
    """
    filtered_properties = []
    for prop in properties:
        protein_id, mw, pI = prop
        if isinstance(pI, tuple):
            original_pI, shifted_pI = pI
            check_pI = shifted_pI if use_shifted else original_pI
        else:
            check_pI = pI
            
        if pI_start <= check_pI < pI_end:
            filtered_properties.append(prop)
    return filtered_properties

def get_pI_range(properties1, properties2):
    """Calculate the overall pI range using original pI values."""
    all_pIs = []
    for props in [properties1, properties2]:
        for prop in props:
            pI = prop[2]
            if isinstance(pI, tuple):
                all_pIs.append(pI[0])  # Use original pI
            else:
                all_pIs.append(pI)
    return min(all_pIs), max(all_pIs)


def scan_available_organisms():
    """
    Scan abundance folder for available organisms and return a dictionary mapping
    organism names to their IDs.
    """
    organisms = {}
    abundance_dir = os.path.join('data', 'abundance')
    
    if not os.path.exists(abundance_dir):
        return organisms

    # Scan all files in abundance directory
    for filename in os.listdir(abundance_dir):
        if filename.endswith('-WHOLE_ORGANISM-integrated.txt'):
            filepath = os.path.join(abundance_dir, filename)
            organism_id = filename.split('-')[0]  # Get ID from filename
            
            # Read the first few lines to get organism name
            try:
                with open(filepath, 'r') as file:
                    for line in file:
                        if line.startswith('#name:'):
                            # Extract organism name from the line
                            full_name = line.split('-')[0].replace('#name:', '').strip()
                            # Add to dictionary
                            organisms[full_name] = organism_id
                            break
            except Exception as e:
                print(f"Error reading {filename}: {str(e)}")
                continue
    
    return organisms

def shift_pI(properties, shift_amount):
    """
    Shift pI values by a specified amount (-1 to +1).
    Returns properties with both original and shifted pI values.
    """
    shifted_properties = []
    for prop in properties:
        protein_id, mw, original_pI = prop
        shifted_pI = original_pI + shift_amount
        # Store both original and shifted pI values
        shifted_properties.append((protein_id, mw, original_pI, shifted_pI))
    return shifted_properties, shift_amount


def scale_MW(properties, scale_factor):
    """
    Scale molecular weights by a factor for visualization.
    Returns properties with both original and scaled MW values.
    """
    scaled_properties = []
    for protein_id, mw, pI in properties:
        scaled_mw = mw * scale_factor
        # Store both original and scaled MW: (original_mw, scaled_mw)
        scaled_properties.append((protein_id, (mw, scaled_mw), pI))
    return scaled_properties

def shift_pI(properties, shift_amount):
    """
    Shift pI values by a specified amount (-1 to +1).
    Returns properties with both original and shifted pI values.
    """
    shifted_properties = []
    for protein_id, mw, original_pI in properties:
        shifted_pI = original_pI + shift_amount
        # Keep the same 3-value tuple structure but store both pIs in the third position
        shifted_properties.append((protein_id, mw, (original_pI, shifted_pI)))
    return shifted_properties, shift_amount