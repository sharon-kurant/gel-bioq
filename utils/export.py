import yaml
import io
import csv
import zipfile
import pandas as pd
from datetime import datetime

def save_parameters_yaml(parameters):
    """Save parameters to YAML format."""
    return yaml.dump(parameters, sort_keys=False)

def save_plot_to_bytes(fig):
    """Convert matplotlib figure to bytes."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    buf.seek(0)
    return buf

def save_csv_to_bytes(data_dict):
    """Convert dictionary to CSV format."""
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerows(data_dict.items())
    return buf.getvalue()

def create_download_package(
    parameters,
    gel_plot_fig,
    gel_plot_data,
    capillary_figs,
    capillary_data,
    x_values,
    organism1,
    organism2
):
    """Create a ZIP file containing all analysis results."""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        # Save parameters YAML
        yaml_content = save_parameters_yaml(parameters)
        zip_file.writestr(f'parameters_{timestamp}.yaml', yaml_content)
        
        # Save 2D gel plot image and data
        gel_plot_bytes = save_plot_to_bytes(gel_plot_fig)
        zip_file.writestr(f'2d_gel_plot_{timestamp}.png', gel_plot_bytes.getvalue())
        
        # Save 2D gel plot data as CSV
        # Fill missing values with NaN for proper alignment
        max_length = max(len(gel_plot_data['protein_ids']),
                        len(gel_plot_data['pI_values']),
                        len(gel_plot_data['mw_values']),
                        len(gel_plot_data['abundance1']),
                        len(gel_plot_data.get('abundance2', [])))
        
        df_gel = pd.DataFrame({
            'Protein_ID': gel_plot_data['protein_ids'] + [None] * (max_length - len(gel_plot_data['protein_ids'])),
            'pI': gel_plot_data['pI_values'] + [None] * (max_length - len(gel_plot_data['pI_values'])),
            'MW_kDa': gel_plot_data['mw_values'] + [None] * (max_length - len(gel_plot_data['mw_values'])),
            f'{organism1}_abundance': gel_plot_data['abundance1'] + [None] * (max_length - len(gel_plot_data['abundance1'])),
            f'{organism2}_abundance': (gel_plot_data.get('abundance2', []) + [None] * 
                                     (max_length - len(gel_plot_data.get('abundance2', []))))
        })
        
        csv_buffer = io.StringIO()
        df_gel.to_csv(csv_buffer, index=False)
        zip_file.writestr(f'2d_gel_data_{timestamp}.csv', csv_buffer.getvalue())
        
        # Save capillary plots and data
        for i, (fig, data) in enumerate(zip(capillary_figs, capillary_data)):
            # Save plot
            cap_bytes = save_plot_to_bytes(fig)
            zip_file.writestr(f'capillary_{i+1}_{timestamp}.png', cap_bytes.getvalue())
            
            # Save CSV data
            df = pd.DataFrame({
                'Molecular_Weight_kDa': x_values,
                f'{organism1}_Volume': data['y1'],
                f'{organism2}_Volume': data['y2'],
                'Sum_Volume': data['y_sum']
            })
            csv_buffer = io.StringIO()
            df.to_csv(csv_buffer, index=False)
            zip_file.writestr(f'capillary_{i+1}_data_{timestamp}.csv', csv_buffer.getvalue())
    
    zip_buffer.seek(0)
    return zip_buffer