# utils/export.py

import io
import json
import zipfile
import numpy as np

def create_download_package(
        parameters: dict,
        heatmap_fig,
        heatmap_data: np.ndarray,
        capillary_figs: list = None,
        capillary_data: list = None,
        x_values: np.ndarray = None
    ) -> io.BytesIO:
    """
    Bundle parameters, figures, and raw data into an in‑memory ZIP.

    Contents:
      - parameters.json
      - heatmap.png
      - heatmap.csv
      - capillary_i.png (if capillary_figs provided)
      - capillary_i.csv (if capillary_data & x_values provided)

    Returns
    -------
    buf : io.BytesIO
        A seekable buffer containing the ZIP file.
    """
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, mode='w', compression=zipfile.ZIP_DEFLATED) as z:
        # 1) parameters.json
        params_bytes = json.dumps(parameters, indent=2).encode('utf-8')
        z.writestr('parameters.json', params_bytes)

        # 2) heatmap figure
        fig_io = io.BytesIO()
        heatmap_fig.savefig(fig_io, format='png', bbox_inches='tight')
        fig_io.seek(0)
        z.writestr('heatmap.png', fig_io.read())

        # 3) heatmap data as CSV
        csv_io = io.StringIO()
        np.savetxt(csv_io, heatmap_data, delimiter=',')
        z.writestr('heatmap.csv', csv_io.getvalue())

        # 4) capillary figures
        if capillary_figs:
            for idx, fig in enumerate(capillary_figs, start=1):
                cap_io = io.BytesIO()
                fig.savefig(cap_io, format='png', bbox_inches='tight')
                cap_io.seek(0)
                z.writestr(f'capillary_{idx}.png', cap_io.read())

        # 5) capillary data CSVs
        if capillary_data is not None and x_values is not None:
            for idx, y in enumerate(capillary_data, start=1):
                cap_csv = io.StringIO()
                data = np.column_stack((x_values, y))
                np.savetxt(
                    cap_csv,
                    data,
                    delimiter=',',
                    header='MW_kDa,Abundance',
                    comments=''
                )
                z.writestr(f'capillary_{idx}.csv', cap_csv.getvalue())

    buf.seek(0)
    return buf

# import yaml
# import io
# import csv
# import zipfile
# import pandas as pd
# from datetime import datetime

# def save_parameters_yaml(parameters):
#     """Save parameters to YAML format."""
#     return yaml.dump(parameters, sort_keys=False)

# def save_plot_to_bytes(fig):
#     """Convert matplotlib figure to bytes."""
#     buf = io.BytesIO()
#     fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
#     buf.seek(0)
#     return buf

# def save_csv_to_bytes(data_dict):
#     """Convert dictionary to CSV format."""
#     buf = io.StringIO()
#     writer = csv.writer(buf)
#     writer.writerows(data_dict.items())
#     return buf.getvalue()

# def create_download_package(
#     parameters,
#     gel_plot_fig,
#     gel_plot_data,
#     capillary_figs,
#     capillary_data,
#     x_values,
#     organism1,
#     organism2
# ):
#     """Create a ZIP file containing all analysis results."""
#     timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
#     zip_buffer = io.BytesIO()
#     with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
#         # Save parameters YAML
#         yaml_content = save_parameters_yaml(parameters)
#         zip_file.writestr(f'parameters_{timestamp}.yaml', yaml_content)
        
#         # Save 2D gel plot image and data
#         gel_plot_bytes = save_plot_to_bytes(gel_plot_fig)
#         zip_file.writestr(f'2d_gel_plot_{timestamp}.png', gel_plot_bytes.getvalue())
        
#         # Save 2D gel plot data as CSV
#         # Fill missing values with NaN for proper alignment
#         max_length = max(len(gel_plot_data['protein_ids']),
#                         len(gel_plot_data['pI_values']),
#                         len(gel_plot_data['mw_values']),
#                         len(gel_plot_data['abundance1']),
#                         len(gel_plot_data.get('abundance2', [])))
        
#         df_gel = pd.DataFrame({
#             'Protein_ID': gel_plot_data['protein_ids'] + [None] * (max_length - len(gel_plot_data['protein_ids'])),
#             'pI': gel_plot_data['pI_values'] + [None] * (max_length - len(gel_plot_data['pI_values'])),
#             'MW_kDa': gel_plot_data['mw_values'] + [None] * (max_length - len(gel_plot_data['mw_values'])),
#             f'{organism1}_abundance': gel_plot_data['abundance1'] + [None] * (max_length - len(gel_plot_data['abundance1'])),
#             f'{organism2}_abundance': (gel_plot_data.get('abundance2', []) + [None] * 
#                                      (max_length - len(gel_plot_data.get('abundance2', []))))
#         })
        
#         csv_buffer = io.StringIO()
#         df_gel.to_csv(csv_buffer, index=False)
#         zip_file.writestr(f'2d_gel_data_{timestamp}.csv', csv_buffer.getvalue())
        
#         # Save capillary plots and data
#         for i, (fig, data) in enumerate(zip(capillary_figs, capillary_data)):
#             # Save plot
#             cap_bytes = save_plot_to_bytes(fig)
#             zip_file.writestr(f'capillary_{i+1}_{timestamp}.png', cap_bytes.getvalue())
            
#             # Save CSV data
#             df = pd.DataFrame({
#                 'Molecular_Weight_kDa': x_values,
#                 f'{organism1}_Volume': data['y1'],
#                 f'{organism2}_Volume': data['y2'],
#                 'Sum_Volume': data['y_sum']
#             })
#             csv_buffer = io.StringIO()
#             df.to_csv(csv_buffer, index=False)
#             zip_file.writestr(f'capillary_{i+1}_data_{timestamp}.csv', csv_buffer.getvalue())
    
#     zip_buffer.seek(0)
#     return zip_buffer