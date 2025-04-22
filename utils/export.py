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
        z.writestr('parameters.json',
                   json.dumps(parameters, indent=2).encode('utf-8'))

        # 2) heatmap figure & data
        fig_io = io.BytesIO()
        heatmap_fig.savefig(fig_io, format='png', bbox_inches='tight')
        fig_io.seek(0)
        z.writestr('heatmap.png', fig_io.read())

        csv_io = io.StringIO()
        np.savetxt(csv_io, heatmap_data, delimiter=',')
        z.writestr('heatmap.csv', csv_io.getvalue())

        # 3) capillary figures
        if capillary_figs:
            # Combined overlay = first fig
            combined_io = io.BytesIO()
            capillary_figs[0].savefig(combined_io, format='png', bbox_inches='tight')
            combined_io.seek(0)
            z.writestr('capillaries_combined.png', combined_io.read())

            # Individual capillaries = the rest
            for idx, fig in enumerate(capillary_figs[1:], start=1):
                cap_io = io.BytesIO()
                fig.savefig(cap_io, format='png', bbox_inches='tight')
                cap_io.seek(0)
                z.writestr(f'capillary_{idx}.png', cap_io.read())

        # 4) capillary data CSVs
        if capillary_data is not None and x_values is not None:
            # Combined profile = first entry
            combined_csv = io.StringIO()
            combined_arr = np.column_stack((x_values, capillary_data[0]))
            np.savetxt(
                combined_csv,
                combined_arr,
                delimiter=',',
                header='MW_kDa,Abundance',
                comments=''
            )
            z.writestr('capillaries_combined.csv', combined_csv.getvalue())

            # Individual capillaries
            for idx, y in enumerate(capillary_data[1:], start=1):
                cap_csv = io.StringIO()
                arr = np.column_stack((x_values, y))
                np.savetxt(
                    cap_csv,
                    arr,
                    delimiter=',',
                    header='MW_kDa,Abundance',
                    comments=''
                )
                z.writestr(f'capillary_{idx}.csv', cap_csv.getvalue())

    buf.seek(0)
    return buf


