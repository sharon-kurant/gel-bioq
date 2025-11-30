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

    # ─── Session state ────────────────────────────────────────
    if 'generated' not in st.session_state:
        st.session_state.generated = False

    # ─── Sidebar ───────────────────────────────────────────────
    st.sidebar.title("2‑D Gel Settings")
    ORGANISMS = scan_available_organisms()
    org_names = list(ORGANISMS.keys())
    org1 = st.sidebar.selectbox("Organism 1", org_names, index=0)
    org2 = st.sidebar.selectbox("Organism 2", org_names, index=1)
    org1_id, org2_id = ORGANISMS[org1], ORGANISMS[org2] # e.g. "882"
    
    if org1 == org2:
        st.sidebar.error("Please pick two different organisms")

    # — Filters —
    st.sidebar.subheader("Filters")
    min_mw = st.sidebar.number_input("Min MW (Da)", min_value=1_000, value=1_000)
    max_mw = st.sidebar.number_input("Max MW (Da)", min_value=min_mw, value=500_000)
    min_abund = st.sidebar.slider("Min normalized abundance", 0.0, 1.0, 0.0)

    # — Sample ratios —
    st.sidebar.subheader("Sample ratios (%)")
    ratio1 = st.sidebar.slider("Organism 1 ratio", 0, 100, 80)
    ratio2 = 100 - ratio1
    st.sidebar.markdown(f"**Organism 2 ratio:** {ratio2}%")

    # — Rendering —
    st.sidebar.subheader("Rendering")
    sigma_x_factor = st.sidebar.number_input(
        label="σ_x factor", min_value=0.001, step=0.0001, max_value=0.02, value=0.006, format="%.4f"
    )
    sigma_y_factor = st.sidebar.number_input(
        label="σ_y factor", min_value=0.001, step=0.0001, max_value=0.02, value=0.004, format="%.4f"
    )
    bleed_sigma = st.sidebar.slider("Stain bleed σ", 0.5, 8.0, 3.5)
    smear_strength = st.sidebar.slider("Smear strength", 0.0, 1.0, 0.45)
    paper_texture = st.sidebar.slider("Paper texture", 0.0, 0.6, 0.28)
    scan_noise = st.sidebar.slider("Scan noise", 0.0, 0.6, 0.2)
    tone_gamma = st.sidebar.slider("Tone curve γ", 0.3, 2.0, 0.8)
    haze_strength = st.sidebar.slider("Bloom strength", 0.0, 1.0, 0.25)
    grid_w = st.sidebar.number_input("Grid width", 200, 2000, 900)
    grid_h = st.sidebar.number_input("Grid height", 200, 2000, 1100)
    random_seed = st.sidebar.number_input("Random seed", value=1337, step=1)

    # — Capillary analysis —
    st.sidebar.subheader("Capillary analysis")
    cap_amount = st.sidebar.number_input("Number of capillaries", 1, 50, 6)
    smoothing = st.sidebar.slider("Smoothing σ", 0.0, 10.0, 1.5)
    
    if not st.session_state.generated:
        st.info(
            """
            1) Choose two proteomes from `data/fasta/` + `data/abundance/` and set their mix.
            2) Filter by molecular weight and normalized abundance.
            3) Use the **Rendering** section to sculpt the gel look: spot widths (σₓ/σᵧ),
               stain bleed, smear, paper texture, scan noise, tone curve, bloom, grid size, seed.
            4) Generate to view the gel and capillary traces; download a ZIP of the outputs.
            """
        )

    # — Generate button —
    if st.sidebar.button("Generate"):
        st.session_state.generated = True

        with st.spinner("Generating heatmap..."):
            # 1) Load & preprocess both organisms
            fasta_path1 = f"data/fasta/fasta.v11.5.{org1_id}.fa"
            rec1 = parse_fasta(fasta_path1)
            props1 = compute_protein_properties(rec1)
            props1 = filter_by_molecular_weight(props1, min_mw, max_mw)
            props1 = filter_by_abundance_threshold(props1, min_abund)
            abund1 = parse_abundance(f"data/abundance/{org1_id}-WHOLE_ORGANISM-integrated.txt")

            fasta_path2 = f"data/fasta/fasta.v11.5.{org2_id}.fa"
            rec2 = parse_fasta(fasta_path2)
            props2 = compute_protein_properties(rec2)
            props2 = filter_by_molecular_weight(props2, min_mw, max_mw)
            props2 = filter_by_abundance_threshold(props2, min_abund)
            abund2 = parse_abundance(f"data/abundance/{org2_id}-WHOLE_ORGANISM-integrated.txt")

            norm1, norm2 = normalize_abundance(abund1, abund2, ratio1, ratio2)
            merged1 = [(*t, norm1[t[0]]) for t in props1 if t[0] in norm1]
            merged2 = [(*t, norm2[t[0]]) for t in props2 if t[0] in norm2]
            combined = merged1 + merged2

            # 2) Generate heatmap
            heat, extent = generate_combined_heatmap(
                combined,
                grid_size=(grid_h, grid_w),
                sigma_x_factor=sigma_x_factor,
                sigma_y_factor=sigma_y_factor,
                bleed_sigma=bleed_sigma,
                smear_strength=smear_strength,
                paper_texture=paper_texture,
                scan_noise=scan_noise,
                tone_gamma=tone_gamma,
                haze_strength=haze_strength,
                random_seed=random_seed,
            )
            heat_fig = plot_heatmap(heat, extent)

        # 3) Tabbed display
        tab_heat, tab_caps = st.tabs(["Heatmap", "Capillaries"])

        with tab_heat:
            st.header("2D Gel Heatmap")
            st.pyplot(heat_fig)

        with tab_caps:
            st.header("Capillary Traces")
            # now we explicitly render each figure
            cap_figs, cap_data, cap_x = plot_capillary_traces(heat, extent, cap_amount, smoothing)
            for fig in cap_figs:
                st.pyplot(fig)

        # 4) Download package
        params = dict(
            org1=org1, org2=org2, min_mw=min_mw, max_mw=max_mw,
            min_abund=min_abund, ratio1=ratio1, ratio2=ratio2,
            rendering=dict(
                sigma_x_factor=sigma_x_factor,
                sigma_y_factor=sigma_y_factor,
                bleed_sigma=bleed_sigma,
                smear_strength=smear_strength,
                paper_texture=paper_texture,
                scan_noise=scan_noise,
                tone_gamma=tone_gamma,
                bloom=haze_strength,
                grid=(grid_h, grid_w),
                seed=random_seed,
            ),
            capillaries=(cap_amount, smoothing)
        )
        buf = create_download_package(
            params,
            heat_fig,
            heat,
            capillary_figs=cap_figs,       # or supply if you capture them
            capillary_data=cap_data,
            x_values=cap_x
        )
        st.download_button("Download Results", data=buf, file_name="results.zip")


if __name__ == "__main__":
    main()
