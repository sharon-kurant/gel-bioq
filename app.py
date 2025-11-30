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
    org1 = st.sidebar.selectbox("Organism 1", org_names, index=0)
    org2 = st.sidebar.selectbox("Organism 2", org_names, index=1)
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

    # 5) Realism tweaks
    st.sidebar.subheader("Realism tweaks")
    spot_irregularity = st.sidebar.checkbox("Spot irregularity", value=True)
    spot_spread_jitter = st.sidebar.slider("Spot σ jitter", 0.5, 2.0, (0.75, 1.35))
    spot_angle_jitter = st.sidebar.checkbox("Angle jitter", value=True)
    background_intensity = st.sidebar.slider("Background stain", 0.0, 0.2, 0.02)
    texture_strength = st.sidebar.slider("Texture strength", 0.0, 0.5, 0.18)
    texture_scale = st.sidebar.slider("Texture scale", 0.2, 2.0, 0.85)
    haze_strength = st.sidebar.slider("Haze strength", 0.0, 1.0, 0.35)
    haze_sigma = st.sidebar.slider("Haze σ", 1.0, 20.0, 6.0)
    dynamic_range_gamma = st.sidebar.slider("Dynamic range γ", 0.3, 1.5, 0.65)

    # 6) Dropout
    apply_dropout = st.sidebar.checkbox("Random dropout", value=False)
    if apply_dropout:
        drop_frac_range = st.sidebar.slider(
            "Drop fraction range", 0.0, 0.5, (0.0, 0.2)
        )
    else:
        drop_frac_range = (0.0, 0.0)

    # 7) Abundance variation
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
    #sigma_x_factor = st.sidebar.slider("σ_x factor", 0.001, 0.02, 0.005)
    #sigma_y_factor = st.sidebar.slider("σ_y factor", 0.001, 0.02, 0.005)
    # st.sidebar.slider("Min normalized abundance", 0.0, 1.0, 0.0)
    sigma_x_factor = st.sidebar.number_input(label="σ_x factor", min_value=0.001, step=0.0001, max_value=0.02, value=0.005, format="%.4f")
    sigma_y_factor = st.sidebar.number_input(label="σ_y factor", min_value=0.001, step=0.0001, max_value=0.02, value=0.005, format="%.4f")
    # — Grid size & seed —
    st.sidebar.subheader("Grid & Seed")
    grid_w = st.sidebar.number_input("Grid width", 100, 2000, 300)
    grid_h = st.sidebar.number_input("Grid height", 100, 2000, 300)
    random_seed = st.sidebar.number_input("Random seed", value=42, step=1)
@@ -243,77 +255,97 @@ def main():

            norm1, norm2 = normalize_abundance(abund1, abund2, ratio1, ratio2)
            merged1 = [(*t, norm1[t[0]]) for t in props1 if t[0] in norm1]
            merged2 = [(*t, norm2[t[0]]) for t in props2 if t[0] in norm2]
            combined = merged1 + merged2

            # 2) Generate heatmap
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
                spot_irregularity=spot_irregularity,
                spot_spread_jitter=spot_spread_jitter,
                spot_angle_jitter=spot_angle_jitter,
                background_intensity=background_intensity,
                texture_strength=texture_strength,
                texture_scale=texture_scale,
                haze_strength=haze_strength,
                haze_sigma=haze_sigma,
                dynamic_range_gamma=dynamic_range_gamma,
                apply_dropout=apply_dropout,
                drop_frac_range=drop_frac_range,
                apply_abundance_variation=apply_abundance_variation,
                abundance_var_range=abundance_var_range,
                abundance_var_sd=abundance_var_sd,
                sigma_x_factor=sigma_x_factor,
                sigma_y_factor=sigma_y_factor,
                random_seed=random_seed
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
            streaks=(apply_streaks, streak_orient, streak_prob),
            trains=(apply_spot_trains, train_N, train_decay, train_offset),
            smile=(apply_smile, smile_rel_amp, smile_curve, smile_pow, smile_s_coef),
            edges=(apply_edges, edge_prob, edge_strength),
            realism=(
                spot_irregularity,
                spot_spread_jitter,
                spot_angle_jitter,
                background_intensity,
                texture_strength,
                texture_scale,
                haze_strength,
                haze_sigma,
                dynamic_range_gamma
            ),
            dropout=(apply_dropout, drop_frac_range),
            abund_var=(apply_abundance_variation, abundance_var_range, abundance_var_sd),
            sigma_factors=(sigma_x_factor, sigma_y_factor),
            grid=(grid_h, grid_w), seed=random_seed,
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
