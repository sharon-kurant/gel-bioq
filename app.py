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
    grid_w = st.sidebar.number_input("Grid width", 100, 2000, 300)
    grid_h = st.sidebar.number_input("Grid height", 100, 2000, 300)
    random_seed = st.sidebar.number_input("Random seed", value=42, step=1)

    # — Capillary analysis —
    st.sidebar.subheader("Capillary analysis")
    cap_amount = st.sidebar.number_input("Number of capillaries", 1, 50, 10)
    smoothing = st.sidebar.slider("Smoothing σ", 0.0, 10.0, 1.0)
    
    if not st.session_state.generated:
        st.info(
            """
            ## 1. Organism selectors  
            - **Organism 1 / Organism 2**  
            - Pick two different proteomes from `data/fasta/` & `data/abundance/`.  
            - Their relative proportions in the mix are set by the “Sample ratios” slider.

            ---

            ## 2. Filters  
            - **Min / Max MW (Da)**  
            - Only proteins between these molecular‐weight bounds are shown.  
            - **Min normalized abundance**  
            - Discard low‐abundance spots below this threshold (0.0–1.0).

            ---

            ## 3. Sample mix ratios  
            - **Organism 1 ratio (%)**  
            - 0–100 % slider; Organism 2 gets the complement.  
            - Controls how strongly each proteome contributes.

            ---

            ## 4. Heatmap artefacts  
            Toggle each on/off; when on, adjust its parameters:

            ### A. Streaks  
            - **Orientation**: `both` / `horizontal` / `vertical`  
            - **Probability**: fraction of spots turned into smears (0.0–1.0)

            ### B. Spot trains  
            - **N parents**: top‑abundance spots spawning satellites (exact or range)  
            - **Decay factor**: 0.0–1.0, geometric drop per satellite  
            - **pI offset**: multiples of σₓ between satellites

            ### C. Smile distortion  
            - **Amplitude**: 0.0–0.1 (fraction of gel height)  
            - **Curve**: `quadratic` or `S`  
            - **Power**: exponent on xₙₒᵣₘ (≥ 1.0)  
            - **S‑coef**: asymmetry control when `S` is chosen

            ### D. Edge fading  
            - **Probability**: chance each side (L/R/T/B) is affected  
            - **Strength**: `(start,end)` fade factors over the outer 8 %

            ### E. Random dropout  
            - **Drop fraction**: uniform range 0.0–0.5 of proteins to remove

            ### F. Abundance variation  
            - **Range**: `(min,max)` truncation for normal multipliers around 1.0  
            - **σ**: standard deviation of that normal distribution

            ---

            ## 5. Blob σ‑multipliers  
            - **σₓ factor** / **σᵧ factor** (0.001–0.02)  
            - Scale each Gaussian spot’s width in pI/MW axes.

            ---

            ## 6. Grid & Seed  
            - **Width / Height** (100–2000)  
            - Pixel resolution (larger = finer detail).  
            - **Random seed**  
            - Fix to reproduce artefacts, dropout & noise.

            ---

            ## 7. Capillary analysis  
            - **# Capillaries** (1–50)  
            - Number of vertical pI slices to sum into line plots.  
            - **Smoothing σ** (0.0–10.0)  
            - Gaussian σ: small = sharp peaks; large = smoother curves.

            ---

            ### Generate  
            1. Click **Generate**.  
            2. A spinner appears (“Generating heatmap…”) while computations run.  
            3. **Heatmap** tab: your in‑silico 2D gel.  
            4. **Capillaries** tab: overlaid + individual line traces.  
            5. **Download Results** gives a ZIP with:  
            - `parameters.json`  
            - `heatmap.png` + `.csv`  
            - `capillaries_combined.png` + `.csv`  
            - `capillary_1.png/.csv`, `capillary_2.png/.csv`, …  

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