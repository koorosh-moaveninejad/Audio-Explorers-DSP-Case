import streamlit as st
import pandas as pd
from src.constants import COSMIC_CSS, APP_TITLE, APP_SUBTITLE
from src.io_utils import unpack_zip_to_temp, load_matlab_export
from src.visualizer import (
    line_fig,
    match_fig,
    fft_fig,
    grid_waveforms_fig,
    grid_fft_fig,
    dataframe_download_button,
)

st.set_page_config(
    page_title="TEOAE Cosmic Analyzer",
    page_icon="🪐",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(COSMIC_CSS, unsafe_allow_html=True)
st.markdown(f'<div class="cosmic-title">{APP_TITLE}</div>', unsafe_allow_html=True)
st.markdown(f'<div class="cosmic-subtitle">{APP_SUBTITLE}</div>', unsafe_allow_html=True)

with st.sidebar:
    st.markdown("### Input portal")

    source_mode = st.radio(
        "Choose data source",
        ["Local folder path", "ZIP upload"],
        help="Use a MATLAB export folder path if Streamlit runs on the same machine as your exported results. Use ZIP upload for browser-based usage.",
    )

    DEFAULT_EXPORT_PATH = "/Users/kourosh/Desktop/University/Self Study/Audio-Explorers20206/Diagnostics DSP/Patient Data/gui_export"

    base_path = DEFAULT_EXPORT_PATH
    uploaded_zip = None

    if source_mode == "Local folder path":
        base_path = st.text_input(
            "MATLAB export folder",
            value=DEFAULT_EXPORT_PATH,
            placeholder="/Users/you/.../gui_export",
        )
        st.caption("The folder should contain final_mapping.csv, template_scores.csv, patient_results.mat, and summary.json.")
    else:
        uploaded_zip = st.file_uploader(
            "Upload ZIP archive of the MATLAB export",
            type=["zip"]
        )
        st.caption("ZIP should contain final_mapping.csv, template_scores.csv, patient_results.mat, and summary.json.")

    run_btn = st.button("Run cosmic analysis", use_container_width=True, type="primary")

if "analysis_bundle" not in st.session_state:
    st.session_state.analysis_bundle = None
    st.session_state.temp_holder = None

if run_btn:
    try:
        if source_mode == "Local folder path":
            if not base_path.strip():
                st.error("Please provide the export folder path.")
            else:
                st.session_state.analysis_bundle = load_matlab_export(base_path.strip())
        else:
            if uploaded_zip is None:
                st.error("Please upload a ZIP archive first.")
            else:
                temp_holder, extracted_root = unpack_zip_to_temp(uploaded_zip)
                st.session_state.temp_holder = temp_holder
                st.session_state.analysis_bundle = load_matlab_export(extracted_root)
    except Exception as e:
        st.session_state.analysis_bundle = None
        st.error(f"Loading failed: {e}")

bundle = st.session_state.analysis_bundle

if bundle is None:
    st.markdown(
        """
        <div class="cosmic-card">
        <h3>Mission briefing</h3>
        <p>This app reproduces your MATLAB workflow in Python and adds:</p>
            <ul>
            <li>folder or ZIP-based input</li>
            <li>repeatability-based PASS / REFER decision</li>
            <li>global competitive template assignment</li>
            <li>patient-by-patient waveform and spectrum comparison</li>
            <li>gallery views for all patients</li>
        </ul>
        <p>Load data from the sidebar, then launch the analysis.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    final_df = bundle["final_df"]
    patient_results = bundle["patient_results"]
    scores_df = bundle["scores_df"]

    c1, c2, c3 = st.columns(3)
    c1.metric("Patients", len(patient_results))
    c2.metric("PASS (RepeatCorr)", int((final_df["Result"] == "PASS").sum()))
    c3.metric("Templates", len(bundle["templates"]))

    tabs = st.tabs([
        "Overview",
        "Patient explorer",
        "Template matches",
        "Spectrum explorer",
        "Waveform gallery",
        "Spectrum gallery",
        "Downloads",
    ])

    with tabs[0]:
         st.markdown("### Final mission mapping")
         st.table(final_df.set_index(pd.Index([""] * len(final_df))))

    with tabs[1]:
        patient_ids = [r["PatientID"] for r in patient_results]
        selected_pid = st.selectbox("Select patient", patient_ids)
        res = next(r for r in patient_results if r["PatientID"] == selected_pid)

        a, b, c = st.columns(3)
        a.metric("Assigned template", res.get("assigned_template", "N/A"))
        b.metric("Repeat corr", f"{res['repeat_corr']:.2f}")
        c.metric("Match score", f"{res.get('assigned_score', 0):.1f}%")
        # st.caption(f"Template-independent OAE exists: {'Yes' if res['oae_exists_rule'] else 'No'}")

        st.plotly_chart(
            line_fig(res["t_ms"], res["oae_clean"], f"Patient {selected_pid} · Estimated OAE", name="Estimated OAE"),
            use_container_width=True,
            key=f"patient_line_{selected_pid}",
        )

        st.plotly_chart(
            match_fig(res),
            use_container_width=True,
            key=f"patient_match_{selected_pid}",
        )

        st.plotly_chart(
            fft_fig(res),
            use_container_width=True,
            key=f"patient_fft_{selected_pid}",
        )

    with tabs[2]:
        st.markdown("### All template scores")
        st.dataframe(scores_df, use_container_width=True, hide_index=True)
        selected_pid = st.selectbox("Inspect template scores for patient", [r["PatientID"] for r in patient_results], key="templ_sel")
        patient_scores = scores_df[scores_df["PatientID"] == selected_pid].sort_values("Score", ascending=False)
        st.dataframe(patient_scores, use_container_width=True, hide_index=True)

    # with tabs[3]:
    #     st.plotly_chart(scatter_quality_fig(existence_df), use_container_width=True)
    #     st.markdown("### Quality table")
    #     st.dataframe(existence_df, use_container_width=True, hide_index=True)

    with tabs[3]:
        patient_ids = [r["PatientID"] for r in patient_results]
        selected_pid_fft = st.selectbox("Select patient for spectrum", patient_ids, key="fft_patient")
        res_fft = next(r for r in patient_results if r["PatientID"] == selected_pid_fft)
        st.plotly_chart(fft_fig(res_fft), use_container_width=True, key="fft_fig_spectrum_tab")

    with tabs[4]:
        st.plotly_chart(
        grid_waveforms_fig(patient_results),
        use_container_width=True,
        key="waveform_gallery",
)

    with tabs[5]:
        st.plotly_chart(grid_fft_fig(patient_results), use_container_width=True, key="spectrum_gallery")

    with tabs[6]:
        col1, col2 = st.columns(2)
        with col1:
            dataframe_download_button(final_df, "final_mapping")
            # dataframe_download_button(existence_df, "oae_existence")
        with col2:
            dataframe_download_button(scores_df, "template_scores")

        import json
        summary_json = json.dumps(
            {
                "final_mapping": final_df.to_dict(orient="records"),
            },
            indent=2,
        )
        st.download_button(
            "Download summary JSON",
            data=summary_json,
            file_name="teoae_summary.json",
            mime="application/json",
            use_container_width=True,
        )

st.markdown("<div style='height: 40px'></div>", unsafe_allow_html=True)