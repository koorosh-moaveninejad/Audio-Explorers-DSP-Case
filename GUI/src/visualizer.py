import json

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st


def line_fig(x, y, title, name="Signal", dashed=False):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines",
            name=name,
            line=dict(dash="dash" if dashed else "solid"),
        )
    )
    fig.update_layout(
        title=title,
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=30, r=20, t=50, b=30),
        height=360,
    )
    return fig


def match_fig(res):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=res["t_crop_ms"], y=res["est_norm"], mode="lines", name="Estimated"))
    fig.add_trace(go.Scatter(x=res["t_crop_ms"], y=res["matched_norm"], mode="lines", name="Matched Template", line=dict(dash="dash")))
    fig.update_layout(
        title=f"Patient {res['PatientID']} · Estimated vs Matched Template",
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis_title="Time (ms)",
        yaxis_title="Normalized amplitude",
        margin=dict(l=30, r=20, t=50, b=40),
        height=380,
    )
    fig.update_xaxes(range=[4, 16])
    return fig


def scatter_quality_fig(df):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=df["SNR_dB"],
            y=df["RepeatCorr"],
            mode="markers+text",
            text=df["PatientID"],
            textposition="top center",
            marker=dict(size=14),
            name="Patients",
        )
    )
    fig.add_vline(x=3, line_dash="dash")
    fig.add_hline(y=0.35, line_dash="dash")
    fig.update_layout(
        title="SNR vs Split-Half Reproducibility",
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis_title="SNR (dB)",
        yaxis_title="Repeatability correlation",
        height=520,
        margin=dict(l=30, r=20, t=50, b=40),
    )
    return fig


def grid_waveforms_fig(results):
    cols = 3
    rows = (len(results) + cols - 1) // cols
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=[f"Patient {r['PatientID']}" for r in results])
    for i, r in enumerate(results):
        rr = i // cols + 1
        cc = i % cols + 1
        fig.add_trace(
            go.Scatter(x=r["t_ms"], y=r["oae_clean"], mode="lines", showlegend=False),
            row=rr,
            col=cc,
        )
        fig.update_xaxes(title_text="ms", range=[0, 20], row=rr, col=cc)
    fig.update_layout(
        title="All Estimated OAEs",
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        height=max(420, 280 * rows),
        margin=dict(l=20, r=20, t=60, b=30),
    )
    return fig


def dataframe_download_button(df, name):
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label=f"Download {name} CSV",
        data=csv,
        file_name=f"{name}.csv",
        mime="text/csv",
        use_container_width=True,
    )