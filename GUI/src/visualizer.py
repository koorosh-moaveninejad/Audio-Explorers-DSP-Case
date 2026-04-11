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
    cols = 1
    rows = (len(results) + cols - 1) // cols

    fig = make_subplots(
        rows=rows,
        cols=cols,
        subplot_titles=[
            f"Patient {r['PatientID']} | "
            f"{r.get('result','N/A')} | "
            f"Score: {r.get('assigned_score', r.get('best_score', 0)*100):.1f}%"
            for r in results
        ]
    )

    for i, r in enumerate(results):
        rr = i // cols + 1
        cc = i % cols + 1

        # Choose color based on PASS / REFER
        if r.get("result") == "PASS":
            color_est = "#00ffcc"
        else:
            color_est = "#ff4d4d"

        # Estimated waveform
        fig.add_trace(
            go.Scatter(
                x=r["t_crop_ms"],
                y=r["est_norm"],
                mode="lines",
                name="Estimated",
                line=dict(color=color_est, width=3),
                showlegend=(i == 0),
            ),
            row=rr,
            col=cc,
        )

        # Matched template waveform
        fig.add_trace(
            go.Scatter(
                x=r["t_crop_ms"],
                y=r["matched_norm"],
                mode="lines",
                name="Matched Template",
                line=dict(dash="dash", width=3),
                showlegend=(i == 0),
            ),
            row=rr,
            col=cc,
        )

        fig.update_xaxes(
            title_text="ms",
            range=[4, 16],
            row=rr,
            col=cc
        )

    fig.update_layout(
        title="Waveform Gallery: Estimated vs Matched Template",
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        height=max(500, 400 * rows),
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


def fft_fig(res):
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=res["fft_freq"],
            y=res["fft_mag"],
            mode="lines",
            name="Estimated FFT",
        )
    )

    fig.add_trace(
        go.Scatter(
            x=res["template_fft_freq"],
            y=res["template_fft_mag"],
            mode="lines",
            name="Template FFT",
            line=dict(dash="dash"),
        )
    )

    fig.update_layout(
        title=f"Patient {res['PatientID']} · Frequency Spectrum",
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Magnitude",
        margin=dict(l=30, r=20, t=50, b=40),
        height=380,
    )

    fig.update_xaxes(range=[500, 5000])
    return fig


def grid_fft_fig(results, cols=2):
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    rows = (len(results) + cols - 1) // cols

    fig = make_subplots(
        rows=rows,
        cols=cols,
        subplot_titles=[
            f"Patient {r['PatientID']} | "
            f"{r.get('result', 'N/A')} | "
            f"Score: {r.get('assigned_score', 0):.1f}%"
            for r in results
        ]
    )

    for i, r in enumerate(results):
        rr = i // cols + 1
        cc = i % cols + 1

        est_mag = np.asarray(r["fft_mag"], dtype=float).copy()
        tmpl_mag = np.asarray(r["template_fft_mag"], dtype=float).copy()

        # Normalize each spectrum for visual comparison
        if np.max(np.abs(est_mag)) > 0:
            est_mag = est_mag / np.max(np.abs(est_mag))

        if np.max(np.abs(tmpl_mag)) > 0:
            tmpl_mag = tmpl_mag / np.max(np.abs(tmpl_mag))

        color_map = {
            "PASS": "#00ffcc",
            "REFER": "#ff4d4d"
        }
        color_est = color_map.get(r.get("result"), "#aaaaaa")

        fig.add_trace(
            go.Scatter(
                x=r["fft_freq"],
                y=est_mag,
                mode="lines",
                name="Estimated FFT",
                line=dict(color=color_est, width=3),
                showlegend=(i == 0),
            ),
            row=rr,
            col=cc,
        )

        fig.add_trace(
            go.Scatter(
                x=r["template_fft_freq"],
                y=tmpl_mag,
                mode="lines",
                name="Template FFT",
                line=dict(dash="dash", width=3),
                showlegend=(i == 0),
            ),
            row=rr,
            col=cc,
        )

        fig.update_xaxes(
            title_text="Frequency (Hz)",
            range=[500, 5000],
            row=rr,
            col=cc
        )

        fig.update_yaxes(
            title_text="Normalized Magnitude",
            range=[0, 1.05],
            row=rr,
            col=cc
        )

    fig.update_layout(
        title="Spectrum Gallery: Estimated vs Matched Template",
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        height=max(500, 420 * rows),
        margin=dict(l=20, r=20, t=70, b=30),
    )

    return fig