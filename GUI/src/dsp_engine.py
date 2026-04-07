import numpy as np
import pandas as pd
import streamlit as st
from scipy.signal import butter, filtfilt, correlate

from src.constants import (
    DEFAULT_BANDPASS_HIGH,
    DEFAULT_BANDPASS_LOW,
    DEFAULT_FILTER_ORDER,
    DEFAULT_HEAD_LEN,
    DEFAULT_MAX_LAG,
    DEFAULT_NOISE_END,
    DEFAULT_NOISE_START,
    DEFAULT_REPEAT_CORR_THRESHOLD,
    DEFAULT_RESPONSE_END,
    DEFAULT_RESPONSE_START,
    DEFAULT_SNR_THRESHOLD,
)
from src.io_utils import read_wav_any, load_json_any


def rms(x):
    x = np.asarray(x).ravel()
    if len(x) == 0:
        return 0.0
    return float(np.sqrt(np.mean(np.square(x))))


def safe_corrcoef(x, y):
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()
    if len(x) == 0 or len(y) == 0:
        return 0.0
    if np.allclose(np.std(x), 0) or np.allclose(np.std(y), 0):
        return 0.0
    c = np.corrcoef(x, y)[0, 1]
    if np.isnan(c):
        return 0.0
    return float(c)


def bandpass_filter(x, fs, low=DEFAULT_BANDPASS_LOW, high=DEFAULT_BANDPASS_HIGH, order=DEFAULT_FILTER_ORDER):
    nyq = 0.5 * fs
    low_n = max(low / nyq, 1e-6)
    high_n = min(high / nyq, 0.999999)
    if high_n <= low_n:
        return x.copy()
    b, a = butter(order, [low_n, high_n], btype="band")
    return filtfilt(b, a, x)


def align_to_reference(seg, ref, max_lag=DEFAULT_MAX_LAG, head_len=DEFAULT_HEAD_LEN):
    seg = np.asarray(seg).ravel()
    ref = np.asarray(ref).ravel()
    head = min(head_len, len(seg), len(ref))
    seg_part = seg[:head]
    ref_part = ref[:head]
    corr = correlate(seg_part, ref_part, mode="full")
    lags = np.arange(-len(ref_part) + 1, len(seg_part))
    mask = (lags >= -max_lag) & (lags <= max_lag)
    corr = corr[mask]
    lags = lags[mask]
    if len(corr) == 0:
        return seg.copy()
    best_lag = int(lags[np.argmax(corr)])
    return np.roll(seg, -best_lag)


def extract_epochs(rec, active, epoch_size, ref=None):
    starts = np.where(np.diff(np.concatenate([[0], (active > 0.5).astype(int)])) == 1)[0]
    epochs = []
    if ref is not None:
        ref = np.asarray(ref).ravel()

    for s in starts:
        e = s + epoch_size
        if e <= len(rec):
            seg = rec[s:e].astype(np.float64)
            seg = seg - np.mean(seg)
            if ref is None:
                ref = seg[: min(DEFAULT_HEAD_LEN, len(seg))].copy()
                shifted = seg
            else:
                shifted = align_to_reference(seg, ref)
            epochs.append(shifted)

    if len(epochs) == 0:
        return np.zeros((epoch_size, 0), dtype=np.float64), ref
    return np.column_stack(epochs), ref


def process_patient_folder(patient_dir, templates):
    p_id = patient_dir.name.replace("patient_", "")

    rec_path = patient_dir / f"Patient_{p_id}_rec.wav"
    info_path = patient_dir / f"Patient_{p_id}_info.json"

    if not rec_path.exists() or not info_path.exists():
        raise FileNotFoundError(f"Missing rec/info files for patient {p_id}")

    fs, rec = read_wav_any(rec_path)
    info = load_json_any(info_path)
    epoch_size = int(info["epochSize"])

    epoch_data = []
    ref = None
    for t in ["A", "B", "C", "D"]:
        active_path = patient_dir / f"Patient_{p_id}_active{t}.wav"
        if not active_path.exists():
            raise FileNotFoundError(f"Missing active{t} file for patient {p_id}")
        fs_active, act = read_wav_any(active_path)
        if fs_active != fs:
            raise ValueError(f"Sample-rate mismatch in patient {p_id}, active{t}")
        epochs, ref = extract_epochs(rec, act, epoch_size, ref=ref)
        epoch_data.append(epochs)

    min_e = min(arr.shape[1] for arr in epoch_data) if epoch_data else 0
    if min_e == 0:
        oae_raw = np.zeros(epoch_size, dtype=np.float64)
        oae_raw1 = np.zeros(epoch_size, dtype=np.float64)
        oae_raw2 = np.zeros(epoch_size, dtype=np.float64)
    else:
        acc = np.zeros(epoch_size, dtype=np.float64)
        acc1 = np.zeros(epoch_size, dtype=np.float64)
        acc2 = np.zeros(epoch_size, dtype=np.float64)
        v = v1 = v2 = 0

        for k in range(min_e):
            quad = epoch_data[0][:, k] + epoch_data[1][:, k] + epoch_data[2][:, k] + epoch_data[3][:, k]
            if np.all(np.isfinite(quad)):
                acc += quad
                v += 1
                if (k % 2) == 0:
                    acc1 += quad
                    v1 += 1
                else:
                    acc2 += quad
                    v2 += 1

        oae_raw = acc / max(v, 1)
        oae_raw1 = acc1 / max(v1, 1)
        oae_raw2 = acc2 / max(v2, 1)

    t = np.arange(epoch_size) / fs
    win_env = np.exp(-((t - 0.008) ** 2) / (2 * 0.004**2))

    oae_clean = bandpass_filter(oae_raw * win_env, fs)
    oae_clean1 = bandpass_filter(oae_raw1 * win_env, fs)
    oae_clean2 = bandpass_filter(oae_raw2 * win_env, fs)

    noise_est = (
        np.mean(epoch_data[0], axis=1) + np.mean(epoch_data[1], axis=1) - np.mean(epoch_data[2], axis=1) - np.mean(epoch_data[3], axis=1)
        if min_e > 0 else np.zeros(epoch_size, dtype=np.float64)
    )
    noise_clean = bandpass_filter(noise_est, fs)

    oae_clean = oae_clean - 0.1 * noise_clean
    oae_clean1 = oae_clean1 - 0.1 * noise_clean
    oae_clean2 = oae_clean2 - 0.1 * noise_clean

    resp_idx = (t > DEFAULT_RESPONSE_START) & (t < DEFAULT_RESPONSE_END)
    # noise_idx = (t > DEFAULT_NOISE_START) & (t < DEFAULT_NOISE_END)

    # resp_rms = rms(oae_clean[resp_idx])
    # noise_rms = rms(oae_clean[noise_idx])
    # snr_db = 20 * np.log10(resp_rms / max(noise_rms, np.finfo(float).eps))
    repeat_corr = safe_corrcoef(oae_clean1[resp_idx], oae_clean2[resp_idx])

    template_scores = []
    best_template_name = "N/A"
    best_score = 0.0
    best_target_crop = np.zeros(np.sum(resp_idx), dtype=np.float64)
    oae_crop = oae_clean[resp_idx]

    for name, target in templates.items():
        target = np.asarray(target, dtype=np.float64).ravel()
        if len(target) < len(oae_clean):
            target_adj = np.pad(target, (0, len(oae_clean) - len(target)))
        else:
            target_adj = target[: len(oae_clean)]
        target_crop = target_adj[resp_idx]

        if len(oae_crop) == 0 or len(target_crop) == 0 or np.std(oae_crop) == 0 or np.std(target_crop) == 0:
            score = 0.0
        else:
            corr = correlate(oae_crop, target_crop, mode="full")
            denom = np.linalg.norm(oae_crop) * np.linalg.norm(target_crop)
            score = float(np.max(corr / denom)) if denom > 0 else 0.0
            if not np.isfinite(score) or score < 0:
                score = 0.0

        template_scores.append({"PatientID": p_id, "Template": name, "Score": score})
        if score > best_score:
            best_score = score
            best_template_name = name
            best_target_crop = target_crop.copy()

    est_norm = oae_crop / np.linalg.norm(oae_crop) if np.linalg.norm(oae_crop) > 0 else oae_crop.copy()
    target_norm = best_target_crop / np.linalg.norm(best_target_crop) if np.linalg.norm(best_target_crop) > 0 else best_target_crop.copy()

    if len(est_norm) > 0 and len(target_norm) > 0:
        corr = correlate(est_norm, target_norm, mode="full")
        lags = np.arange(-len(target_norm) + 1, len(est_norm))
        best_lag = int(lags[np.argmax(corr)]) if len(corr) else 0
        matched_norm = np.roll(target_norm, best_lag)
    else:
        matched_norm = np.zeros_like(est_norm)

    return {
        "PatientID": p_id,
        "fs": fs,
        "t_ms": t * 1000.0,
        "t_crop_ms": t[resp_idx] * 1000.0,
        "oae_clean": oae_clean,
        "oae_crop": oae_crop,
        "est_norm": est_norm,
        "matched_norm": matched_norm,
        "best_template": best_template_name,
        "best_score": best_score,
        "repeat_corr": float(repeat_corr),
        "template_scores": template_scores,
    }


def run_analysis(root_dir, templates, patient_dirs):
    patient_results = []
    all_scores = []

    progress = st.progress(0, text="Launching the scanner through the nebula...")
    status = st.empty()

    for idx, p_dir in enumerate(patient_dirs, start=1):
        status.info(f"Processing {p_dir.name} ({idx}/{len(patient_dirs)})")
        res = process_patient_folder(p_dir, templates)
        patient_results.append(res)
        all_scores.extend(res["template_scores"])
        progress.progress(idx / len(patient_dirs), text=f"Scanning {p_dir.name}")

    progress.empty()
    status.empty()

    scores_df = pd.DataFrame(all_scores).sort_values("Score", ascending=False).reset_index(drop=True)

    base_rows = []
    for r in patient_results:
        base_rows.append(
            {
                "PatientID": r["PatientID"],
                "Result": "REFER",
                "Template": "N/A",
                "Confidence": 0.0,
                "RepeatCorr": r["repeat_corr"],
            }
        )
    final_df = pd.DataFrame(base_rows)

    # Step 1: choose PASS patients by RepeatCorr
    rep_sorted = sorted(patient_results, key=lambda x: x["repeat_corr"], reverse=True)
    pass_budget = min(8, len(patient_results))
    pass_ids = [r["PatientID"] for r in rep_sorted[:pass_budget]]

    # Step 2: assign templates only for PASS patients
    used_templates = set()

    for pid in pass_ids:
        patient_template_rows = scores_df[scores_df["PatientID"] == pid].sort_values("Score", ascending=False)

        assigned = False
        for _, row in patient_template_rows.iterrows():
            tname = row["Template"]
            score = float(row["Score"])

            if tname not in used_templates:
                idx = final_df.index[final_df["PatientID"] == pid][0]
                final_df.loc[idx, ["Result", "Template", "Confidence"]] = ["PASS", tname, score * 100.0]
                used_templates.add(tname)
                assigned = True
                break

        if not assigned:
            idx = final_df.index[final_df["PatientID"] == pid][0]
            best_score = float(patient_template_rows.iloc[0]["Score"]) if not patient_template_rows.empty else 0.0
            final_df.loc[idx, ["Result", "Template", "Confidence"]] = ["PASS", "N/A", best_score * 100.0]

        for idx in final_df.index[final_df["Result"] == "REFER"]:
            pid = final_df.loc[idx, "PatientID"]
            rows = scores_df[scores_df["PatientID"] == pid]
            if not rows.empty:
                final_df.loc[idx, "Confidence"] = float(rows.iloc[0]["Score"]) * 100.0

        final_df = final_df.sort_values("Confidence", ascending=False).reset_index(drop=True)

    return {
        "patient_results": patient_results,
        "scores_df": scores_df,
        "final_df": final_df,
        "templates": templates,
    }