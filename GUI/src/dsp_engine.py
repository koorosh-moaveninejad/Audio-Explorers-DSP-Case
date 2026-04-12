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

    if "centerFreq" in info:
        target_freqs = np.asarray(info["centerFreq"], dtype=np.float64).ravel()
    else:
        target_freqs = np.array([1000.0, 2000.0, 3000.0, 4000.0], dtype=np.float64)

    epoch_data = []
    ref = None

    for tname in ["A", "B", "C", "D"]:
        active_path = patient_dir / f"Patient_{p_id}_active{tname}.wav"
        if not active_path.exists():
            raise FileNotFoundError(f"Missing active{tname} file for patient {p_id}")

        fs_active, act = read_wav_any(active_path)
        if fs_active != fs:
            raise ValueError(f"Sample-rate mismatch in patient {p_id}, active{tname}")

        starts = np.where(np.diff(np.concatenate([[0], (act > 0.5).astype(int)])) == 1)[0]
        tmp = []

        for s in starts:
            e = s + epoch_size
            if e <= len(rec):
                seg = rec[s:e].astype(np.float64)
                seg = seg - np.mean(seg)

                if ref is None:
                    ref = seg[: min(100, len(seg))].copy()
                    shifted = seg
                else:
                    seg_head = seg[: min(100, len(seg))]
                    ref_head = ref[: min(100, len(ref))]
                    corr = correlate(seg_head, ref_head, mode="full")
                    denom = np.linalg.norm(seg_head) * np.linalg.norm(ref_head)
                    if denom > 0:
                        corr = corr / denom
                    lags = np.arange(-len(ref_head) + 1, len(seg_head))
                    keep = (lags >= -20) & (lags <= 20)
                    corr = corr[keep]
                    lags = lags[keep]
                    best_lag = int(lags[np.argmax(corr)]) if len(corr) else 0
                    shifted = np.roll(seg, -best_lag)

                tmp.append(shifted)

        if len(tmp) == 0:
            epoch_data.append(np.zeros((epoch_size, 0), dtype=np.float64))
        else:
            epoch_data.append(np.column_stack(tmp))

    min_e = min(arr.shape[1] for arr in epoch_data) if epoch_data else 0

    if min_e == 0:
        zero = np.zeros(epoch_size, dtype=np.float64)
        template_scores = [{"PatientID": p_id, "Template": name, "Score": 0.0} for name in templates.keys()]
        return {
            "PatientID": p_id,
            "fs": fs,
            "t_ms": (np.arange(epoch_size) / fs) * 1000.0,
            "t_crop_ms": np.array([], dtype=np.float64),
            "oae_clean": zero,
            "oae_crop": np.array([], dtype=np.float64),
            "est_norm": np.array([], dtype=np.float64),
            "matched_norm": np.array([], dtype=np.float64),
            "best_template": "N/A",
            "best_score": 0.0,
            "repeat_corr": 0.0,
            "template_scores": template_scores,
            "fft_freq": np.array([0.0]),
            "fft_mag": np.array([0.0]),
            "template_fft_freq": np.array([0.0]),
            "template_fft_mag": np.array([0.0]),
        }

    acc = np.zeros(epoch_size, dtype=np.float64)
    acc1 = np.zeros(epoch_size, dtype=np.float64)
    acc2 = np.zeros(epoch_size, dtype=np.float64)
    v_sets = 0

    for k in range(min_e):
        quad = epoch_data[0][:, k] + epoch_data[1][:, k] + epoch_data[2][:, k] + epoch_data[3][:, k]
        if np.all(np.isfinite(quad)):
            acc += quad
            v_sets += 1
            if (k % 2) == 0:
                acc1 += quad
            else:
                acc2 += quad

    t = np.arange(epoch_size) / fs
    win_env = np.exp(-((t - 0.008) ** 2) / (2 * 0.004**2))

    oae_clean = bandpass_filter((acc / max(v_sets, 1)) * win_env, fs)
    oae_clean1 = bandpass_filter((acc1 / max(1, int(np.floor(v_sets / 2)))) * win_env, fs)
    oae_clean2 = bandpass_filter((acc2 / max(1, int(np.floor(v_sets / 2)))) * win_env, fs)

    noise_est = (
        np.mean(epoch_data[0], axis=1)
        + np.mean(epoch_data[1], axis=1)
        - np.mean(epoch_data[2], axis=1)
        - np.mean(epoch_data[3], axis=1)
    )
    oae_clean = oae_clean - 0.1 * bandpass_filter(noise_est, fs)

    resp_idx = (t > DEFAULT_RESPONSE_START) & (t < DEFAULT_RESPONSE_END)
    oae_crop = oae_clean[resp_idx]

    x1 = oae_clean1[resp_idx]
    x2 = oae_clean2[resp_idx]
    repeat_corr = max(0.0, safe_corrcoef(x1, x2))

    if len(oae_crop) > 1:
        nfft = len(oae_crop)
        y_fft = np.abs(np.fft.fft(oae_crop) / nfft)
        p1 = y_fft[: nfft // 2 + 1]
        p1 = p1 / (np.max(y_fft) + np.finfo(float).eps)
        f_axis = fs * np.arange(0, nfft // 2 + 1) / nfft
    else:
        p1 = np.array([0.0])
        f_axis = np.array([0.0])

    template_scores = []
    best_template_name = "N/A"
    best_score = 0.0
    best_target_crop = np.zeros(np.sum(resp_idx), dtype=np.float64)
    best_f_axis_t = np.array([0.0])
    best_p1_t = np.array([0.0])

    for name, target in templates.items():
        target = np.asarray(target, dtype=np.float64).ravel()

        if len(target) < len(oae_clean):
            target_adj = np.pad(target, (0, len(oae_clean) - len(target)))
        else:
            target_adj = target[: len(oae_clean)]

        current_target_crop = target_adj[resp_idx]

        if len(oae_crop) == 0 or len(current_target_crop) == 0 or np.std(oae_crop) == 0 or np.std(current_target_crop) == 0:
            time_score = 0.0
        else:
            corr = correlate(oae_crop, current_target_crop, mode="full")
            denom = np.linalg.norm(oae_crop) * np.linalg.norm(current_target_crop)
            time_score = float(np.max(corr / denom)) if denom > 0 else 0.0
            if not np.isfinite(time_score):
                time_score = 0.0

        if len(current_target_crop) > 1:
            y_fft_t = np.abs(np.fft.fft(current_target_crop) / len(current_target_crop))
            p1_t = y_fft_t[: len(current_target_crop) // 2 + 1]
            p1_t = p1_t / (np.max(y_fft_t) + np.finfo(float).eps)
            f_axis_t = fs * np.arange(0, len(current_target_crop) // 2 + 1) / len(current_target_crop)
        else:
            p1_t = np.array([0.0])
            f_axis_t = np.array([0.0])

        spectral_fit = 0.0
        for f_target in target_freqs:
            bin_idx = int(np.argmin(np.abs(f_axis - f_target))) if len(f_axis) else 0
            spectral_fit += float(p1[bin_idx] * p1_t[bin_idx])
        spectral_fit /= max(len(target_freqs), 1)

        score = 0.7 * max(0.0, time_score) + 0.3 * spectral_fit
        template_scores.append({"PatientID": p_id, "Template": name, "Score": score})

        if score > best_score:
            best_score = score
            best_template_name = name
            best_target_crop = current_target_crop.copy()
            best_f_axis_t = f_axis_t
            best_p1_t = p1_t

    est_norm = oae_crop / (np.max(np.abs(oae_crop)) + np.finfo(float).eps) if len(oae_crop) else oae_crop.copy()
    target_norm = best_target_crop / (np.max(np.abs(best_target_crop)) + np.finfo(float).eps) if len(best_target_crop) else best_target_crop.copy()

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
        "fft_freq": f_axis,
        "fft_mag": p1,
        "template_fft_freq": best_f_axis_t,
        "template_fft_mag": best_p1_t,
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

    temp_names = list(templates.keys())
    scores_df = pd.DataFrame(all_scores)

    # MATLAB final logic: top 8 by RepeatCorr are PASS
    rep_sorted = sorted(patient_results, key=lambda x: x["repeat_corr"], reverse=True)
    pass_budget = min(8, len(patient_results))
    pass_results = rep_sorted[:pass_budget]

    final_df = pd.DataFrame({
        "PatientID": [r["PatientID"] for r in patient_results],
        "Result": ["REFER"] * len(patient_results),
        "Template": ["N/A"] * len(patient_results),
        "Confidence": [0.0] * len(patient_results),
        "RepeatCorr": [r["repeat_corr"] for r in patient_results],
    })

    # Build cost matrix [8 PASS x 8 templates]
    cost_matrix = np.zeros((pass_budget, len(temp_names)), dtype=np.float64)
    for i, r in enumerate(pass_results):
        pid = r["PatientID"]
        patient_scores = scores_df[scores_df["PatientID"] == pid].copy()
        patient_scores = patient_scores.set_index("Template").reindex(temp_names).reset_index()
        cost_matrix[i, :] = patient_scores["Score"].fillna(0.0).to_numpy(dtype=np.float64)

    # MATLAB final logic: greedy global lockout assignment
    temp_cost = cost_matrix.copy()
    for _ in range(min(pass_budget, len(temp_names))):
        lin_idx = int(np.argmax(temp_cost))
        best_val = float(temp_cost.flat[lin_idx])

        if not np.isfinite(best_val):
            break

        r_idx, c_idx = np.unravel_index(lin_idx, temp_cost.shape)
        pid = pass_results[r_idx]["PatientID"]

        row_idx = final_df.index[final_df["PatientID"] == pid][0]
        final_df.loc[row_idx, ["Result", "Template", "Confidence"]] = [
            "PASS",
            temp_names[c_idx],
            best_val * 100.0,
        ]

        temp_cost[r_idx, :] = -np.inf
        temp_cost[:, c_idx] = -np.inf

    # REFER rows remain N/A and 0
    final_df = final_df.sort_values("RepeatCorr", ascending=False).reset_index(drop=True)

    # Sync final_df into patient_results so GUI matches Overview
    final_map = {
        row["PatientID"]: {
            "result": row["Result"],
            "assigned_template": row["Template"],
            "assigned_score": float(row["Confidence"]),
        }
        for _, row in final_df.iterrows()
    }

    for r in patient_results:
        pid = r["PatientID"]
        if pid in final_map:
            r["result"] = final_map[pid]["result"]
            r["assigned_template"] = final_map[pid]["assigned_template"]
            r["assigned_score"] = final_map[pid]["assigned_score"]
        else:
            r["result"] = "N/A"
            r["assigned_template"] = "N/A"
            r["assigned_score"] = 0.0

        # REFER patients should not show matched template
        if r["result"] == "REFER":
            r["matched_norm"] = np.zeros_like(r["est_norm"])
            r["template_fft_freq"] = np.array([0.0])
            r["template_fft_mag"] = np.array([0.0])
        else:
            assigned_template = r["assigned_template"]
            if assigned_template in templates:
                target = np.asarray(templates[assigned_template], dtype=np.float64).ravel()

                # rebuild matched waveform from assigned template
                t_sec = r["t_ms"] / 1000.0
                resp_idx = (t_sec > DEFAULT_RESPONSE_START) & (t_sec < DEFAULT_RESPONSE_END)

                if len(target) < len(r["oae_clean"]):
                    target_adj = np.pad(target, (0, len(r["oae_clean"]) - len(target)))
                else:
                    target_adj = target[: len(r["oae_clean"])]

                target_crop = target_adj[resp_idx]
                target_norm = target_crop / (np.max(np.abs(target_crop)) + np.finfo(float).eps)

                if len(r["est_norm"]) > 0 and len(target_norm) > 0:
                    corr = correlate(r["est_norm"], target_norm, mode="full")
                    lags = np.arange(-len(target_norm) + 1, len(r["est_norm"]))
                    best_lag = int(lags[np.argmax(corr)]) if len(corr) else 0
                    r["matched_norm"] = np.roll(target_norm, best_lag)
                else:
                    r["matched_norm"] = np.zeros_like(r["est_norm"])

                # rebuild FFT of assigned template
                if len(target_crop) > 1:
                    nfft_t = len(target_crop)
                    y_fft_t = np.abs(np.fft.fft(target_crop) / nfft_t)
                    p1_t = y_fft_t[: nfft_t // 2 + 1]
                    p1_t = p1_t / (np.max(p1_t) + np.finfo(float).eps)
                    f_axis_t = r["fs"] * np.arange(0, nfft_t // 2 + 1) / nfft_t
                else:
                    p1_t = np.array([0.0])
                    f_axis_t = np.array([0.0])

                r["template_fft_freq"] = f_axis_t
                r["template_fft_mag"] = p1_t

    return {
        "patient_results": patient_results,
        "scores_df": scores_df,
        "final_df": final_df,
        "templates": templates,
    }
