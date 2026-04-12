import json
import tempfile
import zipfile
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.io import wavfile, loadmat

def load_matlab_export(export_dir):
    export_dir = Path(export_dir)

    final_df = pd.read_csv(export_dir / "final_mapping.csv")
    scores_df = pd.read_csv(export_dir / "template_scores.csv")

    mat = loadmat(export_dir / "patient_results.mat", squeeze_me=True, struct_as_record=False)
    raw = mat["patient_results"]

    if raw.ndim == 0:
        raw = [raw]
    else:
        raw = list(raw)

    patient_results = []
    for r in raw:
        patient_results.append({
            "PatientID": str(r.PatientID),
            "fs": float(r.fs),
            "t_ms": np.asarray(r.t_ms, dtype=float).ravel(),
            "t_crop_ms": np.asarray(r.t_crop_ms, dtype=float).ravel(),
            "oae_clean": np.asarray(r.oae_clean, dtype=float).ravel(),
            "oae_crop": np.asarray(r.oae_crop, dtype=float).ravel(),
            "est_norm": np.asarray(r.est_norm, dtype=float).ravel(),
            "matched_norm": np.asarray(r.matched_norm, dtype=float).ravel(),
            "fft_freq": np.asarray(r.fft_freq, dtype=float).ravel(),
            "fft_mag": np.asarray(r.fft_mag, dtype=float).ravel(),
            "template_fft_freq": np.asarray(r.template_fft_freq, dtype=float).ravel(),
            "template_fft_mag": np.asarray(r.template_fft_mag, dtype=float).ravel(),
            "result": str(r.result),
            "assigned_template": str(r.assigned_template),
            "assigned_score": float(r.assigned_score),
            "repeat_corr": float(r.repeat_corr),
        })

    summary = json.loads((export_dir / "summary.json").read_text(encoding="utf-8"))

    return {
        "patient_results": patient_results,
        "scores_df": scores_df,
        "final_df": final_df,
        "summary": summary,
    }


def read_wav_any(path_or_bytes):
    if isinstance(path_or_bytes, (str, Path)):
        fs, data = wavfile.read(str(path_or_bytes))
    else:
        fs, data = wavfile.read(path_or_bytes)

    data = np.asarray(data)
    if data.ndim > 1:
        data = data[:, 0]

    if np.issubdtype(data.dtype, np.integer):
        maxv = np.iinfo(data.dtype).max
        data = data.astype(np.float64) / maxv
    else:
        data = data.astype(np.float64)

    return fs, data


def load_json_any(path_or_bytes):
    if isinstance(path_or_bytes, (str, Path)):
        with open(path_or_bytes, "r", encoding="utf-8") as f:
            return json.load(f)
    return json.loads(path_or_bytes.read().decode("utf-8"))


def parse_patient_structure(root_dir):
    root = Path(root_dir)
    if not root.exists():
        raise FileNotFoundError(f"Folder not found: {root_dir}")

    template_path = root / "lostOaes.json"
    if not template_path.exists():
        raise FileNotFoundError("lostOaes.json was not found in the selected base folder.")

    patient_dirs = sorted([p for p in root.glob("patient_*") if p.is_dir()])
    if not patient_dirs:
        raise FileNotFoundError("No patient_* folders were found in the selected base folder.")

    templates = load_json_any(template_path)
    return root, templates, patient_dirs


def unpack_zip_to_temp(uploaded_zip):
    temp_dir = tempfile.TemporaryDirectory()
    with zipfile.ZipFile(uploaded_zip, "r") as zf:
        zf.extractall(temp_dir.name)

    root = Path(temp_dir.name)
    children = [p for p in root.iterdir()]
    if len(children) == 1 and children[0].is_dir():
        root = children[0]

    return temp_dir, root