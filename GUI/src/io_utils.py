import json
import tempfile
import zipfile
from pathlib import Path

import numpy as np
from scipy.io import wavfile


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