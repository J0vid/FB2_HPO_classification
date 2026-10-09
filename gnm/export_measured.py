"""Export the measured face region of every GNM subject for the classifier.

Writes float32 C-order [n, n_measured, 3] (mm) + meta CSV with QC status.
Only provenance == measured vertices are used: elsewhere the GNM head is
model geometry and carries no subject information.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REF = Path("/Users/jovid/Documents/Hallgrimsson/gnm_reference")
OUT = Path("/Users/jovid/Documents/Hallgrimsson/gnm_classifier")
sys.path.insert(0, str(REF / "code"))
from gnm_compact import load_subjects  # noqa: E402

qc = pd.read_csv(REF / "qc.csv")
prov = np.load(REF / "subjects_gnm_syndromic_params.npz")["provenance"]
measured = np.flatnonzero(prov == 2)
np.save(OUT / "measured_vertex_ids.npy", measured)
metas = []
with open(OUT / "gnm_measured_f32.bin", "wb") as fh:
    for group in ("nonsyndromic", "syndromic"):
        idx = pd.read_csv(REF / f"subjects_gnm_{group}_index.csv")
        for start in range(0, len(idx), 400):
            rows = np.arange(start, min(start + 400, len(idx)))
            v = load_subjects(group, rows=rows)[:, measured, :].astype(np.float32)
            fh.write(np.ascontiguousarray(v).tobytes())
        q = qc[qc["group"] == group].set_index("row").loc[idx["row"], ["qc_status"]].reset_index(drop=True)
        metas.append(pd.concat([idx.reset_index(drop=True), q], axis=1))
meta = pd.concat(metas, ignore_index=True)
meta.to_csv(OUT / "gnm_meta.csv", index=False)
print(len(meta), "subjects,", len(measured), "measured vertices")
