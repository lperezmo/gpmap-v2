"""JSON / CSV / pickle round-trips."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import streamlit as st
from gpmap import read_csv, read_json, read_pickle, to_csv, to_json, to_pickle
from gpmap.simulate import MountFujiSimulation

st.markdown(
    "`gpmap-v2` round-trips cleanly through JSON, CSV + sidecar metadata, "
    "pickle, and Excel. JSON files carry `schema_version` and CSV writes a "
    "matching `<basename>.meta.json` so the wildtype and mutations dict "
    "survive."
)

sim = MountFujiSimulation(
    wildtype="AAAA",
    mutations={i: ["A", "T"] for i in range(4)},
    field_strength=1.0,
    roughness_width=0.1,
    rng=np.random.default_rng(0),
)

format_choice = st.radio("Format", options=["JSON", "CSV + sidecar", "pickle"], horizontal=True)

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)

    if format_choice == "JSON":
        path = root / "fuji.json"
        to_json(sim, path)
        text = path.read_text(encoding="utf-8")
        st.subheader(f"{path.name}  ({len(text):,} bytes)")
        st.code(json.dumps(json.loads(text), indent=2)[:1200] + "\n...", language="json")
        restored = read_json(path)

    elif format_choice == "CSV + sidecar":
        path = root / "fuji.csv"
        to_csv(sim, path)
        st.subheader(path.name)
        st.code(path.read_text(encoding="utf-8")[:600] + "\n...", language="text")
        meta_path = path.with_suffix(path.suffix + ".meta.json")
        st.subheader(meta_path.name)
        st.code(meta_path.read_text(encoding="utf-8"), language="json")
        restored = read_csv(path)

    else:
        path = root / "fuji.pkl"
        to_pickle(sim, path)
        st.subheader(f"{path.name}  ({path.stat().st_size:,} bytes)")
        st.info("Pickle is opaque bytes. Round-tripping is what matters.")
        restored = read_pickle(path)

    st.subheader("Round-trip check")
    c1, c2, c3 = st.columns(3)
    c1.metric("wildtype match", str(restored.wildtype == sim.wildtype))
    c2.metric("n match", str(restored.n == sim.n))
    c3.metric(
        "phenotypes allclose",
        str(bool(np.allclose(restored.phenotypes, sim.phenotypes))),
    )

with st.expander("Code", icon=":material/code:"):
    st.code(
        """from gpmap import to_json, read_json, to_csv, read_csv, to_pickle, read_pickle

to_json(gpm, "map.json")
gpm2 = read_json("map.json")

to_csv(gpm, "map.csv")        # also writes map.csv.meta.json
gpm3 = read_csv("map.csv")

to_pickle(gpm, "map.pkl")
gpm4 = read_pickle("map.pkl")
""",
        language="python",
    )
