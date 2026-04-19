"""CSV / JSON / pickle / Excel round-trip for GenotypePhenotypeMap."""

from __future__ import annotations

import json
import pickle
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

if TYPE_CHECKING:
    from .core import GenotypePhenotypeMap

SCHEMA_VERSION = "1"


# --- helpers ----------------------------------------------------------------
def _mutations_to_json(mutations: dict[int, list[str] | None]) -> dict[str, list[str] | None]:
    return {str(k): v for k, v in mutations.items()}


def _mutations_from_json(obj: dict[str, list[str] | None]) -> dict[int, list[str] | None]:
    return {int(k): v for k, v in obj.items()}


# --- CSV --------------------------------------------------------------------
def read_csv(path: str | Path, **kwargs: Any) -> GenotypePhenotypeMap:
    from .core import GenotypePhenotypeMap

    p = Path(path)
    meta_path = p.with_suffix(p.suffix + ".meta.json")
    if not meta_path.exists():
        meta_path = p.with_suffix(".meta.json")
    if not meta_path.exists():
        raise FileNotFoundError(
            f"CSV metadata sidecar not found next to {p}; expected <basename>.meta.json"
        )
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    df = pd.read_csv(p, **kwargs)
    return GenotypePhenotypeMap(
        wildtype=meta["wildtype"],
        genotypes=df["genotypes"].tolist(),
        phenotypes=df["phenotypes"].to_numpy(),
        stdeviations=(df["stdeviations"].to_numpy() if "stdeviations" in df.columns else None),
        n_replicates=(df["n_replicates"].to_numpy() if "n_replicates" in df.columns else None),
        mutations=_mutations_from_json(meta["mutations"]),
        site_labels=meta.get("site_labels"),
    )


def to_csv(gpm: GenotypePhenotypeMap, path: str | Path, **kwargs: Any) -> None:
    p = Path(path)
    df = pd.DataFrame(
        {
            "genotypes": gpm.genotypes,
            "phenotypes": gpm.phenotypes,
            "stdeviations": gpm.stdeviations,
            "n_replicates": gpm.n_replicates,
        }
    )
    kwargs.setdefault("index", False)
    df.to_csv(p, **kwargs)
    meta = {
        "schema_version": SCHEMA_VERSION,
        "wildtype": gpm.wildtype,
        "mutations": _mutations_to_json(gpm.mutations),
        "site_labels": gpm.site_labels,
    }
    meta_path = p.with_suffix(p.suffix + ".meta.json")
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")


# --- JSON -------------------------------------------------------------------
def read_json(path: str | Path) -> GenotypePhenotypeMap:
    from .core import GenotypePhenotypeMap

    obj = json.loads(Path(path).read_text(encoding="utf-8"))
    sv = obj.get("schema_version")
    if sv is None:
        import warnings

        warnings.warn(
            "JSON file has no schema_version; treating as v1 legacy format",
            UserWarning,
            stacklevel=2,
        )
    return GenotypePhenotypeMap(
        wildtype=obj["wildtype"],
        genotypes=obj["genotypes"],
        phenotypes=obj["phenotypes"],
        stdeviations=obj.get("stdeviations"),
        n_replicates=obj.get("n_replicates"),
        mutations=_mutations_from_json(obj["mutations"]),
        site_labels=obj.get("site_labels"),
    )


def to_json(gpm: GenotypePhenotypeMap, path: str | Path) -> None:
    obj = {
        "schema_version": SCHEMA_VERSION,
        "wildtype": gpm.wildtype,
        "mutations": _mutations_to_json(gpm.mutations),
        "site_labels": gpm.site_labels,
        "genotypes": gpm.genotypes.tolist(),
        "phenotypes": gpm.phenotypes.tolist(),
        "stdeviations": gpm.stdeviations.tolist(),
        "n_replicates": gpm.n_replicates.tolist(),
    }
    Path(path).write_text(json.dumps(obj, indent=2), encoding="utf-8")


# --- pickle -----------------------------------------------------------------
def read_pickle(path: str | Path) -> GenotypePhenotypeMap:
    with open(path, "rb") as f:
        obj = pickle.load(f)
    from .core import GenotypePhenotypeMap

    if not isinstance(obj, GenotypePhenotypeMap):
        raise TypeError(f"unpickled object is {type(obj).__name__}, expected GenotypePhenotypeMap")
    return obj


def to_pickle(gpm: GenotypePhenotypeMap, path: str | Path) -> None:
    with open(path, "wb") as f:
        pickle.dump(gpm, f)


# --- Excel ------------------------------------------------------------------
def read_excel(path: str | Path, **kwargs: Any) -> GenotypePhenotypeMap:
    from .core import GenotypePhenotypeMap

    df = pd.read_excel(path, sheet_name="data", **kwargs)
    meta_df = pd.read_excel(path, sheet_name="meta")
    meta = {str(row["key"]): row["value"] for _, row in meta_df.iterrows()}
    return GenotypePhenotypeMap(
        wildtype=str(meta["wildtype"]),
        genotypes=df["genotypes"].tolist(),
        phenotypes=df["phenotypes"].to_numpy(),
        stdeviations=(df["stdeviations"].to_numpy() if "stdeviations" in df.columns else None),
        n_replicates=(df["n_replicates"].to_numpy() if "n_replicates" in df.columns else None),
        mutations=_mutations_from_json(json.loads(str(meta["mutations"]))),
        site_labels=json.loads(str(meta["site_labels"])),
    )


def to_excel(gpm: GenotypePhenotypeMap, path: str | Path, **kwargs: Any) -> None:
    df = pd.DataFrame(
        {
            "genotypes": gpm.genotypes,
            "phenotypes": gpm.phenotypes,
            "stdeviations": gpm.stdeviations,
            "n_replicates": gpm.n_replicates,
        }
    )
    meta_df = pd.DataFrame(
        [
            {"key": "schema_version", "value": SCHEMA_VERSION},
            {"key": "wildtype", "value": gpm.wildtype},
            {"key": "mutations", "value": json.dumps(_mutations_to_json(gpm.mutations))},
            {"key": "site_labels", "value": json.dumps(gpm.site_labels)},
        ]
    )
    with pd.ExcelWriter(path) as writer:
        df.to_excel(writer, sheet_name="data", index=False)
        meta_df.to_excel(writer, sheet_name="meta", index=False)
