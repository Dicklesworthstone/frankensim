"""Subprocess bridge from instrgrade to the Rust fs-psycho descriptor CLI.

Invokes `$FS_PSYCHO_BIN` (`cargo build -p fs-psycho --example psycho_descriptors`)
as a subprocess and parses its single-line `fs-psycho-descriptors/1` JSON output.
Never ports or re-implements `fs-psycho` math in Python, and never fabricates
`0` or `NaN` when the binary is absent or a metric refuses.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from instrgrade import SR, write_wav

PSYCHO_METRICS: tuple[str, ...] = (
    "loudness_sone",
    "n5_sone",
    "n50_sone",
    "sharpness_acum",
    "roughness_asper",
    "tnr_db",
    "pr_db",
    "fluctuation_vacil",
    "log_attack_time",
)

LEVEL_DEPENDENT_METRICS: frozenset[str] = frozenset(
    {
        "loudness_sone",
        "n5_sone",
        "n50_sone",
        "sharpness_acum",
        "roughness_asper",
        "tnr_db",
        "pr_db",
        "fluctuation_vacil",
    }
)


def resolve_psycho_bin(bin_override: str | Path | None = None) -> Path | None:
    """Return the path to the `psycho_descriptors` binary if available, else None."""
    raw = str(bin_override) if bin_override is not None else os.environ.get("FS_PSYCHO_BIN", "").strip()
    if not raw:
        return None
    p = Path(raw)
    if p.is_file() and os.access(p, os.X_OK):
        return p.resolve()
    return None


def psycho_bin_sha256(bin_override: str | Path | None = None) -> str | None:
    """Return the hex SHA-256 of `$FS_PSYCHO_BIN` if present and executable, else None."""
    p = resolve_psycho_bin(bin_override)
    if p is None:
        return None
    return hashlib.sha256(p.read_bytes()).hexdigest()


def missing_binary_block(
    reason: str,
    *,
    pa_per_full_scale: float = 1.0,
    calibration_source: str = "assumed",
    field: str = "free",
    wav_sha256: str | None = None,
    sample_rate_hz: int = SR,
    samples: int = 0,
) -> dict[str, Any]:
    """Return a structured refusal block when `$FS_PSYCHO_BIN` is unset or missing."""
    refusal_msg = f"FS_PSYCHO_BIN_MISSING: {reason}"
    return {
        "schema": "fs-psycho-descriptors/1",
        "wav_sha256": wav_sha256,
        "sample_rate": int(sample_rate_hz),
        "sample_rate_hz": int(sample_rate_hz),
        "samples": int(samples),
        "pa_per_full_scale": float(pa_per_full_scale),
        "calibration_source": calibration_source,
        "calibration": {
            "pa_per_full_scale": float(pa_per_full_scale),
            "source": calibration_source,
            "level_unverified": calibration_source == "assumed",
        },
        "field": field,
        "level_unverified": calibration_source == "assumed",
        "loudness_sone": None,
        "n5_sone": None,
        "n50_sone": None,
        "sharpness_acum": None,
        "roughness_asper": None,
        "tnr_db": None,
        "pr_db": None,
        "fluctuation_vacil": None,
        "log_attack_time": None,
        "refusals": {m: refusal_msg for m in PSYCHO_METRICS},
    }


def detect_candidate_calibration(candidate: dict[str, Any], wav_path: Path | None = None) -> tuple[float, str]:
    """Determine `(pa_per_full_scale, calibration_source)` for a candidate WAV.

    Checks `candidate["pa_per_full_scale"]`, then sidecar JSON files next to `wav_path`
    (`<stem>.calibration.json` or `<wav>.meta.json`). Falls back to `(1.0, "assumed")`.
    """
    if "pa_per_full_scale" in candidate and candidate["pa_per_full_scale"] is not None:
        val = float(candidate["pa_per_full_scale"])
        if math.isfinite(val) and val > 0.0:
            src = str(candidate.get("calibration_source", "sidecar"))
            if src not in ("sidecar", "assumed"):
                src = "sidecar"
            return val, src

    if wav_path is not None:
        for sidecar in (
            wav_path.with_suffix(".calibration.json"),
            wav_path.with_suffix(".meta.json"),
            Path(str(wav_path) + ".meta.json"),
        ):
            if sidecar.is_file():
                try:
                    meta = json.loads(sidecar.read_text())
                    if "pa_per_full_scale" in meta and meta["pa_per_full_scale"] is not None:
                        val = float(meta["pa_per_full_scale"])
                        if math.isfinite(val) and val > 0.0:
                            return val, "sidecar"
                except Exception:  # noqa: BLE001
                    pass

    return 1.0, "assumed"


def evaluate_psycho_wav(
    wav_path: str | Path,
    *,
    pa_per_full_scale: float = 1.0,
    calibration_source: str = "assumed",
    field: str = "free",
    bin_override: str | Path | None = None,
    timeout_s: float = 30.0,
) -> dict[str, Any]:
    """Invoke `$FS_PSYCHO_BIN` on `wav_path` and return the parsed `fs-psycho-descriptors/1` dict."""
    wav_p = Path(wav_path)
    wav_sha = hashlib.sha256(wav_p.read_bytes()).hexdigest() if wav_p.is_file() else None
    bin_p = resolve_psycho_bin(bin_override)
    if bin_p is None:
        raw = str(bin_override) if bin_override is not None else os.environ.get("FS_PSYCHO_BIN", "")
        reason = f"FS_PSYCHO_BIN={raw!r} is unset or not an executable file"
        return missing_binary_block(
            reason,
            pa_per_full_scale=pa_per_full_scale,
            calibration_source=calibration_source,
            field=field,
            wav_sha256=wav_sha,
        )

    cmd = [
        str(bin_p),
        "--wav",
        str(wav_p),
        "--pa-per-full-scale",
        repr(float(pa_per_full_scale)),
        "--calibration-source",
        str(calibration_source),
        "--field",
        str(field),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
    if proc.returncode != 0:
        err_msg = (proc.stderr or proc.stdout or f"exit {proc.returncode}").strip()
        refusal_msg = f"FS_PSYCHO_CLI_ERROR: {err_msg}"
        return {
            "schema": "fs-psycho-descriptors/1",
            "wav_sha256": wav_sha,
            "sample_rate": SR,
            "sample_rate_hz": SR,
            "samples": 0,
            "pa_per_full_scale": float(pa_per_full_scale),
            "calibration_source": calibration_source,
            "calibration": {
                "pa_per_full_scale": float(pa_per_full_scale),
                "source": calibration_source,
                "level_unverified": calibration_source == "assumed",
            },
            "field": field,
            "level_unverified": calibration_source == "assumed",
            **{m: None for m in PSYCHO_METRICS},
            "refusals": {m: refusal_msg for m in PSYCHO_METRICS},
        }

    line = proc.stdout.strip().splitlines()[-1]
    parsed = json.loads(line)
    cal = parsed.get("calibration") or {}
    cal_src = cal.get("source") or parsed.get("calibration_source") or calibration_source
    cal_pa = cal.get("pa_per_full_scale") if "pa_per_full_scale" in cal else parsed.get("pa_per_full_scale", pa_per_full_scale)
    lvl_unv = bool(cal.get("level_unverified", cal_src == "assumed"))
    parsed["pa_per_full_scale"] = float(cal_pa)
    parsed["calibration_source"] = str(cal_src)
    parsed["level_unverified"] = lvl_unv
    return parsed


def evaluate_psycho_signal(
    x: np.ndarray,
    *,
    pa_per_full_scale: float = 1.0,
    calibration_source: str = "assumed",
    field: str = "free",
    bin_override: str | Path | None = None,
    scratch_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Write mono 48 kHz float32 signal `x` to a temporary WAV and evaluate via `$FS_PSYCHO_BIN`."""
    arr = np.asarray(x, dtype=np.float32).reshape(-1)
    bin_p = resolve_psycho_bin(bin_override)
    if bin_p is None:
        raw = str(bin_override) if bin_override is not None else os.environ.get("FS_PSYCHO_BIN", "")
        reason = f"FS_PSYCHO_BIN={raw!r} is unset or not an executable file"
        return missing_binary_block(
            reason,
            pa_per_full_scale=pa_per_full_scale,
            calibration_source=calibration_source,
            field=field,
            sample_rate_hz=SR,
            samples=int(arr.size),
        )

    base_dir = Path(scratch_dir) if scratch_dir is not None else Path(tempfile.gettempdir())
    base_dir.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(arr.tobytes()).hexdigest()[:16]
    wav_path = base_dir / f"instrgrade_psycho_{os.getpid()}_{digest}.wav"
    write_wav(wav_path, arr, sr=SR)
    return evaluate_psycho_wav(
        wav_path,
        pa_per_full_scale=pa_per_full_scale,
        calibration_source=calibration_source,
        field=field,
        bin_override=bin_p,
    )


def compare_psycho_blocks(
    ours: dict[str, Any],
    real: dict[str, Any],
    ruler_blocks: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Compare candidate (`ours`) and reference (`real`) psychoacoustic blocks.

    Marks any loudness/level-dependent comparison as `level_unverified: True` whenever
    either `ours` or `real` has `calibration_source == "assumed"`.
    When `ruler_blocks` (real-vs-real alternate recordings) are provided, reports
    whether `ours` falls inside the `[min, max]` spread of the real recordings.
    """
    ours_assumed = ours.get("calibration_source") == "assumed"
    real_assumed = real.get("calibration_source") == "assumed"
    any_assumed = bool(ours_assumed or real_assumed)

    metrics_cmp: dict[str, Any] = {}
    for m in PSYCHO_METRICS:
        v_o = ours.get(m)
        v_r = real.get(m)
        level_unverified = any_assumed if m in LEVEL_DEPENDENT_METRICS else False
        delta = (
            float(v_o - v_r)
            if (v_o is not None and v_r is not None and math.isfinite(v_o) and math.isfinite(v_r))
            else None
        )
        entry: dict[str, Any] = {
            "ours": v_o,
            "real": v_r,
            "delta": delta,
            "level_unverified": level_unverified,
        }
        if ruler_blocks:
            r_vals = [
                float(rb[m])
                for rb in [real, *ruler_blocks]
                if rb is not None and rb.get(m) is not None and math.isfinite(float(rb[m]))
            ]
            if len(r_vals) >= 2 and v_o is not None and math.isfinite(float(v_o)):
                lo, hi = min(r_vals), max(r_vals)
                entry["real_spread_min"] = lo
                entry["real_spread_max"] = hi
                entry["inside_real_spread"] = bool(lo <= float(v_o) <= hi)
            else:
                entry["real_spread_min"] = None
                entry["real_spread_max"] = None
                entry["inside_real_spread"] = None
        metrics_cmp[m] = entry

    return {
        "level_unverified": any_assumed,
        "ours_calibration_source": ours.get("calibration_source"),
        "real_calibration_source": real.get("calibration_source"),
        "ours": ours,
        "real": real,
        "metrics": metrics_cmp,
    }
