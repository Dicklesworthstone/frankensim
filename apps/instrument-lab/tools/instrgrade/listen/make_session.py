#!/usr/bin/env python3
"""Build a self-contained offline 2AFC blind listening or spectrogram-pair session (IRP M7.1-W3).

Usage:
  python3 make_session.py --spec SPEC.json --out SESSION_DIR
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
INSTRGRADE_DIR = THIS_DIR.parent
if str(INSTRGRADE_DIR) not in sys.path:
    sys.path.insert(0, str(INSTRGRADE_DIR))
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from equivalence import SessionSpec, build_balanced_trials  # noqa: E402
from instrgrade import SR, highpass, hz_of, load, midi_of, onset, write_wav  # noqa: E402
from psycho import evaluate_psycho_signal, resolve_psycho_bin  # noqa: E402


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rms_first_second(x: np.ndarray) -> float:
    seg = x[: min(len(x), SR)]
    if seg.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(np.square(seg, dtype=np.float64))))


def match_audio_pair(
    real_wav: Path,
    ours_wav: Path,
    out_real: Path,
    out_ours: Path,
    *,
    note: str | None = None,
    matching_cfg: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Apply onset alignment and level matching (loudness or rms-1s) to a WAV pair."""
    cfg = matching_cfg or {"level": "loudness", "onset": True, "pitch": "setup-protocol-only"}
    f0_hint = hz_of(midi_of(note)) if note else None
    x_real = highpass(load(real_wav), f0_hint)
    x_ours = highpass(load(ours_wav), f0_hint)

    onset_aligned = bool(cfg.get("onset", True))
    if onset_aligned:
        o_real = int(onset(x_real))
        o_ours = int(onset(x_ours))
        x_real = x_real[o_real:]
        x_ours = x_ours[o_ours:]
    else:
        o_real = 0
        o_ours = 0

    max_samples = int(4 * SR)
    x_real = np.asarray(x_real[:max_samples], dtype=np.float64)
    x_ours = np.asarray(x_ours[:max_samples], dtype=np.float64)

    requested_level = str(cfg.get("level", "loudness"))
    applied_level = "rms-1s"
    sone_real: float | None = None
    sone_ours: float | None = None

    if requested_level == "loudness" and resolve_psycho_bin() is not None:
        pr = evaluate_psycho_signal(x_real, pa_per_full_scale=1.0, calibration_source="assumed")
        po = evaluate_psycho_signal(x_ours, pa_per_full_scale=1.0, calibration_source="assumed")
        sr_val = pr.get("loudness_sone")
        so_val = po.get("loudness_sone")
        if sr_val is not None and so_val is not None and sr_val > 1e-6 and so_val > 1e-6:
            sone_real = float(sr_val)
            sone_ours = float(so_val)
            # ISO 532-1 sone ~ doubling per 10 dB => amplitude gain ~ sone^(-1 / 0.6)
            g_real = float(sone_real ** (-1.0 / 0.6))
            g_ours = float(sone_ours ** (-1.0 / 0.6))
            applied_level = "loudness"
        else:
            g_real = 1.0 / max(_rms_first_second(x_real), 1e-12)
            g_ours = 1.0 / max(_rms_first_second(x_ours), 1e-12)
    else:
        g_real = 1.0 / max(_rms_first_second(x_real), 1e-12)
        g_ours = 1.0 / max(_rms_first_second(x_ours), 1e-12)

    a = x_real * g_real
    b = x_ours * g_ours
    pk = float(max(np.max(np.abs(a)) if a.size else 1.0, np.max(np.abs(b)) if b.size else 1.0, 1e-12))
    scale = 0.89 / pk
    a = (a * scale).astype(np.float32)
    b = (b * scale).astype(np.float32)

    out_real.parent.mkdir(parents=True, exist_ok=True)
    write_wav(out_real, a, sr=SR)
    write_wav(out_ours, b, sr=SR)

    return {
        "level_method_requested": requested_level,
        "level_method": applied_level,
        "onset_aligned": onset_aligned,
        "onset_samples_real": o_real,
        "onset_samples_ours": o_ours,
        "gain_real": float(g_real * scale),
        "gain_ours": float(g_ours * scale),
        "loudness_sone_real_pre": sone_real,
        "loudness_sone_ours_pre": sone_ours,
        "pitch_policy": str(cfg.get("pitch", "setup-protocol-only")),
        "recording_chain_emulation": str(cfg.get("recording_chain", "dry-matched")),
        "level_unverified": True,
        "source_real_sha256": _sha256_file(real_wav),
        "source_ours_sha256": _sha256_file(ours_wav),
        "matched_real_sha256": _sha256_file(out_real),
        "matched_ours_sha256": _sha256_file(out_ours),
    }


def _render_html(spec: SessionSpec, trials: list[dict[str, Any]], spec_sha256: str) -> str:
    """Render a self-contained static HTML+JS 2AFC session page that works from `file://`."""
    # Redact ground-truth order/real_position from visual DOM while keeping trial metadata for export
    trials_json = json.dumps(trials, indent=2)
    family_esc = html.escape(spec.family)
    modality_esc = html.escape(spec.modality)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>2AFC Equivalence Session — {family_esc} ({modality_esc})</title>
<style>
  body {{ font-family: system-ui, -apple-system, sans-serif; max-width: 860px; margin: 2rem auto; padding: 0 1rem; color: #111; background: #fafafa; }}
  .card {{ background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 1.25rem; margin-bottom: 1rem; }}
  .stimuli {{ display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin: 1rem 0; }}
  button {{ padding: 0.6rem 1.1rem; border-radius: 6px; border: 1px solid #444; background: #fff; cursor: pointer; font-size: 0.95rem; }}
  button.primary {{ background: #111; color: #fff; }}
  img.spec {{ max-width: 100%; border: 1px solid #ccc; border-radius: 4px; }}
  .meta {{ color: #555; font-size: 0.88rem; }}
</style>
</head>
<body>
  <h1>Blind 2AFC Equivalence Session — {family_esc}</h1>
  <p class="meta">Modality: <strong>{modality_esc}</strong> · Seed: <code>{spec.seed}</code> · Spec SHA-256: <code>{spec_sha256[:16]}…</code></p>
  <div class="card">
    <label for="listener-select">Listener ID (pseudonymous): </label>
    <select id="listener-select"></select>
    <span id="progress-label" style="margin-left: 1rem;"></span>
  </div>
  <div id="trial-container" class="card"></div>
  <div class="card">
    <button id="download-btn" class="primary">Download responses JSON</button>
  </div>
<script>
const SPEC_SHA256 = {json.dumps(spec_sha256)};
const ALL_TRIALS = {trials_json};
let responsesByTrial = {{}};

const listeners = [...new Set(ALL_TRIALS.map(t => t.listener_id))];
const selectEl = document.getElementById("listener-select");
listeners.forEach(lid => {{
  const opt = document.createElement("option");
  opt.value = lid;
  opt.textContent = lid;
  selectEl.appendChild(opt);
}});

let currentIdx = 0;
function activeTrials() {{
  const lid = selectEl.value;
  return ALL_TRIALS.filter(t => t.listener_id === lid);
}}

function renderCurrent() {{
  const list = activeTrials();
  const prog = document.getElementById("progress-label");
  const box = document.getElementById("trial-container");
  if (!list.length) {{
    box.textContent = "No trials for listener.";
    return;
  }}
  if (currentIdx >= list.length) currentIdx = list.length - 1;
  const tr = list[currentIdx];
  const answered = Object.keys(responsesByTrial).length;
  prog.textContent = `Trial ${{currentIdx + 1}} of ${{list.length}} (Total recorded: ${{answered}} / ${{ALL_TRIALS.length}})`;
  const mediaA = tr.modality === "spectrogram"
    ? `<img class="spec" src="${{tr.stimulus_a}}" alt="Stimulus A">`
    : `<audio controls preload="auto" src="${{tr.stimulus_a}}"></audio>`;
  const mediaB = tr.modality === "spectrogram"
    ? `<img class="spec" src="${{tr.stimulus_b}}" alt="Stimulus B">`
    : `<audio controls preload="auto" src="${{tr.stimulus_b}}"></audio>`;
  box.innerHTML = `
    <h3>Which one is the real recording? (Item: ${{tr.item_id}})</h3>
    <div class="stimuli">
      <div><h4>Stimulus A</h4>${{mediaA}}<p><button onclick="choose(1)">A is real</button></p></div>
      <div><h4>Stimulus B</h4>${{mediaB}}<p><button onclick="choose(2)">B is real</button></p></div>
    </div>
  `;
}}

window.choose = function(pos) {{
  const list = activeTrials();
  const tr = list[currentIdx];
  responsesByTrial[tr.trial_index] = {{
    ...tr,
    choice_position: pos,
    correct: pos === tr.real_position
  }};
  if (currentIdx + 1 < list.length) {{
    currentIdx += 1;
  }}
  renderCurrent();
}};

selectEl.addEventListener("change", () => {{ currentIdx = 0; renderCurrent(); }});
document.getElementById("download-btn").addEventListener("click", () => {{
  const payload = Object.values(responsesByTrial).sort((a, b) => a.trial_index - b.trial_index);
  const blob = new Blob([JSON.stringify(payload, null, 2)], {{type: "application/json"}});
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = "responses.json";
  a.click();
}});
renderCurrent();
</script>
</body>
</html>
"""


def build_session(
    spec: SessionSpec,
    out_dir: str | Path,
    *,
    base_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Build a self-contained static 2AFC session directory (`index.html`, `stimuli/`, `manifest.jsonl`)."""
    out_p = Path(out_dir)
    stim_dir = out_p / "stimuli"
    stim_dir.mkdir(parents=True, exist_ok=True)
    root_p = Path(base_dir) if base_dir is not None else Path.cwd()

    prepared_items: list[dict[str, Any]] = []
    item_manifests: list[dict[str, Any]] = []

    for item in spec.items:
        iid = str(item["id"])
        real_src = Path(item["real"])
        if not real_src.is_absolute():
            real_src = root_p / real_src
        ours_src = Path(item["ours"])
        if not ours_src.is_absolute():
            ours_src = root_p / ours_src

        modality = str(item.get("modality", spec.modality))
        if modality == "spectrogram":
            ext_r = real_src.suffix or ".png"
            ext_o = ours_src.suffix or ".png"
            rel_real = f"stimuli/{iid}--real{ext_r}"
            rel_ours = f"stimuli/{iid}--ours{ext_o}"
            dst_real = out_p / rel_real
            dst_ours = out_p / rel_ours
            shutil.copyfile(real_src, dst_real)
            shutil.copyfile(ours_src, dst_ours)
            matching_log = {
                "modality": "spectrogram",
                "level_method": "panel-peak-normalized",
                "onset_aligned": True,
                "pitch_policy": str(spec.matching.get("pitch", "setup-protocol-only")),
                "level_unverified": True,
                "source_real_sha256": _sha256_file(real_src),
                "source_ours_sha256": _sha256_file(ours_src),
                "matched_real_sha256": _sha256_file(dst_real),
                "matched_ours_sha256": _sha256_file(dst_ours),
            }
        else:
            rel_real = f"stimuli/{iid}--real.wav"
            rel_ours = f"stimuli/{iid}--ours.wav"
            dst_real = out_p / rel_real
            dst_ours = out_p / rel_ours
            matching_log = match_audio_pair(
                real_src,
                ours_src,
                dst_real,
                dst_ours,
                note=item.get("note"),
                matching_cfg=spec.matching,
            )

        prepared_items.append(
            {
                **item,
                "real": rel_real,
                "ours": rel_ours,
                "modality": modality,
                "matching_log": matching_log,
            }
        )
        item_manifests.append(
            {
                "id": iid,
                "modality": modality,
                "real": rel_real,
                "ours": rel_ours,
                "matching_log": matching_log,
            }
        )

    prepared_spec = SessionSpec(
        family=spec.family,
        items=prepared_items,
        listeners_planned=spec.listeners_planned,
        items_per_listener=spec.items_per_listener,
        design=spec.design,
        band=spec.band,
        alpha=spec.alpha,
        seed=spec.seed,
        matching=spec.matching,
        modality=spec.modality,
    )
    spec_bytes = json.dumps(prepared_spec.to_dict(), sort_keys=True).encode("utf-8")
    spec_sha256 = _sha256_bytes(spec_bytes)

    trials = build_balanced_trials(prepared_spec)
    html_text = _render_html(prepared_spec, trials, spec_sha256)
    (out_p / "index.html").write_text(html_text, encoding="utf-8")
    (out_p / "trials.json").write_text(json.dumps(trials, indent=2) + "\n", encoding="utf-8")

    git_proc = subprocess.run(
        ["git", "-C", str(THIS_DIR), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        timeout=10,
    )
    manifest_entry = {
        "event": "session_built",
        "schema": "instrgrade-listening-session/1",
        "spec_sha256": spec_sha256,
        "family": prepared_spec.family,
        "modality": prepared_spec.modality,
        "seed": prepared_spec.seed,
        "design": prepared_spec.design,
        "band": list(prepared_spec.band),
        "alpha": prepared_spec.alpha,
        "listeners_planned": prepared_spec.listeners_planned,
        "items_per_listener": prepared_spec.items_per_listener,
        "n_trials": len(trials),
        "items": item_manifests,
        "host": platform.node(),
        "git_sha": git_proc.stdout.strip() or None,
        "command": [sys.executable, *sys.argv],
        "utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path = out_p / "manifest.jsonl"
    with manifest_path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(manifest_entry) + "\n")

    return {
        "spec": prepared_spec,
        "spec_sha256": spec_sha256,
        "trials": trials,
        "manifest_path": manifest_path,
        "index_html_path": out_p / "index.html",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build an offline 2AFC equivalence session.")
    parser.add_argument("--spec", required=True, help="Path to SessionSpec JSON")
    parser.add_argument("--out", required=True, help="Output session directory")
    args = parser.parse_args(argv)

    spec_path = Path(args.spec)
    spec_data = json.loads(spec_path.read_text(encoding="utf-8"))
    spec = SessionSpec.from_dict(spec_data)
    build_session(spec, Path(args.out), base_dir=spec_path.parent)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
