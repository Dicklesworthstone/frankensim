#!/usr/bin/env bash
# Bowed-string playable-region sweep + listening receipt
# (bead frankensim-music-v8-root-3ez8g.7.5).
#
# Runs the Schelleng-style (F_n, v_bow) regime sweep through the release
# build and renders the open-string stroke listening receipt. Both artifacts
# are deterministic on one host; the sweep JSONL is the retained receipt.
#
# Usage: scripts/e2e/music/bowed_string_sweep.sh [OUT_DIR]
set -u -o pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
OUT_DIR="${1:-$REPO_ROOT/.e2e-out/bowed-sweep}"
mkdir -p "$OUT_DIR"

# Retain execution context separately from deterministic numerical sweep rows.
# Run the pre-change and corrected trees into distinct directories in the
# central verification lane; never overwrite one with the other.
STARTED="$(python3 -c 'import time; print(time.monotonic())')"
record_execution() {
  local result="$?"
  python3 - "$OUT_DIR" "$REPO_ROOT" "$STARTED" "$result" <<'PY'
import json
import pathlib
import platform
import subprocess
import sys
import time

out, repo, started, code = sys.argv[1:]
head = subprocess.check_output(["git", "-C", repo, "rev-parse", "HEAD"], text=True).strip()
row = {"schema": "frankensim.bowed-sweep.execution.v1", "kind": "execution",
       "host": platform.node(), "git_sha": head,
       "command": ["scripts/e2e/music/bowed_string_sweep.sh", out],
       "elapsed_s": time.monotonic() - float(started), "exit_code": int(code),
       "sweep": "sweep.jsonl", "listening_receipt": "listening-receipt.jsonl",
       "rebaseline_cause": "DEF-06: ODE ramp joins kinetic endpoint; static capacity unchanged"}
with (pathlib.Path(out) / "execution.jsonl").open("a", encoding="utf-8") as sink:
    sink.write(json.dumps(row, sort_keys=True) + "\n")
PY
  return "$result"
}
trap record_execution EXIT

cargo run -r -q -p fs-couple --bin bowed_sweep -- \
  --out "$OUT_DIR/sweep.jsonl" || exit $?

cargo run -r -q -p fs-couple --bin bowed_listen -- \
  --wav "$OUT_DIR/open_string_stroke.wav" \
  --receipt "$OUT_DIR/listening-receipt.jsonl" || exit $?

echo "artifacts: $OUT_DIR/sweep.jsonl $OUT_DIR/open_string_stroke.wav $OUT_DIR/listening-receipt.jsonl"
