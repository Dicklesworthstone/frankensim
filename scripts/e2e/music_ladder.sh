#!/usr/bin/env bash
# Thin retained-evidence wrapper around a natively built convergence_ladder
# test binary. Build/host selection belongs to the central verification lane.
set -euo pipefail
OUT_DIR=""
BINARY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) [[ $# -ge 2 ]] || exit 2; OUT_DIR="$2"; shift 2 ;;
    --binary) [[ $# -ge 2 ]] || exit 2; BINARY="$2"; shift 2 ;;
    -h|--help)
      printf '%s\n' 'usage: music_ladder.sh --binary <convergence_ladder-test-executable> --out-dir <fresh-path>'
      exit 0 ;;
    *) printf 'unknown argument: %s\n' "$1" >&2; exit 2 ;;
  esac
done
[[ -n "$OUT_DIR" && -n "$BINARY" && -x "$BINARY" ]] || {
  printf '%s\n' 'an executable native test binary and a fresh output directory are required' >&2; exit 2;
}
[[ ! -e "$OUT_DIR" && ! -L "$OUT_DIR" ]] || {
  printf 'refusing existing output directory: %s\n' "$OUT_DIR" >&2; exit 2;
}
mkdir -- "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
export FRANKENSIM_LADDER_OUT_DIR="$OUT_DIR"
export FRANKENSIM_LADDER_GIT_SHA="$(git rev-parse HEAD)"
export FRANKENSIM_LADDER_HOST="$(uname -n)"
printf -v FRANKENSIM_LADDER_COMMAND '%q --test-threads=1 --nocapture' "$BINARY"
export FRANKENSIM_LADDER_COMMAND
printf '%s\n' "$FRANKENSIM_LADDER_COMMAND" > "$OUT_DIR/command.txt"
"$BINARY" --test-threads=1 --nocapture 2> "$OUT_DIR/stderr.log" \
  | tee "$OUT_DIR/stdout.log" \
  | awk 'index($0,"IRPLADDER ") { print substr($0,index($0,"IRPLADDER ")+10); fflush(); }' \
  | tee "$OUT_DIR/receipts.jsonl"
python3 -I - "$OUT_DIR" <<'PY'
import json
import pathlib
import sys
root = pathlib.Path(sys.argv[1])
rows = [json.loads(line) for line in (root / "receipts.jsonl").read_text().splitlines()]
terminals = [row for row in rows if row.get("stage") == "terminal"]
if len(terminals) != 1 or terminals[0].get("verdict") != "measured-no-realism-claim":
    raise SystemExit("missing unique successful joint-ladder terminal receipt")
for axis in ("dt", "fft", "tolerance"):
    live = {row["realized"] for row in rows if row.get("stage") == "rung" and row.get("axis") == axis and row.get("live")}
    if len(live) < 3:
        raise SystemExit(f"fewer than three realized live rungs for {axis}")
for path in ("reed.wav", "reed.provenance.json", "reed-ladder.report"):
    if not (root / path).is_file():
        raise SystemExit(f"missing retained render/report artifact: {path}")
PY
