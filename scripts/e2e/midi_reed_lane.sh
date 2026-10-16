#!/usr/bin/env bash
# Called by music_lane.sh --fixture midi-reed. Uses already-built Rust tools;
# never invokes Cargo. Software replay only, not reed realism qualification.
set -euo pipefail
OUT_DIR="" BINARY="" DRY_RUN=0 ENFORCE_BASELINE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --binary) BINARY="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --enforce-baseline) ENFORCE_BASELINE=1; shift ;;
    *) printf 'unknown midi-reed option: %s\n' "$1" >&2; exit 2 ;;
  esac
done
[[ -n "$OUT_DIR" && ! -e "$OUT_DIR" && -x "$BINARY" ]] || { printf 'midi-reed requires a fresh output directory and music_render binary\n' >&2; exit 2; }
IMPORTER="$(dirname "$BINARY")/examples/midi_reed"
[[ -x "$IMPORTER" ]] || { printf 'missing %s; centrally build cargo build -p fs-couple --example midi_reed\n' "$IMPORTER" >&2; exit 2; }
if [[ "$DRY_RUN" == 1 ]]; then
  printf 'midi-reed: %s OUT/import; %s reed OUT/a.wav --seconds 0.5 --schedule OUT/import/score.gesture; replay/refusal/baseline\n' "$IMPORTER" "$BINARY"
  exit 0
fi
mkdir -p "$OUT_DIR"
RECEIPTS="$OUT_DIR/receipts.jsonl"
RUN_HEAD="$(git rev-parse HEAD)"
RUN_HOST="$(uname -n)"
RUN_COMMAND="$IMPORTER OUT/import; $BINARY reed OUT/a.wav --seconds 0.5 --schedule OUT/import/score.gesture"
START_SECONDS=$SECONDS
"$IMPORTER" "$OUT_DIR/import" > "$OUT_DIR/import.stdout" 2> "$OUT_DIR/import.stderr"
META="$OUT_DIR/import/performance-provenance.json"
receipt() {
  python3 -I - "$META" "$1" "$2" "${3:-}" "$RUN_HEAD" "$RUN_HOST" "$RUN_COMMAND" "$((SECONDS-START_SECONDS))" <<'PY' | tee -a "$RECEIPTS"
import json,sys
meta=json.load(open(sys.argv[1]))
meta.update(suite="music-e2e",lane="midi-reed",stage=sys.argv[2],verdict=sys.argv[3],wav_blake3=sys.argv[4] or None,git_sha=sys.argv[5],host=sys.argv[6],command=sys.argv[7],elapsed_s=int(sys.argv[8]))
print(json.dumps(meta,separators=(",",":")))
PY
}
receipt import lowered-and-model-resolved
"$IMPORTER" "$OUT_DIR/replay-import" > "$OUT_DIR/replay-import.stdout" 2> "$OUT_DIR/replay-import.stderr"
for file in score.mid score.performance score.gesture score.assembly performance-provenance.json; do
  cmp "$OUT_DIR/import/$file" "$OUT_DIR/replay-import/$file"
done
render() { "$BINARY" reed "$1" --seconds 0.5 --schedule "$OUT_DIR/import/score.gesture"; }
for run in a b; do render "$OUT_DIR/$run.wav" > "$OUT_DIR/render-$run.stdout" 2> "$OUT_DIR/render-$run.stderr"; done
cmp "$OUT_DIR/a.wav" "$OUT_DIR/b.wav"
cmp "$OUT_DIR/a.provenance.json" "$OUT_DIR/b.provenance.json"
WAV_HASH="$(python3 -I - "$OUT_DIR/a.provenance.json" "$META" <<'PY'
import json,sys
r=json.load(open(sys.argv[1])); p=json.load(open(sys.argv[2]))
assert r["gesture_schedule"]["blake3"]==p["schedule_blake3"]
assert len(r["wav_blake3"])==64
print(r["wav_blake3"])
PY
)"
receipt replay bit-identical "$WAV_HASH"
python3 -I - "$OUT_DIR" "$META" <<'PY'
import json,pathlib,sys
out=pathlib.Path(sys.argv[1]); performance=json.load(open(sys.argv[2]))
for run in ("a","b"):
    render=json.load(open(out/(run+".provenance.json")))
    render["performance"]=performance
    with open(out/(run+".performance-provenance.json"),"x") as f:
        json.dump(render,f,separators=(",",":")); f.write("\n")
PY
if render "$OUT_DIR/a.wav" > "$OUT_DIR/refusal.stdout" 2>&1; then
  printf 'overwrite refusal arm did not fire\n' >&2; exit 2
fi
grep -q 'refuses to overwrite evidence' "$OUT_DIR/refusal.stdout"
cmp "$OUT_DIR/a.wav" "$OUT_DIR/b.wav"
if "$IMPORTER" "$OUT_DIR/unmapped" --refuse-unmapped > "$OUT_DIR/unmapped.stdout" 2>&1; then
  printf 'unmapped-note refusal arm did not fire\n' >&2; exit 2
fi
grep -q 'UnknownNote' "$OUT_DIR/unmapped.stdout"
[[ ! -e "$OUT_DIR/unmapped" ]]
receipt refusal-arm refused-as-required "$WAV_HASH"
BASELINE="$(python3 -I - scripts/e2e/music-baselines.json "$WAV_HASH" <<'PY'
import json,sys
recorded=json.load(open(sys.argv[1])).get("midi-reed")
print("no-baseline" if not recorded else "matched" if recorded==sys.argv[2] else "drifted")
PY
)"
receipt baseline "$BASELINE" "$WAV_HASH"
if [[ "$ENFORCE_BASELINE" == 1 && "$BASELINE" != matched ]]; then
  printf 'midi-reed baseline is %s; retain and register the first qualified native hash before enforcement\n' "$BASELINE" >&2; exit 2
fi
receipt terminal green "$WAV_HASH"
