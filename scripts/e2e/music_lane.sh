#!/usr/bin/env bash
# Music-program e2e lane engine (music bead frankensim-music-v8-root-3ez8g.16).
#
# Renders ONE pinned fixture through the committed music_render binary,
# proves same-host replay (two renders, byte-identical WAV + sidecar),
# checks the never-overwrite refusal arm, and emits one JSON-lines
# receipt per stage plus a terminal receipt chaining them. Artifacts
# (WAV + provenance sidecar) are the retained evidence.
#
# This validates software determinism and encoding law only. It is not
# a physical-validation, realism, or listening claim; those live in the
# claims registry's gates and listening receipts.
#
# Baseline discipline: WAV hashes are ONE-HOST bit-deterministic (the
# music stack claims no cross-ISA replay; zero cross-ISA goldens exist).
# The committed baseline in scripts/e2e/music-baselines.json is therefore
# compared as a NOTE by default (verdict matched|drifted|no-baseline) and
# enforced only under --enforce-baseline (the recording host / CI lane).
# A drift under enforcement is a golden event: explain the semantic cause
# and update the baseline in the same commit, never loosen the check.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

readonly BASELINES="scripts/e2e/music-baselines.json"

FIXTURE=""
OUT_DIR=""
BINARY=""
SECONDS_ARG="0.25"
ENFORCE_BASELINE=0
DRY_RUN=0
QUALIFICATION_LOG=""

usage() {
  printf '%s\n' \
    'usage: scripts/e2e/music_lane.sh --fixture <reed|string|reed-flanged-baffled> --out-dir PATH [options]' \
    '' \
    'Options:' \
    '  --binary PATH        music_render executable (default: auto-discover under' \
    '                       $CARGO_TARGET_DIR or target/, debug then release).' \
    '  --seconds S          Render length in seconds (default 0.25).' \
    '  --qualification-log PATH  Native baffled_pressure test stdout (required for reed-flanged-baffled).' \
    '  --enforce-baseline   Fail (not note) on baseline hash drift for this host lane.' \
    '  --dry-run            Validate inputs, print the plan, create nothing.' \
    '  -h|--help            This help.' \
    '' \
    'The out dir must not exist; the lane creates it and refuses to overwrite' \
    'evidence. Receipts go to <out-dir>/receipts.jsonl.'
}

die() {
  printf '{"suite":"music-e2e","lane":"%s","verdict":"refused","elapsed_s":%s,"reason":"%s"}\n' \
    "${FIXTURE:-unset}" "$SECONDS" "$*" >&2
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fixture) [[ $# -ge 2 ]] || die '--fixture requires reed|string'; FIXTURE="$2"; shift 2 ;;
    --out-dir) [[ $# -ge 2 ]] || die '--out-dir requires a path'; OUT_DIR="$2"; shift 2 ;;
    --binary)  [[ $# -ge 2 ]] || die '--binary requires a path'; BINARY="$2"; shift 2 ;;
    --seconds) [[ $# -ge 2 ]] || die '--seconds requires a number'; SECONDS_ARG="$2"; shift 2 ;;
    --qualification-log) [[ $# -ge 2 ]] || die '--qualification-log requires a path'; QUALIFICATION_LOG="$2"; shift 2 ;;
    --enforce-baseline) ENFORCE_BASELINE=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ "$FIXTURE" == "reed" || "$FIXTURE" == "string" || "$FIXTURE" == "reed-flanged-baffled" ]] || die '--fixture must be reed, string or reed-flanged-baffled'
[[ -n "$OUT_DIR" ]] || die '--out-dir is required'
[[ ! -e "$OUT_DIR" ]] || die "refusing to overwrite existing out dir: $OUT_DIR"
ON_AXIS_ERROR=""
if [[ "$FIXTURE" == "reed-flanged-baffled" ]]; then
  [[ "$SECONDS_ARG" == "0.25" ]] || die 'reed-flanged-baffled duration is pinned to 0.25 s in its performance file'
  [[ -f "$QUALIFICATION_LOG" ]] || die 'reed-flanged-baffled requires --qualification-log from the native baffled_pressure qualification tests'
  ON_AXIS_ERROR="$(python3 -I -c '
import json,sys
rows=[]
for line in open(sys.argv[1]):
    line=line.strip()
    if not line.startswith("{"): continue
    row=json.loads(line)
    if row.get("case")=="baffled_pressure_on_axis_piston": rows.append(row)
assert len(rows)==3 and {r["rung"] for r in rows}=={2,10,100}
assert all(r["verdict"]=="pass" and r["tol"]==1e-3 and 0<=r["err"]<=1e-3 for r in rows)
print(max(r["err"] for r in rows))' "$QUALIFICATION_LOG")" || die 'on-axis qualification is missing or outside its 1e-3 relative band'
fi

if [[ -z "$BINARY" ]]; then
  for candidate in \
    "${CARGO_TARGET_DIR:-target}/debug/music_render" \
    "${CARGO_TARGET_DIR:-target}/release/music_render" \
    "target/debug/music_render" "target/release/music_render"; do
    if [[ -x "$candidate" ]]; then BINARY="$candidate"; break; fi
  done
fi
[[ -n "$BINARY" && -x "$BINARY" ]] || \
  die 'music_render binary not found; build it (cargo build -p fs-couple --bin music_render) or pass --binary'

if [[ $DRY_RUN -eq 1 ]]; then
  printf '{"suite":"music-e2e","lane":"%s","verdict":"dry-run","binary":"%s","seconds":%s,"out_dir":"%s"}\n' \
    "$FIXTURE" "$BINARY" "$SECONDS_ARG" "$OUT_DIR"
  exit 0
fi

mkdir -p "$OUT_DIR"
RECEIPTS="$OUT_DIR/receipts.jsonl"
: > "$RECEIPTS"

RUN_HEAD="$(git rev-parse HEAD)"
RUN_DIRTY="$(git status --porcelain --untracked-files=no | wc -l | tr -d ' ')"
RUN_HOST="$(uname -n)"
RUN_COMMAND="$BINARY $FIXTURE --seconds $SECONDS_ARG"
if [[ "$FIXTURE" == "reed-flanged-baffled" ]]; then
  RUN_COMMAND="$BINARY wind crates/fs-couple/examples/reed-flanged-baffled.performance OUT.wav --acceptance"
fi
receipt() {
  python3 -I -c 'import json,sys; row=json.loads(sys.argv[1]); row.update(git_sha=sys.argv[2],dirty=int(sys.argv[3])!=0,host=sys.argv[4],command=sys.argv[5],profile="supplied-binary",elapsed_s=int(sys.argv[6])); row.setdefault("contact_law_source","none/contact-free-string-control" if row.get("lane")=="string" else "not-retained-in-this-receipt"); [row.setdefault(key,None) for key in ("k_per_measure","alpha","chi","weight_measure","card_hash")]; print(json.dumps(row,separators=(",",":")))' \
    "$1" "$RUN_HEAD" "$RUN_DIRTY" "$RUN_HOST" "$RUN_COMMAND" "$SECONDS" | tee -a "$RECEIPTS"
}

render() {
  if [[ "$FIXTURE" == "reed-flanged-baffled" ]]; then
    "$BINARY" wind crates/fs-couple/examples/reed-flanged-baffled.performance "$1" --acceptance
  else
    "$BINARY" "$FIXTURE" "$1" --seconds "$SECONDS_ARG"
  fi
}

hash_field() { # extract wav_blake3 from a provenance sidecar
  sed -n 's/.*"wav_blake3":"\([0-9a-f]*\)".*/\1/p' "$1"
}

# Stage 1+2: two independent renders (replay proof).
for run in a b; do
  WAV="$OUT_DIR/$FIXTURE-$run.wav"
  if ! render "$WAV" > "$OUT_DIR/render-$run.stdout"; then
    die "render $run refused; see $OUT_DIR/render-$run.stdout"
  fi
  [[ -s "$WAV" ]] || die "render $run produced no WAV"
  [[ -s "$OUT_DIR/$FIXTURE-$run.provenance.json" ]] || die "render $run produced no sidecar"
done
HASH_A="$(hash_field "$OUT_DIR/$FIXTURE-a.provenance.json")"
HASH_B="$(hash_field "$OUT_DIR/$FIXTURE-b.provenance.json")"
[[ -n "$HASH_A" && "$HASH_A" == "$HASH_B" ]] || \
  die "replay hash mismatch: $HASH_A vs $HASH_B (same-host determinism regression)"
cmp -s "$OUT_DIR/$FIXTURE-a.wav" "$OUT_DIR/$FIXTURE-b.wav" || \
  die 'replay byte mismatch: WAV files differ despite matching sidecars'
receipt "{\"stage\":\"replay\",\"lane\":\"$FIXTURE\",\"verdict\":\"bit-identical\",\"wav_blake3\":\"$HASH_A\",\"seconds\":$SECONDS_ARG}"

# Stage 3: the never-overwrite refusal arm must fire.
if render "$OUT_DIR/$FIXTURE-a.wav" \
    > "$OUT_DIR/refusal.stdout" 2>&1; then
  die 'overwrite refusal arm did not fire; the lane cannot trust the artifacts'
fi
REFUSAL_REASON='refuses to overwrite evidence'
if [[ "$FIXTURE" == "reed-flanged-baffled" ]]; then
  REFUSAL_REASON='output/sidecar paths must be distinct and must not already exist'
fi
grep -q "$REFUSAL_REASON" "$OUT_DIR/refusal.stdout" || \
  die 'overwrite refusal fired with the wrong reason'
receipt "{\"stage\":\"refusal-arm\",\"lane\":\"$FIXTURE\",\"verdict\":\"refused-as-required\"}"

# Stage 4: baseline comparison (note by default; gate under enforcement).
BASELINE_VERDICT="no-baseline"
if [[ -f "$BASELINES" ]]; then
  RECORDED="$(sed -n "s/.*\"$FIXTURE\": *\"\([0-9a-f]*\)\".*/\1/p" "$BASELINES" | head -1)"
  if [[ -n "$RECORDED" ]]; then
    if [[ "$RECORDED" == "$HASH_A" ]]; then BASELINE_VERDICT="matched"; else BASELINE_VERDICT="drifted"; fi
  fi
fi
receipt "{\"stage\":\"baseline\",\"lane\":\"$FIXTURE\",\"verdict\":\"$BASELINE_VERDICT\",\"recorded\":\"${RECORDED:-}\",\"observed\":\"$HASH_A\",\"determinism_class\":\"one-host\",\"change_cause\":\"M4.1-W1: physical contact measures and explicit friction parameters; string lane is a contact-free control\"}"
if [[ $ENFORCE_BASELINE -eq 1 && "$BASELINE_VERDICT" == "drifted" ]]; then
  die "baseline drift under enforcement: recorded ${RECORDED:-} vs observed $HASH_A (golden event: explain and update in the same commit)"
fi
if [[ $ENFORCE_BASELINE -eq 1 && "$FIXTURE" == "reed-flanged-baffled" && "$BASELINE_VERDICT" == "no-baseline" ]]; then
  die 'reed-flanged-baffled baseline has no native hash yet; record the first qualified one-host run before enforcement'
fi

if [[ "$FIXTURE" == "reed-flanged-baffled" ]]; then
  OBSERVATION="$(python3 -I -c 'import json,sys; p=json.load(open(sys.argv[1])); r=p["reed_input"]; o=r["outlet_receiver"]; assert r["observation_policy"]=="acceptance"; print(json.dumps(dict(stage="observation",observer="baffled_outlet",termination="flanged",rings=o["radial_rings"],angles=o["angular_points"],band_hz=o["maximum_frequency_hz"],filter_delay_samples=0,exterior_gas=o["gas"],on_axis_check_rel_err=float(sys.argv[2]),qualification=sys.argv[3],verdict="pass"),separators=(",",":")))' "$OUT_DIR/$FIXTURE-a.provenance.json" "$ON_AXIS_ERROR" "$QUALIFICATION_LOG")"
  receipt "$OBSERVATION"
fi

# Terminal receipt: chain the stages.
STAGES='["replay","refusal-arm","baseline"]'
if [[ "$FIXTURE" == "reed-flanged-baffled" ]]; then STAGES='["replay","refusal-arm","baseline","observation"]'; fi
receipt "{\"stage\":\"terminal\",\"lane\":\"$FIXTURE\",\"verdict\":\"green\",\"stages\":$STAGES,\"wav_blake3\":\"$HASH_A\",\"baseline\":\"$BASELINE_VERDICT\",\"elapsed_s\":$SECONDS}"
