#!/usr/bin/env bash
# IRP M3.2-W1: central release qualification, no Mac builds.
# Usage: lbm_qualification.sh --lane d3q27-invariance|dissipation-ladder \
#        --out DIR [--runner native|rch]
# The optional rch runner is for central Linux scheduling, never this Mac.
# Each lane retains the real test output and execution-stamped JSONL.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
LANE="" OUT="" RUNNER="native"
while (($#)); do
  case "$1" in
    --lane) LANE="${2:?missing lane}"; shift 2 ;;
    --out|--out-dir) OUT="${2:?missing output directory}"; shift 2 ;;
    --runner) RUNNER="${2:?missing runner}"; shift 2 ;;
    *) echo "usage: $0 --lane d3q27-invariance|dissipation-ladder --out DIR [--runner native|rch]" >&2; exit 2 ;;
  esac
done
case "$LANE" in
  d3q27-invariance) TEST=rotational_invariance_d3q27_vs_d3q19_control ;;
  dissipation-ladder) TEST=collision_dissipation_ladder ;;
  *) echo "LBM qualification: declare one supported --lane" >&2; exit 2 ;;
esac
[[ -n "$OUT" ]] || { echo "LBM qualification: --out required" >&2; exit 2; }
[[ "$(uname -s)" != Darwin ]] || { echo "LBM qualification: Mac Cargo/RCH refused; use central Linux" >&2; exit 2; }
[[ ! -e "$OUT" ]] || { echo "LBM qualification: existing output refused: $OUT" >&2; exit 3; }
case "$RUNNER" in native|rch) ;; *) echo "unknown runner: $RUNNER" >&2; exit 2 ;; esac
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
cd "$ROOT"
CMD=(env CARGO_BUILD_JOBS=12 cargo test --release --locked -p fs-lbm --test d3q27_battery "$TEST" -- --exact --nocapture)
if [[ "$RUNNER" == rch ]]; then
  export RCH_REQUIRE_REMOTE=1
  CMD=(rch exec -- "${CMD[@]}")
fi
STARTED="$(python3 -I -c 'import time; print(time.monotonic())')"
set +e
"${CMD[@]}" 2>&1 | tee "$OUT/native.log"
CODE=${PIPESTATUS[0]}
set -e
python3 -I - "$ROOT" "$OUT" "$LANE" "$STARTED" "$CODE" "${CMD[@]}" <<'PY'
import json,pathlib,platform,subprocess,sys,time
root,out,lane,started,code,*command=sys.argv[1:]
out=pathlib.Path(out)
sha=subprocess.check_output(['git','-C',root,'rev-parse','HEAD'],text=True).strip()
dirty=bool(subprocess.check_output(['git','-C',root,'status','--porcelain'],text=True).strip())
context=dict(git_sha=sha,host=platform.node(),command=command,profile='release',dirty=dirty,
             native_exit_code=int(code),lane=lane)
rows=[]
for line in (out/'native.log').read_text().splitlines():
    at=line.find('{"suite":"fs-lbm/d3q27"')
    if at<0:continue
    r=json.loads(line[at:])
    if r.get('lane')!=lane:continue
    detail=r.pop('detail',{})
    r.update(detail)
    r.update(context)
    # null means this observable is not part of this lane, not a zero result.
    for key in ('rungs','observed_order','rate_spread_per_rung','control_ratio','tolerance'):
        r.setdefault(key,None)
    rows.append(r)
terminal=[r for r in rows if r.get('stage')=='terminal']
expected=2 if lane=='dissipation-ladder' else 1
passed=(int(code)==0 and len(terminal)==expected and all(r['verdict']=='pass' for r in rows))
rows.append(dict(context,suite='fs-lbm/d3q27',stage='execution-terminal',
                 verdict='pass' if passed else 'fail',elapsed_s=time.monotonic()-float(started),
                 acceptance=False,reason='module qualification only; no jet/edge-tone/NSF instrument claim'))
with (out/'run.jsonl').open('x') as f:
    for r in rows:f.write(json.dumps(r,sort_keys=True,allow_nan=False)+'\n')
if not passed:raise SystemExit(int(code) or 1)
print(json.dumps(dict(lane=lane,verdict='pass',receipt=str(out/'run.jsonl'))))
PY
