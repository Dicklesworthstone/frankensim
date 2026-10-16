#!/usr/bin/env bash
# Bar-S1 native structure/radiation/contact ladder and real Iowa grading.
# Run on central Linux; never run a Cargo workload on the shared Mac.
# Usage: IRP_GRADE_SCRATCH=/path/containing/grade/library.json \
#   scripts/e2e/music/bar_s1.sh --out-dir /tmp/bar-s1
# Grader assets: SCRATCH/grade/library.json and SCRATCH/real/wav48, or
# IRP_UIOWA_WAV48 for the latter. No download/conversion or playback occurs.
# The script writes the requested grader JSON into data/instruments/grades/bars;
# MagentaBadger commits those generated results after reviewing native proof.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
OUT=""
SCRATCH="${IRP_GRADE_SCRATCH:-}"
while (($#)); do
  case "$1" in
    --out-dir) OUT="${2:?missing --out-dir value}"; shift 2 ;;
    --grade-scratch) SCRATCH="${2:?missing --grade-scratch value}"; shift 2 ;;
    *) echo "usage: $0 --out-dir DIR [--grade-scratch DIR]" >&2; exit 2 ;;
  esac
done
[[ -n "$OUT" && -n "$SCRATCH" ]] || {
  echo "BAR-S1: --out-dir and IRP_GRADE_SCRATCH/--grade-scratch are required" >&2; exit 2;
}
[[ "$(uname -s)" != Darwin ]] || {
  echo "BAR-S1: native builds belong to central Linux; Mac build refused" >&2; exit 2;
}
[[ -f "$SCRATCH/grade/library.json" ]] || {
  echo "BAR-S1: missing grader library: $SCRATCH/grade/library.json" >&2; exit 2;
}
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
SCRATCH="$(cd "$SCRATCH" && pwd)"
cd "$ROOT"
export CARGO_BUILD_JOBS=12
export BAR_S1_OUT_DIR="$OUT"
STARTED="$(python3 -I -c 'import time; print(time.monotonic())')"
record_execution() {
  local code="$?"
  python3 -I - "$OUT" "$ROOT" "$SCRATCH" "$STARTED" "$code" <<'PY'
import json, pathlib, platform, subprocess, sys, time
out, root, scratch, started, code = sys.argv[1:]
def git(*args):
    return subprocess.check_output(['git', '-C', root, *args], text=True).strip()
row = dict(stage='execution', host=platform.node(), git_sha=git('rev-parse', 'HEAD'),
           dirty=bool(git('status', '--porcelain')), profile='release',
           command=['scripts/e2e/music/bar_s1.sh', '--out-dir', out, '--grade-scratch', scratch],
           native_command=['cargo','test','--release','--locked','-p','fs-couple',
                           '--test','bar_s1_linear3_bar','--','--nocapture','--test-threads=1'],
           exit_code=int(code), elapsed_s=time.monotonic()-float(started),
           verdict='pass' if int(code)==0 else 'fail', acceptance=False)
with (pathlib.Path(out)/'execution.jsonl').open('a') as f:
    f.write(json.dumps(row, sort_keys=True)+'\n')
PY
  return "$code"
}
trap record_execution EXIT

# All named tests run, including sparse-vs-dense, live windows/work caps,
# the allowed one-cent refusal/TET10 route, and bit-identical render replay.
cargo test --release --locked -p fs-couple --test bar_s1_linear3_bar -- \
  --nocapture --test-threads=1 2>&1 | tee "$OUT/native.log"

python3 -I - "$ROOT" "$OUT" <<'PY'
import json, pathlib, subprocess, sys
root, out = map(pathlib.Path, sys.argv[1:])
rows=[]
for line in (out/'native.log').read_text().splitlines():
    start=line.find('{"suite":"bar-s1"')
    if start >= 0:
        rows.append(json.loads(line[start:]))
required={'uniform-prism-vs-free-free-Timoshenko', 'undercut-modal-uncertainty',
          'panel-pressure-convergence', 'Hertz-rigid-target-contact-time',
          'six-rigid-modes', 'sparse-vs-dense-first-ten-elastic',
          'window-and-work-cap-live', 'render-replay-and-time-convergence'}
assert required <= {r['case'] for r in rows}, 'native receipt missing required case'
assert all(r['verdict'] in ('pass','refused') for r in rows), 'native gate failed'
assert all(r['case'] in ('undercut-modal-uncertainty','undercut-band')
           for r in rows if r['verdict']=='refused'), 'unexpected refusal'
(out/'physics.jsonl').write_text(''.join(json.dumps(r,sort_keys=True)+'\n' for r in rows))
sha=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
candidates=[dict(id='bar-s1-before',label='Appendix A xylophone strike 1: uniform rosewood prism',
    wav=str(root/'data/listening/xylophone-strike.wav'),t0=0.0,t1=1.0,
    inst='xylophone',note='Db5',dyn='mf',art='',kind='decaying',
    expected_db=28.4,tolerance_db=0.1,
    render_command='historical committed xylophone_bar.rs fixture'),
    dict(id='bar-s1-after',label='Bar-S1: 3-D designed undercut, BM receiver, material Hertz; diagnostic',
    wav=str(out/'bar-s1-after.wav'),inst='xylophone',note='Db5',dyn='mf',art='',kind='decaying',
    render_git_sha=sha,render_command='cargo test --release --locked -p fs-couple --test bar_s1_linear3_bar -- --nocapture --test-threads=1')]
(out/'candidates.json').write_text(json.dumps(candidates,indent=2)+'\n')
PY

uv run --with numpy --with matplotlib python3 -I \
  apps/instrument-lab/tools/instrgrade/grade_run.py \
  "$SCRATCH" "$OUT/candidates.json" "$OUT/grade" | tee "$OUT/grade.log"

python3 -I - "$ROOT" "$OUT" <<'PY'
import json, pathlib, sys
root,out=map(pathlib.Path,sys.argv[1:])
results=json.loads((out/'grade/results.json').read_text())
by_id={r['id']:r for r in results}
assert set(by_id)=={'bar-s1-before','bar-s1-after'}
receipts=[json.loads(l) for l in (out/'grade/receipt.jsonl').read_text().splitlines()]
assert receipts[-1]['verdict']=='pass', 'grader terminal failure'
for r in results:
    assert 'error' not in r and r.get('baselines') and r.get('rulers'), 'baseline/ruler missing'
before=by_id['bar-s1-before']['distance']['spec_db']
after=by_id['bar-s1-after']['distance']['spec_db']
assert abs(before-28.4)<=0.1, 'historical 28.4 dB baseline not reproduced'
dest=root/'data/instruments/grades/bars'
dest.mkdir(parents=True,exist_ok=True)
for r in results:
    # Existing committed evidence must reproduce; do not overwrite it silently.
    path=dest/(r['id']+'.json')
    if path.exists():
        assert json.loads(path.read_text())==r, f'grade drift at {path}; review/rebaseline explicitly'
    else:
        with path.open('x') as f:
            json.dump(r,f,indent=1); f.write('\n')
metadata=json.loads((out/'bar-s1-render.json').read_text())
metadata.update(stage='grade-and-baselines',before_db=before,after_db=after,
    grader_json_path='data/instruments/grades/bars/bar-s1-after.json',
    baseline_json_path='data/instruments/grades/bars/bar-s1-before.json',
    native_log=str(out/'native.log'),grader_receipt=str(out/'grade/receipt.jsonl'),
    verdict='pass',acceptance=False)
(out/'receipt.jsonl').write_text(json.dumps(metadata,sort_keys=True)+'\n')
print(json.dumps(dict(stage='complete',before_db=before,after_db=after,
                     receipt=str(out/'receipt.jsonl'),acceptance=False)))
PY
