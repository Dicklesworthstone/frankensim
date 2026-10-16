#!/usr/bin/env bash
# Thin deterministic replay/refusal lane. Requires a supplied binary and native
# qualification stdout; never builds or manufactures a passing numerical row.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
OUT=""; BIN=""; QUAL=""; DRY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT="$2"; shift 2 ;;
    --binary) BIN="$2"; shift 2 ;;
    --qualification-log) QUAL="$2"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    *) printf '%s\n' "unsupported modal-broadband argument: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$OUT" && ! -e "$OUT" && -x "$BIN" && -f "$QUAL" ]] || {
  printf '%s\n' 'need fresh --out-dir, executable --binary, and --qualification-log from broadband_modal_transfer native tests' >&2; exit 2;
}
if [[ "$DRY" == 1 ]]; then printf '%s\n' 'modal-broadband: supplied native qualification -> analytic sphere -> replay -> refusal -> baseline -> terminal'; exit 0; fi
mkdir -p "$OUT"
python3 -I - "$OUT" "$QUAL" <<'PY'
import cmath,json,math,pathlib,sys
out=pathlib.Path(sys.argv[1]); log=pathlib.Path(sys.argv[2]).read_text()
assert 'test result: ok.' in log and 'test result: FAILED' not in log
rows=[json.loads(line) for line in log.splitlines() if line.startswith('{')]
drives=[r for r in rows if r.get('case')=='off-resonance']
assert {r['ratio'] for r in drives}=={.25,.5,2.,4.,8.}
assert all(r['amplitude_rel_error']<=1e-6 and r['phase_error_rad']<=1e-6 and r['held_out_max_rel']<=1e-6 for r in drives)
axes={name:[r for r in rows if r.get('axis')==name] for name in ('poles','fit-grid','fractional-delay')}
for ladder in axes.values():
    assert len(ladder)==3 and all(b['held_out_max_rel']<a['held_out_max_rel'] for a,b in zip(ladder,ladder[1:]))
    assert ladder[-1]['held_out_max_rel']<=1e-6
qualified=dict(rungs=axes['poles'],fit_grid_rungs=axes['fit-grid'],delay_rungs=axes['fractional-delay'],
    off_resonance_max_rel_error=max(r['amplitude_rel_error'] for r in drives),
    off_resonance_max_phase_error_rad=max(r['phase_error_rad'] for r in drives))
(out/'qualification.json').write_text(json.dumps(qualified))
w0=math.tau*500; a=.05; r=1.; rho=1.2; c=343.; delay=(r-a)/c
def h(w): return -1j*w*rho*cmath.exp(1j*w*delay)/(4*math.pi*r*(1-1j*w*a/c))
at_mode=h(w0)
base=f'''frankensim-modal-performance-v1
sample_rate_hz 48000
samples 12000
full_scale_pa 1
limits 0.9 1 100000 1000000000 10000000
compile_limits 4 16
voices 1
voice retain-state 1 1
mode {w0:.17g} 0.01 {at_mode.real:.17g} {at_mode.imag:.17g} 0 0
port 0 1
events 2
force 0 0 0 1
force 240 0 0 0
'''
(out/'legacy.performance').write_text(base)
low=w0*.125; high=w0*9
lines=[base.rstrip(), 'broadband 1',f'transfer analytic:sphere:a=.05:r=1:rho=1.2:c=343 {delay:.17g} 0 1 32 128 1e-6 256 127']
for i in range(256):
    w=low+(high-low)*i/255; z=h(w); lines.append(f'frequency {w:.17g} {z.real:.17g} {z.imag:.17g}')
for i in range(127):
    w=low+(high-low)*(i+.371)/127; z=h(w); lines.append(f'withheld {w:.17g} {z.real:.17g} {z.imag:.17g}')
(out/'sphere.performance').write_text('\n'.join(lines)+'\n')
PY
START=$SECONDS
"$BIN" modal "$OUT/sphere.performance" "$OUT/a.wav" --block 37 > "$OUT/a.stdout"
"$BIN" modal "$OUT/sphere.performance" "$OUT/b.wav" --block 37 > "$OUT/b.stdout"
cmp "$OUT/a.wav" "$OUT/b.wav"
cmp "$OUT/a.provenance.json" "$OUT/b.provenance.json"
if "$BIN" modal "$OUT/sphere.performance" "$OUT/a.wav" --block 37 > "$OUT/refusal.stdout" 2> "$OUT/refusal.stderr"; then
  printf '%s\n' 'overwrite refusal did not fire' >&2; exit 1
fi
cmp "$OUT/a.wav" "$OUT/b.wav"
cmp "$OUT/a.provenance.json" "$OUT/b.provenance.json"
"$BIN" modal "$OUT/legacy.performance" "$OUT/legacy.wav" --block 37 > "$OUT/legacy.stdout"
python3 -I - "$OUT" "$(git rev-parse HEAD)" "$(uname -n)" "$BIN" "$((SECONDS-START))" <<'PY'
import json,pathlib,sys
out=pathlib.Path(sys.argv[1]); broad=json.loads((out/'a.provenance.json').read_text()); old=json.loads((out/'legacy.provenance.json').read_text())
assert broad['clipped_samples']==0 and broad['peak_pa']>0 and old['clipped_samples']==0
r=broad['modal_input']['broadband'][0]
assert r['authority']=='EstimateOnly' and r['held_out_max_rel']<=1e-6 and r['producer_relative_error_estimate']==0
assert broad['wav_blake3']!=old['wav_blake3'], 'off-resonance receiver mutation is not live'
q=json.loads((out/'qualification.json').read_text())
common=dict(suite='music-e2e',lane='modal-broadband',git_sha=sys.argv[2],host=sys.argv[3],command=f'{sys.argv[4]} modal sphere.performance OUT.wav --block 37',wall_s=int(sys.argv[5]),producer_id=r['producer_id'],fit_order=r['fit_order'],delay_samples=r['delay_samples'],tolerance=1e-6,authority='EstimateOnly',**q)
with (out/'receipts.jsonl').open('x') as f:
    for stage in ('replay','refusal-arm','baseline','terminal'):
        row=dict(common,stage=stage,verdict='pass',wav_blake3=broad['wav_blake3'],legacy_wav_blake3=old['wav_blake3'])
        line=json.dumps(row,separators=(',',':')); print(line); f.write(line+'\n')
PY
