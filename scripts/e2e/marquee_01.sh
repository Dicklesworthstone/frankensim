#!/usr/bin/env bash
# Actual study/ledger/report/package workflow (q61wp.20): the scalar thermal
# study, then the free-boundary 2-D elasticity study (plain and projected-area).
# The plain elasticity study must end on its material-area target (q61wp.16.1).
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
BINARY="${FRANKENSIM_BIN:-${CARGO_TARGET_DIR:-${REPO_ROOT}/target}/debug/frankensim}"
case "$COMMAND" in
  --list)
    printf '%s\n' thermal_study report_svg package_checker budget_resume wrong_units wrong_load wrong_area \
      elasticity_study elasticity_package_checker elasticity_budget_resume projected_area_study \
      elasticity_wrong_units elasticity_wrong_area elasticity_asymmetric_band
    exit 0 ;;
  --check|--self-test)
    test -x "$BINARY"
    command -v python3 >/dev/null
    printf '%s\n' 'preflight only: executable and JSON reader available; no study executed'
    exit 0 ;;
  --run|--retain) ;;
  *) printf '%s\n' 'usage: marquee_01.sh [--list|--check|--self-test|--run|--retain]' >&2; exit 2 ;;
esac
if ! test -x "$BINARY"; then
  printf '%s\n' 'Set FRANKENSIM_BIN to a built binary; compile through DSR/RCH before running this lane.' >&2
  exit 1
fi
# A fresh directory preserves earlier evidence and caller-owned databases.
ARTIFACT_BASE="${ARTIFACT_DIR:-${REPO_ROOT}/target/marquee-01}"
mkdir -p "$ARTIFACT_BASE"
RUN_DIR="$(mktemp -d "${ARTIFACT_BASE}/run.XXXXXX")"
STUDY="${REPO_ROOT}/examples/marquee/thermal-2d.fsim"
field() { python3 -c 'import json,sys; v=json.load(open(sys.argv[1])); print(v[sys.argv[2]])' "$1" "$2"; }
"$BINARY" --json study "$STUDY" "$RUN_DIR/study.db" > "$RUN_DIR/study.json"
RUN_ID="$(field "$RUN_DIR/study.json" run_id)"
"$BINARY" --json report "$RUN_ID" "$RUN_DIR/study.db" > "$RUN_DIR/report.json"
"$BINARY" --json package "$RUN_ID" "$RUN_DIR/study.db" > "$RUN_DIR/package.json"
set +e
"$BINARY" --json study "$STUDY" "$RUN_DIR/partial.db" --budget 2 > "$RUN_DIR/partial.json" 2> "$RUN_DIR/partial.stderr"
PARTIAL_EXIT=$?
set -e
test "$PARTIAL_EXIT" -eq 6
PARTIAL_ID="$(field "$RUN_DIR/partial.json" run_id)"
"$BINARY" --json study --resume "$PARTIAL_ID" "$RUN_DIR/partial.db" > "$RUN_DIR/resumed.json"
python3 - "$STUDY" "$RUN_DIR" <<'PY'
import pathlib,sys
source=pathlib.Path(sys.argv[1]).read_text()
root=pathlib.Path(sys.argv[2])
for name,old,new in [
    ('units',':unit "1"',':unit "J"'),
    ('load',':load-region "domain-unit-source"',':load-region "unsupported-region"'),
    ('area',':volume-fraction 0.853',':volume-fraction 1.5'),
]:
    assert old in source
    (root/(name+'.fsim')).write_text(source.replace(old,new))
PY
for name in units load area; do
  set +e
  "$BINARY" --json study "$RUN_DIR/$name.fsim" "$RUN_DIR/$name.db" > "$RUN_DIR/$name.json" 2> "$RUN_DIR/$name.stderr"
  ACTUAL_EXIT=$?
  set -e
  test "$ACTUAL_EXIT" -eq 4
  test ! -e "$RUN_DIR/$name.db"
done
# Elasticity: the canonical plain study, its budget stop and resume, the
# projected-area variant, and three admission refusals before any ledger write.
ELASTIC="${REPO_ROOT}/examples/marquee/bracket-2d.fsim"
PROJECTED="${REPO_ROOT}/examples/marquee/bracket-projected-volume-2d.fsim"
"$BINARY" --json study "$ELASTIC" "$RUN_DIR/elastic.db" > "$RUN_DIR/elastic.json"
ELASTIC_ID="$(field "$RUN_DIR/elastic.json" run_id)"
(cd "$RUN_DIR" && "$BINARY" --json report "$ELASTIC_ID" "$RUN_DIR/elastic.db" > "$RUN_DIR/elastic-report.json")
(cd "$RUN_DIR" && "$BINARY" --json package "$ELASTIC_ID" "$RUN_DIR/elastic.db" > "$RUN_DIR/elastic-package.json")
set +e
"$BINARY" --json study "$ELASTIC" "$RUN_DIR/elastic-partial.db" --budget 2 > "$RUN_DIR/elastic-partial.json" 2> "$RUN_DIR/elastic-partial.stderr"
PARTIAL_EXIT=$?
set -e
test "$PARTIAL_EXIT" -eq 6
ELASTIC_PARTIAL_ID="$(field "$RUN_DIR/elastic-partial.json" run_id)"
"$BINARY" --json study --resume "$ELASTIC_PARTIAL_ID" "$RUN_DIR/elastic-partial.db" > "$RUN_DIR/elastic-resumed.json"
"$BINARY" --json study "$PROJECTED" "$RUN_DIR/projected.db" > "$RUN_DIR/projected.json"
python3 - "$ELASTIC" "$RUN_DIR" <<'PY'
import pathlib,sys
source=pathlib.Path(sys.argv[1]).read_text()
root=pathlib.Path(sys.argv[2])
for name,old,new in [
    ('elastic-units',':unit "J"',':unit "K"'),
    ('elastic-area',':volume-fraction 0.45',':volume-fraction 1.5'),
    ('elastic-band',':load-band (0.375 0.625)',':load-band (0.25 0.625)'),
]:
    assert old in source
    (root/(name+'.fsim')).write_text(source.replace(old,new))
PY
for name in elastic-units elastic-area elastic-band; do
  set +e
  "$BINARY" --json study "$RUN_DIR/$name.fsim" "$RUN_DIR/$name.db" > "$RUN_DIR/$name.json" 2> "$RUN_DIR/$name.stderr"
  ACTUAL_EXIT=$?
  set -e
  test "$ACTUAL_EXIT" -eq 4
  test ! -e "$RUN_DIR/$name.db"
done
python3 - "$RUN_DIR" "$BINARY" <<'PY'
import hashlib,json,pathlib,sys
root=pathlib.Path(sys.argv[1])
read=lambda name:json.loads((root/(name+'.json')).read_text())
full,partial,resumed=map(read,['study','partial','resumed'])
assert full['status']=='completed'
assert partial['status']=='budget-exhausted'
assert partial['receipt']['iterations_completed']==2
assert resumed['status']=='completed'
assert full['receipt']['trace_hash']==resumed['receipt']['trace_hash']
report,package=map(read,['report','package'])
assert report['status']==package['status']=='ok'
assert package['checker']=='pass' and package['checker_authority']=='structural-integrity-only'
html=pathlib.Path(report['report_html']).read_text()
assert '<polyline' in html and 'DWR estimate' in html
summary=json.loads(pathlib.Path(report['report_json']).read_text())
assert summary['iterations_completed']==8
assert summary['final_compliance']>0 and abs(summary['final_area']-0.853)<1e-12
assert pathlib.Path(package['package']).stat().st_size>0
for name,code in [('units','study-objective-dimension-mismatch'),('load','study-load-non-boundary'),('area','study-volume-fraction-out-of-bounds')]:
    assert code in (root/(name+'.stderr')).read_text()
elastic,epartial,eresumed,projected=map(read,['elastic','elastic-partial','elastic-resumed','projected'])
assert elastic['status']=='completed' and eresumed['status']=='completed'
assert epartial['status']=='budget-exhausted' and epartial['receipt']['iterations_completed']==2
assert elastic['receipt']['trace_hash']==eresumed['receipt']['trace_hash']
assert projected['status']=='completed'
ereport,epackage=map(read,['elastic-report','elastic-package'])
assert ereport['status']==epackage['status']=='ok'
assert epackage['checker']=='pass' and epackage['checker_authority']=='structural-integrity-only'
esummary=json.loads((root/ereport['report_json']).read_text())
assert esummary['iterations_completed']==32 and esummary['final_compliance_j']>0
# The scheduled area projection (q61wp.16.1) must land the canonical bracket
# on its 0.45 target; before it, the 32-step endpoint fell inside a swing.
assert esummary['volume_constraint_satisfied'] is True, esummary
assert abs(esummary['final_material_area_m2']-0.45) <= 0.001, esummary
assert (root/epackage['package']).stat().st_size>0
for name,code in [('elastic-units','study-objective-dimension-mismatch'),('elastic-area','study-volume-fraction-out-of-bounds'),('elastic-band','cli-study-elasticity-load')]:
    assert code in (root/(name+'.stderr')).read_text()
receipt={'schema':'frankensim.ci.thermal-study-e2e.v1',
         'binary_sha256':hashlib.sha256(pathlib.Path(sys.argv[2]).read_bytes()).hexdigest(),
         'run':full['run_id'],'trace_hash':full['receipt']['trace_hash'],
         'elasticity_run':elastic['run_id'],'elasticity_trace_hash':elastic['receipt']['trace_hash'],
         'checks':['real-study','retained-report-svg','package-checker-export','budget-partial','resume-same-trace','wrong-units','wrong-load','wrong-area',
                   'elasticity-study','elasticity-package-checker-export','elasticity-budget-partial','elasticity-resume-same-trace','projected-area-study',
                   'elasticity-wrong-units','elasticity-wrong-area','elasticity-asymmetric-band'],
         'no_claim':'Scalar normalized thermal radius study and the free-boundary 2-D elasticity study (plain and projected-area). Package export invokes the structural checker. The plain elasticity study is asserted to end on its area target (q61wp.16.1). No physical validation or guaranteed PDE bound.'}
(root/'marquee-e2e-summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(root/'marquee-e2e-summary.json')
PY
# --retain deliberately keeps the same immutable run directory. It does not
# rewrite a tracked summary with a claim about a dirty or different source HEAD.
printf 'study artifacts retained at %s\n' "$RUN_DIR"
