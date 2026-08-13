#!/usr/bin/env bash
#
# e2e_extreal_arithmetic_audit.sh — arithmetic-root mutation campaign
# (bead frankensim-extreal-program-f85xj.3.7, terminal owner).
#
# Proves, with a kill matrix, that a defect planted at the certified-
# arithmetic trust root is DETECTED by the committed batteries — the
# scariest failure is a wrong interval every consumer happily accepts.
#
# Design constraints this script is built around:
#   * The shared working tree is NEVER mutated. Mutants are applied to a
#     scratch COPY of the 7-crate trust-root cone (fs-ivl + fs-math +
#     leaf deps), assembled tools/oracle-style as its own tiny workspace.
#     Agents and the repo sweeper commit the live tree at any moment, so
#     in-tree mutation would risk committing a planted defect; the
#     integrity stage hashes the real crates before and after and refuses
#     on any drift.
#   * A mutant that fails to BUILD is `inconclusive-build`, never
#     "killed". A red positive control aborts the whole campaign as
#     inconclusive: kills are only meaningful against a green baseline.
#   * Every mutation is an exact-match single-site text replacement whose
#     match count is asserted (=1) before the run; the applied file must
#     differ from the original (application integrity).
#
# Usage:
#   scripts/ci/e2e_extreal_arithmetic_audit.sh [--artifact-dir PATH]
#       [--mutant NAME]        # run one mutant (default: all)
#       [--keep-scratch]       # retain the scratch tree for debugging
#
# Kill layers, cheapest first (a mutant is charged to the FIRST layer
# that kills it):
#   L1  fs-ivl committed battery (unit + conformance + corpus + fuzz +
#       directed-rounding audit + eft bridge + predicates + taylor)
# A mutant surviving every layer is a SURVIVOR: this script exits nonzero
# and the finding must become a test-gap bead.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CONE=(fs-ivl fs-math fs-evidence fs-blake3 fs-obs fs-casebook fs-propcheck)

ARTIFACT_DIR=""
ONLY_MUTANT=""
KEEP_SCRATCH=0
SCRATCH_OVERRIDE=""
CONTROL_ONLY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --artifact-dir) ARTIFACT_DIR="${2:-}"; shift 2 ;;
    --mutant)       ONLY_MUTANT="${2:-}"; shift 2 ;;
    --scratch)      SCRATCH_OVERRIDE="${2:-}"; KEEP_SCRATCH=1; shift 2 ;;
    --control-only) CONTROL_ONLY=1; shift ;;
    --keep-scratch) KEEP_SCRATCH=1; shift ;;
    -h|--help)      sed -n '3,40p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 2 ;;
    *) printf 'FATAL: unknown argument %s\n' "$1" >&2; exit 2 ;;
  esac
done

if [[ -z "${ARTIFACT_DIR}" ]]; then
  ARTIFACT_DIR="$(mktemp -d "${TMPDIR:-/tmp}/extreal-audit-XXXXXX")"
fi
mkdir -p "${ARTIFACT_DIR}"
if [[ -n "${SCRATCH_OVERRIDE}" ]]; then
  SCRATCH="${SCRATCH_OVERRIDE}"
  mkdir -p "${SCRATCH}"
else
  SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/extreal-audit-scratch-XXXXXX")"
fi
TARGET_BASE="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_mutation"
MATRIX="${ARTIFACT_DIR}/kill-matrix.jsonl"
touch "${MATRIX}"

log() {
  printf '{"suite":"fs-extreal-arithmetic-audit","kind":"%s","message":"%s"}\n' "$1" "$2" \
    | tee -a "${MATRIX}" >&2
}

cleanup() {
  if [[ "${KEEP_SCRATCH}" -eq 0 ]]; then rm -rf "${SCRATCH}"; fi
}
trap cleanup EXIT

# ------------------------------------------------------------ integrity pre
integrity_hash() {
  # Content hash of the REAL trust-root crates (tracked files only, via the
  # git index so other agents' unrelated edits elsewhere don't perturb it).
  (cd "${REPO_ROOT}" && git ls-files -s crates/fs-ivl crates/fs-math | shasum -a 256 | cut -d' ' -f1)
}
PRE_HASH="$(integrity_hash)"
log integrity "pre-campaign trust-root index hash ${PRE_HASH}"

# ------------------------------------------------------- scratch assembly
assemble_scratch() {
  local dest="$1"
  mkdir -p "${dest}/crates"
  for crate in "${CONE[@]}"; do
    cp -R "${REPO_ROOT}/crates/${crate}" "${dest}/crates/${crate}"
  done
  cp "${REPO_ROOT}/rust-toolchain.toml" "${dest}/"
  # fs-math's dev-only frankenscipy oracle is a PATH reference whose
  # manifest cargo insists on reading even though a non-member's dev-deps
  # never build; the sibling does not exist under the scratch root, so the
  # copied manifest drops that single dev-dep line. Code is untouched and
  # the positive control validates the modified cone.
  python3 - "${dest}/crates/fs-math/Cargo.toml" <<'PYEOF'
import sys
path = sys.argv[1]
lines = [l for l in open(path).read().splitlines(keepends=True)
         if 'fsci-special' not in l]
open(path, 'w').writelines(lines)
PYEOF
  # Minimal root: fs-ivl is the ONLY member, so fs-math's dev-only
  # frankenscipy oracle reference is never resolved (the tools/oracle
  # precedent). The real workspace.lints tables are carried verbatim so
  # `[lints] workspace = true` in every cone crate keeps meaning.
  {
    printf '[workspace]\nmembers = ["crates/fs-ivl"]\nresolver = "2"\n\n'
    python3 - "$REPO_ROOT/Cargo.toml" <<'PYEOF'
import sys
text = open(sys.argv[1]).read()
keep, active = [], False
for line in text.splitlines():
    if line.startswith('['):
        active = line.startswith('[workspace.lints') or line.startswith('[workspace.package')
    if active:
        keep.append(line)
print('\n'.join(keep))
PYEOF
  } > "${dest}/Cargo.toml"
}

# apply_mutation <scratch> <name>  — exact single-site replacement.
apply_mutation() {
  local dest="$1" name="$2"
  python3 - "$dest" "$name" <<'PYEOF'
import sys, pathlib
dest, name = sys.argv[1], sys.argv[2]

# operator -> (relative file, exact old, new). Each `old` must occur
# EXACTLY once; anything else aborts the mutant as inapplicable (source
# drift), which the campaign reports rather than mis-scoring.
MUTATIONS = {
    "flip-next-up": (
        "crates/fs-math/src/lib.rs",
        "pub fn next_up(x: f64) -> f64 {\n    x.next_up()\n}",
        "pub fn next_up(x: f64) -> f64 {\n    x.next_down()\n}",
    ),
    "flip-next-down": (
        "crates/fs-math/src/lib.rs",
        "pub fn next_down(x: f64) -> f64 {\n    x.next_down()\n}",
        "pub fn next_down(x: f64) -> f64 {\n    x.next_up()\n}",
    ),
    "drop-next-up": (
        "crates/fs-math/src/lib.rs",
        "pub fn next_up(x: f64) -> f64 {\n    x.next_up()\n}",
        "pub fn next_up(x: f64) -> f64 {\n    x\n}",
    ),
    "budget-off-by-one": (
        "crates/fs-math/src/det.rs",
        "pub const EXP_ULP_BUDGET: u64 = 3;",
        "pub const EXP_ULP_BUDGET: u64 = 2;",
    ),
    "swap-add-outward": (
        "crates/fs-ivl/src/interval.rs",
        "    /// Addition, outward-rounded.\n    fn add(self, o: Interval) -> Interval {\n        let lo = enclose_rounded_binary(self.lo + o.lo, self.lo.is_finite() && o.lo.is_finite());\n        let hi = enclose_rounded_binary(self.hi + o.hi, self.hi.is_finite() && o.hi.is_finite());",
        "    /// Addition, outward-rounded.\n    fn add(self, o: Interval) -> Interval {\n        let hi = enclose_rounded_binary(self.lo + o.lo, self.lo.is_finite() && o.lo.is_finite());\n        let lo = enclose_rounded_binary(self.hi + o.hi, self.hi.is_finite() && o.hi.is_finite());",
    ),
    "weaken-orient-filter": (
        "crates/fs-ivl/src/predicates.rs",
        "const CCWERRBOUND_A: f64 = (3.0 + 16.0 * EPS) * EPS;",
        "const CCWERRBOUND_A: f64 = (3.0 + 16.0 * EPS) * EPS * 0.015625;",
    ),
}

if name not in MUTATIONS:
    print(f"UNKNOWN-MUTANT {name}", file=sys.stderr)
    sys.exit(3)
rel, old, new = MUTATIONS[name]
path = pathlib.Path(dest) / rel
text = path.read_text()
count = text.count(old)
if count != 1:
    print(f"INAPPLICABLE {name}: {count} matches in {rel}", file=sys.stderr)
    sys.exit(4)
path.write_text(text.replace(old, new))
# Application integrity: the mutated file must now differ.
assert path.read_text() != text
print(rel)
PYEOF
}

MUTANTS=(flip-next-up flip-next-down drop-next-up budget-off-by-one swap-add-outward weaken-orient-filter)
if [[ -n "${ONLY_MUTANT}" ]]; then MUTANTS=("${ONLY_MUTANT}"); fi

# run_layer1 <label> -> exit code; log retained per run. ONE workspace and
# ONE target dir: mutants apply in place and are restored from the real
# repo file afterwards, so successive invocations build incrementally
# instead of paying a cold 7-crate build per mutant.
run_layer1() {
  local label="$1"
  local out="${ARTIFACT_DIR}/${label}.log"
  local rc=0
  (
    cd "${BASE}"
    CARGO_TARGET_DIR="${TARGET_BASE}/shared" cargo test -q -p fs-ivl
  ) > "${out}" 2>&1 || rc=$?
  return ${rc}
}

# restore_file <repo-relative path>: revert the mutated file from the real
# repo (whose integrity hash brackets the campaign).
restore_file() {
  cp "${REPO_ROOT}/$1" "${BASE}/$1"
}

# --------------------------------------------------- positive control
BASE="${SCRATCH}/base"
if [[ ! -d "${BASE}" ]]; then
  assemble_scratch "${BASE}"
fi
if [[ -n "${ONLY_MUTANT}" && -f "${ARTIFACT_DIR}/control.log" ]] \
  && grep -q 'positive control green' "${MATRIX}"; then
  log control "reusing the recorded green positive control"
else
  log phase "positive control: unmutated scratch cone must be green"
  if ! run_layer1 control; then
    log control "POSITIVE CONTROL RED - campaign inconclusive (see control.log)"
    printf '{"mutant":"__control__","outcome":"inconclusive-control-red"}\n' >> "${MATRIX}"
    exit 2
  fi
  log control "positive control green"
fi
if [[ "${CONTROL_ONLY}" -eq 1 ]]; then
  log summary "control-only invocation complete"
  exit 0
fi

# --------------------------------------------------------- the campaign
SURVIVORS=0
INCONCLUSIVE=0
for mutant in "${MUTANTS[@]}"; do
  log phase "mutant ${mutant}"
  MUTATED_FILE="$(apply_mutation "${BASE}" "${mutant}" 2>> "${ARTIFACT_DIR}/apply.log")" || {
    printf '{"mutant":"%s","outcome":"inapplicable-source-drift"}\n' "${mutant}" >> "${MATRIX}"
    INCONCLUSIVE=$((INCONCLUSIVE + 1))
    continue
  }
  start_s=${SECONDS}
  if run_layer1 "${mutant}"; then
    restore_file "${MUTATED_FILE}"
    printf '{"mutant":"%s","outcome":"SURVIVED","layer":"L1-fs-ivl","seconds":%s}\n' \
      "${mutant}" "$((SECONDS - start_s))" >> "${MATRIX}"
    log survivor "${mutant} SURVIVED the fs-ivl battery - file a test-gap bead"
    SURVIVORS=$((SURVIVORS + 1))
  else
    rc=$?
    restore_file "${MUTATED_FILE}"
    # Distinguish red tests (kill) from a broken build (inconclusive):
    # cargo test exits 101 for test failures AND build errors, so classify
    # by the presence of test-harness output.
    if grep -qE 'test result: FAILED|panicked at' "${ARTIFACT_DIR}/${mutant}.log"; then
      # -q cargo output lists failing test names under a `failures:`
      # heading rather than as `test ... FAILED` lines; harvest from both
      # shapes and never let an empty harvest kill the script (pipefail).
      killer="$( (grep -A20 '^failures:$' "${ARTIFACT_DIR}/${mutant}.log" \
                  | grep -oE '^    [a-z0-9_:]+' | head -3 \
                  | sed 's/^    //' | paste -sd, -) 2>/dev/null || true)"
      printf '{"mutant":"%s","outcome":"killed","layer":"L1-fs-ivl","killed_by":"%s","seconds":%s}\n' \
        "${mutant}" "${killer:-panic}" "$((SECONDS - start_s))" >> "${MATRIX}"
      log kill "${mutant} killed by L1 (${killer:-panic})"
    else
      printf '{"mutant":"%s","outcome":"inconclusive-build","seconds":%s}\n' \
        "${mutant}" "$((SECONDS - start_s))" >> "${MATRIX}"
      log inconclusive "${mutant} did not build - inconclusive, never counted as killed"
      INCONCLUSIVE=$((INCONCLUSIVE + 1))
    fi
  fi
done

# ----------------------------------------------------- integrity post
POST_HASH="$(integrity_hash)"
if [[ "${PRE_HASH}" != "${POST_HASH}" ]]; then
  log integrity "REAL TREE CHANGED DURING THE CAMPAIGN (${PRE_HASH} -> ${POST_HASH}); results void"
  printf '{"mutant":"__integrity__","outcome":"void-tree-drift"}\n' >> "${MATRIX}"
  exit 2
fi
log integrity "post-campaign trust-root index hash unchanged"

# ------------------------------------------------------------- summary
TOTAL=${#MUTANTS[@]}
KILLED=$((TOTAL - SURVIVORS - INCONCLUSIVE))
printf '{"mutant":"__summary__","total":%d,"killed":%d,"survivors":%d,"inconclusive":%d}\n' \
  "${TOTAL}" "${KILLED}" "${SURVIVORS}" "${INCONCLUSIVE}" >> "${MATRIX}"
log summary "total=${TOTAL} killed=${KILLED} survivors=${SURVIVORS} inconclusive=${INCONCLUSIVE}"
log summary "kill matrix: ${MATRIX}"

if [[ "${SURVIVORS}" -ne 0 || "${INCONCLUSIVE}" -ne 0 ]]; then
  printf 'FAILED: %d survivor(s), %d inconclusive of %d mutants\n' \
    "${SURVIVORS}" "${INCONCLUSIVE}" "${TOTAL}" >&2
  exit 1
fi
printf 'OK: all %d arithmetic-root mutants killed by the committed batteries\n' "${TOTAL}"
