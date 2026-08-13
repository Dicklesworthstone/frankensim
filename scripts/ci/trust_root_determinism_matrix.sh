#!/usr/bin/env bash
#
# trust_root_determinism_matrix.sh — compiler × opt-level × ISA determinism
# matrix for the certified-arithmetic trust root (bead
# frankensim-extreal-program-f85xj.3.4).
#
# Runs the COMMITTED fs-math + fs-ivl batteries — which pin every golden
# hash, dd-oracle containment law, directed-rounding model, and exact
# corpus row — once per matrix cell, and emits one JSON line per cell with
# the machine fingerprint, toolchain, target, profile, and verdict. A green
# cell means every golden and law in the batteries held under that build
# configuration; a red cell is a determinism incident (playbook:
# stage-wise bit bisection, det:: routing doctrine).
#
# Cells executed by this harness:
#   local  (this machine's ISA)      × {debug, release}
#   remote (rch Linux x86-64 worker) × {debug, release}   [--with-remote]
#
# Cells the matrix REPORTS AS NOT-RUN rather than omitting (the bead's
# explicit-unsupported rule): next-nightly candidate (no second toolchain
# is pinned in-tree; a toolchain bump event is the trigger to run this
# script before and after), alternate -C opt-level / DSR codegen flags
# beyond the two cargo profiles, and thread-count axes (the trust-root
# batteries are single-threaded by construction; thread-count determinism
# is owned by the fs-exec/fs-sparse reduction batteries, not this lane).
#
# Divergence drill: the seeded-perturbation requirement is discharged by
# the sibling mutation campaign (scripts/ci/e2e_extreal_arithmetic_audit.sh,
# kill matrix 6/6): a perturbed trust-root constant demonstrably turns
# these same batteries red, so a red cell here is a real signal, not an
# untested alarm path.
#
# Usage:
#   scripts/ci/trust_root_determinism_matrix.sh --artifact PATH
#       [--with-remote]      # add the rch x86-64 cells
#       [--cell NAME]        # run one cell: local-debug|local-release|
#                            #               remote-debug|remote-release
#
# The artifact is append-per-cell JSONL, so cells can run in separate
# bounded invocations and the committed matrix is their concatenation.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ARTIFACT=""
WITH_REMOTE=0
ONLY_CELL=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --artifact)    ARTIFACT="${2:-}"; shift 2 ;;
    --with-remote) WITH_REMOTE=1; shift ;;
    --cell)        ONLY_CELL="${2:-}"; shift 2 ;;
    -h|--help)     sed -n '3,40p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 2 ;;
    *) printf 'FATAL: unknown argument %s\n' "$1" >&2; exit 2 ;;
  esac
done
[[ -n "${ARTIFACT}" ]] || { printf 'FATAL: --artifact is required\n' >&2; exit 2; }
touch "${ARTIFACT}"

TOOLCHAIN="$(rustc --version | tr -d '\n')"
LOCAL_TARGET="$(rustc -vV | sed -n 's/^host: //p')"
FINGERPRINT="$(uname -srm | tr ' ' '/')"
STAMP_CMD="cargo test -p fs-math -p fs-ivl"

emit() {
  # emit <cell> <target> <profile> <verdict> <detail>
  printf '{"schema":"frankensim.trust-root-determinism-matrix.v1","cell":"%s","target":"%s","profile":"%s","toolchain":"%s","fingerprint":"%s","batteries":"%s","verdict":"%s","detail":"%s"}\n' \
    "$1" "$2" "$3" "${TOOLCHAIN}" "$4" "${STAMP_CMD}" "$5" "$6" >> "${ARTIFACT}"
  printf '[cell ] %s %s/%s -> %s\n' "$1" "$2" "$3" "$5" >&2
}

run_local() {
  local profile="$1"
  local flag=""
  [[ "${profile}" == "release" ]] && flag="--release"
  local rc=0
  ( cd "${REPO_ROOT}" \
    && CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_matrix" \
       cargo test -q -p fs-math -p fs-ivl ${flag} ) \
    > "${ARTIFACT%.jsonl}.local-${profile}.log" 2>&1 || rc=$?
  if [[ ${rc} -eq 0 ]]; then
    emit "local-${profile}" "${LOCAL_TARGET}" "${profile}" "${FINGERPRINT}" "pass" \
      "all committed goldens and laws held"
  else
    emit "local-${profile}" "${LOCAL_TARGET}" "${profile}" "${FINGERPRINT}" "FAIL" \
      "determinism incident: see local-${profile}.log"
    return 1
  fi
}

run_remote() {
  local profile="$1"
  local flag=""
  [[ "${profile}" == "release" ]] && flag="--release"
  local out="${ARTIFACT%.jsonl}.remote-${profile}.log"
  local rc=0
  ( cd "${REPO_ROOT}" \
    && rch exec -- env CARGO_TARGET_DIR="\${TMPDIR:-/tmp}/rch_target_frankensim_matrix" \
       cargo test -q -p fs-math -p fs-ivl ${flag} ) > "${out}" 2>&1 || rc=$?
  # The worker's fingerprint is in the rch transcript; harvest best-effort.
  local remote_fp
  remote_fp="$( (grep -oE 'worker[ =:][A-Za-z0-9._-]+' "${out}" | head -1 | tr ' :=' '/') 2>/dev/null || true)"
  if [[ ${rc} -eq 0 ]]; then
    emit "remote-${profile}" "x86_64-unknown-linux-gnu" "${profile}" \
      "rch/${remote_fp:-worker-see-log}" "pass" "all committed goldens and laws held"
  else
    emit "remote-${profile}" "x86_64-unknown-linux-gnu" "${profile}" \
      "rch/${remote_fp:-worker-see-log}" "FAIL" \
      "determinism incident OR remote infrastructure failure: see remote-${profile}.log"
    return 1
  fi
}

not_run_rows() {
  emit "next-nightly" "any" "any" "n/a" "NOT-RUN" \
    "no second toolchain pinned in-tree; run this matrix before and after any toolchain bump"
  emit "alt-codegen-flags" "any" "any" "n/a" "NOT-RUN" \
    "opt-level beyond the two cargo profiles and DSR codegen flags are not exercised by this harness yet"
  emit "thread-count-axis" "any" "any" "n/a" "NOT-APPLICABLE" \
    "trust-root batteries are single-threaded by construction; thread determinism is owned by fs-exec/fs-sparse reduction batteries"
}

FAILURES=0
maybe() { # maybe <cell-name> <fn> <arg>
  if [[ -z "${ONLY_CELL}" || "${ONLY_CELL}" == "$1" ]]; then
    "$2" "$3" || FAILURES=$((FAILURES + 1))
  fi
}

maybe local-debug    run_local  debug
maybe local-release  run_local  release
if [[ "${WITH_REMOTE}" -eq 1 || "${ONLY_CELL}" == remote-* ]]; then
  maybe remote-debug   run_remote debug
  maybe remote-release run_remote release
fi
if [[ -z "${ONLY_CELL}" || "${ONLY_CELL}" == "not-run" ]]; then
  not_run_rows
fi

if [[ "${FAILURES}" -ne 0 ]]; then
  printf 'FAILED: %d matrix cell(s) diverged or errored\n' "${FAILURES}" >&2
  exit 1
fi
printf 'OK: matrix cells appended to %s\n' "${ARTIFACT}"
