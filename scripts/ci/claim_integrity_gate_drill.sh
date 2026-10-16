#!/usr/bin/env bash
# claim_integrity_gate_drill.sh (bead frankensim-extreal-program-f85xj.2.3) —
# the promotion gate, exercised end to end against the REAL repository state.
#
# The gate's rule: a capability may not be promoted while an open
# `severity:default-path` claim-integrity defect is in scope for it
# (docs/CLAIM_INTEGRITY.md, docs/MATURITY_LEVELS.md).
#
# Unit tests pin the decision algebra with fixtures. This lane proves the
# assembled thing works on the live tree, which is a different claim: it drives
# `capability-maturity.json` through real level changes and shows the real
# `.beads/issues.jsonl` inventory blocking and admitting real promotions.
#
# Three drills:
#   A  BLOCK   — promote a capability whose crate scope overlaps an open
#                severity:default-path defect. Must FAIL, naming the bead.
#   B  ADMIT   — promote a capability with no defect in scope. Must PASS.
#   C  DEMOTE  — lower a capability's level. Must PASS even with defects open,
#                because demotion is how the registry stays honest and must
#                never be procedurally harder than promotion.
#
# The registry is mutated in place and ALWAYS put back from the committed copy
# by the exit trap, including on failure or interrupt. Nothing here writes to
# beads; the inventory is read-only input.
#
# Usage:  scripts/ci/claim_integrity_gate_drill.sh
# Exit:   0 all drills behaved as specified, 1 a drill misbehaved, 2 setup refused.
set -euo pipefail

cd "$(dirname "$0")/../.."

REGISTRY="capability-maturity.json"
CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_check"
export CARGO_TARGET_DIR

if ! git diff --quiet -- "$REGISTRY"; then
  printf '{"lane":"claim-integrity-gate-drill","verdict":"REFUSED","reason":"%s has uncommitted changes; the drill rewrites it and would clobber them"}\n' "$REGISTRY"
  exit 2
fi

# Always put the committed registry back, whatever happens.
cleanup() { git show "HEAD:$REGISTRY" > "$REGISTRY" 2>/dev/null || true; }
trap cleanup EXIT INT TERM

# The gate exits nonzero when it BLOCKS, which is drill A's expected outcome,
# so the exit status is swallowed here and the verdict is read from the rows.
gate() { sh -c 'cargo run -q -p xtask -- check-claim-integrity' 2>/dev/null || true; }

# Set one capability's level, and append an evidence ref when promoting so the
# maturity check's own level bars do not mask the gate's verdict.
set_level() {
  python3 - "$1" "$2" "${3:-}" <<'PY'
import json, sys
capability, level, extra = sys.argv[1], sys.argv[2], sys.argv[3]
doc = json.load(open("capability-maturity.json"))
for entry in doc["capabilities"]:
    if entry["id"] == capability:
        entry["level"] = level
        if extra:
            kind, ref = extra.split("=", 1)
            entry["evidence"].append({"kind": kind, "ref": ref})
        break
else:
    sys.exit(f"no such capability: {capability}")
open("capability-maturity.json", "w").write(json.dumps(doc, indent=2) + "\n")
PY
}

fail=0
report() { printf '  %-8s %-26s %s\n' "$1" "$2" "$3" >&2; }

printf 'claim-integrity promotion gate drill\n\n' >&2

# ---------------------------------------------------------------- drill A ---
# fs-wasm is in scope for an open severity:default-path defect, so promoting
# the browser-flagships capability must be refused.
cleanup
set_level "wasm.browser-flagships" "L2" "test=crates/fs-wasm/src/campaigns.rs::flowcert_gates_accuracy_on_convergence"
out=$(gate)
if printf '%s' "$out" | grep -q '"check":"claim-integrity-gate".*wasm.browser-flagships.*BLOCKED'; then
  bead=$(printf '%s' "$out" | grep -o 'defect frankensim-[a-z0-9.-]*' | head -1)
  report "PASS" "A block" "promotion refused, naming ${bead:-<none>}"
else
  report "FAIL" "A block" "an in-scope open P0 did NOT block the promotion"
  printf '%s\n' "$out" >&2
  fail=1
fi

# ---------------------------------------------------------------- drill B ---
# fs-sparse is in scope for no open severity:default-path defect.
cleanup
set_level "numerics.sparse-assembly" "L3" "lane=scripts/ci/quality_lanes.sh"
out=$(gate)
if printf '%s' "$out" | grep -q 'numerics.sparse-assembly.*"verdict":"promotion-admitted"'; then
  report "PASS" "B admit" "promotion with no defect in scope was admitted"
else
  report "FAIL" "B admit" "a clean promotion was not admitted"
  printf '%s\n' "$out" >&2
  fail=1
fi

# ---------------------------------------------------------------- drill C ---
# Demotion must pass even though defects are open.
cleanup
set_level "evidence.colour-algebra" "L1"
out=$(gate)
if printf '%s' "$out" | grep -q '"verdict":"no-promotion"'; then
  report "PASS" "C demote" "demotion passed the gate untouched"
else
  report "FAIL" "C demote" "a demotion was gated; it must never be"
  printf '%s\n' "$out" >&2
  fail=1
fi

cleanup
printf '\n' >&2
if [[ "$fail" -ne 0 ]]; then
  printf '{"lane":"claim-integrity-gate-drill","verdict":"FAIL"}\n'
  printf 'A drill misbehaved. The gate is the control that makes E02 more than\n' >&2
  printf 'bookkeeping; do not land a change that weakens it.\n' >&2
  exit 1
fi
printf '{"lane":"claim-integrity-gate-drill","verdict":"PASS","drills":3}\n'
printf 'All three drills behaved as specified; %s restored to its committed state.\n' "$REGISTRY" >&2
