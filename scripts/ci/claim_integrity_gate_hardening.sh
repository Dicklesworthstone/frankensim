#!/usr/bin/env bash
# claim_integrity_gate_hardening.sh (bead frankensim-extreal-program-f85xj.2.4)
# — seeded-fault drills proving the claim-integrity promotion gate FAILS CLOSED.
#
# A broken gate is worse than a missing one: it manufactures institutional
# confidence. `claim_integrity_gate_drill.sh` shows the gate works on healthy
# inputs; this suite corrupts the inputs on purpose and demands the gate still
# refuse. In the spirit of the Gauntlet's G4 chaos tier.
#
# Every drill asserts BOTH directions — with the fault present the gate must
# misbehave in the specified safe way, and with the fault removed it must
# recover. A drill that only checked the fault could pass against a gate that
# refuses everything, which is not a working gate either.
#
#   A  MID-FLUSH     a truncated issues.jsonl (another agent flushing) must
#                    REFUSE with an actionable diagnostic, never pass.
#   B  TYPO'D LABEL  `claimintegrity` makes a defect invisible to the gate; it
#                    must be surfaced as suspicious rather than silently lost.
#   C  STALE SCOPE   a defect scoped to a crate that no longer exists matches
#                    nothing, so it must be globalized, not disarmed.
#   D  ORDER DRIFT   a promotion recorded before the P0 was filed must still be
#                    blocked; the gate reads state, not chronology.
#   E  PERFORMANCE   the gate must finish within a bounded time on the REAL
#                    ~2000-row inventory, so nobody is tempted to skip it.
#
# Faults are injected into FIXTURE inventories via FSIM_CLAIM_INTEGRITY_BEADS.
# The real .beads store is never written. The registry IS mutated (a promotion
# has to be pending for the gate to weigh anything) and is always put back from
# the committed copy by the exit trap.
#
# Writes a deterministic summary table to
# .doctor/claim-integrity-gate-hardening.txt — no timestamps, stable order, so
# two clean runs produce byte-identical artifacts.
#
# Usage:  scripts/ci/claim_integrity_gate_hardening.sh
# Exit:   0 every drill behaved as specified, 1 a drill misbehaved, 2 setup refused.
set -euo pipefail

cd "$(dirname "$0")/../.."

REGISTRY="capability-maturity.json"
ARTIFACT_DIR=".doctor"
ARTIFACT="$ARTIFACT_DIR/claim-integrity-gate-hardening.txt"
PERF_BUDGET_S=30
CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_check"
export CARGO_TARGET_DIR

if ! git diff --quiet -- "$REGISTRY"; then
  printf '{"lane":"claim-integrity-gate-hardening","verdict":"REFUSED","reason":"%s has uncommitted changes; these drills rewrite it"}\n' "$REGISTRY"
  exit 2
fi

WORK=$(mktemp -d "${TMPDIR:-/tmp}/ci-gate-hardening.XXXXXX")
# Between drills only the registry is reset; the fixture dir must survive until
# the end, so teardown is a separate function bound to the trap.
restore_registry() { git show "HEAD:$REGISTRY" > "$REGISTRY" 2>/dev/null || true; }
teardown() { restore_registry; rm -r -- "$WORK" 2>/dev/null || true; }
trap teardown EXIT INT TERM

# The gate exits nonzero when it blocks or refuses — both are expected outcomes
# here, so the status is swallowed and verdicts are read from the rows.
gate() { sh -c 'cargo run -q -p xtask -- check-claim-integrity' 2>/dev/null || true; }

# Promote a capability so the gate always has something to weigh.
promote() {
  python3 - "$1" "$2" <<'PY'
import json, sys
capability, ref = sys.argv[1], sys.argv[2]
doc = json.load(open("capability-maturity.json"))
for entry in doc["capabilities"]:
    if entry["id"] == capability:
        entry["level"] = "L3"
        entry["evidence"].append({"kind": "lane", "ref": ref})
        break
else:
    sys.exit(f"no such capability: {capability}")
open("capability-maturity.json", "w").write(json.dumps(doc, indent=2) + "\n")
PY
}

# A one-row fixture inventory naming a gating defect scoped to $2.
fixture_row() {
  printf '{"id":"%s","issue_type":"bug","status":"open","labels":["%s","severity:default-path","crate:%s"]}\n' \
    "$1" "${3:-claim-integrity}" "$2"
}

rows=()
fail=0

record() { # name | expected | actual | verdict
  rows+=("$(printf '%-14s | %-46s | %-46s | %s' "$1" "$2" "$3" "$4")")
  if [[ "$4" != "PASS" ]]; then fail=1; fi
  printf '  %-6s %-14s expected=%s actual=%s\n' "$4" "$1" "$2" "$3" >&2
}

printf 'claim-integrity gate hardening drills (seeded faults)\n\n' >&2

# --------------------------------------------------------------- drill A ----
restore_registry; promote "numerics.sparse-assembly" "scripts/ci/quality_lanes.sh"
fixture_row "fx-a" "fs-sparse" > "$WORK/a-clean.jsonl"
# Fault: a row cut off mid-write, as a concurrent flush would leave it.
{ cat "$WORK/a-clean.jsonl"; printf '{"id":"fx-a2","issue_type":"bu'; } > "$WORK/a-torn.jsonl"

out=$(FSIM_CLAIM_INTEGRITY_BEADS="$WORK/a-torn.jsonl" gate)
if printf '%s' "$out" | grep -q 'mid-flush'; then a_fault="refused (mid-flush)"; else a_fault="DID NOT REFUSE"; fi
out=$(FSIM_CLAIM_INTEGRITY_BEADS="$WORK/a-clean.jsonl" gate)
if printf '%s' "$out" | grep -q 'BLOCKED'; then a_ok="blocked on the intact row"; else a_ok="NO VERDICT"; fi
if [[ "$a_fault" == "refused (mid-flush)" && "$a_ok" == "blocked on the intact row" ]]; then v=PASS; else v=FAIL; fi
record "A mid-flush" "refuse on torn row; block when intact" "$a_fault; $a_ok" "$v"

# --------------------------------------------------------------- drill B ----
restore_registry; promote "numerics.sparse-assembly" "scripts/ci/quality_lanes.sh"
fixture_row "fx-b" "fs-sparse" "claimintegrity" > "$WORK/b-typo.jsonl"
fixture_row "fx-b" "fs-sparse" > "$WORK/b-fixed.jsonl"

out=$(FSIM_CLAIM_INTEGRITY_BEADS="$WORK/b-typo.jsonl" gate)
b_seen=$(printf '%s' "$out" | grep -c 'suspicious-label' || true)
b_blocked=$(printf '%s' "$out" | grep -c 'BLOCKED' || true)
# The typo'd defect is genuinely invisible to blocking — that is the hazard —
# so the requirement is that it be REPORTED, not that it block.
if [[ "$b_seen" -ge 1 && "$b_blocked" -eq 0 ]]; then b_fault="surfaced as suspicious, not silent"; else b_fault="TYPO WAS SILENT"; fi
out=$(FSIM_CLAIM_INTEGRITY_BEADS="$WORK/b-fixed.jsonl" gate)
if printf '%s' "$out" | grep -q 'BLOCKED'; then b_ok="blocks once spelled correctly"; else b_ok="NO VERDICT"; fi
if [[ "$b_fault" == "surfaced as suspicious, not silent" && "$b_ok" == "blocks once spelled correctly" ]]; then v=PASS; else v=FAIL; fi
record "B typo label" "warn on near-miss; block when correct" "$b_fault; $b_ok" "$v"

# --------------------------------------------------------------- drill C ----
restore_registry; promote "numerics.sparse-assembly" "scripts/ci/quality_lanes.sh"
# Fault: scoped to a crate that does not exist, so a naive matcher finds nothing.
fixture_row "fx-c" "fs-crate-that-was-renamed-away" > "$WORK/c-stale.jsonl"
fixture_row "fx-c" "fs-geom" > "$WORK/c-real.jsonl"

out=$(FSIM_CLAIM_INTEGRITY_BEADS="$WORK/c-stale.jsonl" gate)
if printf '%s' "$out" | grep -q 'scope-globalized' && printf '%s' "$out" | grep -q 'BLOCKED'; then
  c_fault="globalized and blocked"
else
  c_fault="STALE SCOPE DISARMED THE DEFECT"
fi
# An out-of-scope but REAL crate must not block this capability.
out=$(FSIM_CLAIM_INTEGRITY_BEADS="$WORK/c-real.jsonl" gate)
if printf '%s' "$out" | grep -q 'promotion-admitted'; then c_ok="real out-of-scope crate does not block"; else c_ok="OVER-BLOCKED"; fi
if [[ "$c_fault" == "globalized and blocked" && "$c_ok" == "real out-of-scope crate does not block" ]]; then v=PASS; else v=FAIL; fi
record "C stale scope" "globalize unknown; respect known scope" "$c_fault; $c_ok" "$v"

# --------------------------------------------------------------- drill D ----
restore_registry; promote "numerics.sparse-assembly" "scripts/ci/quality_lanes.sh"
# Fault: the P0 is filed "after" the promotion — it sits last in the file and
# carries a later created_at. The gate must read STATE, not chronology.
{
  printf '{"id":"fx-d-old","issue_type":"bug","status":"closed","created_at":"2020-01-01T00:00:00Z","labels":["claim-integrity","severity:default-path","crate:fs-sparse"]}\n'
  printf '{"id":"fx-d-new","issue_type":"bug","status":"open","created_at":"2099-12-31T23:59:59Z","labels":["claim-integrity","severity:default-path","crate:fs-sparse"]}\n'
} > "$WORK/d-drift.jsonl"
# Recovery: the same late-filed defect, closed.
printf '{"id":"fx-d-new","issue_type":"bug","status":"closed","created_at":"2099-12-31T23:59:59Z","labels":["claim-integrity","severity:default-path","crate:fs-sparse"]}\n' > "$WORK/d-closed.jsonl"

out=$(FSIM_CLAIM_INTEGRITY_BEADS="$WORK/d-drift.jsonl" gate)
if printf '%s' "$out" | grep -q 'fx-d-new'; then d_fault="blocked by the late-filed P0"; else d_fault="CHRONOLOGY DEFEATED THE GATE"; fi
out=$(FSIM_CLAIM_INTEGRITY_BEADS="$WORK/d-closed.jsonl" gate)
if printf '%s' "$out" | grep -q 'promotion-admitted'; then d_ok="admits once it is closed"; else d_ok="STILL BLOCKED"; fi
if [[ "$d_fault" == "blocked by the late-filed P0" && "$d_ok" == "admits once it is closed" ]]; then v=PASS; else v=FAIL; fi
record "D order drift" "block regardless of filing order" "$d_fault; $d_ok" "$v"

# --------------------------------------------------------------- drill E ----
restore_registry; promote "numerics.sparse-assembly" "scripts/ci/quality_lanes.sh"
real_rows=$(wc -l < .beads/issues.jsonl | tr -d ' ')
start=$(date +%s)
out=$(gate)                       # no override: the REAL inventory
elapsed=$(( $(date +%s) - start ))
printf '  (drill E measured %ss over %s rows, budget %ss)\n' "$elapsed" "$real_rows" "$PERF_BUDGET_S" >&2
# The measured duration goes to the log, NOT the artifact: a wall-clock number
# in the table would make the artifact differ run to run, which would falsify
# the determinism claim in its own header.
if [[ "$elapsed" -le "$PERF_BUDGET_S" ]] && printf '%s' "$out" | grep -q '"check":"claim-integrity-gate"'; then
  e_fault="within budget on the real inventory"
  v=PASS
else
  e_fault="EXCEEDED BUDGET (see log for the measured duration)"
  v=FAIL
fi
record "E performance" "<= ${PERF_BUDGET_S}s on the real inventory" "$e_fault" "$v"

# ------------------------------------------------------------------ report --
restore_registry
mkdir -p "$ARTIFACT_DIR"
{
  printf 'claim-integrity gate hardening drills (bead f85xj.2.4)\n'
  printf 'Seeded faults against the promotion gate. Deterministic: no timestamps.\n\n'
  printf '%-14s | %-46s | %-46s | %s\n' "DRILL" "EXPECTED" "ACTUAL" "VERDICT"
  printf -- '---------------+------------------------------------------------+------------------------------------------------+--------\n'
  printf '%s\n' "${rows[@]}"
} > "$ARTIFACT"

printf '\n' >&2
printf 'summary table written to %s\n' "$ARTIFACT" >&2
if [[ "$fail" -ne 0 ]]; then
  printf '{"lane":"claim-integrity-gate-hardening","verdict":"FAIL","drills":5}\n'
  printf 'A hardening drill misbehaved. A gate that does not fail closed is worse\n' >&2
  printf 'than no gate: it manufactures confidence. Do not land this.\n' >&2
  exit 1
fi
printf '{"lane":"claim-integrity-gate-hardening","verdict":"PASS","drills":5}\n'
